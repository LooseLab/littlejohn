"""
Ray Core implementation of the LittleJohn workflow engine
- Specialized per-queue actors (preprocessing, bed_conversion, mgmt, cnv, target, fusion, classification, slow)
- Central Coordinator actor for dedup (1 running + 1 pending per (sample_id, job_type)), triggers, stats
- Handlers registered for real LittleJohn job types, with resource hints per job
- tqdm-based live monitor similar to the original implementation

Usage (example):

    pip install "ray>=2.30" tqdm watchdog

    python ray_littlejohn_core.py \
        --plan preprocessing:preprocessing \
        --plan bed_conversion:bed_conversion \
        --plan mgmt:mgmt \
        --plan sturgeon:sturgeon \
        --paths /data/incoming

Replace/extend the plan and paths as needed. If you already have a Watchdog-based
file watcher, call submit_jobs() on the Coordinator in your event handlers instead.
"""

from __future__ import annotations

import os
import time
import asyncio
import argparse
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple, Callable
import inspect
from pathlib import Path
import itertools

import ray
from tqdm import tqdm
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer

# Optional GUI hook integration
try:
    from littlejohn.gui_launcher import (
        send_gui_update as _gui_send_update,
        UpdateType as _GUIUpdateType,
        launch_gui as _gui_launch,
    )
except Exception as e:
    raise Exception(f"GUI not available: {e}")

    def _gui_send_update(*args, **kwargs):
        return None

    _GUIUpdateType = None

    def _gui_launch(*args, **kwargs):
        return None


# Import handlers
try:
    from littlejohn.analysis.bam_preprocessor import (
        bam_preprocessing_handler as _preprocessing_handler,
    )
except Exception:
    _preprocessing_handler = None

try:
    from littlejohn.analysis.bed_conversion import (
        bed_conversion_handler as _bed_conversion_handler,
    )
except Exception:
    _bed_conversion_handler = None

try:
    from littlejohn.analysis.mgmt_analysis import mgmt_handler as _mgmt_handler
except Exception:
    _mgmt_handler = None

try:
    from littlejohn.analysis.cnv_analysis import cnv_handler as _cnv_handler
except Exception:
    _cnv_handler = None

try:
    from littlejohn.analysis.target_analysis import target_handler as _target_handler
except Exception:
    _target_handler = None

try:
    from littlejohn.analysis.fusion_analysis import fusion_handler as _fusion_handler
except Exception:
    _fusion_handler = None

try:
    from littlejohn.analysis.sturgeon_analysis import (
        sturgeon_handler as _sturgeon_handler,
    )
except Exception:
    _sturgeon_handler = None

try:
    from littlejohn.analysis.nanodx_analysis import nanodx_handler as _nanodx_handler
except Exception:
    _nanodx_handler = None

try:
    from littlejohn.analysis.nanodx_analysis import (
        pannanodx_handler as _pannanodx_handler,
    )
except Exception:
    _pannanodx_handler = None

try:
    from littlejohn.analysis.random_forest_analysis import (
        random_forest_handler as _random_forest_handler,
    )
except Exception:
    _random_forest_handler = None

# Optional logging helper
try:
    from littlejohn.logging_config import (
        get_job_logger as _get_job_logger,
        configure_logging as _configure_logging,
    )
except Exception:

    def _get_job_logger(job_id: str, job_type: str, filepath: str):
        class _L:
            def info(self, *a, **k):
                pass

            def debug(self, *a, **k):
                pass

            def error(self, *a, **k):
                pass

        return _L()

    def _configure_logging(
        global_level: str = "INFO", job_levels: Optional[Dict[str, str]] = None
    ):
        return None


GLOBAL_LOG_LEVEL: str = "INFO"


# ---------- Shared DTOs ----------
@dataclass
class WorkflowContext:
    filepath: str
    metadata: Dict[str, Any] = field(default_factory=dict)
    results: Dict[str, Any] = field(default_factory=dict)
    history: List[str] = field(default_factory=list)
    errors: List[Dict[str, Any]] = field(default_factory=list)

    def add_metadata(self, key: str, value: Any) -> None:
        self.metadata[key] = value

    def add_result(self, job_type: str, result: Any) -> None:
        self.results[job_type] = result
        self.history.append(job_type)

    def add_error(self, job_type: str, error: str) -> None:
        self.errors.append(
            {"job_type": job_type, "error": error, "timestamp": time.time()}
        )

    def get_sample_id(self) -> str:
        bam_md = self.metadata.get("bam_metadata", {})
        return bam_md.get("sample_id", "unknown")


@dataclass
class Job:
    job_id: int
    job_type: str
    origin: str
    workflow: List[str]
    step: int
    context: WorkflowContext

    def next_job(self) -> Optional["Job"]:
        if self.step + 1 >= len(self.workflow):
            return None
        next_step = self.workflow[self.step + 1]
        if ":" not in next_step:
            raise ValueError(f"Invalid workflow step: {next_step}")
        queue_type, next_type = next_step.split(":", 1)
        return Job(
            self.job_id,
            next_type,
            queue_type,
            self.workflow,
            self.step + 1,
            self.context,
        )


# ---------- Queue mapping & triggers ----------
QUEUE_TO_TYPES: Dict[str, Set[str]] = {
    "preprocessing": {"preprocessing"},
    "bed_conversion": {"bed_conversion"},
    "mgmt": {"mgmt"},
    "cnv": {"cnv"},
    "target": {"target"},
    "fusion": {"fusion"},
    "classification": {"sturgeon", "nanodx", "pannanodx"},
    "slow": {"random_forest"},
}

TRIGGERS: Dict[str, List[str]] = {
    # preprocessing -> analyses
    "preprocessing": ["bed_conversion", "mgmt", "cnv", "target", "fusion"],
    # bed_conversion -> classifiers
    "bed_conversion": ["sturgeon", "nanodx", "pannanodx", "random_forest"],
}

# Per-sample de-duplication to avoid output races. Ensure only one job of these types
# runs concurrently per sample (max 1 running + 1 pending).
# Only deduplicate classifiers/slow per sample; analysis/bed_conversion should queue, not skip
DEDUP_TYPES: Set[str] = {"sturgeon", "nanodx", "pannanodx", "random_forest"}

# Per-sample serialization by type to avoid races on shared per-sample outputs
# Disabled: treat CNV like other analysis jobs (mgmt/target/fusion)
SERIALIZE_BY_TYPE_PER_SAMPLE: Set[str] = set()

# Classification job types (single global pipeline per type)
CLASSIFICATION_TYPES: Set[str] = {"sturgeon", "nanodx", "pannanodx", "random_forest"}

# Job types that must be serialized per sample (no overlap across these types)
## Removed per-sample cross-type serialization to allow one job per type globally

RESOURCE_HINTS: Dict[str, Dict[str, Any]] = {
    # Tune these to your cluster
    "preprocessing": {"num_cpus": 1},
    "bed_conversion": {"num_cpus": 1},
    "mgmt": {"num_cpus": 1},
    "cnv": {"num_cpus": 1},
    "target": {"num_cpus": 1},
    "fusion": {"num_cpus": 1},
    # Classifiers do not require GPU by default (CPU-only)
    "sturgeon": {"num_cpus": 1},
    "nanodx": {"num_cpus": 1},
    "pannanodx": {"num_cpus": 1},
    "random_forest": {"num_cpus": 1},
}

# ---------- Utilities ----------
_job_id_counter = itertools.count(1000)


def job_queue_of(job_type: str) -> str:
    for q, types in QUEUE_TO_TYPES.items():
        if job_type in types:
            return q
    return "slow"


# ---------- Real handler wrappers as Ray tasks ----------
# Each wrapper returns an updated WorkflowContext.

NEEDS_WORK_DIR: Set[str] = {
    "bed_conversion",
    "mgmt",
    "cnv",
    "target",
    "fusion",
    "sturgeon",
    "nanodx",
    "pannanodx",
    "random_forest",
}


def _wrap_real_handler(
    py_handler: Optional[Callable[[Job], None]], job_type: str
) -> Callable[[Job], WorkflowContext]:
    def _impl(job: Job) -> WorkflowContext:
        # Ensure per-process logging honors global level
        try:
            _configure_logging(global_level=GLOBAL_LOG_LEVEL)
        except Exception:
            pass
        logger = _get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
        try:
            if py_handler is None:
                raise RuntimeError(
                    f"No implementation available for job_type='{job_type}'"
                )
            # Call user's handler; if they mutate job.context in-place, we just return it
            work_dir = job.context.metadata.get("work_dir")
            # Ensure output directories exist when a work_dir is provided
            if work_dir:
                try:
                    os.makedirs(work_dir, exist_ok=True)
                except Exception:
                    pass
                try:
                    sid = (
                        job.context.get_sample_id()
                        if hasattr(job.context, "get_sample_id")
                        else None
                    )
                    if sid and sid != "unknown":
                        os.makedirs(os.path.join(work_dir, sid), exist_ok=True)
                except Exception:
                    pass
            try:
                sig = inspect.signature(py_handler)
                if (
                    "work_dir" in sig.parameters
                    and work_dir
                    and (job_type in NEEDS_WORK_DIR)
                ):
                    py_handler(job, work_dir=work_dir)
                else:
                    py_handler(job)
            except Exception:
                # Fallback to legacy call
                py_handler(job)
            # Log completion or skip based on handler-set context
            skipped = False
            try:
                skipped = bool(job.context.metadata.get("skip_downstream", False))
            except Exception:
                skipped = False
            if skipped and job_type == "preprocessing":
                logger.warning(
                    f"{job_type} skipped (deliberate): {os.path.basename(job.context.filepath)}"
                )
                try:
                    # Inform GUI/CLI via GUI hook if available
                    _gui_send_update(
                        _GUIUpdateType.LOG_MESSAGE,
                        {
                            "message": f"Skipped (too many reads): {os.path.basename(job.context.filepath)}",
                            "level": "WARNING",
                        },
                        priority=0,
                    )
                except Exception:
                    pass
            else:
                logger.info(f"{job_type} completed")
        except Exception as e:
            job.context.add_error(job_type, str(e))
            logger.error(f"{job_type} failed: {e}")
            # still return context with error recorded
        finally:
            # Mark result if user handler didn't
            if job_type not in job.context.results:
                job.context.add_result(job_type, f"{job_type}_ok")
        return job.context

    return _impl


# Build remote versions with no fixed options; Pool will apply .options(**RESOURCE_HINTS[job_type])


# Real wrappers
preprocessing_handler_remote = ray.remote(
    _wrap_real_handler(_preprocessing_handler, "preprocessing")
)
bed_conversion_handler_remote = ray.remote(
    _wrap_real_handler(_bed_conversion_handler, "bed_conversion")
)
mgmt_handler_remote = ray.remote(_wrap_real_handler(_mgmt_handler, "mgmt"))
cnv_handler_remote = ray.remote(_wrap_real_handler(_cnv_handler, "cnv"))
target_handler_remote = ray.remote(_wrap_real_handler(_target_handler, "target"))
fusion_handler_remote = ray.remote(_wrap_real_handler(_fusion_handler, "fusion"))
sturgeon_handler_remote = ray.remote(_wrap_real_handler(_sturgeon_handler, "sturgeon"))
nanodx_handler_remote = ray.remote(_wrap_real_handler(_nanodx_handler, "nanodx"))
pannanodx_handler_remote = ray.remote(
    _wrap_real_handler(_pannanodx_handler, "pannanodx")
)
random_forest_handler_remote = ray.remote(
    _wrap_real_handler(_random_forest_handler, "random_forest")
)


# ---------- TypeProcessor & Coordinator Actors ----------
@ray.remote
class TypeProcessor:
    def __init__(self, job_type: str, remote_func, resource_options: Dict[str, Any]):
        self.job_type = job_type
        self.remote_func = remote_func
        self.resource_options = resource_options or {}

    def update_handler(self, remote_func, resource_options: Dict[str, Any]):
        self.remote_func = remote_func
        self.resource_options = resource_options or {}

    async def process(self, job: Job):
        # Submit the real handler task and await the result to avoid nested ObjectRefs
        ref = self.remote_func.options(**(self.resource_options or {})).remote(job)
        try:
            ctx = await ref
        except Exception:
            # Propagate error to coordinator; it will handle ok=False
            raise
        return ctx


@ray.remote
class Coordinator:
    def __init__(self, analysis_workers: int = 1, preset: Optional[str] = None):
        # dedup maps
        self.pending: Dict[Tuple[str, str], int] = {}
        self.running: Dict[Tuple[str, str], int] = {}
        # stats
        self.completed: List[int] = []
        self.failed: List[int] = []
        self.completed_by_type: Dict[str, int] = {}
        self.failed_by_type: Dict[str, int] = {}
        self.submitted_by_type: Dict[str, int] = {}
        self.total_enqueued = 0
        self.total_skipped = 0
        self.active: Dict[int, Dict[str, Any]] = {}

        self.analysis_workers = analysis_workers

        # Per-sample per-type serialization tracking
        self.running_by_type_sample: Dict[Tuple[str, str], int] = {}
        self.pending_by_type_sample: Dict[Tuple[str, str], int] = {}
        self.waiting_by_type_sample: Dict[Tuple[str, str], List[Job]] = {}

        # Per-sample one-shot scheduling for classifiers/slow
        self.scheduled_by_sample: Dict[str, Set[str]] = {}

        # Global one-in-flight (plus at most one waiting) per classification type
        self.classif_pending_by_type: Dict[str, int] = {}
        self.classif_waiting_by_type: Dict[str, Optional[Job]] = {}

        # create one TypeProcessor actor per job type
        self.processors: Dict[str, Any] = {}
        # inflight tracking: ObjectRef -> job
        self._inflight: Dict[Any, Job] = {}
        # backpressure: limit outstanding submissions per job type to avoid
        # massive pending tasks and worker explosion
        self.max_inflight_per_type: int = 64
        self.inflight_by_type: Dict[str, int] = {}
        self.waiting_by_type_global: Dict[str, List[Job]] = {}
        # per-sample tracking for GUI
        self.samples_by_id: Dict[str, Dict[str, Any]] = {}
        # defer setup/registration to async setup()

        # Preset controls how processing actors are created and their concurrency
        # None -> legacy per-type actors; otherwise use grouped Pool actors
        self.preset: Optional[str] = preset
        self.using_pools: bool = False
        # Keep references to any CPU limiter actors
        self._cpu_limiters: List[Any] = []

    # ----- handler registration API -----
    async def register_handler(self, job_type: str, remote_func) -> None:
        opts = RESOURCE_HINTS.get(job_type, {})
        proc = self.processors.get(job_type)
        if proc is not None:
            await proc.update_handler.remote(remote_func, opts)

    async def setup(self) -> None:
        """Async setup to register default handlers and start drain loop."""
        # remote function handles (job_type -> remote function)
        registrations = [
            ("preprocessing", preprocessing_handler_remote),
            ("bed_conversion", bed_conversion_handler_remote),
            ("mgmt", mgmt_handler_remote),
            ("cnv", cnv_handler_remote),
            ("target", target_handler_remote),
            ("fusion", fusion_handler_remote),
            ("sturgeon", sturgeon_handler_remote),
            ("nanodx", nanodx_handler_remote),
            ("pannanodx", pannanodx_handler_remote),
            ("random_forest", random_forest_handler_remote),
        ]

        preset = (self.preset or "").lower().strip()
        if preset in {"p2i", "standard"}:
            # Optionally reserve CPUs to cap global availability
            try:
                total_cpus = float((ray.cluster_resources() or {}).get("CPU", 0))
            except Exception:
                total_cpus = 0.0
            desired_cap = 2.0 if preset == "p2i" else 4.0
            reserve = max(0.0, total_cpus - desired_cap)
            try:
                if reserve >= 1.0:
                    # Reserve floor(reserve) CPUs via a dummy actor
                    reserve_int = int(reserve)

                    @ray.remote
                    class _CpuLimiter:
                        async def ping(self):
                            return True

                    try:
                        limiter = _CpuLimiter.options(num_cpus=reserve_int).remote()
                    except Exception:
                        limiter = _CpuLimiter.remote()
                    # Keep reference so it's not GC'd
                    self._cpu_limiters.append(limiter)
            except Exception:
                pass

            # Grouped Pool actors
            groups = {
                "prep": ["preprocessing", "bed_conversion"],
                "analysis": ["mgmt", "cnv", "target", "fusion"],
                "classif": ["sturgeon", "nanodx", "pannanodx"],
                "rf": ["random_forest"],
            }

            # Determine concurrency per pool based on preset
            def pool_parallel(name: str) -> int:
                if preset == "p2i":
                    # Concurrency 1 everywhere
                    return 1
                # standard preset
                if name == "analysis":
                    return max(1, int(self.analysis_workers))
                return 1

            # Create pools and register handlers
            pools: Dict[str, Any] = {}
            for pool_name, job_types in groups.items():
                par = pool_parallel(pool_name)
                actor_name = f"pool_{pool_name}"
                try:
                    pool = Pool.options(
                        name=actor_name, max_concurrency=max(1, int(par)), num_cpus=0
                    ).remote(pool_name, par)
                except Exception:
                    try:
                        # Some Ray versions don't allow num_cpus=0; use a tiny fractional CPU to avoid reservation deadlock
                        pool = Pool.options(
                            name=actor_name,
                            max_concurrency=max(1, int(par)),
                            num_cpus=0.001,
                        ).remote(pool_name, par)
                    except Exception:
                        pool = Pool.options(max_concurrency=max(1, int(par))).remote(
                            pool_name, par
                        )
                pools[pool_name] = pool
                # Wire coordinator callback
                try:
                    # Fire-and-forget to avoid blocking if actor cannot start immediately
                    pool.set_coordinator_name.remote("littlejohn_coordinator")
                except Exception:
                    pass
                # Register handlers on the pool and map processors
                for jt, rf in registrations:
                    if jt in job_types:
                        opts = RESOURCE_HINTS.get(jt, {})
                        try:
                            await pool.register_handler.remote(jt, rf, opts)
                        except Exception:
                            pass
                        self.processors[jt] = pool
            self.using_pools = True
        elif preset == "high":
            # High-performance: per-type actors (legacy model), actors reserve 1 CPU by default
            for jt, rf in registrations:
                opts = RESOURCE_HINTS.get(jt, {})
                max_conc = 1
                if jt in {"mgmt", "cnv", "target", "fusion"}:
                    max_conc = max(1, int(self.analysis_workers))
                name = f"typeproc_{jt}"
                try:
                    proc = TypeProcessor.options(
                        name=name, max_concurrency=max_conc
                    ).remote(jt, rf, opts)
                except Exception:
                    proc = TypeProcessor.options(max_concurrency=max_conc).remote(
                        jt, rf, opts
                    )
                self.processors[jt] = proc
            self.using_pools = False
        else:
            # Legacy: one TypeProcessor per job type
            for jt, rf in registrations:
                opts = RESOURCE_HINTS.get(jt, {})
                # Set actor concurrency via .options on the actor itself
                max_conc = 1
                if jt in {"mgmt", "cnv", "target", "fusion"}:
                    max_conc = max(1, int(self.analysis_workers))
                name = f"typeproc_{jt}"
                try:
                    proc = TypeProcessor.options(
                        name=name, max_concurrency=max_conc
                    ).remote(jt, rf, opts)
                except Exception:
                    proc = TypeProcessor.options(max_concurrency=max_conc).remote(
                        jt, rf, opts
                    )
                self.processors[jt] = proc
            self.using_pools = False
        # start drain loop to observe completions
        try:
            asyncio.create_task(self._drain_loop())
        except Exception:
            pass

    # ----- submission & lifecycle -----
    async def submit_jobs(self, jobs: List[Job]) -> None:
        for job in jobs:
            sample_id = job.context.get_sample_id()
            if job.job_type in DEDUP_TYPES:
                key = (sample_id, job.job_type)
                run = self.running.get(key, 0)
                pend = self.pending.get(key, 0)
                if run + pend >= 2 or pend >= 1:
                    self.total_skipped += 1
                    continue
                self.pending[key] = pend + 1
            # Per-sample per-type serialization; queue instead of overlapping
            if job.job_type in SERIALIZE_BY_TYPE_PER_SAMPLE:
                sid = sample_id or "unknown"
                key_ts = (job.job_type, sid)
                ra = self.running_by_type_sample.get(key_ts, 0)
                pa = self.pending_by_type_sample.get(key_ts, 0)
                if (ra + pa) > 0:
                    self.waiting_by_type_sample.setdefault(key_ts, []).append(job)
                    self.pending_by_type_sample[key_ts] = pa + 1
                    continue
                else:
                    self.running_by_type_sample[key_ts] = ra + 1

            # submit to dedicated job-type processor
            proc = self.processors.get(job.job_type)
            if proc is None:
                continue
            # Backpressure: if too many outstanding for this type, queue locally
            inflight_for_type = int(self.inflight_by_type.get(job.job_type, 0))
            if inflight_for_type >= self.max_inflight_per_type:
                self.waiting_by_type_global.setdefault(job.job_type, []).append(job)
                # do not count as submitted yet; GUI totals reflect actual submissions
                continue
            # Mark submission only when actually dispatching to a processor
            q = job_queue_of(job.job_type)
            self.total_enqueued += 1
            try:
                self.submitted_by_type[job.job_type] = (
                    self.submitted_by_type.get(job.job_type, 0) + 1
                )
            except Exception:
                pass
            self.active[job.job_id] = {
                "job_type": job.job_type,
                "filepath": job.context.filepath,
                "queue": q,
                "start_time": time.time(),
            }
            # Update per-sample totals/active only when a real sample_id is known
            try:
                sid = sample_id or "unknown"
                if sid != "unknown":
                    ent = self.samples_by_id.get(sid)
                    if ent is None:
                        ent = {
                            "sample_id": sid,
                            "active_jobs": 0,
                            "total_jobs": 0,
                            "completed_jobs": 0,
                            "failed_jobs": 0,
                            "job_types": set(),
                            "last_seen": time.time(),
                        }
                        self.samples_by_id[sid] = ent
                    ent["total_jobs"] += 1
                    ent["active_jobs"] += 1
                    try:
                        ent["job_types"].add(job.job_type)
                    except Exception:
                        pass
                    ent["last_seen"] = time.time()
            except Exception:
                pass
            if getattr(self, "using_pools", False):
                try:
                    proc.enqueue.remote(job)
                except Exception:
                    pass
            else:
                ref = proc.process.remote(job)
                self._inflight[ref] = job
            self.inflight_by_type[job.job_type] = inflight_for_type + 1

    async def _on_finish(
        self, job: Job, ok: bool, ctx: WorkflowContext, err: Optional[str]
    ):
        # capture start info for duration logging, then update dedup maps
        active_info = self.active.pop(job.job_id, None)
        start_time = None
        queue_name = None
        filepath = None
        try:
            if isinstance(active_info, dict):
                start_time = active_info.get("start_time")
                queue_name = active_info.get("queue")
                filepath = active_info.get("filepath")
        except Exception:
            pass
        # update dedup maps
        if job.job_type in DEDUP_TYPES:
            key = (job.context.get_sample_id(), job.job_type)
            self.running[key] = max(0, self.running.get(key, 0) - 1)
            if self.running[key] == 0:
                self.running.pop(key, None)

        # update stats
        # Detect deliberate skip (large BAM)
        is_skipped = False
        try:
            is_skipped = (job.job_type == "preprocessing") and bool(
                ctx.metadata.get("skip_downstream", False)
            )
        except Exception:
            is_skipped = False
        status = (
            "completed"
            if (ok and not is_skipped)
            else ("skipped" if is_skipped else "failed")
        )

        # Append per-job duration entry
        try:
            end_time = time.time()
            duration = (end_time - float(start_time)) if start_time else None
            sample_id = None
            try:
                sample_id = (
                    ctx.get_sample_id() if hasattr(ctx, "get_sample_id") else "unknown"
                )
            except Exception:
                sample_id = "unknown"
            work_dir = None
            try:
                work_dir = ctx.metadata.get("work_dir") or job.context.metadata.get(
                    "work_dir"
                )
            except Exception:
                work_dir = None
            out_path = (
                os.path.join(work_dir, "job_durations.csv")
                if work_dir
                else "job_durations.csv"
            )
            try:
                parent = os.path.dirname(out_path)
                if parent:
                    os.makedirs(parent, exist_ok=True)
            except Exception:
                pass
            new_file = not os.path.exists(out_path)
            with open(out_path, "a", encoding="utf-8") as fh:
                if new_file:
                    fh.write(
                        "timestamp_iso,sample_id,job_id,job_type,queue,filepath,duration_seconds,status,error\n"
                    )
                ts_iso = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(end_time))
                dur_str = (
                    f"{duration:.3f}" if isinstance(duration, (int, float)) else ""
                )
                fp = filepath or (ctx.filepath if hasattr(ctx, "filepath") else "")
                err_short = (err or "").replace("\n", " ").replace(",", " ")[:500]
                fh.write(
                    f"{ts_iso},{sample_id},{job.job_id},{job.job_type},{queue_name or ''},{os.path.basename(fp)},{dur_str},{status},{err_short}\n"
                )
        except Exception:
            pass
        if ok:
            if not is_skipped:
                self.completed.append(job.job_id)
                self.completed_by_type[job.job_type] = (
                    self.completed_by_type.get(job.job_type, 0) + 1
                )
            else:
                # Track skip for stats
                self.total_skipped += 1
            # trigger downstream (preferred). For bed_conversion, schedule
            # classifiers/slow once per sample instead of per-file.
            requested = job.context.metadata.get("original_workflow", [])
            req_types = [s.split(":", 1)[1] for s in requested] if requested else None
            if job.job_type == "bed_conversion":
                sample_id = (
                    ctx.get_sample_id() if hasattr(ctx, "get_sample_id") else "unknown"
                )
                # For each requested classification/slow type, allow at most
                # one submitted at a time globally; keep at most one waiting.
                for t in TRIGGERS.get("bed_conversion", []):
                    if req_types is not None and t not in req_types:
                        continue
                    if t not in CLASSIFICATION_TYPES:
                        continue
                    busy = self.classif_pending_by_type.get(t, 0) >= 1
                    if not busy:
                        q = job_queue_of(t)
                        j = Job(next(_job_id_counter), t, q, [f"{q}:{t}"], 0, ctx)
                        # submit directly and mark pending
                        proc = self.processors.get(t)
                        if proc is not None:
                            if getattr(self, "using_pools", False):
                                try:
                                    proc.enqueue.remote(j)
                                except Exception:
                                    pass
                            else:
                                ref = proc.process.remote(j)
                                self._inflight[ref] = j
                            self.classif_pending_by_type[t] = (
                                self.classif_pending_by_type.get(t, 0) + 1
                            )
                            # Record submission for totals used by GUI
                            try:
                                self.submitted_by_type[j.job_type] = (
                                    self.submitted_by_type.get(j.job_type, 0) + 1
                                )
                            except Exception:
                                pass
                            self.total_enqueued += 1
                            # Update per-sample aggregate immediately for GUI samples view
                            try:
                                sid_local = (
                                    ctx.get_sample_id()
                                    if hasattr(ctx, "get_sample_id")
                                    else "unknown"
                                )
                                if sid_local != "unknown":
                                    ent_local = self.samples_by_id.get(sid_local)
                                    if ent_local is None:
                                        ent_local = {
                                            "sample_id": sid_local,
                                            "active_jobs": 0,
                                            "total_jobs": 0,
                                            "completed_jobs": 0,
                                            "failed_jobs": 0,
                                            "job_types": set(),
                                            "last_seen": time.time(),
                                        }
                                        self.samples_by_id[sid_local] = ent_local
                                    ent_local["total_jobs"] += 1
                                    ent_local["active_jobs"] += 1
                                    try:
                                        ent_local["job_types"].add(j.job_type)
                                    except Exception:
                                        pass
                                    ent_local["last_seen"] = time.time()
                            except Exception:
                                pass
                            self.active[j.job_id] = {
                                "job_type": j.job_type,
                                "filepath": j.context.filepath,
                                "queue": q,
                                "start_time": time.time(),
                            }
                    else:
                        # store single waiting job if none recorded yet
                        if self.classif_waiting_by_type.get(t) is None:
                            q = job_queue_of(t)
                            self.classif_waiting_by_type[t] = Job(
                                next(_job_id_counter), t, q, [f"{q}:{t}"], 0, ctx
                            )
            else:
                triggered_jobs: List[Job] = []
                for t in TRIGGERS.get(job.job_type, []):
                    if (req_types is None) or (t in req_types):
                        q = job_queue_of(t)
                        triggered_jobs.append(
                            Job(next(_job_id_counter), t, q, [f"{q}:{t}"], 0, ctx)
                        )
                if (not is_skipped) and triggered_jobs:
                    await self.submit_jobs(triggered_jobs)
                else:
                    # advance linear chain if any
                    nxt = job.next_job()
                    if (not is_skipped) and nxt:
                        await self.submit_jobs([nxt])
        else:
            # Only record as failed if not a deliberate skip
            if not is_skipped:
                self.failed.append(job.job_id)
                self.failed_by_type[job.job_type] = (
                    self.failed_by_type.get(job.job_type, 0) + 1
                )

        # per-sample finish updates
        try:
            sid2 = (
                job.context.get_sample_id()
                if hasattr(job.context, "get_sample_id")
                else "unknown"
            )
            ent2 = self.samples_by_id.get(sid2)
            if ent2 is not None:
                if ent2.get("active_jobs", 0) > 0:
                    ent2["active_jobs"] -= 1
                # Do not count deliberately skipped preprocessing as completed or failed
                is_skipped = False
                try:
                    is_skipped = (job.job_type == "preprocessing") and bool(
                        job.context.metadata.get("skip_downstream", False)
                    )
                except Exception:
                    is_skipped = False
                if not is_skipped:
                    if ok:
                        ent2["completed_jobs"] = ent2.get("completed_jobs", 0) + 1
                    else:
                        ent2["failed_jobs"] = ent2.get("failed_jobs", 0) + 1
                ent2["last_seen"] = time.time()
        except Exception:
            pass

        # Promote next waiting job for this (type, sample) if any (CNV-only if enabled)
        if job.job_type in SERIALIZE_BY_TYPE_PER_SAMPLE:
            sid = job.context.get_sample_id() or "unknown"
            key_ts = (job.job_type, sid)
            if key_ts in self.running_by_type_sample:
                self.running_by_type_sample[key_ts] = max(
                    0, self.running_by_type_sample.get(key_ts, 0) - 1
                )
                if self.running_by_type_sample[key_ts] == 0:
                    self.running_by_type_sample.pop(key_ts, None)
            queue_list = self.waiting_by_type_sample.get(key_ts, [])
            if queue_list:
                nxt_job = queue_list.pop(0)
                if not queue_list:
                    self.waiting_by_type_sample.pop(key_ts, None)
                self.pending_by_type_sample[key_ts] = max(
                    0, self.pending_by_type_sample.get(key_ts, 0) - 1
                )
                if self.pending_by_type_sample[key_ts] == 0:
                    self.pending_by_type_sample.pop(key_ts, None)
                self.running_by_type_sample[key_ts] = (
                    self.running_by_type_sample.get(key_ts, 0) + 1
                )
                proc2 = self.processors.get(nxt_job.job_type)
                if proc2 is not None:
                    ref2 = proc2.process.remote(nxt_job)
                    self._inflight[ref2] = nxt_job
                    # Record submission for totals used by GUI
                    try:
                        self.submitted_by_type[nxt_job.job_type] = (
                            self.submitted_by_type.get(nxt_job.job_type, 0) + 1
                        )
                    except Exception:
                        pass
                    self.total_enqueued += 1
                    # Update per-sample aggregate immediately for GUI samples view
                    try:
                        sid_local2 = (
                            nxt_job.context.get_sample_id()
                            if hasattr(nxt_job.context, "get_sample_id")
                            else "unknown"
                        )
                        if sid_local2 != "unknown":
                            ent_local2 = self.samples_by_id.get(sid_local2)
                            if ent_local2 is None:
                                ent_local2 = {
                                    "sample_id": sid_local2,
                                    "active_jobs": 0,
                                    "total_jobs": 0,
                                    "completed_jobs": 0,
                                    "failed_jobs": 0,
                                    "job_types": set(),
                                    "last_seen": time.time(),
                                }
                                self.samples_by_id[sid_local2] = ent_local2
                            ent_local2["total_jobs"] += 1
                            ent_local2["active_jobs"] += 1
                            try:
                                ent_local2["job_types"].add(nxt_job.job_type)
                            except Exception:
                                pass
                            ent_local2["last_seen"] = time.time()
                    except Exception:
                        pass
                    self.active[nxt_job.job_id] = {
                        "job_type": nxt_job.job_type,
                        "filepath": nxt_job.context.filepath,
                        "queue": job_queue_of(nxt_job.job_type),
                        "start_time": time.time(),
                    }

        # Promote next classification job for this type if any (single global pipeline)
        # Global backpressure release: for the finished job type, submit one waiting
        # job if we are below the max_inflight threshold.
        try:
            jt = job.job_type
            # decrement inflight count for this type
            if self.inflight_by_type.get(jt, 0) > 0:
                self.inflight_by_type[jt] -= 1
                if self.inflight_by_type[jt] == 0:
                    self.inflight_by_type.pop(jt, None)
            queue_global = self.waiting_by_type_global.get(jt, [])
            if (
                queue_global
                and int(self.inflight_by_type.get(jt, 0)) < self.max_inflight_per_type
            ):
                nxt_job_g = queue_global.pop(0)
                if not queue_global:
                    self.waiting_by_type_global.pop(jt, None)
                proc_g = self.processors.get(jt)
                if proc_g is not None:
                    if getattr(self, "using_pools", False):
                        try:
                            proc_g.enqueue.remote(nxt_job_g)
                        except Exception:
                            pass
                    else:
                        ref_g = proc_g.process.remote(nxt_job_g)
                        self._inflight[ref_g] = nxt_job_g
                    self.inflight_by_type[jt] = (
                        int(self.inflight_by_type.get(jt, 0)) + 1
                    )
                    self.submitted_by_type[jt] = self.submitted_by_type.get(jt, 0) + 1
                    self.total_enqueued += 1
                    # Update per-sample aggregate for GUI samples view
                    try:
                        sid_local_g = (
                            nxt_job_g.context.get_sample_id()
                            if hasattr(nxt_job_g.context, "get_sample_id")
                            else "unknown"
                        )
                        if sid_local_g != "unknown":
                            ent_local_g = self.samples_by_id.get(sid_local_g)
                            if ent_local_g is None:
                                ent_local_g = {
                                    "sample_id": sid_local_g,
                                    "active_jobs": 0,
                                    "total_jobs": 0,
                                    "completed_jobs": 0,
                                    "failed_jobs": 0,
                                    "job_types": set(),
                                    "last_seen": time.time(),
                                }
                                self.samples_by_id[sid_local_g] = ent_local_g
                            ent_local_g["total_jobs"] += 1
                            ent_local_g["active_jobs"] += 1
                            try:
                                ent_local_g["job_types"].add(nxt_job_g.job_type)
                            except Exception:
                                pass
                            ent_local_g["last_seen"] = time.time()
                    except Exception:
                        pass
                    self.active[nxt_job_g.job_id] = {
                        "job_type": nxt_job_g.job_type,
                        "filepath": nxt_job_g.context.filepath,
                        "queue": job_queue_of(nxt_job_g.job_type),
                        "start_time": time.time(),
                    }
        except Exception:
            pass
        if job.job_type in CLASSIFICATION_TYPES:
            # decrement pending for this type
            if self.classif_pending_by_type.get(job.job_type, 0) > 0:
                self.classif_pending_by_type[job.job_type] -= 1
                if self.classif_pending_by_type[job.job_type] == 0:
                    self.classif_pending_by_type.pop(job.job_type, None)
            nxt = self.classif_waiting_by_type.get(job.job_type)
            if nxt is not None:
                self.classif_waiting_by_type[job.job_type] = None
                proc3 = self.processors.get(job.job_type)
                if proc3 is not None:
                    if getattr(self, "using_pools", False):
                        try:
                            proc3.enqueue.remote(nxt)
                        except Exception:
                            pass
                    else:
                        ref3 = proc3.process.remote(nxt)
                        self._inflight[ref3] = nxt
                    # Record submission for totals used by GUI
                    try:
                        self.submitted_by_type[job.job_type] = (
                            self.submitted_by_type.get(job.job_type, 0) + 1
                        )
                    except Exception:
                        pass
                    self.classif_pending_by_type[job.job_type] = 1
                    self.total_enqueued += 1
                    # Update per-sample aggregate immediately for GUI samples view
                    try:
                        sid_local3 = (
                            nxt.context.get_sample_id()
                            if hasattr(nxt, "context")
                            and hasattr(nxt.context, "get_sample_id")
                            else "unknown"
                        )
                        if sid_local3 != "unknown":
                            ent_local3 = self.samples_by_id.get(sid_local3)
                            if ent_local3 is None:
                                ent_local3 = {
                                    "sample_id": sid_local3,
                                    "active_jobs": 0,
                                    "total_jobs": 0,
                                    "completed_jobs": 0,
                                    "failed_jobs": 0,
                                    "job_types": set(),
                                    "last_seen": time.time(),
                                }
                                self.samples_by_id[sid_local3] = ent_local3
                            ent_local3["total_jobs"] += 1
                            ent_local3["active_jobs"] += 1
                            try:
                                ent_local3["job_types"].add(nxt.job_type)
                            except Exception:
                                pass
                            ent_local3["last_seen"] = time.time()
                    except Exception:
                        pass
                    self.active[nxt.job_id] = {
                        "job_type": nxt.job_type,
                        "filepath": nxt.context.filepath,
                        "queue": job_queue_of(nxt.job_type),
                        "start_time": time.time(),
                    }

    async def _drain_loop(self):
        while True:
            try:
                if not self._inflight:
                    await asyncio.sleep(0.1)
                    continue
                refs = list(self._inflight.keys())
                ready, _ = ray.wait(refs, num_returns=1, timeout=0.1)
                for r in ready:
                    job = self._inflight.pop(r, None)
                    ok, ctx, err = True, None, None
                    try:
                        ctx = await r
                    except Exception as e:
                        ok, err = False, str(e)
                        ctx = job.context if job else WorkflowContext(filepath="")
                    if job:
                        await self._on_finish(job, ok, ctx, err)
            except Exception:
                await asyncio.sleep(0.1)

    async def mark_running(self, job: Job):
        if job.job_type in DEDUP_TYPES:
            key = (job.context.get_sample_id(), job.job_type)
            self.pending[key] = max(0, self.pending.get(key, 0) - 1)
            if self.pending[key] == 0:
                self.pending.pop(key, None)
            self.running[key] = self.running.get(key, 0) + 1

    async def stats(self) -> Dict[str, Any]:
        # per-queue active summary (also expose as 'active_by_worker' for GUI compat)
        active_by_queue: Dict[str, List[Dict[str, Any]]] = {}
        now = time.time()
        for jid, info in self.active.items():
            q = info["queue"]
            active_by_queue.setdefault(q, []).append(
                {
                    "job_id": jid,
                    "job_type": info["job_type"],
                    "filepath": info["filepath"],
                    "duration": int(now - info["start_time"]),
                }
            )

        # build category counts
        def _cat_of(jt: str) -> str:
            if jt == "preprocessing":
                return "preprocessing"
            if jt in {"mgmt", "cnv", "target", "fusion"}:
                return jt
            if jt in CLASSIFICATION_TYPES:
                return "classification"
            return "other"

        running_counts: Dict[str, int] = {
            "preprocessing": 0,
            "mgmt": 0,
            "cnv": 0,
            "target": 0,
            "fusion": 0,
            "classification": 0,
            "other": 0,
        }
        totals_counts: Dict[str, int] = {
            "preprocessing": 0,
            "mgmt": 0,
            "cnv": 0,
            "target": 0,
            "fusion": 0,
            "classification": 0,
            "other": 0,
        }
        for info in self.active.values():
            running_counts[_cat_of(info["job_type"])] += 1
        # totals are the submitted (ever) per type/category
        for jt, c in self.submitted_by_type.items():
            totals_counts[_cat_of(jt)] += int(c or 0)
        # Count waiting (serialized) jobs so the monitor can account for them
        try:
            waiting_serialized = sum(
                len(v) for v in getattr(self, "waiting_by_type_sample", {}).values()
            )
        except Exception:
            waiting_serialized = 0
        # Count global backpressure queues (per-type)
        try:
            waiting_global = sum(
                len(v) for v in getattr(self, "waiting_by_type_global", {}).values()
            )
        except Exception:
            waiting_global = 0
        # Samples payload for GUI
        samples_payload: List[Dict[str, Any]] = []
        try:
            for sid, ent in self.samples_by_id.items():
                samples_payload.append(
                    {
                        "sample_id": sid,
                        "active_jobs": ent.get("active_jobs", 0),
                        "total_jobs": ent.get("total_jobs", 0),
                        "completed_jobs": ent.get("completed_jobs", 0),
                        "failed_jobs": ent.get("failed_jobs", 0),
                        "job_types": list(ent.get("job_types", set())),
                        "last_seen": ent.get("last_seen", now),
                    }
                )
        except Exception:
            samples_payload = []

        return {
            "completed": len(self.completed),
            "failed": len(self.failed),
            "total_enqueued": self.total_enqueued,
            "total_skipped": self.total_skipped,
            "completed_by_type": dict(self.completed_by_type),
            "failed_by_type": dict(self.failed_by_type),
            "active_by_queue": active_by_queue,
            "active_by_worker": active_by_queue,  # compatibility for GUI table
            "active_count": len(self.active),
            "waiting_serialized": waiting_serialized,
            "waiting_global": waiting_global,
            "running_by_category": running_counts,
            "totals_by_category": totals_counts,
            "samples": samples_payload,
        }


@ray.remote
class Pool:
    def __init__(self, queue_name: str, max_parallel: int = 1):
        self.queue_name = queue_name
        self.max_parallel = max(1, int(max_parallel))
        # job_type -> (remote_func, resource_options)
        self.handlers: Dict[str, Tuple[Any, Dict[str, Any]]] = {}
        # callback to Coordinator actor
        self._coordinator = None
        self._coordinator_name: Optional[str] = None
        self._inflight: Set[ray.ObjectRef] = set()
        self._running_count: int = 0

    async def set_coordinator_name(self, coord_name: str):
        # store Coordinator actor name; resolve handle lazily
        self._coordinator_name = coord_name
        try:
            self._coordinator = ray.get_actor(coord_name)
        except Exception:
            self._coordinator = None
        return True

    def _get_coordinator(self):
        if self._coordinator is None and self._coordinator_name:
            try:
                self._coordinator = ray.get_actor(self._coordinator_name)
            except Exception:
                self._coordinator = None
        return self._coordinator

    async def register_handler(
        self, job_type: str, remote_func, resource_options: Dict[str, Any]
    ):
        self.handlers[job_type] = (remote_func, resource_options)

    async def enqueue(self, job: Job):
        # tell coordinator that we're starting this job (updates dedup maps)
        # We cannot call back into Coordinator directly from here without a reference;
        # pass coordinator handle via callback's bound actor method (Ray passes actor ref implicitly).
        # Instead, we expose a small trampoline: the driver calls mark_running before enqueue.
        # For simplicity, Coordinator calls mark_running() before Pool.enqueue().
        tup = self.handlers.get(job.job_type)
        if tup is None:
            ctx = job.context
            ctx.add_error(
                job.job_type,
                f"No handler for {job.job_type} in queue {self.queue_name}",
            )
            coord = self._get_coordinator()
            if coord is not None:
                await coord._on_finish.remote(
                    job, False, ctx, f"no_handler:{job.job_type}"
                )
            return
        remote_func, opts = tup
        # Apply resource hints at call site
        # Bound concurrency to max_parallel to emulate thread runner
        while self._running_count >= self.max_parallel:
            await asyncio.sleep(0.005)
        ref = remote_func.options(**opts).remote(job)
        self._inflight.add(ref)
        self._running_count += 1

        async def _wait(r):
            ok, ctx, err = True, None, None
            try:
                ctx = await r
            except Exception as e:
                ok, err = False, str(e)
                ctx = job.context
                ctx.add_error(job.job_type, err)
            coord2 = self._get_coordinator()
            if coord2 is not None:
                await coord2._on_finish.remote(job, ok, ctx, err)
            self._inflight.discard(r)
            if self._running_count > 0:
                self._running_count -= 1

        asyncio.create_task(_wait(ref))

    # Adapter so callers using TypeProcessor-style .process() keep working
    async def process(self, job: Job):
        return await self.enqueue(job)


# ---------- Classifier & Runner ----------


def default_file_classifier(filepath: str, plan: List[str]) -> List[Job]:
    ctx = WorkflowContext(filepath)
    ctx.add_metadata("filename", os.path.basename(filepath))
    ctx.add_metadata("created", time.time())
    try:
        ctx.add_metadata("file_size", os.path.getsize(filepath))
    except OSError:
        ctx.add_metadata("file_size", 0)
    ctx.add_metadata("original_workflow", plan)

    job_id = next(_job_id_counter)
    first = plan[0]
    if ":" not in first:
        raise ValueError(f"Invalid workflow step: {first}")
    q, jt = first.split(":", 1)
    return [Job(job_id, jt, q, plan, 0, ctx)]


def _matches_any_pattern(p: Path, patterns: Optional[List[str]]) -> bool:
    if not patterns:
        return True
    try:
        return any(p.match(pat) for pat in patterns)
    except Exception:
        return True


def _matches_no_ignores(p: Path, ignore_patterns: Optional[List[str]]) -> bool:
    if not ignore_patterns:
        return True
    try:
        return not any(p.match(pat) for pat in ignore_patterns)
    except Exception:
        return True


async def submit_existing_paths(
    coord,
    paths: List[str],
    plan: List[str],
    patterns: Optional[List[str]] = None,
    ignore_patterns: Optional[List[str]] = None,
    recursive: bool = True,
    work_dir: Optional[str] = None,
) -> None:
    seed_jobs: List[Job] = []
    for p in paths:
        pth = Path(p)
        if pth.is_file():
            if _matches_any_pattern(pth, patterns) and _matches_no_ignores(
                pth, ignore_patterns
            ):
                jobs = default_file_classifier(str(pth), plan)
                if work_dir:
                    for j in jobs:
                        j.context.add_metadata("work_dir", work_dir)
                seed_jobs += jobs
        elif pth.is_dir():
            # Stream directory walk in small batches; avoid building huge lists
            walker = pth.rglob("*") if recursive else pth.glob("*")
            batch: List[Job] = []
            for f in walker:
                if not f.is_file():
                    continue
                if not _matches_any_pattern(f, patterns) or not _matches_no_ignores(
                    f, ignore_patterns
                ):
                    continue
                jobs = default_file_classifier(str(f), plan)
                if work_dir:
                    for j in jobs:
                        j.context.add_metadata("work_dir", work_dir)
                batch += jobs
                if len(batch) >= 256:
                    await coord.submit_jobs.remote(batch)
                    batch = []
            if batch:
                await coord.submit_jobs.remote(batch)
    if seed_jobs:
        # Submit in bounded batches to avoid flooding the coordinator/actors
        BATCH = 256
        for i in range(0, len(seed_jobs), BATCH):
            await coord.submit_jobs.remote(seed_jobs[i : i + BATCH])


async def tqdm_monitor(coord, continuous: bool = False) -> None:
    # Create queue-specific progress bars (similar to your original). If
    # continuous=True, do not exit when queues drain; keep running for
    # watchdog-driven new work until interrupted.
    order = [
        "preprocessing",
        "bed_conversion",
        "mgmt",
        "cnv",
        "target",
        "fusion",
        "classification",
        "slow",
    ]
    bars = {
        q: tqdm(desc=q.title(), unit="jobs", position=i, leave=True)
        for i, q in enumerate(order)
    }
    overall = tqdm(
        desc="Overall Progress", unit="jobs", position=len(order), leave=True
    )

    try:
        last_completed = 0
        while True:
            s = await coord.stats.remote()
            # Build completed per queue (completed + failed)
            completed_by_type = s.get("completed_by_type", {})
            failed_by_type = s.get("failed_by_type", {})

            completed_queue_counts: Dict[str, int] = {q: 0 for q in order}
            for jt, c in completed_by_type.items():
                q = job_queue_of(jt)
                completed_queue_counts[q] += c
            for jt, c in failed_by_type.items():
                q = job_queue_of(jt)
                completed_queue_counts[q] += c

            # Active per queue
            active_by_queue = s.get("active_by_queue", {})

            # Update per-queue bars
            for q in order:
                completed_in_q = completed_queue_counts.get(q, 0)
                active_in_q = len(active_by_queue.get(q, []))
                total_for_q = completed_in_q + active_in_q
                b = bars[q]
                b.total = max(b.total or 0, total_for_q)
                b.n = completed_in_q
                # show up to 2 active jobs
                active_jobs_info = active_by_queue.get(q, [])
                if active_jobs_info:
                    parts = []
                    for j in active_jobs_info[:2]:
                        fn = os.path.basename(j.get("filepath", ""))
                        parts.append(f"{j['job_type']}:{fn}({j['duration']}s)")
                    b.set_postfix_str(" | ".join(parts))
                else:
                    b.set_postfix_str("")
                b.refresh()

            # Overall
            total_with_waiting = (
                s.get("total_enqueued", 0) - s.get("total_skipped", 0)
            ) + int(s.get("waiting_serialized", 0) or 0)
            overall.total = max(overall.total or 0, total_with_waiting) or 0
            overall.n = s["completed"] + s["failed"]
            overall.set_postfix_str(
                f"Active:{s['active_count']} | Completed:{s['completed']} | "
                f"Failed:{s['failed']} | Skipped:{s['total_skipped']} | Total:{s['total_enqueued']}"
            )
            overall.refresh()

            # exit condition: only in non-continuous mode (batch). In continuous
            # mode, keep monitoring to allow watchdog to enqueue new files.
            if (
                (not continuous)
                and s["active_count"] == 0
                and overall.n >= overall.total
            ):
                break

            await asyncio.sleep(0.5)
    finally:
        for b in bars.values():
            b.close()
        overall.close()


class RayFileWatcher(FileSystemEventHandler):
    def __init__(
        self,
        coord,
        plan: List[str],
        patterns: Optional[List[str]] = None,
        ignore_patterns: Optional[List[str]] = None,
        recursive: bool = True,
        work_dir: Optional[str] = None,
    ):
        self.coord = coord
        self.plan = plan
        self.patterns = patterns or ["*"]
        self.ignore_patterns = ignore_patterns or []
        self.recursive = recursive
        self.work_dir = work_dir
        self.processed: set[str] = set()
        # Rate limiting / batching
        self._pending_jobs: List[Job] = []
        self._last_flush: float = time.time()
        self._flush_interval_s: float = 0.5
        self._batch_size: int = 128

    def _should_process(self, fp: str) -> bool:
        p = Path(fp)
        if any(p.match(ip) for ip in self.ignore_patterns):
            return False
        if self.patterns:
            return any(p.match(pt) for pt in self.patterns)
        return True

    def _flush_if_needed(self, force: bool = False) -> None:
        try:
            if not self._pending_jobs:
                return
            now = time.time()
            if (
                force
                or len(self._pending_jobs) >= self._batch_size
                or (now - self._last_flush) >= self._flush_interval_s
            ):
                batch = self._pending_jobs[: self._batch_size]
                del self._pending_jobs[: self._batch_size]
                try:
                    ray.get(self.coord.submit_jobs.remote(batch))
                except Exception:
                    pass
                self._last_flush = now
        except Exception:
            pass

    def _handle(self, fp: str):
        if fp in self.processed:
            return
        if not self._should_process(fp):
            return
        self.processed.add(fp)
        jobs = default_file_classifier(fp, self.plan)
        if self.work_dir:
            for j in jobs:
                j.context.add_metadata("work_dir", self.work_dir)
        # enqueue and flush under rate limiter
        self._pending_jobs.extend(jobs)
        self._flush_if_needed()

    def on_created(self, event):
        if not event.is_directory:
            self._handle(event.src_path)

    def on_modified(self, event):
        if not event.is_directory:
            self._handle(event.src_path)

    # Ensure remaining jobs are flushed on shutdown
    def flush_remaining(self):
        while self._pending_jobs:
            self._flush_if_needed(force=True)


async def run(
    plan: List[str],
    paths: List[str],
    analysis_workers: int = 1,
    process_existing: bool = True,
    monitor: bool = True,
    watch: bool = False,
    patterns: Optional[List[str]] = None,
    ignore_patterns: Optional[List[str]] = None,
    recursive: bool = True,
    work_dir: Optional[str] = None,
    log_level: str = "INFO",
    preset: Optional[str] = None,
):
    global GLOBAL_LOG_LEVEL
    GLOBAL_LOG_LEVEL = (log_level or "INFO").upper()
    # Ensure any previous coordinator is terminated to avoid stale-code actors
    try:
        old = ray.get_actor("littlejohn_coordinator")
        try:
            ray.kill(old)
        except Exception:
            pass
    except Exception:
        pass
    # Start a fresh coordinator (not detached)
    try:
        coord = Coordinator.options(name="littlejohn_coordinator", num_cpus=0).remote(
            analysis_workers=analysis_workers, preset=preset
        )
    except Exception:
        coord = Coordinator.options(name="littlejohn_coordinator").remote(
            analysis_workers=analysis_workers, preset=preset
        )
    # run async setup without blocking the event loop
    try:
        await coord.setup.remote()
    except Exception:
        pass

    if process_existing and paths:
        await submit_existing_paths(
            coord,
            paths,
            plan,
            patterns=patterns,
            ignore_patterns=ignore_patterns,
            recursive=recursive,
            work_dir=work_dir,
        )

    observer = None
    watcher = None
    if watch and paths:
        observer = Observer()
        watcher = RayFileWatcher(
            coord,
            plan,
            patterns=patterns,
            ignore_patterns=ignore_patterns,
            recursive=recursive,
            work_dir=work_dir,
        )
        for p in paths:
            if Path(p).is_dir():
                observer.schedule(watcher, p, recursive=recursive)
        observer.start()

    tasks = []
    if monitor:
        tasks.append(asyncio.create_task(tqdm_monitor(coord, continuous=watch)))

    # Optional GUI integration
    try:
        if _GUIUpdateType is None:
            print("GUI not launched: GUI modules unavailable on this environment.")
        elif not work_dir:
            print("GUI not launched: --work-dir not provided.")
        else:
            launcher = _gui_launch(
                workflow_runner=None, workflow_steps=plan, monitored_directory=work_dir
            )
            try:
                url = (
                    launcher.get_gui_url()
                    if hasattr(launcher, "get_gui_url")
                    else "http://localhost:8081"
                )
                print(f"NiceGUI ready to go on {url}")
            except Exception:
                print("NiceGUI launched.")

            async def _publish_gui():
                while True:
                    try:
                        s = await coord.stats.remote()
                        # Basic workflow status & progress
                        total = int(
                            (s.get("total_enqueued", 0) - s.get("total_skipped", 0))
                            or 0
                        )
                        completed = int(s.get("completed", 0) or 0)
                        failed = int(s.get("failed", 0) or 0)
                        active = int(s.get("active_count", 0) or 0)
                        progress = (completed + failed) / total if total > 0 else 0.0
                        _gui_send_update(
                            _GUIUpdateType.WORKFLOW_STATUS,
                            {"is_running": True, "start_time": time.time()},
                            priority=1,
                        )
                        _gui_send_update(
                            _GUIUpdateType.PROGRESS_UPDATE,
                            {
                                "progress": progress,
                                "completed": completed,
                                "failed": failed,
                                "total": total,
                            },
                            priority=1,
                        )

                        # Queue status (running/total by category). If the coordinator
                        # doesn't provide category totals/running, derive them from the
                        # available fields so the GUI can still display correct values.
                        ru = s.get("running_by_category", {}) or {}
                        tu = s.get("totals_by_category", {}) or {}

                        # Derive running by category when missing using active_by_queue
                        if not ru or sum(int(v or 0) for v in ru.values()) == 0:
                            ru = {
                                "preprocessing": 0,
                                "mgmt": 0,
                                "cnv": 0,
                                "target": 0,
                                "fusion": 0,
                                "classification": 0,
                                "other": 0,
                            }
                            abq = s.get("active_by_queue", {}) or {}
                            for qname, jobs in abq.items():
                                n = len(jobs) if isinstance(jobs, list) else 0
                                if qname in {
                                    "preprocessing",
                                    "mgmt",
                                    "cnv",
                                    "target",
                                    "fusion",
                                    "classification",
                                }:
                                    ru[qname] += n
                                elif qname in {"bed_conversion", "slow", "other"}:
                                    ru["other"] += n

                        # Derive totals by category when missing using completed/failed + active
                        if not tu or sum(int(v or 0) for v in tu.values()) == 0:
                            tu = {
                                "preprocessing": 0,
                                "mgmt": 0,
                                "cnv": 0,
                                "target": 0,
                                "fusion": 0,
                                "classification": 0,
                                "other": 0,
                            }
                            cbt = s.get("completed_by_type", {}) or {}
                            fbt = s.get("failed_by_type", {}) or {}

                            def _cat_of_local(jt: str) -> str:
                                if jt == "preprocessing":
                                    return "preprocessing"
                                if jt in {"mgmt", "cnv", "target", "fusion"}:
                                    return jt
                                if jt in CLASSIFICATION_TYPES:
                                    return "classification"
                                return "other"

                            for jt, c in cbt.items():
                                tu[_cat_of_local(jt)] += int(c or 0)
                            for jt, c in fbt.items():
                                tu[_cat_of_local(jt)] += int(c or 0)
                            # Add currently active jobs
                            abq2 = s.get("active_by_queue", {}) or {}
                            for qname, jobs in abq2.items():
                                n = len(jobs) if isinstance(jobs, list) else 0
                                if qname in {
                                    "preprocessing",
                                    "mgmt",
                                    "cnv",
                                    "target",
                                    "fusion",
                                    "classification",
                                }:
                                    tu[qname] += n
                                elif qname in {"bed_conversion", "slow", "other"}:
                                    tu["other"] += n

                        queue_payload = {
                            "preprocessing": {
                                "running": int(ru.get("preprocessing", 0) or 0),
                                "total": int(tu.get("preprocessing", 0) or 0),
                            },
                            "mgmt": {
                                "running": int(ru.get("mgmt", 0) or 0),
                                "total": int(tu.get("mgmt", 0) or 0),
                            },
                            "cnv": {
                                "running": int(ru.get("cnv", 0) or 0),
                                "total": int(tu.get("cnv", 0) or 0),
                            },
                            "target": {
                                "running": int(ru.get("target", 0) or 0),
                                "total": int(tu.get("target", 0) or 0),
                            },
                            "fusion": {
                                "running": int(ru.get("fusion", 0) or 0),
                                "total": int(tu.get("fusion", 0) or 0),
                            },
                            "classification": {
                                "running": int(ru.get("classification", 0) or 0),
                                "total": int(tu.get("classification", 0) or 0),
                            },
                            "other": {
                                "running": int(ru.get("other", 0) or 0),
                                "total": int(tu.get("other", 0) or 0),
                            },
                        }
                        _gui_send_update(
                            _GUIUpdateType.QUEUE_UPDATE, queue_payload, priority=2
                        )
                        # Also emit a human-readable log line so users see queue summaries in Live Logs
                        try:
                            pre = queue_payload["preprocessing"]
                            an_run = sum(
                                int(queue_payload[q]["running"])
                                for q in ["mgmt", "cnv", "target", "fusion"]
                            )
                            an_tot = sum(
                                int(queue_payload[q]["total"])
                                for q in ["mgmt", "cnv", "target", "fusion"]
                            )
                            cl = queue_payload["classification"]
                            ot = queue_payload["other"]
                            summary = (
                                f"Queues | Pre:{pre['running']}/{pre['total']} "
                                f"| An:{an_run}/{an_tot} "
                                f"| Cl:{cl['running']}/{cl['total']} "
                                f"| Ot:{ot['running']}/{ot['total']}"
                            )
                            _gui_send_update(
                                _GUIUpdateType.LOG_MESSAGE,
                                {"message": summary, "level": "INFO"},
                                priority=0,
                            )
                        except Exception:
                            pass

                        # Active jobs table
                        rows = []
                        for qname, jobs in (s.get("active_by_queue", {}) or {}).items():
                            for j in jobs:
                                rows.append(
                                    {
                                        "job_id": str(j.get("job_id", "")),
                                        "job_type": j.get("job_type", ""),
                                        "filepath": j.get("filepath", ""),
                                        "worker_name": qname,
                                        "duration": int(j.get("duration", 0) or 0),
                                        "progress": 0.0,
                                    }
                                )
                        _gui_send_update(
                            _GUIUpdateType.JOB_UPDATE, {"active_jobs": rows}, priority=1
                        )

                        # Samples overview
                        samples = s.get("samples", []) or []
                        _gui_send_update(
                            _GUIUpdateType.SAMPLES_UPDATE,
                            {"samples": samples},
                            priority=1,
                        )

                        await asyncio.sleep(1.0)
                    except asyncio.CancelledError:
                        break
                    except Exception:
                        await asyncio.sleep(1.0)

            tasks.append(asyncio.create_task(_publish_gui()))
    except Exception as e:
        print(f"Warning: GUI failed to launch: {e}")

    try:
        if tasks:
            await asyncio.gather(*tasks)
    finally:
        if watcher is not None:
            try:
                watcher.flush_remaining()
            except Exception:
                pass
        if observer is not None:
            observer.stop()
            observer.join()


def parse_args():
    p = argparse.ArgumentParser(description="Ray Core LittleJohn workflow")
    p.add_argument(
        "--plan",
        action="append",
        default=[],
        help="Workflow plan entries like 'queue:job_type' in order",
    )
    p.add_argument(
        "--paths",
        nargs="*",
        default=[],
        help="Files/dirs to process immediately (optional)",
    )
    p.add_argument(
        "--analysis-workers",
        type=int,
        default=1,
        help="Parallel workers per analysis queue",
    )
    p.add_argument("--no-monitor", action="store_true", help="Disable tqdm monitor")
    p.add_argument(
        "--no-process-existing",
        action="store_true",
        help="Do not process existing files/dirs",
    )
    p.add_argument(
        "--watch",
        action="store_true",
        help="Watch directories for new files (Watchdog)",
    )
    p.add_argument(
        "--patterns",
        action="append",
        default=None,
        help="Glob patterns to include (repeatable)",
    )
    p.add_argument(
        "--ignore",
        action="append",
        default=None,
        help="Glob patterns to ignore (repeatable)",
    )
    p.add_argument(
        "--no-recursive",
        action="store_true",
        help="Disable recursive directory watching",
    )
    p.add_argument(
        "--ray-address", default=None, help="Ray cluster address (None = local)"
    )
    p.add_argument(
        "--preset",
        default=None,
        choices=["p2i", "standard", "high"],
        help="Execution preset controlling actor grouping and concurrency",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    try:
        # Expose dashboard on all interfaces when supported
        ray.init(
            address=args.ray_address,
            include_dashboard=True,
            dashboard_host=os.environ.get("RAY_DASHBOARD_HOST", "0.0.0.0"),
        )
    except TypeError:
        # Older Ray versions may not support dashboard args
        ray.init(address=args.ray_address)

    if not args.plan:
        # default example
        args.plan = [
            "preprocessing:preprocessing",
            "bed_conversion:bed_conversion",
            "mgmt:mgmt",
            "sturgeon:sturgeon",
        ]

    try:
        asyncio.run(
            run(
                plan=args.plan,
                paths=args.paths,
                analysis_workers=args.analysis_workers,
                process_existing=not args.no_process_existing,
                monitor=not args.no_monitor,
                watch=args.watch,
                patterns=args.patterns,
                ignore_patterns=args.ignore,
                recursive=not args.no_recursive,
                preset=args.preset,
            )
        )
    except KeyboardInterrupt:
        pass
