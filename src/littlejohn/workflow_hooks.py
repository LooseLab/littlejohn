"""
Workflow hooks for integrating with the GUI message queue system.

This module provides hooks that send real-time updates to the GUI
without blocking workflow execution.
"""

import logging
import time
from typing import Any, List, Dict

from .gui_launcher import send_gui_update, UpdateType


def install_workflow_hooks(workflow_runner: Any, workflow_steps: List[str], monitored_directory: str) -> None:
    """Install hooks into the workflow system for GUI updates.

    The hooks are intentionally non-invasive: we wrap public methods and start a
    lightweight polling thread that reads `manager.get_stats()` to drive the UI.
    """
    try:
        logging.info("Installing workflow hooks for GUI integration")

        # Prime the UI with an initial status
        try:
            send_gui_update(
                UpdateType.WORKFLOW_STATUS,
                {
                    "is_running": True,
                    "start_time": time.time(),
                    "workflow_steps": workflow_steps,
                    "monitored_directory": monitored_directory,
                },
                priority=10,
            )
        except Exception:
            pass

        _install_runner_hooks(workflow_runner, workflow_steps, monitored_directory)
        _install_manager_hooks(workflow_runner)

        logging.info("✅ Workflow hooks installed for real-time GUI monitoring")
    except Exception as exc:
        logging.error(f"Failed to install workflow hooks: {exc}")
        logging.info("GUI will run without real-time updates")


def _infer_sample_id_from_path(filepath: str) -> str:
    """Best-effort sample id inference from a BAM path.

    Prefer the nearest ancestor directory whose name is not a BAM/FASTQ bucket
    (e.g., bam_fail, bam_pass, bam) and is not simply 'pass'/'fail'.
    Falls back to the immediate parent directory name.
    """
    try:
        from pathlib import Path
        p = Path(filepath)
        disallowed = {"bam_fail", "bam_pass", "bam", "fastq", "fastq_fail", "fastq_pass", "pass", "fail"}
        # Look up to 5 ancestors for a reasonable sample folder
        current = p.parent
        for _ in range(5):
            name = current.name
            if not name:
                break
            lname = name.lower()
            if lname not in disallowed and "bam" not in lname and "fastq" not in lname:
                return name
            current = current.parent
        return p.parent.name or "unknown"
    except Exception:
        return "unknown"


def _install_runner_hooks(workflow_runner: Any, workflow_steps: List[str], monitored_directory: str) -> None:
    """Wrap `run_workflow` to emit start/stop UI events."""
    try:
        if not hasattr(workflow_runner, "run_workflow"):
            return

        original_run = workflow_runner.run_workflow

        def hooked_run_workflow(*args, **kwargs):
            logging.info("[GUI] run_workflow wrapper: start event will be sent")
            send_gui_update(
                UpdateType.WORKFLOW_STATUS,
                {
                    "is_running": True,
                    "start_time": time.time(),
                    "workflow_steps": workflow_steps,
                    "monitored_directory": monitored_directory,
                },
                priority=10,
            )
            try:
                return original_run(*args, **kwargs)
            finally:
                logging.info("[GUI] run_workflow wrapper: stop event will be sent")
                send_gui_update(
                    UpdateType.WORKFLOW_STATUS,
                    {"is_running": False, "stop_time": time.time()},
                    priority=10,
                )

        workflow_runner.run_workflow = hooked_run_workflow
    except Exception as exc:
        logging.debug(f"Failed to install runner hooks: {exc}")


def _install_manager_hooks(workflow_runner: Any) -> None:
    """Start a polling thread based on the runner's manager."""
    try:
        manager = getattr(workflow_runner, "manager", None) or getattr(workflow_runner, "workflow_manager", None)
        if manager is None:
            logging.debug("No workflow manager found on runner; skipping manager hooks")
            return
        logging.info("[GUI] Starting polling thread for GUI updates")
        _start_polling_updates(manager)
    except Exception as exc:
        logging.debug(f"Failed to install manager hooks: {exc}")


def _start_polling_updates(manager: Any, interval_seconds: float = 1.5) -> None:
    """Start a background thread that polls `manager.get_stats()` and sends GUI updates."""
    import threading

    def poll() -> None:
        last_log_sent = 0.0
        while True:
            try:
                stats = manager.get_stats()
                now = time.time()
                logging.info(
                    f"[GUI] Poll: active={stats.get('active_jobs', 0)}, completed={stats.get('completed', 0)}, "
                    f"failed={stats.get('failed', 0)}, total={stats.get('total_actual_jobs', 0)}"
                )

                # Queue status
                try:
                    queue_sizes = stats.get("queue_sizes", {})
                    queue_update: Dict[str, Dict[str, int]] = {
                        "preprocessing": {"running": 0, "total": queue_sizes.get("preprocessing", 0)},
                        "mgmt": {"running": 0, "total": queue_sizes.get("mgmt", 0)},
                        "cnv": {"running": 0, "total": queue_sizes.get("cnv", 0)},
                        "target": {"running": 0, "total": queue_sizes.get("target", 0)},
                        "fusion": {"running": 0, "total": queue_sizes.get("fusion", 0)},
                        "classification": {"running": 0, "total": queue_sizes.get("classification", 0)},
                        "other": {"running": 0, "total": queue_sizes.get("slow", 0)},
                    }

                    active_by_worker: Dict[str, List[Dict[str, Any]]] = stats.get("active_by_worker", {})  # type: ignore
                    for worker_name, jobs in active_by_worker.items():
                        if worker_name.startswith("PreprocessingWorker"):
                            queue_update["preprocessing"]["running"] += len(jobs)
                        elif worker_name.startswith("MGMTWorker"):
                            queue_update["mgmt"]["running"] += len(jobs)
                        elif worker_name.startswith("CNVWorker"):
                            queue_update["cnv"]["running"] += len(jobs)
                        elif worker_name.startswith("TargetWorker"):
                            queue_update["target"]["running"] += len(jobs)
                        elif worker_name.startswith("FusionWorker"):
                            queue_update["fusion"]["running"] += len(jobs)
                        elif worker_name.startswith("ClassificationWorker"):
                            queue_update["classification"]["running"] += len(jobs)
                        elif worker_name.startswith("SlowWorker"):
                            queue_update["other"]["running"] += len(jobs)

                    send_gui_update(UpdateType.QUEUE_UPDATE, queue_update, priority=4)
                except Exception:
                    pass

                # Active jobs table
                try:
                    job_rows: List[Dict[str, Any]] = []
                    for worker_name, jobs in stats.get("active_by_worker", {}).items():  # type: ignore
                        for job in jobs:
                            job_rows.append(
                                {
                                    "job_id": "",
                                    "job_type": job.get("job_type", ""),
                                    "filepath": job.get("filepath", ""),
                                    "worker_name": worker_name,
                                    "duration": int(job.get("duration", 0)),
                                    "progress": 0.0,
                                }
                            )
                    send_gui_update(UpdateType.JOB_UPDATE, {"active_jobs": job_rows}, priority=5)
                except Exception:
                    pass

                # Overall progress (include counts for the UI summary)
                try:
                    total_processed = stats.get('total_processed', 0)
                    total_actual_jobs = stats.get('total_actual_jobs', 0)
                    progress = (
                        min(0.999, float(total_processed) / float(total_actual_jobs))
                        if total_actual_jobs > 0
                        else 0.0
                    )
                    send_gui_update(
                        UpdateType.PROGRESS_UPDATE,
                        {
                            "progress": progress,
                            "completed": stats.get('completed', 0),
                            "failed": stats.get('failed', 0),
                            "total": total_actual_jobs,
                        },
                        priority=2,
                    )
                except Exception:
                    pass

                # Periodic log heartbeat (every ~20s) so users see activity in GUI
                if now - last_log_sent > 20:
                    try:
                        msg = (
                            f"Heartbeat: active={stats.get('active_jobs', 0)}, "
                            f"completed={stats.get('completed', 0)}, failed={stats.get('failed', 0)}"
                        )
                        send_gui_update(UpdateType.LOG_MESSAGE, {"level": "INFO", "message": msg}, priority=1)
                        last_log_sent = now
                    except Exception:
                        pass

                # Samples overview using manager-provided per-sample stats when available
                try:
                    if 'samples' in stats and isinstance(stats['samples'], list):
                        send_gui_update(UpdateType.SAMPLES_UPDATE, {"samples": stats['samples']}, priority=1)
                    else:
                        # Fallback to aggregation from active_by_worker (legacy)
                        samples_map: Dict[str, Dict[str, Any]] = {}
                        for _worker, jobs in stats.get('active_by_worker', {}).items():  # type: ignore
                            for job in jobs:
                                fp = job.get('filepath', '')
                                sid = job.get('sample_id') or _infer_sample_id_from_path(fp)
                                if sid.lower() in {"bam_fail", "bam_pass", "bam", "pass", "fail"}:
                                    sid = _infer_sample_id_from_path(fp)
                                entry = samples_map.setdefault(
                                    sid,
                                    {
                                        "sample_id": sid,
                                        "active_jobs": 0,
                                        "total_jobs": 0,
                                        "completed_jobs": 0,
                                        "failed_jobs": 0,
                                        "job_types": set(),
                                        "last_seen": now,
                                    },
                                )
                                entry["active_jobs"] += 1
                                jt = job.get("job_type", "")
                                if jt:
                                    entry["job_types"].add(jt)
                                entry["last_seen"] = now
                        global_completed = stats.get("completed", 0)
                        global_failed = stats.get("failed", 0)
                        global_total = stats.get("total_actual_jobs", 0)
                        samples_payload = []
                        for entry in samples_map.values():
                            samples_payload.append(
                                {
                                    "sample_id": entry["sample_id"],
                                    "active_jobs": entry["active_jobs"],
                                    "total_jobs": global_total,
                                    "completed_jobs": global_completed,
                                    "failed_jobs": global_failed,
                                    "job_types": sorted(entry["job_types"]),
                                    "last_seen": entry["last_seen"],
                                }
                            )
                        send_gui_update(UpdateType.SAMPLES_UPDATE, {"samples": samples_payload}, priority=1)
                except Exception:
                    pass

            except Exception as exc:
                logging.debug(f"Polling error: {exc}")
            finally:
                time.sleep(interval_seconds)

    thread = threading.Thread(target=poll, daemon=True, name="LittleJohn-GUI-Polling")
    thread.start()
