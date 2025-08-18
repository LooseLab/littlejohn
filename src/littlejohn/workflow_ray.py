"""
Ray-powered workflow engine for LittleJohn.

This module provides a distributed, priority-aware job orchestration layer built on
top of Ray. It mirrors the semantics of the threaded workflow manager used elsewhere
in LittleJohn while adding the ability to run multiple jobs concurrently across
dedicated queues.

High-level architecture:

- Jobs and context
  - A `Job` encapsulates a unit of work (e.g., "preprocessing", "mgmt", "cnv").
  - Each job carries a shared `WorkflowContext` that accumulates per-file metadata,
    results, execution history, and errors as the workflow progresses.

- Queues and priorities
  - Jobs are mapped to queue types (e.g., `preprocessing`, `bed_conversion`,
    `mgmt`, `cnv`, `target`, `fusion`, `classification`, `slow`).
  - The `RayWorkflowManager` maintains per-priority in-memory queues and fairly
    drains them from highest to lowest priority (submitting at most one job per
    priority "tick" to prevent starvation and provide smooth progress across
    queues).

- Execution model (Ray)
  - Each queue type has one or more Ray actors (`JobProcessor`) responsible for
    executing jobs of that queue. Actors are configured with CPU and concurrency
    limits to spread work while avoiding oversubscription.
  - Handlers are simple callables that implement the actual domain logic for a
    given `job_type`. The manager is transport/coordination only.

- Triggering downstream work
  - Only the first job from a user-provided `workflow_plan` is enqueued initially
    (typically `preprocessing`).
  - As jobs complete, `Job.can_trigger_jobs()` checks declared dependencies and the
    original plan to spawn downstream jobs (e.g., `bed_conversion` after
    `preprocessing`, then `sturgeon`/`nanodx`/`pannanodx`/`random_forest` after
    `bed_conversion`). If no parallel triggers are applicable, a linear next step
    fallback (`Job.next_job()`) keeps the workflow moving.

- Operational safeguards and UX
  - A watchdog can mark jobs as timed-out per queue type and restart actors for
    recovery.
  - Per-sample de-duplication limits how many concurrent classification-like jobs
    (e.g., `sturgeon`, `nanodx`) can run for the same sample.
  - Rich, pull-based statistics power CLI/GUI progress bars without coupling the
    orchestration layer to any specific UI.

Usage pattern:

1) Instantiate `RayWorkflowRunner` (or `RayWorkflowManager` directly).
2) Register job handlers per queue and job type using `register_handler(...)`.
3) Provide a `workflow_plan` and input directory to `run_workflow(...)`.
4) The file watcher discovers inputs, creates the first job per file, and the
   manager handles the rest.

This module intentionally avoids embedding domain-specific logic inside handlers;
it focuses on scheduling, coordination, and observability. See inline docstrings
throughout for details on data flow and the scheduling rules.
"""

import os
import logging
import time
import threading
import queue
import itertools
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Any, Set, Union, Tuple
import pickle

import ray
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
from tqdm import tqdm

from littlejohn.logging_config import get_job_logger


# === Shared Context for Each File ===
@dataclass
class WorkflowContext:
    """
    Shared, mutable state carried across all jobs for a single input file.

    The context serves as the "glue" of the workflow:
    - `metadata` stores arbitrary properties (e.g., `bam_metadata`,
      `original_workflow`).
    - `results` holds outputs keyed by `job_type` for later cross-step usage.
    - `history` captures the sequence of job types that have recorded a result,
      which is also used to evaluate downstream dependency satisfaction.
    - `errors` collects structured error entries for post-mortem inspection.

    Attributes:
    - filepath: Absolute path to the input file (e.g., a BAM) this context refers to.
    - metadata: Free-form key/value store enriched by handlers and orchestration.
    - results: Outputs keyed by job type.
    - history: Ordered list of job types that have produced results.
    - errors: List of error dicts with `job_type`, `error`, and `timestamp`.
    """
    filepath: str
    metadata: dict = field(default_factory=dict)
    results: dict = field(default_factory=dict)
    history: list = field(default_factory=list)
    errors: list = field(default_factory=list)

    def add_metadata(self, key: str, value: Any) -> None:
        """Add metadata to the context."""
        self.metadata[key] = value

    def add_result(self, job_type: str, result: Any) -> None:
        """Add a result from a job to the context."""
        self.results[job_type] = result
        self.history.append(job_type)

    def add_error(self, job_type: str, error: str) -> None:
        """Add an error to the context."""
        self.errors.append({"job_type": job_type, "error": error, "timestamp": time.time()})

    def get_summary(self) -> dict:
        """Get a summary of the workflow execution."""
        return {
            "filepath": self.filepath,
            "metadata": self.metadata,
            "results": self.results,
            "history": self.history,
            "errors": self.errors,
            "success": len(self.errors) == 0
        }

    def get_sample_id(self) -> str:
        """Get the sample ID from metadata."""
        bam_metadata = self.metadata.get('bam_metadata', {})
        return bam_metadata.get('sample_id', 'unknown')


# === Job Object ===
@dataclass
class Job:
    """
    A unit of work that runs on a specific queue and advances the workflow.

    A job contains just enough information for the orchestration layer to route it
    to the right executor and for handlers to perform their work.

    Attributes:
    - job_id: Monotonic identifier for tracking and logging.
    - job_type: Logical type (e.g., "preprocessing", "mgmt").
    - context: The shared `WorkflowContext` for this file/sample.
    - origin: Queue type where this job was first enqueued (informational).
    - workflow: The original user plan as `"queue:job_type"` strings.
    - step: Index into the `workflow` (used by `next_job()` linear fallback).
    - dependencies: Which job types must have completed before this job can run.
    - triggers: Which job types this job can unlock upon completion.
    """
    job_id: int
    job_type: str
    context: WorkflowContext
    origin: str  # queue type
    workflow: List[str]
    step: int = 0
    dependencies: Set[str] = field(default_factory=set)  # Job types this job depends on
    triggers: Set[str] = field(default_factory=set)  # Job types this job can trigger

    def next_job(self) -> Optional['Job']:
        """
        Create the linear next job defined in the `workflow` plan, if any.

        This is a conservative fallback used only when there are no applicable
        parallel triggers. It preserves a simple linear progression and guards
        against malformed workflow steps.
        """
        if self.step + 1 < len(self.workflow):
            next_step = self.workflow[self.step + 1]
            if not isinstance(next_step, str):
                raise ValueError(f"Expected string in workflow step, got {type(next_step)}: {next_step}")
            
            if ':' not in next_step:
                raise ValueError(f"Invalid workflow step format (missing ':'): {next_step}")
                
            queue_type, next_type = next_step.split(":", 1)
            return Job(
                job_id=next(_job_id_counter),
                job_type=next_type,
                context=self.context,
                origin=queue_type,
                workflow=self.workflow,
                step=self.step + 1
            )
        return None

    def get_sample_id(self) -> str:
        """Get the sample ID for this job."""
        return self.context.get_sample_id()

    def can_trigger_jobs(self) -> List['Job']:
        """
        Compute parallel downstream jobs unlocked by this job's completion.

        - Dependency graph is defined locally here to keep orchestration simple.
        - Only jobs listed in the original workflow plan are eligible to trigger.
        - A dependency is considered satisfied if the dep is present in
          `context.history` (or is the current job).
        """
        triggered_jobs = []
        
        # Define job dependencies and triggers
        job_dependencies = {
            # Jobs that depend on preprocessing
            "bed_conversion": {"preprocessing"},
            "mgmt": {"preprocessing"},
            "cnv": {"preprocessing"},
            "target": {"preprocessing"},
            "fusion": {"preprocessing"},
            
            # Jobs that depend on bed_conversion
            "sturgeon": {"bed_conversion"},
            "nanodx": {"bed_conversion"},
            "pannanodx": {"bed_conversion"},
            "random_forest": {"bed_conversion"},
        }
        
        # Get the original workflow plan to check what jobs are actually requested
        original_workflow = self.context.metadata.get("original_workflow", [])
        if original_workflow:
            workflow_job_types = []
            for step in original_workflow:
                if isinstance(step, str) and ':' in step:
                    queue_type, job_type = step.split(':', 1)
                    workflow_job_types.append(job_type)
            original_workflow = workflow_job_types
        else:
            original_workflow = []
            for step in self.workflow:
                if isinstance(step, str) and ':' in step:
                    queue_type, job_type = step.split(':', 1)
                    original_workflow.append(job_type)
        
        # Check if this job can trigger other jobs
        for job_type, dependencies in job_dependencies.items():
            if self.job_type in dependencies:
                # Check if all dependencies are met
                all_deps_met = True
                for dep in dependencies:
                    if dep not in self.context.history and dep != self.job_type:
                        all_deps_met = False
                        break
                
                # Only trigger if the job is in the original workflow plan
                if all_deps_met and job_type in original_workflow:
                    queue_type = self._get_queue_type_for_job(job_type)
                    triggered_job = Job(
                        job_id=next(_job_id_counter),
                        job_type=job_type,
                        context=self.context,
                        origin=queue_type,
                        workflow=[f"{queue_type}:{job_type}"],
                        step=0,
                        dependencies=dependencies,
                        triggers=set()
                    )
                    triggered_jobs.append(triggered_job)
        
        return triggered_jobs

    def _get_queue_type_for_job(self, job_type: str) -> str:
        """Get the queue type for a given job type."""
        queue_mapping = {
            "preprocessing": "preprocessing",
            "bed_conversion": "bed_conversion",
            "mgmt": "mgmt",
            "cnv": "cnv",
            "target": "target",
            "fusion": "fusion",
            "sturgeon": "classification",
            "nanodx": "classification",
            "pannanodx": "classification",
            "random_forest": "slow",
        }
        return queue_mapping.get(job_type, "slow")

    def _get_display_queue_for_job(self, job_type: str) -> str:
        """Return the display queue key for a job type.

        In legacy analysis mode, surface per-type queues (mgmt/cnv/target/fusion)
        even though they run on a unified actor.
        """
        q_map = self._get_queue_type_for_job(job_type)
        if not self.use_separate_analysis_queues and job_type in {'mgmt', 'cnv', 'target', 'fusion'}:
            return job_type
        return q_map


# Global job ID counter
_job_id_counter = itertools.count(1)


# === Ray Actors for Job Processing ===
@ray.remote
class JobProcessor:
    """
    Ray actor responsible for executing jobs for a specific queue.

    The manager submits pickled `Job` objects to these actors. Each actor
    owns counters for processed/failed jobs and is configured for CPU and logical
    concurrency. Handlers are injected by the manager and looked up by `job_type`.
    """
    
    def __init__(self, queue_type: str, handlers: Dict[str, Callable], log_level: str, verbose: bool):
        self.queue_type = queue_type
        self.handlers = handlers
        self.processed_jobs = 0
        self.failed_jobs = 0
        self.log_level = log_level
        self.verbose = verbose
        
        # Configure logging for this actor
        self._configure_logging()
        
    def _configure_logging(self):
        """
        Configure per-actor logging.

        Logging in Ray actors can get noisy and interfere with progress bars;
        we gate user-facing prints behind `verbose` while still setting the
        desired global level using the project's logging config.
        """
        try:
            from littlejohn.logging_config import configure_logging
            configure_logging(global_level=self.log_level)
            # Only print initialization message if verbose mode is enabled
            # This prevents interference with progress bars
            if self.verbose:
                print(f"JobProcessor for {self.queue_type} initialized with log level: {self.log_level}")
        except Exception as e:
            # Only print warnings if verbose mode is enabled
            if self.verbose:
                print(f"Warning: Could not configure logging for {self.queue_type} actor: {e}")
        
    def process_job(self, job_data: bytes) -> Dict[str, Any]:
        """
        Execute one job and return a structured result payload.

        Returns a dict with keys such as:
        - status: 'success' or 'error'
        - job_id, job_type
        - context: the (possibly updated) `WorkflowContext`
        - triggered_jobs: list of parallel `Job`s to enqueue next
        - next_job: linear fallback `Job` if no triggers are applicable
        - processing_start_time: when the handler began (used for accurate durations)
        """
        try:
            # Deserialize the job
            job = pickle.loads(job_data)
        except Exception as e:
            logger.error(f"Failed to deserialize job data: {e}")
            return {
                'status': 'error',
                'error': f"Deserialization error: {e}"
            }

        # Get job-specific logger
        logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)

        # Set log level for this job if needed
        if hasattr(logger, 'setLevel'):
            import logging
            level_map = {
                "DEBUG": logging.DEBUG,
                "INFO": logging.INFO,
                "WARNING": logging.WARNING,
                "ERROR": logging.ERROR
            }
            logger.setLevel(level_map.get(self.log_level, logging.INFO))

        logger.info(f"Starting job {job.job_id} ({job.job_type}) for {job.context.filepath}")
        logger.info(f"Job {job.job_id} ({job.job_type}) state transition: QUEUED -> PROCESSING")

        # Record the actual processing start time (excludes queue waiting time)
        processing_start_time = time.time()

        # Process the job
        if job.job_type in self.handlers:
            try:
                self.handlers[job.job_type](job)
                logger.info(f"Job {job.job_id} ({job.job_type}) completed successfully")
                logger.info(f"Job {job.job_id} ({job.job_type}) state transition: PROCESSING -> COMPLETED")
                self.processed_jobs += 1
                # Post-handler check: ensure preprocessing determined a sample_id
                try:
                    if job.job_type == 'preprocessing':
                        sid = job.context.get_sample_id() if hasattr(job.context, 'get_sample_id') else None
                        if not sid or sid == 'unknown':
                            msg = (
                                f"Unknown sample_id after preprocessing. File: {job.context.filepath}. "
                                f"This will impact downstream per-sample dedup/tracking."
                            )
                            job.context.add_error('preprocessing', msg)
                            logger.error(msg)
                except Exception as e:
                    logger.error(f"Error during post-handler check: {e}")

                # Compute triggers explicitly so we can log them for diagnostics
                triggered_jobs = job.can_trigger_jobs()
                downstream = [tj.job_type for tj in triggered_jobs]
                logger.info(f"Job {job.job_id} triggered downstream jobs: {downstream}")

                # Return updated context and triggered jobs
                return {
                    'status': 'success',
                    'job_id': job.job_id,
                    'job_type': job.job_type,
                    'context': job.context,
                    # Return only job types to avoid cross-process job_id collisions
                    'triggered_job_types': downstream,
                    'processing_start_time': processing_start_time
                }
            except Exception as handler_error:
                error_msg = f"Handler execution failed: {str(handler_error)}"
                job.context.add_error(job.job_type, error_msg)
                logger.error(f"Job {job.job_id} failed: {error_msg}")
                logger.info(f"Job {job.job_id} ({job.job_type}) state transition: PROCESSING -> FAILED")
                self.failed_jobs += 1
                return {
                    'status': 'error',
                    'job_id': job.job_id,
                    'job_type': job.job_type,
                    'context': job.context,
                    'error': error_msg,
                    'processing_start_time': processing_start_time
                }
        else:
            available_handlers = list(self.handlers.keys())
            error_msg = f"No handler found for job type: {job.job_type}. Available handlers: {available_handlers}"
            job.context.add_error(job.job_type, error_msg)
            logger.error(f"Job {job.job_id} failed: {error_msg}")
            logger.info(f"Job {job.job_id} ({job.job_type}) state transition: QUEUED -> NO_HANDLER")
            self.failed_jobs += 1
            return {
                'status': 'error',
                'job_id': job.job_id,
                'job_type': job.job_type,
                'context': job.context,
                'error': error_msg,
                'processing_start_time': processing_start_time
            }
    
    def get_stats(self) -> Dict[str, int]:
        """Get processing statistics."""
        return {
            'processed_jobs': self.processed_jobs,
            'failed_jobs': self.failed_jobs
        }


# === Ray-based Workflow Manager ===
class RayWorkflowManager:
    """
    Priority-based scheduler that coordinates Ray actors and workflow state.

    Responsibilities:
    - Maintain per-priority queues and submit work fairly from highest to lowest.
    - Track active/running jobs, durations, and completions/failures by type.
    - Enforce per-sample de-duplication for classification-like workloads.
    - Surface aggregated stats for CLI/GUI while remaining UI-agnostic.
    - Optionally restart actors when a watchdog timeout is exceeded.

    Configuration knobs include the number of workers per queue and whether to
    use separate analysis queues (`mgmt`, `cnv`, `target`, `fusion`) or a unified
    legacy `analysis` queue.
    """

    def __init__(self, verbose: bool = False, analysis_workers: int = 1, use_separate_analysis_queues: bool = True, log_level: str = "INFO", preprocessing_workers: int = 1, bed_workers: int = 1, fair_submit: bool = True):
        # Initialize Ray if not already initialized
        if not ray.is_initialized():
            # Suppress Ray logging to prevent interference with progress bars
            self._suppress_ray_logging()
            # Start local Ray with dashboard enabled when possible
            try:
                ray.init(
                    ignore_reinit_error=True,
                    include_dashboard=True,
                    dashboard_host="0.0.0.0",
                    dashboard_port=8265,
                )
            except TypeError:
                # Fallback for Ray versions without dashboard args
                ray.init(ignore_reinit_error=True)
            try:
                url = getattr(ray, "get_dashboard_url", lambda: None)()
                if url and self.verbose:
                    print(f"Ray dashboard: {url}")
            except Exception:
                pass
        
        self.verbose = verbose
        self.use_separate_analysis_queues = use_separate_analysis_queues
        self.analysis_workers_count = analysis_workers
        self.log_level = log_level
        # Allow disabling fairness throttling via arg or env var
        try:
            env_fair = os.environ.get("LJ_FAIR_SUBMIT")
            if env_fair is not None:
                fair_submit = (env_fair != "0")
        except Exception:
            pass
        self.fair_submit = bool(fair_submit)
        # Explicitly announce fairness mode to the CLI
        try:
            src = "env LJ_FAIR_SUBMIT" if env_fair is not None else "argument/default"
            status = "enabled" if self.fair_submit else "disabled (unthrottled)"
            print(f"Scheduling fairness throttling: {status} [{src}]")
        except Exception:
            pass

        # Delta diagnostics (env LJ_DELTA_DEBUG=N; LJ_DELTA_DEBUG_INTERVAL=secs)
        try:
            self.delta_debug_n = int(os.environ.get("LJ_DELTA_DEBUG", "0"))
        except Exception:
            self.delta_debug_n = 0
        try:
            self.delta_debug_interval = float(os.environ.get("LJ_DELTA_DEBUG_INTERVAL", "10"))
        except Exception:
            self.delta_debug_interval = 10.0
        self._last_delta_debug_print_ts = 0.0
        
        # Job handlers for each queue type
        self.job_handlers_preprocessing: Dict[str, Callable] = {}
        self.job_handlers_bed_conversion: Dict[str, Callable] = {}
        
        if use_separate_analysis_queues:
            self.job_handlers_mgmt: Dict[str, Callable] = {}
            self.job_handlers_cnv: Dict[str, Callable] = {}
            self.job_handlers_target: Dict[str, Callable] = {}
            self.job_handlers_fusion: Dict[str, Callable] = {}
        else:
            self.job_handlers_analysis: Dict[str, Callable] = {}
        
        self.job_handlers_classification: Dict[str, Callable] = {}
        self.job_handlers_slow: Dict[str, Callable] = {}
        
        # Worker counts
        self.preprocessing_workers_count = preprocessing_workers if isinstance(preprocessing_workers, int) and preprocessing_workers > 0 else 1
        self.bed_workers_count = bed_workers if isinstance(bed_workers, int) and bed_workers > 0 else 1

        # Ray actors for job processing
        self.processors = {}
        self.running = True
        
        # Job tracking
        self.completed_jobs = []
        self.failed_jobs = []
        self.active_jobs = {}
        self.job_start_times = {}
        self.completed_jobs_by_type = {}
        self.failed_jobs_by_type = {}
        # Track skipped jobs (e.g., deduplicated) by type
        self.skipped_jobs_by_type: Dict[str, int] = {}
        # Track submissions and results by job_id for reconciliation
        self._submissions_by_job_id: Dict[int, Dict[str, Any]] = {}
        self._results_by_job_id: Dict[int, Dict[str, Any]] = {}
        # Track futures to ensure completion is observed even if active map is overwritten
        self._pending_futures: Dict[Any, int] = {}
        # File-centric registry for stable, monotonic counting per display queue
        self._files_registry: Dict[str, Dict[str, Any]] = {}
        self._registry_submitted_by_queue: Dict[str, int] = {k: 0 for k in ['preprocessing','bed_conversion','mgmt','cnv','target','fusion','analysis','classification','slow']}
        self._registry_completed_by_queue: Dict[str, int] = {k: 0 for k in ['preprocessing','bed_conversion','mgmt','cnv','target','fusion','analysis','classification','slow']}
        # Monotonic submitted counters per display queue for accurate totals
        self.submitted_jobs_by_queue: Dict[str, int] = {
            'preprocessing': 0,
            'bed_conversion': 0,
            'mgmt': 0,
            'cnv': 0,
            'target': 0,
            'fusion': 0,
            'analysis': 0,
            'classification': 0,
            'slow': 0,
        }
        self.progress_lock = threading.Lock()
        self.total_jobs_enqueued = 0
        self.total_jobs_skipped = 0
        
        # Job deduplication
        self.deduplicate_job_types: Set[str] = {"sturgeon", "nanodx", "pannanodx", "random_forest"}
        self.running_jobs_by_sample: Dict[str, Dict[str, int]] = {}
        self.pending_jobs_by_sample: Dict[str, Dict[str, int]] = {}
        self.job_tracking_lock = threading.Lock()
        
        # Watchdog: maximum runtime per queue before considering a job stuck (seconds)
        # Keep conservative defaults; preprocessing can stall on malformed BAMs
        self.max_runtime_seconds_by_queue: Dict[str, int] = {
            'preprocessing': 1800,
            'bed_conversion': 1800,
            # Legacy unified analysis queue
            'analysis': 1800,
            # Per-type entries used in legacy-mode display mapping
            'mgmt': 1800,
            'cnv': 1800,
            'target': 1800,
            'fusion': 1800,
            # Classification and slow
            'classification': 1800,
            'slow': 1800,
        }
        
        # Per-sample tracking for GUI/monitoring (parity with threaded manager)
        # sample_id -> {
        #   'sample_id', 'active_jobs', 'total_jobs', 'completed_jobs', 'failed_jobs', 'job_types' (set), 'last_seen'
        # }
        self.samples_by_id: Dict[str, Dict[str, Any]] = {}
        
        # Priority system
        self.queue_priorities = {
            # High priority - critical path jobs
            "preprocessing": 10,      # Highest priority - must complete first
            "bed_conversion": 9,      # High priority - needed for downstream jobs
            
            # Medium priority - analysis jobs that can run in parallel
            "mgmt": 5,
            "cnv": 5,
            "target": 5,
            "fusion": 5,
            "analysis": 5,  # Legacy mode
            
            # Lower priority - classification and slow jobs
            "classification": 3,
            "slow": 1,
        }
        
        # Priority-based job queues (lists hold job_info dicts)
        self.priority_queues = {
            10: [],  # preprocessing
            9: [],   # bed_conversion
            5: [],   # analysis jobs (mgmt/cnv/target/fusion or unified analysis)
            3: [],   # classification
            1: [],   # slow jobs
        }
        
        # Initialize processors
        self._initialize_processors()
    
    def _suppress_ray_logging(self):
        """Suppress Ray logging to prevent interference with progress bars."""
        import logging
        
        # Suppress all Ray-related loggers
        ray_loggers = [
            "ray", "ray.worker", "ray.remote", "ray.actor", "ray.util",
            "ray.raylet", "ray.gcs", "ray.core", "ray.serve"
        ]
        
        for logger_name in ray_loggers:
            logger = logging.getLogger(logger_name)
            logger.setLevel(logging.ERROR)
            # Also disable propagation to prevent messages from bubbling up
            logger.propagate = False
    
    def _initialize_processors(self):
        """
        Create Ray actors per queue with appropriate CPU and concurrency limits.

        - `preprocessing` and `bed_conversion` dedicate a CPU per worker and run
          one job at a time to ensure predictable throughput for critical-path
          steps.
        - Analysis queues can be split (one actor per queue) or unified under the
          legacy `analysis` actor with adjustable logical concurrency.
        - `classification` and `slow` queues each get their own actor; logical
          concurrency is configurable but CPU is kept at 1 to avoid oversubscription.
        """
        # Define concurrency targets
        analysis_concurrency = max(1, self.analysis_workers_count)

        # Preprocessing processors (configurable workers) - dedicate 1 core each
        self.processors['preprocessing'] = [
            JobProcessor.options(name=f"JobProcessor-preprocessing-{i+1}", num_cpus=1, max_concurrency=1).remote(
                "preprocessing", self.job_handlers_preprocessing, self.log_level, self.verbose
            )
            for i in range(self.preprocessing_workers_count)
        ]
        
        # Bed conversion processors (configurable workers) - dedicate 1 core each
        self.processors['bed_conversion'] = [
            JobProcessor.options(name=f"JobProcessor-bed_conversion-{i+1}", num_cpus=1, max_concurrency=1).remote(
                "bed_conversion", self.job_handlers_bed_conversion, self.log_level, self.verbose
            )
            for i in range(self.bed_workers_count)
        ]
        
        # Analysis processors
        if self.use_separate_analysis_queues:
            # Separate processors for each analysis job type
            self.processors['mgmt'] = [
                JobProcessor.options(name=f"JobProcessor-mgmt-{i+1}", num_cpus=1, max_concurrency=1).remote(
                    "mgmt", self.job_handlers_mgmt, self.log_level, self.verbose
                )
                for i in range(self.analysis_workers_count)
            ]
            self.processors['cnv'] = [
                JobProcessor.options(name=f"JobProcessor-cnv-{i+1}", num_cpus=1, max_concurrency=1).remote(
                    "cnv", self.job_handlers_cnv, self.log_level, self.verbose
                )
                for i in range(self.analysis_workers_count)
            ]
            self.processors['target'] = [
                JobProcessor.options(name=f"JobProcessor-target-{i+1}", num_cpus=1, max_concurrency=1).remote(
                    "target", self.job_handlers_target, self.log_level, self.verbose
                )
                for i in range(self.analysis_workers_count)
            ]
            self.processors['fusion'] = [
                JobProcessor.options(name=f"JobProcessor-fusion-{i+1}", num_cpus=1, max_concurrency=1).remote(
                    "fusion", self.job_handlers_fusion, self.log_level, self.verbose
                )
                for i in range(self.analysis_workers_count)
            ]
        else:
            # Single analysis processor
            # Combine mgmt/cnv/target/fusion onto ONE core with concurrent threads
            # Use analysis_workers_count to control concurrency while leaving CPU at 1
            self.processors['analysis'] = [
                JobProcessor.options(name=f"JobProcessor-analysis-1", num_cpus=1, max_concurrency=analysis_concurrency).remote(
                    "analysis", self.job_handlers_analysis, self.log_level, self.verbose
                )
            ]
        
        # Classification processor (1 worker)
        # Give classification its own core but allow concurrent classification jobs
        self.processors['classification'] = [
            JobProcessor.options(name=f"JobProcessor-classification-1", num_cpus=1, max_concurrency=analysis_concurrency).remote(
                "classification", self.job_handlers_classification, self.log_level, self.verbose
            )
        ]
        
        # Slow processor (1 worker)
        # Keep slow tasks isolated on their own core; allow concurrent jobs if desired
        self.processors['slow'] = [
            JobProcessor.options(name=f"JobProcessor-slow-1", num_cpus=1, max_concurrency=analysis_concurrency).remote(
                "slow", self.job_handlers_slow, self.log_level, self.verbose
            )
        ]
    
    def register_handler(self, queue_type: str, job_type: str, handler: Callable[[Job], None]) -> None:
        """
        Register a function to handle a specific `job_type` on a given queue.

        Handlers are injected into the relevant actor(s). The manager itself
        does not implement domain work; it merely routes jobs to the right
        callable.
        """
        if queue_type == "preprocessing":
            self.job_handlers_preprocessing[job_type] = handler
        elif queue_type == "bed_conversion":
            self.job_handlers_bed_conversion[job_type] = handler
        elif queue_type == "analysis":
            if not self.use_separate_analysis_queues:
                self.job_handlers_analysis[job_type] = handler
            else:
                raise ValueError(f"Cannot register 'analysis' queue handler in separate analysis queues mode")
        elif queue_type == "mgmt":
            if self.use_separate_analysis_queues:
                self.job_handlers_mgmt[job_type] = handler
            else:
                raise ValueError(f"Cannot register 'mgmt' queue handler in legacy analysis queue mode")
        elif queue_type == "cnv":
            if self.use_separate_analysis_queues:
                self.job_handlers_cnv[job_type] = handler
            else:
                raise ValueError(f"Cannot register 'cnv' queue handler in legacy analysis queue mode")
        elif queue_type == "target":
            if self.use_separate_analysis_queues:
                self.job_handlers_target[job_type] = handler
            else:
                raise ValueError(f"Cannot register 'target' queue handler in legacy analysis queue mode")
        elif queue_type == "fusion":
            if self.use_separate_analysis_queues:
                self.job_handlers_fusion[job_type] = handler
            else:
                raise ValueError(f"Cannot register 'fusion' queue handler in legacy analysis queue mode")
        elif queue_type == "classification":
            self.job_handlers_classification[job_type] = handler
        elif queue_type == "slow":
            self.job_handlers_slow[job_type] = handler
        else:
            raise ValueError(f"Invalid queue type: {queue_type}")
        
        # Update existing processors with new handler
        self._update_processors_with_handlers()
    
    def _update_processors_with_handlers(self):
        """Update existing Ray actors with new handlers."""
        # This is a simplified approach - in a production system, you might want to
        # recreate the actors or use a more sophisticated handler update mechanism
        if self.verbose:
            print(f"Handlers updated - processors will use new handlers for next jobs")
    
    def _reinitialize_processors(self):
        """Reinitialize Ray actors with current handlers."""
        if self.verbose:
            print("Reinitializing Ray processors with registered handlers...")
        
        # Shutdown existing processors
        for queue_type, processors in self.processors.items():
            for processor in processors:
                try:
                    ray.kill(processor)
                except:
                    pass  # Ignore errors if processor is already dead
        
        # Reinitialize processors with current handlers
        self._initialize_processors()

    # === Per-sample tracking helpers (parity with threaded manager) ===
    def _ensure_sample_entry(self, sample_id: str) -> Dict[str, Any]:
        with self.progress_lock:
            if sample_id not in self.samples_by_id:
                self.samples_by_id[sample_id] = {
                    'sample_id': sample_id,
                    'active_jobs': 0,
                    'total_jobs': 0,
                    'completed_jobs': 0,
                    'failed_jobs': 0,
                    'job_types': set(),
                    'last_seen': time.time(),
                }
            return self.samples_by_id[sample_id]

    def _on_sample_job_started(self, sample_id: str, job_type: str) -> None:
        entry = self._ensure_sample_entry(sample_id)
        with self.progress_lock:
            entry['active_jobs'] += 1
            entry['total_jobs'] += 1
            if job_type:
                entry['job_types'].add(job_type)
            entry['last_seen'] = time.time()
        logging.info(f"Sample {sample_id} started job {job_type}. Active jobs: {entry['active_jobs']}, Total jobs: {entry['total_jobs']}")

    def _on_sample_job_finished(self, sample_id: str, job_type: str, success: bool) -> None:
        entry = self._ensure_sample_entry(sample_id)
        with self.progress_lock:
            if entry['active_jobs'] > 0:
                entry['active_jobs'] -= 1
            if success:
                entry['completed_jobs'] += 1
            else:
                entry['failed_jobs'] += 1
            if job_type:
                entry['job_types'].add(job_type)
            entry['last_seen'] = time.time()
        logging.info(f"Sample {sample_id} finished job {job_type}. Success: {success}. Active jobs: {entry['active_jobs']}, Completed jobs: {entry['completed_jobs']}, Failed jobs: {entry['failed_jobs']}")

    # === File registry helpers (stable counts) ===
    def note_file_seen(self, filepath: str) -> None:
        try:
            if filepath not in self._files_registry:
                self._files_registry[filepath] = {
                    'submitted': set(),
                    'completed': set(),
                }
        except Exception:
            pass

    def _registry_mark_submitted(self, filepath: str, display_queue: str) -> None:
        try:
            self.note_file_seen(filepath)
            entry = self._files_registry[filepath]
            if display_queue not in entry['submitted']:
                entry['submitted'].add(display_queue)
                if display_queue not in self._registry_submitted_by_queue:
                    self._registry_submitted_by_queue[display_queue] = 0
                self._registry_submitted_by_queue[display_queue] += 1
        except Exception:
            pass

    def _registry_mark_completed(self, filepath: str, display_queue: str) -> None:
        try:
            self.note_file_seen(filepath)
            entry = self._files_registry[filepath]
            if display_queue not in entry['completed']:
                entry['completed'].add(display_queue)
                if display_queue not in self._registry_completed_by_queue:
                    self._registry_completed_by_queue[display_queue] = 0
                self._registry_completed_by_queue[display_queue] += 1
        except Exception:
            pass
    
    def add_deduplication_job_type(self, job_type: str) -> None:
        """Add a job type to the deduplication list."""
        self.deduplicate_job_types.add(job_type)
    
    def set_queue_priority(self, queue_type: str, priority: int) -> None:
        """Set the priority for a specific queue type."""
        if priority < 1 or priority > 10:
            raise ValueError("Priority must be between 1 and 10")
        
        self.queue_priorities[queue_type] = priority
        
        # Ensure the priority queue exists
        if priority not in self.priority_queues:
            self.priority_queues[priority] = []
        
        if self.verbose:
            print(f"Set priority for queue '{queue_type}' to {priority}")
    
    def get_queue_priority(self, queue_type: str) -> int:
        """Get the priority for a specific queue type."""
        return self.queue_priorities.get(queue_type, 1)
    
    def get_priority_info(self) -> Dict[str, Any]:
        """Get information about queue priorities and job counts."""
        priority_info = {
            'queue_priorities': self.queue_priorities.copy(),
            'priority_queue_counts': {},
            'active_jobs_by_priority': {}
        }
        
        # Count jobs in each priority queue
        for priority, queue in self.priority_queues.items():
            priority_info['priority_queue_counts'][priority] = len(queue)
        
        # Count active jobs by priority
        for job_id, job_info in self.active_jobs.items():
            priority = job_info.get('priority', 1)
            if priority not in priority_info['active_jobs_by_priority']:
                priority_info['active_jobs_by_priority'][priority] = 0
            priority_info['active_jobs_by_priority'][priority] += 1
        
        return priority_info
    
    def print_priority_status(self) -> None:
        """Print the current status of priority queues."""
        print("\n=== Priority Queue Status ===")
        
        # Print queue priorities
        print("Queue Priorities:")
        for queue_type, priority in sorted(self.queue_priorities.items(), key=lambda x: x[1], reverse=True):
            print(f"  {queue_type}: {priority}")
        
        # Print job counts in each priority queue
        print("\nJobs in Priority Queues:")
        for priority in sorted(self.priority_queues.keys(), reverse=True):
            count = len(self.priority_queues[priority])
            if count > 0:
                print(f"  Priority {priority}: {count} jobs")
        
        # Print active jobs by priority
        active_by_priority = {}
        for job_id, job_info in self.active_jobs.items():
            priority = job_info.get('priority', 1)
            if priority not in active_by_priority:
                active_by_priority[priority] = 0
            active_by_priority[priority] += 1
        
        if active_by_priority:
            print("\nActive Jobs by Priority:")
            for priority in sorted(active_by_priority.keys(), reverse=True):
                print(f"  Priority {priority}: {active_by_priority[priority]} jobs")
        
        print()
    
    def _can_enqueue_job_for_sample(self, job_type: str, sample_id: str) -> bool:
        """Check if a job can be enqueued for the sample (max 2 total with at most 1 pending).

        Mirrors the threaded manager behavior: allow at most one pending job
        per sample/job_type and a combined cap of one running + one pending.
        """
        with self.job_tracking_lock:
            running_count = self.running_jobs_by_sample.get(sample_id, {}).get(job_type, 0)
            pending_count = self.pending_jobs_by_sample.get(sample_id, {}).get(job_type, 0)
            return (running_count + pending_count) < 2 and pending_count < 1
    
    def _mark_job_pending_for_sample(self, job_type: str, sample_id: str) -> None:
        """Mark a job as pending for a sample."""
        with self.job_tracking_lock:
            if sample_id not in self.pending_jobs_by_sample:
                self.pending_jobs_by_sample[sample_id] = {}
            if job_type not in self.pending_jobs_by_sample[sample_id]:
                self.pending_jobs_by_sample[sample_id][job_type] = 0
            self.pending_jobs_by_sample[sample_id][job_type] += 1
    
    def _mark_job_running_for_sample(self, job_type: str, sample_id: str) -> None:
        """Mark a job as running for a sample."""
        with self.job_tracking_lock:
            if sample_id not in self.running_jobs_by_sample:
                self.running_jobs_by_sample[sample_id] = {}
            if job_type not in self.running_jobs_by_sample[sample_id]:
                self.running_jobs_by_sample[sample_id][job_type] = 0
            self.running_jobs_by_sample[sample_id][job_type] += 1
            
            # Decrease pending count
            if sample_id in self.pending_jobs_by_sample and job_type in self.pending_jobs_by_sample[sample_id]:
                self.pending_jobs_by_sample[sample_id][job_type] = max(0, self.pending_jobs_by_sample[sample_id][job_type] - 1)
    
    def _unmark_job_for_sample(self, job_type: str, sample_id: str) -> None:
        """Unmark a job as running for a sample."""
        with self.job_tracking_lock:
            if sample_id in self.running_jobs_by_sample and job_type in self.running_jobs_by_sample[sample_id]:
                self.running_jobs_by_sample[sample_id][job_type] = max(0, self.running_jobs_by_sample[sample_id][job_type] - 1)
    
    def enqueue_jobs(self, jobs: List[Job]) -> None:
        """
        Enqueue one or more jobs into the appropriate priority buckets.

        Additional behaviors:
        - Enforces per-sample de-duplication for certain job types (e.g.,
          classification family) to avoid redundant concurrent work on the same
          sample.
        - Applies a small fairness boost for the classification queue to ensure it
          makes progress even under heavy analysis load.
        """
        for job in jobs:
            # Check deduplication
            if job.job_type in self.deduplicate_job_types:
                sample_id = job.get_sample_id()
                if not self._can_enqueue_job_for_sample(job.job_type, sample_id):
                    # Already have max jobs (1 running + 1 pending) for this sample, skip it
                    try:
                        logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
                        logger.info(f"Skipping {job.job_type} job for sample {sample_id} (max jobs reached: 1 running + 1 pending)")
                    except Exception:
                        pass
                    with self.progress_lock:
                        self.total_jobs_skipped += 1
                        # Track skipped by type
                        if job.job_type not in self.skipped_jobs_by_type:
                            self.skipped_jobs_by_type[job.job_type] = 0
                        self.skipped_jobs_by_type[job.job_type] += 1
                    continue
                self._mark_job_pending_for_sample(job.job_type, sample_id)
            
            # Determine queue type
            queue_type = self._get_queue_type_for_job(job.job_type)
            
            # Get priority for this queue and apply small fairness boost for lower priority queues
            priority = self.queue_priorities.get(queue_type, 1)
            # Encourage some classification throughput even when lots of analysis jobs exist
            if queue_type == 'classification':
                priority = max(priority, 4)
            # Ensure priority bucket exists
            if priority not in self.priority_queues:
                self.priority_queues[priority] = []
            
            # Add job to priority queue
            self.priority_queues[priority].append({
                'job': job,
                'queue_type': queue_type,
                'priority': priority,
                'timestamp': time.time()
            })
            # Diagnostics removed for clean CLI output
            
            # Track the job
            with self.progress_lock:
                self.total_jobs_enqueued += 1
            
            if self.verbose:
                print(f"Enqueued job {job.job_id} ({job.job_type}) with priority {priority}")
        
        # Kick the scheduler to start submitting work immediately
        self._process_priority_queues()
    
    def _process_priority_queues(self) -> None:
        """Process jobs across priorities.
        - Fair mode: one submission per priority per tick.
        - Unthrottled mode: drain all pending jobs high→low in this tick.
        """
        priorities = sorted(self.priority_queues.keys(), reverse=True)
        if self.fair_submit:
            for priority in priorities:
                queue = self.priority_queues[priority]
                if not queue:
                    continue
                job_info = queue[0]
                job = job_info['job']
                queue_type = job_info['queue_type']
                if self._can_process_job(job, queue_type):
                    queue.pop(0)
                    self._submit_job_to_ray(job, queue_type)
        else:
            for priority in priorities:
                queue = self.priority_queues[priority]
                while queue:
                    job_info = queue[0]
                    job = job_info['job']
                    queue_type = job_info['queue_type']
                    if not self._can_process_job(job, queue_type):
                        break
                    queue.pop(0)
                    self._submit_job_to_ray(job, queue_type)
    
    def _can_process_job(self, job: Job, queue_type: str) -> bool:
        """Check if a job can be processed based on available processors."""
        if queue_type in self.processors:
            processors = self.processors[queue_type]
            if processors:
                # Check if any processor is available (not too busy)
                # For now, we'll use a simple check - in a more sophisticated system,
                # we could track processor load
                return True
        return False
    
    def _submit_job_to_ray(self, job: Job, queue_type: str) -> None:
        """
        Submit a serialized `Job` to one of the queue's actors and track it.

        We perform lightweight round-robin across the queue's actors and capture
        bookkeeping necessary for progress bars and watchdog logic.
        """
        processors = self.processors[queue_type]
        if not processors:
            if self.verbose:
                print(f"Warning: No processors available for queue type: {queue_type}")
            return

        # Round-robin selection for load balancing
        processor = processors[len(self.completed_jobs) % len(processors)]

        # Serialize job for Ray
        job_data = pickle.dumps(job)

        # Submit job to Ray; if this fails, requeue the job to avoid dropping it
        try:
            future = processor.process_job.remote(job_data)
        except Exception as submit_exc:
            # Requeue at the head of its priority bucket and log best-effort
            try:
                prio = self.queue_priorities.get(queue_type, 1)
                self.priority_queues.setdefault(prio, []).insert(0, {
                    'job': job,
                    'queue_type': queue_type,
                    'priority': prio,
                    'timestamp': time.time()
                })
            except Exception:
                pass
            if self.verbose:
                print(f"Rescheduled job {job.job_id} ({job.job_type}) after submit error: {submit_exc}")
            return

        # Track the job
        with self.progress_lock:
            sample_id = job.get_sample_id() if hasattr(job, 'get_sample_id') else 'unknown'
            self.active_jobs[job.job_id] = {
                'job_type': job.job_type,
                'filepath': job.context.filepath,
                'sample_id': sample_id,
                'queue_type': queue_type,
                'priority': self.queue_priorities.get(queue_type, 1),
                'start_time': time.time(),
                'processing_start_time': None,  # Will be set when actual processing starts
                'future': future
            }
            self.job_start_times[job.job_id] = time.time()
            # Track future globally for safety
            try:
                self._pending_futures[future] = job.job_id
            except Exception:
                pass
            # Increment submitted counters per display queue for accurate totals
            try:
                display_q = self._get_display_queue_for_job(job.job_type)
                if display_q not in self.submitted_jobs_by_queue:
                    self.submitted_jobs_by_queue[display_q] = 0
                self.submitted_jobs_by_queue[display_q] += 1
                # Stable file-centric submission tracking
                try:
                    self._registry_mark_submitted(job.context.filepath, display_q)
                except Exception:
                    pass
            except Exception:
                pass
            # Track submission by job_id for reconciliation
            try:
                self._submissions_by_job_id[job.job_id] = {
                    'queue_type': queue_type,
                    'display_queue': self._get_display_queue_for_job(job.job_type),
                    'job_type': job.job_type,
                    'filepath': job.context.filepath,
                    'submitted_at': time.time()
                }
            except Exception:
                pass
        # For deduplicated jobs, move from pending -> running now that it's submitted
        if job.job_type in self.deduplicate_job_types and sample_id and sample_id != 'unknown':
            self._mark_job_running_for_sample(job.job_type, sample_id)
        # Per-sample start tracking
        if sample_id and sample_id != 'unknown':
            self._on_sample_job_started(sample_id, job.job_type)

        if self.verbose:
            print(f"Submitted job {job.job_id} ({job.job_type}) to Ray with priority {self.queue_priorities.get(queue_type, 1)}")
    
    def _get_queue_type_for_job(self, job_type: str) -> str:
        """Get the queue type for a given job type."""
        if self.use_separate_analysis_queues:
            queue_mapping = {
                "preprocessing": "preprocessing",
                "bed_conversion": "bed_conversion",
                "mgmt": "mgmt",
                "cnv": "cnv",
                "target": "target",
                "fusion": "fusion",
                "sturgeon": "classification",
                "nanodx": "classification",
                "pannanodx": "classification",
                "random_forest": "slow",
            }
        else:
            queue_mapping = {
                "preprocessing": "preprocessing",
                "bed_conversion": "bed_conversion",
                "mgmt": "analysis",
                "cnv": "analysis",
                "target": "analysis",
                "fusion": "analysis",
                "sturgeon": "classification",
                "nanodx": "classification",
                "pannanodx": "classification",
                "random_forest": "slow",
            }
        return queue_mapping.get(job_type, "slow")

    def _get_display_queue_for_job(self, job_type: str) -> str:
        """Return the display queue key for a job type.

        When using the legacy unified analysis queue, we still want to display
        individual progress bars for mgmt/cnv/target/fusion, so map those job
        types back to their own display queues.
        """
        if (not self.use_separate_analysis_queues) and job_type in {"mgmt", "cnv", "target", "fusion"}:
            return job_type
        return self._get_queue_type_for_job(job_type)
    
    def process_completed_jobs(self) -> None:
        """
        Drain completed Ray futures, propagate triggers, and maintain stats.

        Steps:
        1) Snapshot the active jobs and non-blockingly probe each future.
        2) Apply watchdog timeouts per queue and request queue restarts as needed.
        3) For each completed result:
           - On success: enqueue trigger jobs (preferred) or the linear next job.
           - On error: record failure counts and mark per-sample completion.
        4) Clean up active tracking and update aggregate counters for the UI.
        """
        completed_futures: List[Tuple[int, Dict[str, Any]]] = []

        # Take a stable snapshot of active jobs under lock to avoid concurrent mutation
        with self.progress_lock:
            active_items = list(self.active_jobs.items())
            tracked_futures = list(getattr(self, '_pending_futures', {}).items())

        # Check for completed jobs without mutating dicts during iteration
        for job_id, job_info in active_items:
            future = job_info['future']
            try:
                ready_futures, _ = ray.wait([future], timeout=0)
                if ready_futures:
                    result = ray.get(future)
                    completed_futures.append((job_id, result))
            except Exception as e:
                if self.verbose:
                    print(f"Error checking future for job {job_id}: {e}")
                completed_futures.append((job_id, {
                    'status': 'error',
                    'job_id': job_id,
                    'job_type': job_info.get('job_type', 'unknown'),
                    'error': str(e)
                }))

        # Also check globally tracked futures that might have been overwritten in active map
        for fut, jid in tracked_futures:
            if any(jid == j and fut == info.get('future') for j, info in active_items):
                continue
            try:
                ready_futures, _ = ray.wait([fut], timeout=0)
                if ready_futures:
                    result = ray.get(fut)
                    completed_futures.append((jid, result))
            except Exception as e:
                if self.verbose:
                    print(f"Error checking global future for job {jid}: {e}")
                completed_futures.append((jid, {
                    'status': 'error',
                    'job_id': jid,
                    'job_type': 'unknown',
                    'error': str(e)
                }))

        # Watchdog: detect and fail jobs that exceed max runtime for their queue
        ids_already_completed = {jid for jid, _ in completed_futures}
        need_restart_queues: Set[str] = set()
        now_ts = time.time()
        for job_id, job_info in active_items:
            if job_id in ids_already_completed:
                continue
            jt = job_info.get('job_type', 'unknown')
            q = self._get_queue_type_for_job(jt)
            if not self.use_separate_analysis_queues and jt in {'mgmt','cnv','target','fusion'}:
                q = jt
            max_s = self.max_runtime_seconds_by_queue.get(q)
            if not max_s:
                continue
            start_ts = job_info.get('start_time', now_ts)
            elapsed = now_ts - start_ts
            if elapsed > max_s:
                # Consider this job stuck; mark as failed and request queue restart
                completed_futures.append((job_id, {
                    'status': 'error',
                    'job_id': job_id,
                    'job_type': jt,
                    'error': f'timeout: exceeded {int(max_s)}s without completion',
                }))
                need_restart_queues.add(q)
        
        # Process completed jobs
        for job_id, result in completed_futures:
            # Move fields and remove from active maps under lock
            with self.progress_lock:
                if job_id in self.active_jobs:
                    if 'processing_start_time' in result:
                        self.active_jobs[job_id]['processing_start_time'] = result['processing_start_time']
                    # Remove now that it's completed
                    del self.active_jobs[job_id]
                # Remove start time if present
                if job_id in self.job_start_times:
                    del self.job_start_times[job_id]
                # Remove from global future tracking
                try:
                    for f, jid in list(getattr(self, '_pending_futures', {}).items()):
                        if jid == job_id:
                            self._pending_futures.pop(f, None)
                except Exception:
                    pass
                # Record a result footprint for reconciliation
                try:
                    self._results_by_job_id[job_id] = {
                        'job_type': result.get('job_type', 'unknown'),
                        'status': result.get('status', 'unknown'),
                        'finished_at': time.time()
                    }
                except Exception:
                    pass

            if result['status'] == 'success':
                self.completed_jobs.append(job_id)
                
                # Update context if available
                if 'context' in result:
                    # In a real implementation, you'd want to merge the context
                    # For now, we'll just track the completion
                    pass
                
                # Handle successors. Prefer explicit triggered jobs returned by the actor.
                # If only job types were returned, reconstruct Jobs here to ensure
                # globally unique job IDs (avoid collisions across processes).
                triggered_jobs_to_enqueue: List[Job] = []
                if 'triggered_jobs' in result and result['triggered_jobs']:
                    triggered_jobs_to_enqueue = result['triggered_jobs']
                elif 'triggered_job_types' in result and result['triggered_job_types']:
                    try:
                        ctx = result['context']
                        for jt in result['triggered_job_types']:
                            q = self._get_queue_type_for_job(jt)
                            triggered_jobs_to_enqueue.append(Job(
                                job_id=next(_job_id_counter),
                                job_type=jt,
                                context=ctx,
                                origin=q,
                                workflow=[f"{q}:{jt}"],
                                step=0
                            ))
                    except Exception:
                        triggered_jobs_to_enqueue = []

                if triggered_jobs_to_enqueue:
                    # Log cross-job transition: COMPLETED -> QUEUED for each triggered job
                    try:
                        ctx_for_fp = result.get('context') if isinstance(result, dict) else None
                        fp_name = os.path.basename(ctx_for_fp.filepath) if ctx_for_fp and hasattr(ctx_for_fp, 'filepath') else 'unknown'
                    except Exception:
                        fp_name = 'unknown'
                    for tj in triggered_jobs_to_enqueue:
                        try:
                            logging.info(f"File {fp_name}: {result['job_type']} COMPLETED -> {tj.job_type} QUEUED")
                        except Exception:
                            pass
                    self.enqueue_jobs(triggered_jobs_to_enqueue)
                else:
                    if 'next_job' in result and result['next_job']:
                        maybe_next = result['next_job']
                        # Ensure any provided next job has a fresh id
                        if isinstance(maybe_next, Job):
                            if not isinstance(getattr(maybe_next, 'job_id', None), int):
                                maybe_next.job_id = next(_job_id_counter)
                            # Log cross-job transition for linear progression
                            try:
                                ctx_for_fp2 = result.get('context') if isinstance(result, dict) else None
                                fp_name2 = os.path.basename(ctx_for_fp2.filepath) if ctx_for_fp2 and hasattr(ctx_for_fp2, 'filepath') else 'unknown'
                            except Exception:
                                fp_name2 = 'unknown'
                            try:
                                logging.info(f"File {fp_name2}: {result['job_type']} COMPLETED -> {maybe_next.job_type} QUEUED")
                            except Exception:
                                pass
                            self.enqueue_jobs([maybe_next])
                
                # Track completion by job type
                with self.progress_lock:
                    job_type = result['job_type']
                    if job_type not in self.completed_jobs_by_type:
                        self.completed_jobs_by_type[job_type] = 0
                    self.completed_jobs_by_type[job_type] += 1
                    # Stable file-centric completion tracking
                    try:
                        display_q_res = self._get_display_queue_for_job(job_type)
                        ctx = result.get('context')
                        fp = ctx.filepath if ctx and hasattr(ctx, 'filepath') else None
                        if fp:
                            self._registry_mark_completed(fp, display_q_res)
                    except Exception:
                        pass
                
                # Unmark job for deduplication
                if job_type in self.deduplicate_job_types and 'context' in result:
                    sample_id = result['context'].get_sample_id()
                    self._unmark_job_for_sample(job_type, sample_id)
                # Per-sample finish tracking
                if 'context' in result:
                    sample_id = result['context'].get_sample_id()
                    if sample_id and sample_id != 'unknown':
                        self._on_sample_job_finished(sample_id, result['job_type'], True)
            
            else:
                self.failed_jobs.append(job_id)
                with self.progress_lock:
                    job_type = result.get('job_type', 'unknown')
                    if job_type not in self.failed_jobs_by_type:
                        self.failed_jobs_by_type[job_type] = 0
                    self.failed_jobs_by_type[job_type] += 1
                # Per-sample finish tracking (failure)
                ctx = result.get('context')
                if ctx is not None:
                    sample_id = ctx.get_sample_id()
                    if sample_id and sample_id != 'unknown':
                        self._on_sample_job_finished(sample_id, result.get('job_type', 'unknown'), False)

        # Restart any queues that likely have a wedged actor (best effort)
        for q in need_restart_queues:
            try:
                self._restart_queue_processors(q)
            except Exception as e:
                if self.verbose:
                    print(f"Warning: failed to restart processors for queue '{q}': {e}")

    def _restart_queue_processors(self, queue_type: str) -> None:
        """Restart Ray actors for a specific queue to recover from a stuck handler."""
        procs = self.processors.get(queue_type, [])
        for p in procs:
            try:
                ray.kill(p)
            except Exception:
                pass
        # Recreate actors for this queue using current handlers
        if queue_type == 'preprocessing':
            self.processors['preprocessing'] = [
                JobProcessor.options(num_cpus=1, max_concurrency=1).remote(
                    'preprocessing', self.job_handlers_preprocessing, self.log_level, self.verbose
                )
            ]
    
    def _worker_loop(self) -> None:
        """Main worker loop for processing jobs."""
        while self.running:
            # Process priority queues using the in-memory priority_queues dict
            self._process_priority_queues()
            
            # Process completed jobs
            self.process_completed_jobs()
            
            time.sleep(0.1)  # Small delay to prevent busy waiting
    
    def run(self) -> None:
        """Run the workflow manager."""
        if self.verbose:
            print("Starting Ray-based workflow manager with priority scheduling...")
        
        self.running = True
        
        # Start worker threads
        self.worker_thread = threading.Thread(target=self._worker_loop, daemon=True)
        self.worker_thread.start()
        
        # Wait for worker thread to complete
        self.worker_thread.join()
    
    def stop(self, timeout: float = 30.0) -> bool:
        """Stop the workflow manager."""
        if self.verbose:
            print("Stopping Ray-based workflow manager...")
        
        self.running = False
        
        # Wait for worker thread to finish
        if hasattr(self, 'worker_thread') and self.worker_thread.is_alive():
            self.worker_thread.join(timeout=timeout)
            return not self.worker_thread.is_alive()
        
        return True
    
    def get_stats(self) -> dict:
        """
        Aggregate live metrics for CLI/GUI consumption.

        Returns a nested dict capturing queue sizes, active jobs (both by worker
        and logical queue), and cumulative completions/failures by type and by
        queue. Also includes per-sample summaries to drive list views in the GUI.
        """
        with self.progress_lock:
            # Get active jobs by worker type and also aggregate by queue for CLI progress
            active_by_worker: Dict[str, List[Dict[str, Any]]] = {}
            active_by_queue: Dict[str, List[Dict[str, Any]]] = {
                'preprocessing': [], 'bed_conversion': [], 'mgmt': [], 'cnv': [], 'target': [], 'fusion': [], 'analysis': [], 'classification': [], 'slow': []
            }
            def _worker_name_for_job_type(job_type: str) -> str:
                if job_type == 'preprocessing':
                    return 'PreprocessingWorker-1'
                if job_type == 'bed_conversion':
                    return 'BedConversionWorker-1'
                if job_type == 'mgmt':
                    return 'MGMTWorker-1'
                if job_type == 'cnv':
                    return 'CNVWorker-1'
                if job_type == 'target':
                    return 'TargetWorker-1'
                if job_type == 'fusion':
                    return 'FusionWorker-1'
                if job_type in {'sturgeon', 'nanodx', 'pannanodx'}:
                    return 'ClassificationWorker-1'
                return 'SlowWorker-1'

            for job_id, job_info in self.active_jobs.items():
                job_type = job_info['job_type']
                worker_name = _worker_name_for_job_type(job_type)
                if worker_name not in active_by_worker:
                    active_by_worker[worker_name] = []
                # Use processing_start_time if available (excludes queue waiting time)
                duration = (time.time() - job_info['processing_start_time']) if job_info.get('processing_start_time') else (time.time() - job_info['start_time'])
                sample_id = job_info.get('sample_id') or (os.path.basename(job_info['filepath']).replace('.bam', '') if job_info['filepath'] else 'unknown')
                active_by_worker[worker_name].append({
                    'job_type': job_type,
                    'filepath': job_info['filepath'],
                    'sample_id': sample_id,
                    'duration': duration
                })
                # Aggregate by display queue consistently
                q = self._get_display_queue_for_job(job_type)
                if q not in active_by_queue:
                    active_by_queue[q] = []
                active_by_queue[q].append({
                    'job_type': job_type,
                    'filepath': job_info['filepath'],
                    'sample_id': sample_id,
                    'duration': duration
                })

            # Calculate queue sizes by inspecting priority queues and counting jobs by job_type
            queue_sizes = {
                'preprocessing': 0,
                'bed_conversion': 0,
                'mgmt': 0,
                'cnv': 0,
                'target': 0,
                'fusion': 0,
                'classification': 0,
                'slow': 0,
            }
            for _priority, pending in self.priority_queues.items():
                for job_info in pending:
                    j: Job = job_info['job']
                    jt = j.job_type
                    q = self._get_display_queue_for_job(jt)
                    if q not in queue_sizes:
                        queue_sizes[q] = 0
                    queue_sizes[q] += 1
            
            # Aggregate completions by queue (for CLI progress parity with GUI)
            completed_by_queue: Dict[str, int] = {k: 0 for k in ['preprocessing','bed_conversion','mgmt','cnv','target','fusion','analysis','classification','slow']}
            for jt, count in self.completed_jobs_by_type.items():
                q = self._get_display_queue_for_job(jt)
                completed_by_queue[q] = completed_by_queue.get(q, 0) + (count or 0)

            failed_by_queue: Dict[str, int] = {k: 0 for k in ['preprocessing','bed_conversion','mgmt','cnv','target','fusion','analysis','classification','slow']}
            for jt, count in self.failed_jobs_by_type.items():
                q = self._get_display_queue_for_job(jt)
                failed_by_queue[q] = failed_by_queue.get(q, 0) + (count or 0)

            # Aggregate skipped by queue (dedup skips)
            skipped_by_queue: Dict[str, int] = {k: 0 for k in ['preprocessing','bed_conversion','mgmt','cnv','target','fusion','analysis','classification','slow']}
            for jt, count in self.skipped_jobs_by_type.items():
                q = self._get_display_queue_for_job(jt)
                skipped_by_queue[q] = skipped_by_queue.get(q, 0) + (count or 0)

            # Registry-based submitted/completed totals (file-centric, monotonic)
            submitted_by_queue_reg: Dict[str, int] = self._registry_submitted_by_queue.copy()
            completed_by_queue_reg: Dict[str, int] = self._registry_completed_by_queue.copy()
            # Stabilize submitted for core queues using monotonic counters as an upper bound
            try:
                monotonic_submitted = self.submitted_jobs_by_queue.copy()
            except Exception:
                monotonic_submitted = {}
            for qname in ['preprocessing', 'bed_conversion', 'classification', 'slow']:
                try:
                    submitted_by_queue_reg[qname] = max(
                        int(submitted_by_queue_reg.get(qname, 0) or 0),
                        int(monotonic_submitted.get(qname, 0) or 0)
                    )
                except Exception:
                    pass

            # Reconcile by queue using registry-based submitted/completed counts
            unaccounted_by_queue: Dict[str, int] = {}
            for qname, submitted_q in submitted_by_queue_reg.items():
                completed_q_reg = int(completed_by_queue_reg.get(qname, 0) or 0)
                active_q = len(active_by_queue.get(qname, [])) if isinstance(active_by_queue.get(qname, []), list) else 0
                pending_q = queue_sizes.get(qname, 0) or 0
                unacc = submitted_q - (completed_q_reg + active_q + pending_q)
                if unacc:
                    unaccounted_by_queue[qname] = unacc

            # Prepare samples payload for GUI
            samples_payload: List[Dict[str, Any]] = []
            for sid, info in self.samples_by_id.items():
                samples_payload.append({
                    'sample_id': info['sample_id'],
                    'active_jobs': info['active_jobs'],
                    'total_jobs': info['total_jobs'],
                    'completed_jobs': info['completed_jobs'],
                    'failed_jobs': info['failed_jobs'],
                    'job_types': list(info['job_types']),
                    'last_seen': info['last_seen'],
                })
            
            # Overall processed should include both completed and failed
            overall_processed = len(self.completed_jobs) + len(self.failed_jobs)
            # Derive total from submission registry for consistency across views
            overall_total = 0
            try:
                overall_total = sum(int(v or 0) for v in submitted_by_queue_reg.values())
            except Exception:
                overall_total = self.total_jobs_enqueued

            # Optionally print delta diagnostics for unaccounted job_ids
            try:
                if self.delta_debug_n and (time.time() - self._last_delta_debug_print_ts) >= self.delta_debug_interval:
                    # Build sets for accounted job_ids
                    active_jids_dbg = set(self.active_jobs.keys())
                    pending_jids_dbg: Set[int] = set()
                    for _prio, pending in self.priority_queues.items():
                        for ji in pending:
                            j_dbg: Job = ji['job']
                            pending_jids_dbg.add(j_dbg.job_id)
                    result_jids_dbg = set(self._results_by_job_id.keys())
                    # For each queue with Δ, print up to N examples
                    for qname, delta_val in unaccounted_by_queue.items():
                        if delta_val <= 0:
                            continue
                        examples = []
                        for jid, sub in self._submissions_by_job_id.items():
                            dq = sub.get('display_queue') or sub.get('queue_type') or 'slow'
                            if dq != qname:
                                continue
                            if jid in active_jids_dbg or jid in pending_jids_dbg or jid in result_jids_dbg:
                                continue
                            fp = sub.get('filepath', '')
                            jt = sub.get('job_type', '')
                            examples.append(f"{jid}:{jt}:{os.path.basename(fp)}")
                            if len(examples) >= self.delta_debug_n:
                                break
                        if examples:
                            try:
                                logging.warning(f"Δ debug [{qname}]: {delta_val} unaccounted. Examples: {', '.join(examples)}")
                            except Exception:
                                print(f"Δ debug [{qname}]: {delta_val} unaccounted. Examples: {', '.join(examples)}", flush=True)
                    self._last_delta_debug_print_ts = time.time()
            except Exception:
                pass

            return {
                'total_jobs_enqueued': self.total_jobs_enqueued,
                'total_jobs_skipped': getattr(self, 'total_jobs_skipped', 0),
                'total_processed': overall_processed,
                'total_actual_jobs': overall_total,  # Derived from submitted_by_queue for stability
                'active_jobs': len(self.active_jobs),
                'use_separate_analysis_queues': self.use_separate_analysis_queues,
                # Provide both legacy and explicit keys used by GUI hooks
                'completed': len(self.completed_jobs),
                'failed': len(self.failed_jobs),
                'completed_jobs': len(self.completed_jobs),
                'failed_jobs': len(self.failed_jobs),
                'skipped': getattr(self, 'total_jobs_skipped', 0),
                'completed_by_type': self.completed_jobs_by_type.copy(),
                'failed_by_type': self.failed_jobs_by_type.copy(),
                'skipped_by_type': self.skipped_jobs_by_type.copy(),
                'completed_by_queue': completed_by_queue,
                'failed_by_queue': failed_by_queue,
                'skipped_by_queue': skipped_by_queue,
                'unaccounted_by_queue': unaccounted_by_queue,
                # Submitted totals derived from job_id registry
                'submitted_by_queue': submitted_by_queue_reg,
                'active_by_worker': active_by_worker,
                'active_by_queue': active_by_queue,
                'queue_sizes': queue_sizes,
                'samples': samples_payload
            }

    def is_running(self) -> bool:
        """Check if the workflow manager is running."""
        return self.running


# === Ray-based File Watcher ===
class RayFileWatcher(FileSystemEventHandler):
    """
    Filesystem watcher that discovers inputs and seeds the workflow.

    The watcher calls a `preprocessor_func` (typically a classifier) for each
    relevant file to create the initial job(s). It can process existing files
    at startup and then continue to watch for new/modified files.
    """
    
    def __init__(
        self,
        watch_dir: str,
        preprocessor_func: Callable[[str], List[Job]],
        manager: RayWorkflowManager,
        recursive: bool = True,
        patterns: Optional[List[str]] = None,
        ignore_patterns: Optional[List[str]] = None,
        verbose: bool = False,
        show_progress: bool = True
    ):
        self.watch_dir = watch_dir
        self.preprocessor_func = preprocessor_func
        self.manager = manager
        self.recursive = recursive
        self.patterns = patterns or ["*.bam"]
        self.ignore_patterns = ignore_patterns or []
        self.verbose = verbose
        self.show_progress = show_progress
        self.observer = Observer()
        self.processed_files = set()
    
    def _should_process_file(self, filepath: str) -> bool:
        """Check if a file should be processed based on patterns (Path.match)."""
        if filepath in self.processed_files:
            return False

        path_obj = Path(filepath)

        # Ignore patterns (use Path.match for glob semantics)
        for pattern in self.ignore_patterns:
            try:
                if path_obj.match(pattern):
                    return False
            except Exception:
                # If a pattern is malformed, skip it without blocking processing
                continue

        # Include patterns (match any)
        if self.patterns:
            for pattern in self.patterns:
                try:
                    if path_obj.match(pattern):
                        return True
                except Exception:
                    continue
            return False

        # If no patterns specified, allow all
        return True
    
    def handle_file(self, filepath: str) -> None:
        """Handle a file by creating and enqueuing jobs."""
        if self._should_process_file(filepath):
            try:
                jobs = self.preprocessor_func(filepath)
                if jobs:
                    # Diagnostics removed for clean CLI output
                    self.manager.enqueue_jobs(jobs)
                    self.processed_files.add(filepath)
                    if self.verbose:
                        print(f"Enqueued {len(jobs)} jobs for {filepath}")
            except Exception as e:
                if self.verbose:
                    print(f"Error processing file {filepath}: {e}")
    
    def on_created(self, event) -> None:
        """Handle file creation events."""
        if not event.is_directory:
            self.handle_file(event.src_path)
    
    def on_modified(self, event) -> None:
        """Handle file modification events."""
        if not event.is_directory:
            self.handle_file(event.src_path)
    
    def start(self, process_existing: bool = True) -> None:
        """Start watching the directory."""
        if process_existing:
            self._process_existing_files()
        
        self.observer.schedule(self, self.watch_dir, recursive=self.recursive)
        self.observer.start()
    
    def _process_existing_files(self) -> None:
        """Process existing files in the watch directory."""
        if self.verbose:
            print(f"Processing existing files in {self.watch_dir}...")
        
        for root, dirs, files in os.walk(self.watch_dir):
            for file in files:
                if self._should_process_file(os.path.join(root, file)):
                    self.handle_file(os.path.join(root, file))
    
    def stop(self, timeout: float = 30.0) -> bool:
        """Stop watching the directory."""
        self.observer.stop()
        return self.observer.join(timeout=timeout)


# === Ray-based Workflow Runner ===
class RayWorkflowRunner:
    """
    High-level façade that ties together the watcher, manager, and progress UI.

    Typical usage:
    - Register handlers with `register_handler()`.
    - Call `run_workflow(...)` with a directory and a `workflow_plan`.
    - Optionally enable progress bars and provide a custom classifier function.
    """
    
    def __init__(self, verbose: bool = False, analysis_workers: int = 1, use_separate_analysis_queues: bool = True, log_level: str = "INFO", preprocessing_workers: int = 1, bed_workers: int = 1):
        self.verbose = verbose
        self.log_level = log_level
        self.manager = RayWorkflowManager(
            verbose=verbose,
            analysis_workers=analysis_workers,
            use_separate_analysis_queues=use_separate_analysis_queues,
            log_level=log_level,
            preprocessing_workers=preprocessing_workers,
            bed_workers=bed_workers
        )
        self._watcher: Optional[RayFileWatcher] = None
    
    def register_handler(self, queue_type: str, job_type: str, handler: Callable[[Job], None]) -> None:
        """Register a handler for a specific job type and queue."""
        self.manager.register_handler(queue_type, job_type, handler)
    
    def register_command_handler(self, queue_type: str, job_type: str, command_template: str) -> None:
        """Register a command handler."""
        def handler(job: Job) -> None:
            # Implementation for command handling
            pass
        self.manager.register_handler(queue_type, job_type, handler)
    
    def run_workflow(
        self,
        watch_dir: str,
        workflow_plan: List[str],
        recursive: bool = True,
        patterns: Optional[List[str]] = None,
        ignore_patterns: Optional[List[str]] = None,
        classifier_func: Optional[Callable] = None,
        process_existing: bool = True,
        show_progress: bool = True
    ) -> None:
        """
        Run the workflow end-to-end for files in `watch_dir` using Ray.

        The first step for each discovered file is created by the provided
        `classifier_func` (defaults to `default_file_classifier`). All downstream
        progression is handled automatically by job triggers and the linear
        fallback.
        """
        # Use the provided classifier or default
        if classifier_func is None:
            classifier_func = lambda filepath: default_file_classifier(filepath, workflow_plan)
        
        # Create file watcher
        watcher = RayFileWatcher(
            watch_dir=watch_dir,
            preprocessor_func=classifier_func,
            manager=self.manager,
            recursive=recursive,
            patterns=patterns,
            ignore_patterns=ignore_patterns,
            verbose=self.verbose,
            show_progress=show_progress
        )
        self._watcher = watcher
        
        # Start the workflow
        try:
            watcher.start(process_existing=process_existing)
            
            # Start progress monitoring in a separate thread if show_progress is enabled
            if show_progress:
                progress_thread = threading.Thread(target=self._monitor_progress, args=(watcher,), daemon=True)
                progress_thread.start()
            
            self.manager.run()
        except KeyboardInterrupt:
            if self.verbose:
                print("Stopping workflow...")
            self.manager.stop(timeout=5.0)
            try:
                watcher.stop(timeout=5.0)
            except Exception:
                pass
            if self.verbose:
                print("Workflow stopped by user")
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
            self.manager.stop(timeout=5.0)
            try:
                watcher.stop(timeout=5.0)
            except Exception:
                pass
            raise

    def _monitor_progress(self, watcher) -> None:
        """Monitor and display worker progress in real-time."""
        import time
        
        # Suppress Ray logging to prevent interference with tqdm
        import logging
        ray_logger = logging.getLogger("ray")
        ray_logger.setLevel(logging.ERROR)
        
        # Also suppress other potentially noisy loggers
        logging.getLogger("ray.worker").setLevel(logging.ERROR)
        logging.getLogger("ray.remote").setLevel(logging.ERROR)
        logging.getLogger("ray.actor").setLevel(logging.ERROR)
        
        # Create progress bars for each worker type
        preprocessing_pbar = tqdm(
            desc="Preprocessing",
            unit="jobs",
            position=0,
            leave=True
        )
        bed_conversion_pbar = tqdm(
            desc="Bed Conversion",
            unit="jobs",
            position=1,
            leave=True
        )
        mgmt_pbar = tqdm(
            desc="MGMT",
            unit="jobs",
            position=2,
            leave=True
        )
        cnv_pbar = tqdm(
            desc="CNV",
            unit="jobs",
            position=3,
            leave=True
        )
        target_pbar = tqdm(
            desc="Target",
            unit="jobs",
            position=4,
            leave=True
        )
        fusion_pbar = tqdm(
            desc="Fusion",
            unit="jobs",
            position=5,
            leave=True
        )
        classification_pbar = tqdm(
            desc="Classification",
            unit="jobs",
            position=6,
            leave=True
        )
        slow_pbar = tqdm(
            desc="Slow",
            unit="jobs",
            position=7,
            leave=True
        )
        
        # Create overall progress bar
        overall_pbar = tqdm(
            desc="Overall Progress",
            unit="jobs",
            position=8,
            leave=True
        )
        
        progress_bars = {
            "preprocessing": preprocessing_pbar,
            "bed_conversion": bed_conversion_pbar,
            "mgmt": mgmt_pbar,
            "cnv": cnv_pbar,
            "target": target_pbar,
            "fusion": fusion_pbar,
            "classification": classification_pbar,
            "slow": slow_pbar
        }
        
        # Track job type completion counts
        job_type_completion = {
            "preprocessing": 0,
            "bed_conversion": 0,
            "mgmt": 0,
            "cnv": 0,
            "target": 0,
            "fusion": 0,
            "classification": 0,
            "slow": 0
        }
        
        last_stats = None
        
        try:
            while self.manager.is_running():
                stats = self.manager.get_stats()
                
                # Update overall progress using total actual jobs (excluding skipped)
                total_processed = stats.get("total_processed", 0)
                total_actual_jobs = stats.get("total_actual_jobs", 0)
                
                # Ensure we have valid numbers
                if total_processed is None:
                    total_processed = 0
                if total_actual_jobs is None:
                    total_actual_jobs = 0
                
                if total_actual_jobs > 0:
                    overall_pbar.total = total_actual_jobs
                    overall_pbar.n = total_processed
                    overall_failed = stats.get('failed', 0) or 0
                    overall_skipped = stats.get('total_jobs_skipped', 0) or 0
                    unaccounted_map = stats.get('unaccounted_by_queue', {}) or {}
                    unacc_total = sum(int(v or 0) for v in unaccounted_map.values())
                    overall_pbar.set_postfix_str(
                        f"Active: {stats.get('active_jobs', 0)} | "
                        f"Done: {total_processed}/{total_actual_jobs} | "
                        f"Failed: {overall_failed} | Skipped: {overall_skipped} | "
                        f"Δ: {unacc_total}"
                    )
                
                # Update individual progress bars using per-queue aggregates
                for queue_type, pbar in progress_bars.items():
                    completed_q = 0
                    if "completed_by_queue" in stats:
                        completed_q += stats["completed_by_queue"].get(queue_type, 0) or 0
                    if "failed_by_queue" in stats:
                        completed_q += stats["failed_by_queue"].get(queue_type, 0) or 0

                    qsize = stats.get("queue_sizes", {}).get(queue_type, 0) or 0
                    active_in_q = len(stats.get("active_by_queue", {}).get(queue_type, [])) if "active_by_queue" in stats else 0

                    # Prefer monotonic submitted totals when available for stable progress bars
                    submitted_totals = stats.get("submitted_by_queue", {}) if isinstance(stats.get("submitted_by_queue", {}), dict) else {}
                    if submitted_totals:
                        total_for_queue = submitted_totals.get(queue_type, 0) or (completed_q + active_in_q + qsize)
                    else:
                        total_for_queue = completed_q + active_in_q + qsize
                    pbar.total = total_for_queue
                    pbar.n = min(completed_q, total_for_queue)

                    # Show first active job summary if any, else clear
                    if active_in_q > 0:
                        first_job = stats["active_by_queue"][queue_type][0]
                        sid = first_job.get("sample_id", "unknown")
                        dur = first_job.get("duration", 0)
                        failed_q = stats.get("failed_by_queue", {}).get(queue_type, 0) or 0
                        skipped_q = stats.get("skipped_by_queue", {}).get(queue_type, 0) or 0
                        unacc_q = stats.get("unaccounted_by_queue", {}).get(queue_type, 0) or 0
                        pbar.set_postfix_str(f"{sid}({dur:.0f}s) F:{failed_q} S:{skipped_q} Δ:{unacc_q}")
                    else:
                        failed_q = stats.get("failed_by_queue", {}).get(queue_type, 0) or 0
                        skipped_q = stats.get("skipped_by_queue", {}).get(queue_type, 0) or 0
                        unacc_q = stats.get("unaccounted_by_queue", {}).get(queue_type, 0) or 0
                        if failed_q or skipped_q or unacc_q:
                            pbar.set_postfix_str(f"F:{failed_q} S:{skipped_q} Δ:{unacc_q}")
                        else:
                            pbar.set_postfix_str("")
                    job_type_completion[queue_type] = completed_q
                
                # Update queue sizes display
                if stats.get("queue_sizes"):
                    queue_info = []
                    # In legacy mode, expand analysis queue into per-type sizes when possible
                    q_sizes = stats["queue_sizes"].copy()
                    if not stats.get("use_separate_analysis_queues"):
                        # Prefer active_by_queue counts as a proxy for per-type visibility
                        for qname in ['mgmt','cnv','target','fusion']:
                            size_est = len(stats.get('active_by_queue', {}).get(qname, []))
                            if size_est:
                                q_sizes[qname] = q_sizes.get(qname, 0) + size_est
                    for queue_type, size in q_sizes.items():
                        if size and size > 0:
                            queue_info.append(f"{queue_type}:{size}")
                    
                    if queue_info:
                        overall_pbar.set_postfix_str(
                            f"Active: {stats.get('active_jobs', 0)} | "
                            f"Completed: {total_processed}/{total_actual_jobs} | "
                            f"Queues: {', '.join(queue_info)}"
                        )
                
                time.sleep(0.5)  # Update every 500ms
                
        except Exception as e:
            print(f"Progress monitoring error: {e}")
        finally:
            # Close all progress bars
            for pbar in progress_bars.values():
                pbar.close()
            overall_pbar.close()

    def stop(self, timeout: float = 5.0) -> bool:
        """Gracefully stop the workflow (manager + watcher)."""
        ok = True
        try:
            ok = self.manager.stop(timeout=timeout) and ok
        except Exception:
            ok = False
        try:
            if self._watcher is not None:
                self._watcher.stop(timeout=timeout)
        except Exception:
            ok = False
        return ok


def default_file_classifier(filepath: str, workflow_plan: List[str]) -> List[Job]:
    """
    Create the initial job for a file and stash the full workflow plan.

    This mirrors the behavior of the threaded manager: we enqueue only the first
    step and rely on in-band triggers to progress the workflow in-order.

    Example:
        workflow_plan = [
            "preprocessing:preprocessing",
            "bed_conversion:bed_conversion",
            "mgmt:mgmt",
            "cnv:cnv",
        ]
        jobs = default_file_classifier("/data/sample.bam", workflow_plan)
        # -> returns one Job for "preprocessing" with context that embeds the plan
    """
    jobs: List[Job] = []
    context = WorkflowContext(filepath=filepath)
    
    # Persist the full plan for trigger logic
    context.add_metadata("original_workflow", workflow_plan)
    
    if not workflow_plan:
        return jobs
    
    first_step = workflow_plan[0]
    if ':' not in first_step:
        # Defensive: ignore malformed steps
        return jobs
    
    queue_type, job_type = first_step.split(':', 1)
    jobs.append(Job(
        job_id=next(_job_id_counter),
        job_type=job_type,
        context=context,
        origin=queue_type,
        workflow=workflow_plan,
        step=0,
    ))
    
    return jobs

logging.basicConfig(
    filename='workflow_ray.log',  # Log file name
    level=logging.INFO,           # Log level
    format='%(asctime)s - %(levelname)s - %(message)s'  # Log format
)
