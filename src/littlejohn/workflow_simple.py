"""Simple workflow management system for LittleJohn using threading."""

import os
import time
import threading
import queue
import itertools
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Any, Set

from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
from tqdm import tqdm
from littlejohn.logging_config import get_job_logger


# === Shared Context for Each File ===
@dataclass
class WorkflowContext:
    """Context object that tracks metadata and results for each file being processed."""

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
        self.errors.append(
            {"job_type": job_type, "error": error, "timestamp": time.time()}
        )

    def get_summary(self) -> dict:
        """Get a summary of the workflow execution."""
        return {
            "filepath": self.filepath,
            "metadata": self.metadata,
            "results": self.results,
            "history": self.history,
            "errors": self.errors,
            "success": len(self.errors) == 0,
        }

    def get_sample_id(self) -> str:
        """Get the sample ID from metadata."""
        bam_metadata = self.metadata.get("bam_metadata", {})
        return bam_metadata.get("sample_id", "unknown")

    def clear_all(self) -> None:
        """Clear all metadata/results/errors/history for this context.
        Use after all jobs for this file are finished to allow GC to reclaim memory.
        """
        self.metadata.clear()
        self.results.clear()
        self.history.clear()
        self.errors.clear()


# === Job Object ===
@dataclass
class Job:
    """Represents a single job in the workflow."""

    job_id: int
    job_type: str
    context: WorkflowContext
    origin: str  # "fast" or "slow"
    workflow: List[str]
    step: int = 0
    dependencies: Set[str] = field(default_factory=set)  # Job types this job depends on
    triggers: Set[str] = field(default_factory=set)  # Job types this job can trigger

    def next_job(self) -> Optional["Job"]:
        """Create the next job in the workflow if available."""
        if self.step + 1 < len(self.workflow):
            next_step = self.workflow[self.step + 1]
            # Ensure we have a string to work with
            if not isinstance(next_step, str):
                raise ValueError(
                    f"Expected string in workflow step, got {type(next_step)}: {next_step}"
                )

            if ":" not in next_step:
                raise ValueError(
                    f"Invalid workflow step format (missing ':'): {next_step}"
                )

            queue_type, next_type = next_step.split(":", 1)
            return Job(
                job_id=self.job_id,
                job_type=next_type,
                context=self.context,
                origin=queue_type,
                workflow=self.workflow,
                step=self.step + 1,
            )
        return None

    def get_sample_id(self) -> str:
        """Get the sample ID for this job."""
        return self.context.get_sample_id()

    def can_trigger_jobs(self) -> List["Job"]:
        """Create all jobs that can be triggered by this job's completion."""
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
            # Extract job types from the workflow plan
            workflow_job_types = []
            for step in original_workflow:
                if isinstance(step, str) and ":" in step:
                    queue_type, job_type = step.split(":", 1)
                    workflow_job_types.append(job_type)
            original_workflow = workflow_job_types
        else:
            # Fallback: extract from current workflow
            original_workflow = []
            for step in self.workflow:
                if isinstance(step, str) and ":" in step:
                    queue_type, job_type = step.split(":", 1)
                    original_workflow.append(job_type)

        # Check if this job can trigger other jobs
        for job_type, dependencies in job_dependencies.items():
            if self.job_type in dependencies:
                # This job completion can trigger the dependent job
                # Check if all dependencies are met
                all_deps_met = True
                for dep in dependencies:
                    if dep not in self.context.history and dep != self.job_type:
                        all_deps_met = False
                        break

                # Only trigger if the job is in the original workflow plan
                if all_deps_met and job_type in original_workflow:
                    # Create the triggered job
                    queue_type = self._get_queue_type_for_job(job_type)
                    triggered_job = Job(
                        job_id=next(_job_id_counter),
                        job_type=job_type,
                        context=self.context,
                        origin=queue_type,
                        workflow=[f"{queue_type}:{job_type}"],  # Single step workflow
                        step=0,
                        dependencies=dependencies,
                        triggers=set(),  # No further triggers for now
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


# === Workflow Manager ===
class WorkflowManager:
    """Manages job queues and execution using threading with specialized workers."""

    def __init__(
        self,
        verbose: bool = False,
        analysis_workers: int = 1,
        use_separate_analysis_queues: bool = True,
        preprocessing_workers: int = 1,
        bed_conversion_workers: int = 1,
    ):
        # Specialized queues for different task categories
        self.preprocessing_queue = queue.Queue()  # bam_preprocessing only
        self.bed_conversion_queue = queue.Queue()  # bed_conversion only

        # Analysis queue mode
        self.use_separate_analysis_queues = use_separate_analysis_queues

        if use_separate_analysis_queues:
            # Separate queues for each analysis job type (ensures only one of each type runs at a time)
            self.mgmt_queue = queue.Queue()  # MGMT analysis only
            self.cnv_queue = queue.Queue()  # CNV analysis only
            self.target_queue = queue.Queue()  # Target analysis only
            self.fusion_queue = queue.Queue()  # Fusion analysis only
            self.analysis_queue = None  # Not used in separate mode
        else:
            # Legacy single analysis queue (all analysis types share one queue)
            self.analysis_queue = queue.Queue()  # mgmt, cnv, target, fusion
            self.mgmt_queue = None
            self.cnv_queue = None
            self.target_queue = None
            self.fusion_queue = None

        self.classification_queue = queue.Queue()  # sturgeon, nanodx, pannanodx
        self.slow_queue = queue.Queue()  # slow jobs (legacy)

        self.completed_jobs = []
        self.failed_jobs = []

        # Analysis workers per queue (each analysis type gets its own worker)
        self.analysis_workers_count = analysis_workers
        # Configurable workers for preprocessing and bed conversion (defaults to 1)
        self.preprocessing_workers_count = (
            preprocessing_workers
            if isinstance(preprocessing_workers, int) and preprocessing_workers > 0
            else 1
        )
        self.bed_conversion_workers_count = (
            bed_conversion_workers
            if isinstance(bed_conversion_workers, int) and bed_conversion_workers > 0
            else 1
        )
        self.classification_workers_count = 1
        self.slow_workers_count = 1

        # Specialized job handlers for each queue
        self.job_handlers_preprocessing: Dict[str, Callable[[Job], None]] = {}
        self.job_handlers_bed_conversion: Dict[str, Callable[[Job], None]] = {}

        # Analysis job handlers (mode-dependent)
        if use_separate_analysis_queues:
            # Separate handlers for each analysis job type
            self.job_handlers_mgmt: Dict[str, Callable[[Job], None]] = {}
            self.job_handlers_cnv: Dict[str, Callable[[Job], None]] = {}
            self.job_handlers_target: Dict[str, Callable[[Job], None]] = {}
            self.job_handlers_fusion: Dict[str, Callable[[Job], None]] = {}
            self.job_handlers_analysis = None  # Not used in separate mode
        else:
            # Legacy single analysis handler
            self.job_handlers_analysis: Dict[str, Callable[[Job], None]] = {}
            self.job_handlers_mgmt = None
            self.job_handlers_cnv = None
            self.job_handlers_target = None
            self.job_handlers_fusion = None

        self.job_handlers_classification: Dict[str, Callable[[Job], None]] = {}
        self.job_handlers_slow: Dict[str, Callable[[Job], None]] = {}

        self.running = True
        self.verbose = verbose

        # Job deduplication tracking
        # Jobs that should be deduplicated by sample ID (e.g., sturgeon analysis)
        self.deduplicate_job_types: Set[str] = {
            "sturgeon",
            "nanodx",
            "pannanodx",
            "random_forest",
        }
        # Track running and pending jobs by sample ID for deduplication
        # Allow max 2 jobs per sample: 1 running + 1 pending
        self.running_jobs_by_sample: Dict[str, Dict[str, int]] = (
            {}
        )  # sample_id -> {job_type: count}
        self.pending_jobs_by_sample: Dict[str, Dict[str, int]] = (
            {}
        )  # sample_id -> {job_type: count}

        self.job_tracking_lock = threading.Lock()

        # Track worker threads for graceful shutdown
        self.preprocessing_workers = []
        self.bed_conversion_workers = []

        # Analysis workers (mode-dependent)
        if use_separate_analysis_queues:
            # Separate workers for each analysis job type
            self.mgmt_workers = []
            self.cnv_workers = []
            self.target_workers = []
            self.fusion_workers = []
            self.analysis_workers = []  # Not used in separate mode
        else:
            # Legacy single analysis workers
            self.analysis_workers = []
            self.mgmt_workers = []
            self.cnv_workers = []
            self.target_workers = []
            self.fusion_workers = []

        self.classification_workers = []
        self.slow_workers = []

        # Progress tracking
        self.total_jobs_enqueued = 0
        self.total_jobs_skipped = 0  # Track jobs skipped due to deduplication
        self.active_jobs = {}  # job_id -> job_info
        self.job_start_times = {}  # job_id -> start_time
        self.completed_jobs_by_type = (
            {}
        )  # job_type -> count (successful completions only)
        self.failed_jobs_by_type = {}  # job_type -> count (failed jobs only)
        self.progress_lock = threading.Lock()  # Lock for progress tracking

        # Per-sample tracking for GUI/monitoring
        # sample_id -> {
        #   'sample_id', 'active_jobs', 'total_jobs', 'completed_jobs', 'failed_jobs', 'job_types' (set), 'last_seen'
        # }
        self.samples_by_id: Dict[str, Dict[str, Any]] = {}

        # Define job type mappings to queues
        if use_separate_analysis_queues:
            self.job_queue_mapping = {
                # Preprocessing queue
                "preprocessing": "preprocessing",
                # Bed conversion queue
                "bed_conversion": "bed_conversion",
                # Separate analysis queues (ensures only one of each type runs at a time)
                "mgmt": "mgmt",
                "cnv": "cnv",
                "target": "target",
                "fusion": "fusion",
                "test": "mgmt",  # For testing purposes - map to mgmt queue
                "long": "cnv",  # For testing purposes - map to cnv queue
                "quick": "target",  # For testing purposes - map to target queue
                # Classification queue
                "sturgeon": "classification",
                "nanodx": "classification",
                "pannanodx": "classification",
                # Slow queue
                "random_forest": "slow",
                "sleep": "slow",
                "echo": "slow",
            }
        else:
            self.job_queue_mapping = {
                # Preprocessing queue
                "preprocessing": "preprocessing",
                # Bed conversion queue
                "bed_conversion": "bed_conversion",
                # Legacy single analysis queue
                "mgmt": "analysis",
                "cnv": "analysis",
                "target": "analysis",
                "fusion": "analysis",
                "test": "analysis",  # For testing purposes
                "long": "analysis",  # For testing purposes
                "quick": "analysis",  # For testing purposes
                # Classification queue
                "sturgeon": "classification",
                "nanodx": "classification",
                "pannanodx": "classification",
                # Slow queue
                "random_forest": "slow",
                "sleep": "slow",
                "echo": "slow",
            }

    def register_handler(
        self, queue_type: str, job_type: str, handler: Callable[[Job], None]
    ) -> None:
        """Register a handler for a specific job type and queue."""
        if queue_type == "preprocessing":
            self.job_handlers_preprocessing[job_type] = handler
        elif queue_type == "bed_conversion":
            self.job_handlers_bed_conversion[job_type] = handler
        elif queue_type == "analysis":
            if not self.use_separate_analysis_queues:
                self.job_handlers_analysis[job_type] = handler
            else:
                raise ValueError(
                    f"Cannot register 'analysis' queue handler in separate analysis queues mode"
                )
        elif queue_type == "mgmt":
            if self.use_separate_analysis_queues:
                self.job_handlers_mgmt[job_type] = handler
            else:
                raise ValueError(
                    f"Cannot register 'mgmt' queue handler in legacy analysis queue mode"
                )
        elif queue_type == "cnv":
            if self.use_separate_analysis_queues:
                self.job_handlers_cnv[job_type] = handler
            else:
                raise ValueError(
                    f"Cannot register 'cnv' queue handler in legacy analysis queue mode"
                )
        elif queue_type == "target":
            if self.use_separate_analysis_queues:
                self.job_handlers_target[job_type] = handler
            else:
                raise ValueError(
                    f"Cannot register 'target' queue handler in legacy analysis queue mode"
                )
        elif queue_type == "fusion":
            if self.use_separate_analysis_queues:
                self.job_handlers_fusion[job_type] = handler
            else:
                raise ValueError(
                    f"Cannot register 'fusion' queue handler in legacy analysis queue mode"
                )
        elif queue_type == "classification":
            self.job_handlers_classification[job_type] = handler
        elif queue_type == "slow":
            self.job_handlers_slow[job_type] = handler
        else:
            raise ValueError(f"Invalid queue type: {queue_type}")

    def add_deduplication_job_type(self, job_type: str) -> None:
        """Add a job type to the deduplication list."""
        self.deduplicate_job_types.add(job_type)

    def remove_deduplication_job_type(self, job_type: str) -> None:
        """Remove a job type from the deduplication list."""
        self.deduplicate_job_types.discard(job_type)

    def _can_enqueue_job_for_sample(self, job_type: str, sample_id: str) -> bool:
        """Check if a job can be enqueued for the sample (max 2: 1 running + 1 pending)."""
        with self.job_tracking_lock:
            running_count = 0
            pending_count = 0

            if (
                sample_id in self.running_jobs_by_sample
                and job_type in self.running_jobs_by_sample[sample_id]
            ):
                running_count = self.running_jobs_by_sample[sample_id][job_type]

            if (
                sample_id in self.pending_jobs_by_sample
                and job_type in self.pending_jobs_by_sample[sample_id]
            ):
                pending_count = self.pending_jobs_by_sample[sample_id][job_type]

            # Allow if we have less than 2 total jobs (running + pending)
            # AND we don't exceed 1 pending job
            return (running_count + pending_count) < 2 and pending_count < 1

    def _mark_job_pending_for_sample(self, job_type: str, sample_id: str) -> None:
        """Mark a job as pending for the sample."""
        with self.job_tracking_lock:
            if sample_id not in self.pending_jobs_by_sample:
                self.pending_jobs_by_sample[sample_id] = {}
            if job_type not in self.pending_jobs_by_sample[sample_id]:
                self.pending_jobs_by_sample[sample_id][job_type] = 0
            self.pending_jobs_by_sample[sample_id][job_type] += 1

    def _mark_job_running_for_sample(self, job_type: str, sample_id: str) -> None:
        """Mark a job as running for the sample (move from pending to running)."""
        with self.job_tracking_lock:
            # Remove from pending
            if (
                sample_id in self.pending_jobs_by_sample
                and job_type in self.pending_jobs_by_sample[sample_id]
            ):
                self.pending_jobs_by_sample[sample_id][job_type] -= 1
                if self.pending_jobs_by_sample[sample_id][job_type] <= 0:
                    del self.pending_jobs_by_sample[sample_id][job_type]
                if not self.pending_jobs_by_sample[sample_id]:
                    del self.pending_jobs_by_sample[sample_id]

            # Add to running
            if sample_id not in self.running_jobs_by_sample:
                self.running_jobs_by_sample[sample_id] = {}
            if job_type not in self.running_jobs_by_sample[sample_id]:
                self.running_jobs_by_sample[sample_id][job_type] = 0
            self.running_jobs_by_sample[sample_id][job_type] += 1

    def _unmark_job_for_sample(self, job_type: str, sample_id: str) -> None:
        """Unmark a job as running for the sample (job completed)."""
        with self.job_tracking_lock:
            # Remove from running
            if (
                sample_id in self.running_jobs_by_sample
                and job_type in self.running_jobs_by_sample[sample_id]
            ):
                self.running_jobs_by_sample[sample_id][job_type] -= 1
                if self.running_jobs_by_sample[sample_id][job_type] <= 0:
                    del self.running_jobs_by_sample[sample_id][job_type]
                if not self.running_jobs_by_sample[sample_id]:
                    del self.running_jobs_by_sample[sample_id]

    def enqueue_jobs(self, jobs: List[Job]) -> None:
        """Add jobs to their respective queues with bounded deduplication (max 2: 1 running + 1 pending)."""
        jobs_to_enqueue = []

        for job in jobs:
            sample_id = job.get_sample_id()

            # Check if this job type should be deduplicated
            if job.job_type in self.deduplicate_job_types:
                if not self._can_enqueue_job_for_sample(job.job_type, sample_id):
                    # Already have max jobs (1 running + 1 pending) for this sample, skip it
                    logger = get_job_logger(
                        str(job.job_id), job.job_type, job.context.filepath
                    )
                    logger.info(
                        f"Skipping {job.job_type} job for sample {sample_id} (max jobs reached: 1 running + 1 pending)"
                    )
                    self.total_jobs_skipped += 1  # Track skipped jobs
                    continue
                else:
                    # Mark this job as pending for the sample
                    self._mark_job_pending_for_sample(job.job_type, sample_id)
                    logger = get_job_logger(
                        str(job.job_id), job.job_type, job.context.filepath
                    )
                    logger.debug(
                        f"Marked {job.job_type} job as pending for sample {sample_id}"
                    )

            jobs_to_enqueue.append(job)

        # Enqueue the filtered jobs
        for job in jobs_to_enqueue:
            # Track total jobs enqueued
            self.total_jobs_enqueued += 1

            # Determine which queue to use based on job type
            queue_type = self.job_queue_mapping.get(job.job_type, "slow")

            # Log job enqueuing for debugging
            logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
            logger.info(
                f"Enqueuing job {job.job_id} ({job.job_type}) to {queue_type} queue"
            )

            if queue_type == "preprocessing":
                self.preprocessing_queue.put(job)
            elif queue_type == "bed_conversion":
                self.bed_conversion_queue.put(job)
            elif queue_type == "analysis":
                if not self.use_separate_analysis_queues:
                    self.analysis_queue.put(job)
                else:
                    raise ValueError(
                        f"Cannot enqueue to 'analysis' queue in separate analysis queues mode"
                    )
            elif queue_type == "mgmt":
                if self.use_separate_analysis_queues:
                    self.mgmt_queue.put(job)
                else:
                    raise ValueError(
                        f"Cannot enqueue to 'mgmt' queue in legacy analysis queue mode"
                    )
            elif queue_type == "cnv":
                if self.use_separate_analysis_queues:
                    self.cnv_queue.put(job)
                else:
                    raise ValueError(
                        f"Cannot enqueue to 'cnv' queue in legacy analysis queue mode"
                    )
            elif queue_type == "target":
                if self.use_separate_analysis_queues:
                    self.target_queue.put(job)
                else:
                    raise ValueError(
                        f"Cannot enqueue to 'target' queue in legacy analysis queue mode"
                    )
            elif queue_type == "fusion":
                if self.use_separate_analysis_queues:
                    self.fusion_queue.put(job)
                else:
                    raise ValueError(
                        f"Cannot enqueue to 'fusion' queue in legacy analysis queue mode"
                    )
            elif queue_type == "classification":
                self.classification_queue.put(job)
            elif queue_type == "slow":
                self.slow_queue.put(job)
            else:
                raise ValueError(
                    f"Invalid queue type for job {job.job_type}: {queue_type}"
                )

    def worker(
        self, name: str, job_queue: queue.Queue, handlers: Dict[str, Callable]
    ) -> None:
        """Worker that processes jobs from a queue."""
        while self.running:
            try:
                job = job_queue.get(timeout=1.0)

                # Track job start time and active status (this includes queue waiting time)
                self.job_start_times[job.job_id] = time.time()
                self.active_jobs[job.job_id] = {
                    "job_type": job.job_type,
                    "filepath": job.context.filepath,
                    "sample_id": (
                        job.get_sample_id()
                        if hasattr(job, "get_sample_id")
                        else "unknown"
                    ),
                    "worker": name,
                    "start_time": time.time(),
                    "processing_start_time": None,  # Will be set when actual processing starts
                }

                # Get job-specific logger
                logger = get_job_logger(
                    str(job.job_id), job.job_type, job.context.filepath
                )
                logger.info(
                    f"Starting job {job.job_id} ({job.job_type}) for {job.context.filepath}"
                )

                # Mark job as running for deduplicated job types
                if job.job_type in self.deduplicate_job_types:
                    sample_id = job.get_sample_id()
                    self._mark_job_running_for_sample(job.job_type, sample_id)
                    logger.debug(
                        f"Marked {job.job_type} job as running for sample {sample_id}"
                    )

                # Per-sample start tracking (when sample_id is known)
                sid_start = job.get_sample_id()
                if sid_start and sid_start != "unknown":
                    self._on_sample_job_started(sid_start, job.job_type)

                try:
                    # Set processing start time right before calling the handler (excludes queue waiting time)
                    self.active_jobs[job.job_id]["processing_start_time"] = time.time()

                    handlers[job.job_type](job)

                    # Track this job as completed regardless of whether it has a next job
                    logger.info(
                        f"Job {job.job_id} ({job.job_type}) completed successfully"
                    )
                    self.completed_jobs.append(job.job_id)

                    # Track completion by job type
                    with self.progress_lock:
                        if job.job_type not in self.completed_jobs_by_type:
                            self.completed_jobs_by_type[job.job_type] = 0
                        self.completed_jobs_by_type[job.job_type] += 1

                    # Check if this job can trigger other jobs (parallel execution)
                    triggered_jobs = job.can_trigger_jobs()
                    if triggered_jobs:
                        logger.info(
                            f"Job {job.job_type} completed, triggering {len(triggered_jobs)} parallel jobs: {[j.job_type for j in triggered_jobs]}"
                        )
                        self.enqueue_jobs(triggered_jobs)
                    else:
                        # Fallback to linear workflow if no parallel triggers
                        next_job = job.next_job()
                        if next_job:
                            logger.info(f"Creating next job: {next_job.job_type}")
                            self.enqueue_jobs([next_job])
                        else:
                            logger.info(
                                f"No more jobs in workflow for {job.context.filepath}"
                            )
                            logger.debug(f"Final results: {job.context.results}")
                            # Context clearing temporarily disabled to avoid interfering with downstream queues

                except Exception as e:
                    error_msg = f"Job {job.job_type} failed: {str(e)}"
                    job.context.add_error(job.job_type, error_msg)
                    self.failed_jobs.append(job.job_id)

                    # Track failed jobs by type
                    with self.progress_lock:
                        if job.job_type not in self.failed_jobs_by_type:
                            self.failed_jobs_by_type[job.job_type] = 0
                        self.failed_jobs_by_type[job.job_type] += 1

                    logger.error(f"Job {job.job_id} failed: {error_msg}")
                finally:
                    # Per-sample finish tracking (prefer sample_id after handler updated context)
                    sid_finish = job.get_sample_id()
                    if sid_finish and sid_finish != "unknown":
                        self._on_sample_job_finished(
                            sid_finish, job.job_type, job.job_id in self.completed_jobs
                        )

                    # Remove from active jobs
                    if job.job_id in self.active_jobs:
                        del self.active_jobs[job.job_id]
                    if job.job_id in self.job_start_times:
                        del self.job_start_times[job.job_id]

                    # Always unmark the job as running when it's done (success or failure)
                    if job.job_type in self.deduplicate_job_types:
                        sample_id = job.get_sample_id()
                        self._unmark_job_for_sample(job.job_type, sample_id)
                        logger.debug(
                            f"Unmarked {job.job_type} job as running for sample {sample_id}"
                        )

                job_queue.task_done()

            except queue.Empty:
                continue

    def run(self) -> None:
        """Start the workflow manager with specialized workers for different task categories."""
        # Create preprocessing workers (bam_preprocessing only)
        for i in range(self.preprocessing_workers_count):
            preprocessing_worker = threading.Thread(
                target=self.worker,
                args=(
                    f"PreprocessingWorker-{i+1}",
                    self.preprocessing_queue,
                    self.job_handlers_preprocessing,
                ),
            )
            preprocessing_worker.daemon = True
            self.preprocessing_workers.append(preprocessing_worker)
            preprocessing_worker.start()

        # Create bed conversion workers (bed_conversion only)
        for i in range(self.bed_conversion_workers_count):
            bed_conversion_worker = threading.Thread(
                target=self.worker,
                args=(
                    f"BedConversionWorker-{i+1}",
                    self.bed_conversion_queue,
                    self.job_handlers_bed_conversion,
                ),
            )
            bed_conversion_worker.daemon = True
            self.bed_conversion_workers.append(bed_conversion_worker)
            bed_conversion_worker.start()

        # Create analysis workers based on mode
        if self.use_separate_analysis_queues:
            # Create separate workers for each analysis job type (ensures only one of each type runs at a time)
            for i in range(self.analysis_workers_count):
                # MGMT worker
                mgmt_worker = threading.Thread(
                    target=self.worker,
                    args=(f"MGMTWorker-{i+1}", self.mgmt_queue, self.job_handlers_mgmt),
                )
                mgmt_worker.daemon = True
                self.mgmt_workers.append(mgmt_worker)
                mgmt_worker.start()

                # CNV worker
                cnv_worker = threading.Thread(
                    target=self.worker,
                    args=(f"CNVWorker-{i+1}", self.cnv_queue, self.job_handlers_cnv),
                )
                cnv_worker.daemon = True
                self.cnv_workers.append(cnv_worker)
                cnv_worker.start()

                # Target worker
                target_worker = threading.Thread(
                    target=self.worker,
                    args=(
                        f"TargetWorker-{i+1}",
                        self.target_queue,
                        self.job_handlers_target,
                    ),
                )
                target_worker.daemon = True
                self.target_workers.append(target_worker)
                target_worker.start()

                # Fusion worker
                fusion_worker = threading.Thread(
                    target=self.worker,
                    args=(
                        f"FusionWorker-{i+1}",
                        self.fusion_queue,
                        self.job_handlers_fusion,
                    ),
                )
                fusion_worker.daemon = True
                self.fusion_workers.append(fusion_worker)
                fusion_worker.start()
        else:
            # Create legacy single analysis workers
            for i in range(self.analysis_workers_count):
                analysis_worker = threading.Thread(
                    target=self.worker,
                    args=(
                        f"AnalysisWorker-{i+1}",
                        self.analysis_queue,
                        self.job_handlers_analysis,
                    ),
                )
                analysis_worker.daemon = True
                self.analysis_workers.append(analysis_worker)
                analysis_worker.start()

        # Create classification workers (sturgeon, nanodx, pannanodx)
        for i in range(self.classification_workers_count):
            classification_worker = threading.Thread(
                target=self.worker,
                args=(
                    f"ClassificationWorker-{i+1}",
                    self.classification_queue,
                    self.job_handlers_classification,
                ),
            )
            classification_worker.daemon = True
            self.classification_workers.append(classification_worker)
            classification_worker.start()

        # Create slow workers (legacy)
        for i in range(self.slow_workers_count):
            slow_worker = threading.Thread(
                target=self.worker,
                args=(f"SlowWorker-{i+1}", self.slow_queue, self.job_handlers_slow),
            )
            slow_worker.daemon = True
            self.slow_workers.append(slow_worker)
            slow_worker.start()

    def _try_clear_context(self, context: "WorkflowContext") -> None:
        """Clear a context only when all jobs in the original workflow have produced results.
        This prevents premature clearing that would block downstream triggers.
        """
        original_workflow = context.metadata.get("original_workflow", [])
        required_job_types: Set[str] = set()
        for step in original_workflow:
            if isinstance(step, str) and ":" in step:
                _, job_type = step.split(":", 1)
                required_job_types.add(job_type)
        if not required_job_types:
            return
        completed_job_types = set(context.history)
        # Only clear when every planned job type has at least one result entry
        if required_job_types.issubset(completed_job_types):
            context.clear_all()

        # Wait for all workers to finish with periodic checks for shutdown
        all_workers = (
            self.preprocessing_workers
            + self.bed_conversion_workers
            + self.mgmt_workers
            + self.cnv_workers
            + self.target_workers
            + self.fusion_workers
            + self.classification_workers
            + self.slow_workers
        )

        try:
            while self.running:
                # Check if any workers are still alive
                alive_workers = [w for w in all_workers if w.is_alive()]
                if not alive_workers:
                    break

                # Sleep briefly to allow for interrupt handling
                time.sleep(0.1)

        except KeyboardInterrupt:
            # Signal workers to stop
            self.running = False
            raise

    def stop(self, timeout: float = 30.0) -> bool:
        """
        Stop the workflow manager and clean up running threads.

        Args:
            timeout: Maximum time to wait for workers to finish (seconds)

        Returns:
            True if all workers stopped gracefully, False if timeout occurred
        """
        if not self.running:
            return True

        if self.verbose:
            print(f"[WorkflowManager] Stopping workflow manager (timeout: {timeout}s)")

        # Signal workers to stop
        self.running = False

        # Wait for workers to finish with timeout
        if self.use_separate_analysis_queues:
            all_workers = (
                self.preprocessing_workers
                + self.bed_conversion_workers
                + self.mgmt_workers
                + self.cnv_workers
                + self.target_workers
                + self.fusion_workers
                + self.classification_workers
                + self.slow_workers
            )
        else:
            all_workers = (
                self.preprocessing_workers
                + self.bed_conversion_workers
                + self.analysis_workers
                + self.classification_workers
                + self.slow_workers
            )

        # Wait for each worker with timeout
        for worker in all_workers:
            if worker.is_alive():
                worker.join(timeout=timeout)
                if worker.is_alive():
                    if self.verbose:
                        print(
                            f"[WorkflowManager] Warning: Worker {worker.name} did not stop within timeout"
                        )
                    return False

        if self.verbose:
            print("[WorkflowManager] All workers stopped successfully")

        return True

    def shutdown(self, timeout: float = 30.0) -> bool:
        """
        Alias for stop() method for backward compatibility.

        Args:
            timeout: Maximum time to wait for workers to finish (seconds)

        Returns:
            True if all workers stopped gracefully, False if timeout occurred
        """
        return self.stop(timeout)

    def is_running(self) -> bool:
        """Check if the workflow manager is currently running."""
        return self.running

    def get_stats(self) -> dict:
        """Get statistics about the workflow execution."""
        # Get queue sizes
        preprocessing_queue_size = self.preprocessing_queue.qsize()
        bed_conversion_queue_size = self.bed_conversion_queue.qsize()

        if self.use_separate_analysis_queues:
            mgmt_queue_size = self.mgmt_queue.qsize()
            cnv_queue_size = self.cnv_queue.qsize()
            target_queue_size = self.target_queue.qsize()
            fusion_queue_size = self.fusion_queue.qsize()
            analysis_queue_size = 0  # Not used in separate mode
        else:
            analysis_queue_size = self.analysis_queue.qsize()
            mgmt_queue_size = 0
            cnv_queue_size = 0
            target_queue_size = 0
            fusion_queue_size = 0

        classification_queue_size = self.classification_queue.qsize()
        slow_queue_size = self.slow_queue.qsize()

        # Get active jobs by worker type
        active_by_worker = {}
        for job_info in self.active_jobs.values():
            worker = job_info["worker"]
            if worker not in active_by_worker:
                active_by_worker[worker] = []

            # Use processing_start_time if available (excludes queue waiting time), otherwise fall back to start_time
            if job_info.get("processing_start_time") is not None:
                duration = time.time() - job_info["processing_start_time"]
            else:
                duration = time.time() - job_info["start_time"]

            active_by_worker[worker].append(
                {
                    "job_type": job_info["job_type"],
                    "filepath": job_info["filepath"],
                    "sample_id": job_info.get("sample_id", "unknown"),
                    "duration": duration,
                }
            )

        # Calculate total expected jobs based on completed + failed + active + queued
        # This gives us the current total of all jobs that have been processed or are in progress
        if self.use_separate_analysis_queues:
            total_expected = (
                len(self.completed_jobs)
                + len(self.failed_jobs)
                + len(self.active_jobs)
                + preprocessing_queue_size
                + bed_conversion_queue_size
                + mgmt_queue_size
                + cnv_queue_size
                + target_queue_size
                + fusion_queue_size
                + classification_queue_size
                + slow_queue_size
            )
        else:
            total_expected = (
                len(self.completed_jobs)
                + len(self.failed_jobs)
                + len(self.active_jobs)
                + preprocessing_queue_size
                + bed_conversion_queue_size
                + analysis_queue_size
                + classification_queue_size
                + slow_queue_size
            )

        # If we have no jobs in progress but have completed jobs, use the completed count as the total
        if total_expected == 0 and len(self.completed_jobs) > 0:
            total_expected = len(self.completed_jobs)

        # Ensure total_expected is at least as large as the number of completed jobs
        if total_expected < len(self.completed_jobs):
            total_expected = len(self.completed_jobs)

        # Calculate total jobs that were actually enqueued (excluding skipped jobs)
        total_actual_jobs = self.total_jobs_enqueued

        # Prepare samples payload (copy to avoid race conditions)
        samples_payload: List[Dict[str, Any]] = []
        with self.progress_lock:
            for sid, info in self.samples_by_id.items():
                samples_payload.append(
                    {
                        "sample_id": info["sample_id"],
                        "active_jobs": info["active_jobs"],
                        "total_jobs": info["total_jobs"],
                        "completed_jobs": info["completed_jobs"],
                        "failed_jobs": info["failed_jobs"],
                        "job_types": list(info["job_types"]),
                        "last_seen": info["last_seen"],
                    }
                )

        return {
            "completed": len(self.completed_jobs),
            "failed": len(self.failed_jobs),
            "total_enqueued": self.total_jobs_enqueued,
            "total_skipped": self.total_jobs_skipped,
            "total_actual_jobs": total_actual_jobs,
            "total_processed": len(self.completed_jobs) + len(self.failed_jobs),
            "total_expected": total_expected,
            "active_jobs": len(self.active_jobs),
            "completed_by_type": self.completed_jobs_by_type.copy(),  # Make a copy to avoid threading issues
            "failed_by_type": self.failed_jobs_by_type.copy(),  # Make a copy to avoid threading issues
            "samples": samples_payload,
            "queue_sizes": {
                "preprocessing": preprocessing_queue_size,
                "bed_conversion": bed_conversion_queue_size,
                "mgmt": mgmt_queue_size,
                "cnv": cnv_queue_size,
                "target": target_queue_size,
                "fusion": fusion_queue_size,
                "analysis": analysis_queue_size,
                "classification": classification_queue_size,
                "slow": slow_queue_size,
            },
            "active_by_worker": active_by_worker,
        }

    # === Per-sample tracking helpers ===
    def _ensure_sample_entry(self, sample_id: str) -> Dict[str, Any]:
        with self.progress_lock:
            if sample_id not in self.samples_by_id:
                self.samples_by_id[sample_id] = {
                    "sample_id": sample_id,
                    "active_jobs": 0,
                    "total_jobs": 0,
                    "completed_jobs": 0,
                    "failed_jobs": 0,
                    "job_types": set(),
                    "last_seen": time.time(),
                }
            return self.samples_by_id[sample_id]

    def _on_sample_job_started(self, sample_id: str, job_type: str) -> None:
        entry = self._ensure_sample_entry(sample_id)
        with self.progress_lock:
            entry["active_jobs"] += 1
            entry["total_jobs"] += 1
            if job_type:
                entry["job_types"].add(job_type)
            entry["last_seen"] = time.time()

    def _on_sample_job_finished(
        self, sample_id: str, job_type: str, success: bool
    ) -> None:
        entry = self._ensure_sample_entry(sample_id)
        with self.progress_lock:
            if entry["active_jobs"] > 0:
                entry["active_jobs"] -= 1
            if success:
                entry["completed_jobs"] += 1
            else:
                entry["failed_jobs"] += 1
            if job_type:
                entry["job_types"].add(job_type)
            entry["last_seen"] = time.time()


# === File Watcher ===
class FileWatcher(FileSystemEventHandler):
    """File watcher that integrates with the workflow manager."""

    def __init__(
        self,
        watch_dir: str,
        preprocessor_func: Callable[[str], List[Job]],
        manager: WorkflowManager,
        recursive: bool = True,
        patterns: Optional[List[str]] = None,
        ignore_patterns: Optional[List[str]] = None,
        verbose: bool = False,
        show_progress: bool = True,
    ):
        self.watch_dir = watch_dir
        self.preprocessor_func = preprocessor_func
        self.manager = manager
        self.recursive = recursive
        self.patterns = patterns or ["*"]
        self.ignore_patterns = ignore_patterns or []
        self.verbose = verbose
        self.show_progress = show_progress
        self.observer = Observer()
        self.processed_files = set()

    def _should_process_file(self, filepath: str) -> bool:
        """Check if a file should be processed based on patterns."""
        path = Path(filepath)

        # Check if file matches any ignore patterns
        for pattern in self.ignore_patterns:
            if path.match(pattern):
                return False

        # Check if file matches any watch patterns
        for pattern in self.patterns:
            if path.match(pattern):
                return True

        return False

    def handle_file(self, filepath: str) -> None:
        """Process a new file through the workflow."""
        if filepath in self.processed_files:
            return

        self.processed_files.add(filepath)

        # Use logger for file watcher events
        import logging

        logger = logging.getLogger("littlejohn.filewatcher")
        logger.info(f"Detected new file: {filepath}")

        try:
            jobs = self.preprocessor_func(filepath)
            self.manager.enqueue_jobs(jobs)
            logger.debug(f"Queued {len(jobs)} jobs for {filepath}")
        except Exception as e:
            logger.error(f"Error processing {filepath}: {e}")

    def on_created(self, event) -> None:
        """Handle file creation events."""
        if not event.is_directory and self._should_process_file(event.src_path):
            self.handle_file(event.src_path)

    def on_modified(self, event) -> None:
        """Handle file modification events."""
        if not event.is_directory and self._should_process_file(event.src_path):
            self.handle_file(event.src_path)

    def start(self, process_existing: bool = True) -> None:
        """Start watching the directory."""
        if process_existing:
            self._process_existing_files()

        self.observer.schedule(self, self.watch_dir, recursive=self.recursive)
        self.observer.start()

        # Use logger for file watcher events
        import logging

        logger = logging.getLogger("littlejohn.filewatcher")
        logger.info(f"Watching directory: {self.watch_dir}")
        logger.debug(f"Patterns: {self.patterns}")
        logger.debug(f"Ignore patterns: {self.ignore_patterns}")

    def _process_existing_files(self) -> None:
        """Process existing files that match the patterns."""
        import logging

        logger = logging.getLogger("littlejohn.filewatcher")
        logger.info(f"Processing existing files in: {self.watch_dir}")

        existing_files = []

        # Find all existing files that match patterns
        for pattern in self.patterns:
            if pattern == "*":
                # Default pattern - all files
                if self.recursive:
                    pattern_files = list(Path(self.watch_dir).rglob("*"))
                else:
                    pattern_files = list(Path(self.watch_dir).glob("*"))
            else:
                # Simple pattern like *.bam
                if self.recursive:
                    pattern_files = list(Path(self.watch_dir).rglob(pattern))
                else:
                    pattern_files = list(Path(self.watch_dir).glob(pattern))

            # Filter to only files (not directories) and apply ignore patterns
            for file_path in pattern_files:
                if file_path.is_file():
                    # Check if file should be ignored
                    should_ignore = False
                    if self.ignore_patterns:
                        for ignore_pattern in self.ignore_patterns:
                            if file_path.match(ignore_pattern):
                                should_ignore = True
                                break

                    if not should_ignore:
                        existing_files.append(file_path)

        # Remove duplicates and sort
        existing_files = sorted(set(existing_files))

        if existing_files:
            logger.info(f"Found {len(existing_files)} existing file(s) to process:")
            for file_path in existing_files:
                logger.debug(f"  - {file_path}")

            # Process each existing file with progress bar
            with tqdm(
                total=len(existing_files),
                desc="Processing existing files",
                unit="file",
                disable=not self.show_progress,
            ) as pbar:
                for file_path in existing_files:
                    if self.show_progress:
                        pbar.set_postfix_str(f"Processing: {file_path.name}")
                    logger.debug(f"Processing existing file: {file_path}")
                    self.handle_file(str(file_path))
                    pbar.update(1)
        else:
            logger.info("No existing files found matching the patterns.")

    def stop(self, timeout: float = 30.0) -> bool:
        """Stop watching and shutdown the workflow manager."""
        if self.verbose:
            print("[Watcher] Stopping file watcher and workflow manager")

        # Stop the file observer
        self.observer.stop()
        self.observer.join()

        # Stop the workflow manager
        manager_stopped = self.manager.stop(timeout=timeout)

        if self.verbose:
            print("[Watcher] Stopped")

        return manager_stopped


# === Job Classifier ===
_job_id_counter = itertools.count(1000)


def default_file_classifier(filepath: str, workflow_plan: List[str]) -> List[Job]:
    """Default classifier that creates jobs for a file based on a workflow plan."""
    job_id = next(_job_id_counter)
    ctx = WorkflowContext(filepath)
    ctx.add_metadata("filename", os.path.basename(filepath))
    ctx.add_metadata("created", time.time())

    # Try to get file size, but don't fail if file doesn't exist
    try:
        ctx.add_metadata("file_size", os.path.getsize(filepath))
    except (OSError, FileNotFoundError):
        ctx.add_metadata("file_size", 0)

    # Define job dependencies and triggers for parallel execution
    job_dependencies = {
        "bed_conversion": {"preprocessing"},
        "mgmt": {"preprocessing"},
        "cnv": {"preprocessing"},
        "target": {"preprocessing"},
        "fusion": {"preprocessing"},
        "sturgeon": {"bed_conversion"},
        "nanodx": {"bed_conversion"},
        "pannanodx": {"bed_conversion"},
        "random_forest": {"bed_conversion"},
    }

    # Create preprocessing job if it's the first step
    if workflow_plan and workflow_plan[0].endswith(":preprocessing"):
        # Determine what jobs this preprocessing can trigger
        triggers = set()
        for job_type, deps in job_dependencies.items():
            if "preprocessing" in deps:
                triggers.add(job_type)

        # Store the original workflow plan in context for reference
        ctx.add_metadata("original_workflow", workflow_plan)

        return [
            Job(
                job_id=job_id,
                job_type="preprocessing",
                context=ctx,
                origin="preprocessing",
                workflow=workflow_plan,
                step=0,
                dependencies=set(),
                triggers=triggers,
            )
        ]

    # Create first job in plan (skip preprocessing if not present)
    if len(workflow_plan) > 0:
        first_step = workflow_plan[0]
        if not isinstance(first_step, str):
            raise ValueError(
                f"Expected string in workflow step, got {type(first_step)}: {first_step}"
            )

        if ":" not in first_step:
            raise ValueError(
                f"Invalid workflow step format (missing ':'): {first_step}"
            )

        queue_type, job_type = first_step.split(":", 1)
        dependencies = job_dependencies.get(job_type, set())
        triggers = set()
        for trigger_job, trigger_deps in job_dependencies.items():
            if job_type in trigger_deps:
                triggers.add(trigger_job)

        # Store the original workflow plan in context for reference
        ctx.add_metadata("original_workflow", workflow_plan)

        return [
            Job(
                job_id=job_id,
                job_type=job_type,
                context=ctx,
                origin=queue_type,
                workflow=workflow_plan,
                step=0,
                dependencies=dependencies,
                triggers=triggers,
            )
        ]

    return []


# === Built-in Job Handlers ===
def echo_handler(job: Job) -> None:
    """Simple echo handler for testing."""
    time.sleep(0.1)
    ctx = job.context
    print(f"  [Echo] Processing {ctx.filepath}")
    ctx.add_metadata("echo_processed", True)
    ctx.add_result("echo", "echo_result")


def sleep_handler(job: Job) -> None:
    """Sleep handler for testing slow operations."""
    time.sleep(1.0)
    ctx = job.context
    print(f"  [Sleep] Processing {ctx.filepath}")
    ctx.add_metadata("sleep_processed", True)
    ctx.add_result("sleep", "sleep_result")


def command_handler(job: Job, command_template: str) -> None:
    """Handler that runs shell commands."""
    import subprocess

    ctx = job.context
    command = command_template.replace("{file}", ctx.filepath)

    try:
        result = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
            timeout=30,
        )

        ctx.add_result(
            job.job_type,
            {
                "command": command,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "returncode": result.returncode,
            },
        )

        if result.returncode != 0:
            ctx.add_error(
                job.job_type, f"Command failed with return code {result.returncode}"
            )

    except subprocess.TimeoutExpired:
        ctx.add_error(job.job_type, "Command timed out")
    except Exception as e:
        ctx.add_error(job.job_type, f"Command error: {str(e)}")


# === Workflow Runner ===
class WorkflowRunner:
    """High-level interface for running workflows."""

    def __init__(
        self,
        verbose: bool = False,
        analysis_workers: int = 1,
        use_separate_analysis_queues: bool = True,
        preprocessing_workers: int = 1,
        bed_workers: int = 1,
    ):
        self.manager = WorkflowManager(
            verbose=verbose,
            analysis_workers=analysis_workers,
            use_separate_analysis_queues=use_separate_analysis_queues,
            preprocessing_workers=preprocessing_workers,
            bed_conversion_workers=bed_workers,
        )
        self.verbose = verbose

        # Register default handlers
        self.manager.register_handler("preprocessing", "echo", echo_handler)
        self.manager.register_handler("slow", "sleep", sleep_handler)

    def register_handler(
        self, queue_type: str, job_type: str, handler: Callable[[Job], None]
    ) -> None:
        """Register a custom job handler."""
        self.manager.register_handler(queue_type, job_type, handler)

    def register_command_handler(
        self, queue_type: str, job_type: str, command_template: str
    ) -> None:
        """Register a command handler that runs shell commands."""

        def handler(job: Job) -> None:
            command_handler(job, command_template)

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
        show_progress: bool = True,
    ) -> None:
        """Run a complete workflow."""

        if classifier_func is None:
            classifier_func = lambda filepath: default_file_classifier(
                filepath, workflow_plan
            )

        watcher = FileWatcher(
            watch_dir=watch_dir,
            preprocessor_func=classifier_func,
            manager=self.manager,
            recursive=recursive,
            patterns=patterns,
            ignore_patterns=ignore_patterns,
            verbose=self.verbose,
            show_progress=show_progress,
        )

        # Start worker threads first to ensure queues are consumed immediately
        self.manager.run()

        # Now start watching files (which enqueues jobs)
        watcher.start(process_existing=process_existing)

        # Start progress monitoring if enabled
        progress_thread = None
        if show_progress:
            progress_thread = threading.Thread(
                target=self._monitor_progress, args=(watcher,), daemon=True
            )
            progress_thread.start()

        try:
            # Keep the main thread alive while workers and watcher run
            while self.manager.is_running():
                time.sleep(0.5)
        except KeyboardInterrupt:
            if self.verbose:
                print("[WorkflowRunner] Shutdown requested.")
            print("[WorkflowRunner] Shutdown requested.")
            # Stop the watcher and workflow manager
            graceful_shutdown = watcher.stop(timeout=1.0)

            if self.verbose:
                if graceful_shutdown:
                    print("[WorkflowRunner] Graceful shutdown completed.")
                else:
                    print(
                        "[WorkflowRunner] Warning: Some workers may not have stopped gracefully."
                    )
                stats = self.manager.get_stats()
                print(f"[WorkflowRunner] Final stats: {stats}")
                print(
                    f"[WorkflowRunner] Final completion by type: {stats.get('completed_by_type', {})}"
                )
                print(
                    f"[WorkflowRunner] Final failed by type: {stats.get('failed_by_type', {})}"
                )

    def _monitor_progress(self, watcher) -> None:
        """Monitor and display worker progress in real-time."""
        import time

        # Create progress bars for each worker type
        preprocessing_pbar = tqdm(
            desc="Preprocessing", unit="jobs", position=0, leave=True
        )
        bed_conversion_pbar = tqdm(
            desc="Bed Conversion", unit="jobs", position=1, leave=True
        )
        mgmt_pbar = tqdm(desc="MGMT", unit="jobs", position=2, leave=True)
        cnv_pbar = tqdm(desc="CNV", unit="jobs", position=3, leave=True)
        target_pbar = tqdm(desc="Target", unit="jobs", position=4, leave=True)
        fusion_pbar = tqdm(desc="Fusion", unit="jobs", position=5, leave=True)
        classification_pbar = tqdm(
            desc="Classification", unit="jobs", position=6, leave=True
        )
        slow_pbar = tqdm(desc="Slow", unit="jobs", position=7, leave=True)

        # Create overall progress bar
        overall_pbar = tqdm(
            desc="Overall Progress", unit="jobs", position=8, leave=True
        )

        progress_bars = {
            "preprocessing": preprocessing_pbar,
            "bed_conversion": bed_conversion_pbar,
            "mgmt": mgmt_pbar,
            "cnv": cnv_pbar,
            "target": target_pbar,
            "fusion": fusion_pbar,
            "classification": classification_pbar,
            "slow": slow_pbar,
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
            "slow": 0,
        }

        last_stats = None

        try:
            while self.manager.is_running():
                stats = self.manager.get_stats()

                # Update overall progress using total actual jobs (excluding skipped)
                total_processed = stats["total_processed"]
                total_actual_jobs = stats["total_actual_jobs"]

                if total_actual_jobs > 0:
                    overall_pbar.total = total_actual_jobs
                    overall_pbar.n = total_processed
                    overall_pbar.set_postfix_str(
                        f"Active: {stats['active_jobs']} | "
                        f"Completed: {stats['completed']} | "
                        f"Failed: {stats['failed']} | "
                        f"Skipped: {stats['total_skipped']} | "
                        f"Total: {total_actual_jobs}"
                    )

                # Update queue-specific progress bars
                for queue_name, pbar in progress_bars.items():
                    queue_size = stats["queue_sizes"][queue_name]

                    # Count active jobs for this queue type
                    active_in_queue = 0
                    active_jobs_info = []

                    for worker_name, jobs in stats["active_by_worker"].items():
                        for job in jobs:
                            # Check if this worker belongs to the current queue
                            if (
                                queue_name == "preprocessing"
                                and worker_name.startswith("PreprocessingWorker")
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif (
                                queue_name == "bed_conversion"
                                and worker_name.startswith("BedConversionWorker")
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif queue_name == "mgmt" and worker_name.startswith(
                                "MGMTWorker"
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif queue_name == "cnv" and worker_name.startswith(
                                "CNVWorker"
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif queue_name == "target" and worker_name.startswith(
                                "TargetWorker"
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif queue_name == "fusion" and worker_name.startswith(
                                "FusionWorker"
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif (
                                queue_name == "classification"
                                and worker_name.startswith("ClassificationWorker")
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )
                            elif queue_name == "slow" and worker_name.startswith(
                                "SlowWorker"
                            ):
                                active_in_queue += 1
                                filename = job["filepath"].split("/")[-1]
                                duration = int(job["duration"])
                                active_jobs_info.append(
                                    f"{job['job_type']}:{filename}({duration}s)"
                                )

                    # Calculate completed jobs for this queue type (both successful and failed)
                    completed_in_queue = 0
                    for job_type, count in stats["completed_by_type"].items():
                        # Get queue type from the manager's job queue mapping
                        queue_type = self.manager.job_queue_mapping.get(
                            job_type, "slow"
                        )
                        if queue_type == queue_name:
                            completed_in_queue += count

                    # Add failed jobs to the completion count
                    for job_type, count in stats.get("failed_by_type", {}).items():
                        # Get queue type from the manager's job queue mapping
                        queue_type = self.manager.job_queue_mapping.get(
                            job_type, "slow"
                        )
                        if queue_type == queue_name:
                            completed_in_queue += count

                    # Calculate total jobs for this queue (completed + active + queued)
                    total_for_queue = completed_in_queue + active_in_queue + queue_size

                    # Update the progress bar with total and current values
                    if total_for_queue > 0:
                        pbar.total = total_for_queue
                        pbar.n = completed_in_queue

                    # Update the progress bar description with queue info
                    pbar.set_description(
                        f"{queue_name.title()} (Q:{queue_size} A:{active_in_queue} C:{completed_in_queue})"
                    )

                    # Update active jobs info
                    if active_jobs_info:
                        pbar.set_postfix_str(
                            " | ".join(active_jobs_info[:2])
                        )  # Show first 2 active jobs

                # Check if we should stop monitoring
                if not self.manager.is_running():
                    break

                time.sleep(1.0)  # Update every second

        except KeyboardInterrupt:
            pass
        finally:
            # Close all progress bars
            for pbar in progress_bars.values():
                pbar.close()
            overall_pbar.close()
