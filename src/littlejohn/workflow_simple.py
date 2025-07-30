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
    """Represents a single job in the workflow."""
    job_id: int
    job_type: str
    context: WorkflowContext
    origin: str  # "fast" or "slow"
    workflow: List[str]
    step: int = 0

    def next_job(self) -> Optional['Job']:
        """Create the next job in the workflow if available."""
        if self.step + 1 < len(self.workflow):
            queue_type, next_type = self.workflow[self.step + 1].split(":", 1)
            return Job(
                job_id=self.job_id,
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


# === Workflow Manager ===
class WorkflowManager:
    """Manages job queues and execution using threading with specialized workers."""

    def __init__(self, verbose: bool = False):
        # Specialized queues for different task categories
        self.preprocessing_queue = queue.Queue()  # bam_preprocessing, bed_conversion
        self.analysis_queue = queue.Queue()       # mgmt, cnv
        self.classification_queue = queue.Queue() # sturgeon, nanodx, pannanodx
        self.slow_queue = queue.Queue()           # slow jobs (legacy)
        
        self.completed_jobs = []
        self.failed_jobs = []
        
        # Specialized job handlers for each queue
        self.job_handlers_preprocessing: Dict[str, Callable[[Job], None]] = {}
        self.job_handlers_analysis: Dict[str, Callable[[Job], None]] = {}
        self.job_handlers_classification: Dict[str, Callable[[Job], None]] = {}
        self.job_handlers_slow: Dict[str, Callable[[Job], None]] = {}
        
        self.running = True
        self.verbose = verbose
        
        # Job deduplication tracking
        # Jobs that should be deduplicated by sample ID (e.g., sturgeon analysis)
        self.deduplicate_job_types: Set[str] = {"sturgeon", "nanodx", "pannanodx", "random_forest"}
        # Track running and pending jobs by sample ID for deduplication
        # Allow max 2 jobs per sample: 1 running + 1 pending
        self.running_jobs_by_sample: Dict[str, Dict[str, int]] = {}  # sample_id -> {job_type: count}
        self.pending_jobs_by_sample: Dict[str, Dict[str, int]] = {}  # sample_id -> {job_type: count}
        self.job_tracking_lock = threading.Lock()
        
        # Track worker threads for graceful shutdown
        self.preprocessing_workers = []
        self.analysis_workers = []
        self.classification_workers = []
        self.slow_workers = []
        
        # Define job type mappings to queues
        self.job_queue_mapping = {
            # Preprocessing queue
            "preprocessing": "preprocessing",
            "bed_conversion": "preprocessing",
            
            # Analysis queue  
            "mgmt": "analysis",
            "cnv": "analysis",
            "target": "analysis",
            "fusion": "analysis",
            
            # Classification queue
            "sturgeon": "classification", 
            "nanodx": "classification",
            "pannanodx": "classification",
            
            # Slow queue (legacy)
            "sleep": "slow",
            "echo": "slow"
        }

    def register_handler(self, queue_type: str, job_type: str, handler: Callable[[Job], None]) -> None:
        """Register a handler for a specific job type and queue."""
        if queue_type == "preprocessing":
            self.job_handlers_preprocessing[job_type] = handler
        elif queue_type == "analysis":
            self.job_handlers_analysis[job_type] = handler
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
            
            if sample_id in self.running_jobs_by_sample and job_type in self.running_jobs_by_sample[sample_id]:
                running_count = self.running_jobs_by_sample[sample_id][job_type]
            
            if sample_id in self.pending_jobs_by_sample and job_type in self.pending_jobs_by_sample[sample_id]:
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
            if sample_id in self.pending_jobs_by_sample and job_type in self.pending_jobs_by_sample[sample_id]:
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
            if sample_id in self.running_jobs_by_sample and job_type in self.running_jobs_by_sample[sample_id]:
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
                    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
                    logger.info(f"Skipping {job.job_type} job for sample {sample_id} (max jobs reached: 1 running + 1 pending)")
                    continue
                else:
                    # Mark this job as pending for the sample
                    self._mark_job_pending_for_sample(job.job_type, sample_id)
                    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
                    logger.debug(f"Marked {job.job_type} job as pending for sample {sample_id}")
            
            jobs_to_enqueue.append(job)
        
        # Enqueue the filtered jobs
        for job in jobs_to_enqueue:
            # Determine which queue to use based on job type
            queue_type = self.job_queue_mapping.get(job.job_type, "slow")
            
            if queue_type == "preprocessing":
                self.preprocessing_queue.put(job)
            elif queue_type == "analysis":
                self.analysis_queue.put(job)
            elif queue_type == "classification":
                self.classification_queue.put(job)
            elif queue_type == "slow":
                self.slow_queue.put(job)
            else:
                raise ValueError(f"Invalid queue type for job {job.job_type}: {queue_type}")

    def worker(self, name: str, job_queue: queue.Queue, handlers: Dict[str, Callable]) -> None:
        """Worker that processes jobs from a queue."""
        while self.running:
            try:
                job = job_queue.get(timeout=1.0)
                
                # Get job-specific logger
                logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
                logger.info(f"Starting job {job.job_id} ({job.job_type}) for {job.context.filepath}")
                
                # Mark job as running for deduplicated job types
                if job.job_type in self.deduplicate_job_types:
                    sample_id = job.get_sample_id()
                    self._mark_job_running_for_sample(job.job_type, sample_id)
                    logger.debug(f"Marked {job.job_type} job as running for sample {sample_id}")
                
                try:
                    handlers[job.job_type](job)
                    
                    next_job = job.next_job()
                    if next_job:
                        self.enqueue_jobs([next_job])
                    else:
                        logger.info(f"Job {job.job_id} finished successfully")
                        logger.debug(f"Results: {job.context.results}")
                        self.completed_jobs.append(job.job_id)
                        
                except Exception as e:
                    error_msg = f"Job {job.job_type} failed: {str(e)}"
                    job.context.add_error(job.job_type, error_msg)
                    self.failed_jobs.append(job.job_id)
                    
                    logger.error(f"Job {job.job_id} failed: {error_msg}")
                finally:
                    # Always unmark the job as running when it's done (success or failure)
                    if job.job_type in self.deduplicate_job_types:
                        sample_id = job.get_sample_id()
                        self._unmark_job_for_sample(job.job_type, sample_id)
                        logger.debug(f"Unmarked {job.job_type} job as running for sample {sample_id}")

                job_queue.task_done()
                
            except queue.Empty:
                continue

    def run(self) -> None:
        """Start the workflow manager with specialized workers for different task categories."""
        # Create preprocessing worker (bam_preprocessing, bed_conversion)
        preprocessing_worker = threading.Thread(
            target=self.worker, 
            args=("PreprocessingWorker", self.preprocessing_queue, self.job_handlers_preprocessing)
        )
        preprocessing_worker.daemon = True
        self.preprocessing_workers.append(preprocessing_worker)
        preprocessing_worker.start()
        
        # Create analysis worker (mgmt, cnv)
        analysis_worker = threading.Thread(
            target=self.worker, 
            args=("AnalysisWorker", self.analysis_queue, self.job_handlers_analysis)
        )
        analysis_worker.daemon = True
        self.analysis_workers.append(analysis_worker)
        analysis_worker.start()
        
        # Create classification worker (sturgeon, nanodx, pannanodx)
        classification_worker = threading.Thread(
            target=self.worker, 
            args=("ClassificationWorker", self.classification_queue, self.job_handlers_classification)
        )
        classification_worker.daemon = True
        self.classification_workers.append(classification_worker)
        classification_worker.start()
        
        # Create slow worker (legacy)
        slow_worker = threading.Thread(
            target=self.worker, 
            args=("SlowWorker", self.slow_queue, self.job_handlers_slow)
        )
        slow_worker.daemon = True
        self.slow_workers.append(slow_worker)
        slow_worker.start()
        
        # Wait for all workers to finish
        all_workers = (self.preprocessing_workers + self.analysis_workers + 
                      self.classification_workers + self.slow_workers)
        for worker in all_workers:
            worker.join()

    def get_stats(self) -> dict:
        """Get statistics about the workflow execution."""
        return {
            "completed": len(self.completed_jobs),
            "failed": len(self.failed_jobs),
            "total": len(self.completed_jobs) + len(self.failed_jobs)
        }


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
        verbose: bool = False
    ):
        self.watch_dir = watch_dir
        self.preprocessor_func = preprocessor_func
        self.manager = manager
        self.recursive = recursive
        self.patterns = patterns or ["*"]
        self.ignore_patterns = ignore_patterns or []
        self.verbose = verbose
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
            
            # Process each existing file
            for file_path in existing_files:
                logger.debug(f"Processing existing file: {file_path}")
                self.handle_file(str(file_path))
        else:
            logger.info("No existing files found matching the patterns.")

    def stop(self) -> None:
        """Stop watching and shutdown the workflow manager."""
        self.manager.running = False
        self.observer.stop()
        self.observer.join()
        
        if self.verbose:
            print("[Watcher] Stopped")


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

    # Create preprocessing job if it's the first step
    if workflow_plan and workflow_plan[0].endswith(":preprocessing"):
        return [Job(job_id=job_id, job_type="preprocessing", context=ctx, origin="preprocessing", workflow=workflow_plan, step=0)]
    
    # Create first job in plan (skip preprocessing if not present)
    if len(workflow_plan) > 0:
        queue_type, job_type = workflow_plan[0].split(":", 1)
        return [Job(job_id=job_id, job_type=job_type, context=ctx, origin=queue_type, workflow=workflow_plan, step=0)]
    
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
        
        ctx.add_result(job.job_type, {
            "command": command,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "returncode": result.returncode
        })
        
        if result.returncode != 0:
            ctx.add_error(job.job_type, f"Command failed with return code {result.returncode}")
            
    except subprocess.TimeoutExpired:
        ctx.add_error(job.job_type, "Command timed out")
    except Exception as e:
        ctx.add_error(job.job_type, f"Command error: {str(e)}")


# === Workflow Runner ===
class WorkflowRunner:
    """High-level interface for running workflows."""
    
    def __init__(self, verbose: bool = False):
        self.manager = WorkflowManager(verbose=verbose)
        self.verbose = verbose
        
        # Register default handlers
        self.manager.register_handler("preprocessing", "echo", echo_handler)
        self.manager.register_handler("slow", "sleep", sleep_handler)

    def register_handler(self, queue_type: str, job_type: str, handler: Callable[[Job], None]) -> None:
        """Register a custom job handler."""
        self.manager.register_handler(queue_type, job_type, handler)

    def register_command_handler(self, queue_type: str, job_type: str, command_template: str) -> None:
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
        process_existing: bool = True
    ) -> None:
        """Run a complete workflow."""
        
        if classifier_func is None:
            classifier_func = lambda filepath: default_file_classifier(filepath, workflow_plan)

        watcher = FileWatcher(
            watch_dir=watch_dir,
            preprocessor_func=classifier_func,
            manager=self.manager,
            recursive=recursive,
            patterns=patterns,
            ignore_patterns=ignore_patterns,
            verbose=self.verbose
        )

        watcher.start(process_existing=process_existing)

        try:
            self.manager.run()
        except KeyboardInterrupt:
            watcher.stop()
            if self.verbose:
                print("[WorkflowRunner] Shutdown requested.")
                stats = self.manager.get_stats()
                print(f"[WorkflowRunner] Final stats: {stats}") 