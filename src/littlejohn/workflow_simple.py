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
        self.preprocessing_queue = queue.Queue()  # bam_preprocessing only
        self.bed_conversion_queue = queue.Queue() # bed_conversion only
        self.analysis_queue = queue.Queue()       # mgmt, cnv, target, fusion
        self.classification_queue = queue.Queue() # sturgeon, nanodx, pannanodx
        self.slow_queue = queue.Queue()           # slow jobs (legacy)
        
        self.completed_jobs = []
        self.failed_jobs = []
        
        # Specialized job handlers for each queue
        self.job_handlers_preprocessing: Dict[str, Callable[[Job], None]] = {}
        self.job_handlers_bed_conversion: Dict[str, Callable[[Job], None]] = {}
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
        self.bed_conversion_workers = []
        self.analysis_workers = []
        self.classification_workers = []
        self.slow_workers = []
        
        # Progress tracking
        self.total_jobs_enqueued = 0
        self.total_jobs_skipped = 0  # Track jobs skipped due to deduplication
        self.active_jobs = {}  # job_id -> job_info
        self.job_start_times = {}  # job_id -> start_time
        self.completed_jobs_by_type = {}  # job_type -> count (successful completions only)
        self.failed_jobs_by_type = {}  # job_type -> count (failed jobs only)
        self.progress_lock = threading.Lock()  # Lock for progress tracking
        
        # Define job type mappings to queues
        self.job_queue_mapping = {
            # Preprocessing queue
            "preprocessing": "preprocessing",
            
            # Bed conversion queue
            "bed_conversion": "bed_conversion",
            
            # Analysis queue  
            "mgmt": "analysis",
            "cnv": "analysis",
            "target": "analysis",
            "fusion": "analysis",
            "test": "analysis",  # For testing purposes
            "long": "analysis",   # For testing purposes
            "quick": "analysis",  # For testing purposes
            
            # Classification queue
            "sturgeon": "classification", 
            "nanodx": "classification",
            "pannanodx": "classification",
            
            # Slow queue
            "random_forest": "slow",
            "sleep": "slow",
            "echo": "slow"
        }

    def register_handler(self, queue_type: str, job_type: str, handler: Callable[[Job], None]) -> None:
        """Register a handler for a specific job type and queue."""
        if queue_type == "preprocessing":
            self.job_handlers_preprocessing[job_type] = handler
        elif queue_type == "bed_conversion":
            self.job_handlers_bed_conversion[job_type] = handler
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
                    self.total_jobs_skipped += 1  # Track skipped jobs
                    continue
                else:
                    # Mark this job as pending for the sample
                    self._mark_job_pending_for_sample(job.job_type, sample_id)
                    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
                    logger.debug(f"Marked {job.job_type} job as pending for sample {sample_id}")
            
            jobs_to_enqueue.append(job)
        
        # Enqueue the filtered jobs
        for job in jobs_to_enqueue:
            # Track total jobs enqueued
            self.total_jobs_enqueued += 1
            
            # Determine which queue to use based on job type
            queue_type = self.job_queue_mapping.get(job.job_type, "slow")
            
            # Log job enqueuing for debugging
            logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
            logger.info(f"Enqueuing job {job.job_id} ({job.job_type}) to {queue_type} queue")
            
            if queue_type == "preprocessing":
                self.preprocessing_queue.put(job)
            elif queue_type == "bed_conversion":
                self.bed_conversion_queue.put(job)
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
                
                # Track job start time and active status
                self.job_start_times[job.job_id] = time.time()
                self.active_jobs[job.job_id] = {
                    'job_type': job.job_type,
                    'filepath': job.context.filepath,
                    'worker': name,
                    'start_time': time.time()
                }
                
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
                    
                    # Track this job as completed regardless of whether it has a next job
                    logger.info(f"Job {job.job_id} ({job.job_type}) completed successfully")
                    self.completed_jobs.append(job.job_id)
                    
                    # Track completion by job type
                    with self.progress_lock:
                        if job.job_type not in self.completed_jobs_by_type:
                            self.completed_jobs_by_type[job.job_type] = 0
                        self.completed_jobs_by_type[job.job_type] += 1
                    
                    # Check if there's a next job in the workflow
                    next_job = job.next_job()
                    if next_job:
                        logger.info(f"Creating next job: {next_job.job_type}")
                        self.enqueue_jobs([next_job])
                    else:
                        logger.info(f"No more jobs in workflow for {job.context.filepath}")
                        logger.debug(f"Final results: {job.context.results}")
                        
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
                    # Remove from active jobs
                    if job.job_id in self.active_jobs:
                        del self.active_jobs[job.job_id]
                    if job.job_id in self.job_start_times:
                        del self.job_start_times[job.job_id]
                    
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
        # Create preprocessing worker (bam_preprocessing only)
        preprocessing_worker = threading.Thread(
            target=self.worker, 
            args=("PreprocessingWorker", self.preprocessing_queue, self.job_handlers_preprocessing)
        )
        preprocessing_worker.daemon = True
        self.preprocessing_workers.append(preprocessing_worker)
        preprocessing_worker.start()
        
        # Create bed conversion worker (bed_conversion only)
        bed_conversion_worker = threading.Thread(
            target=self.worker, 
            args=("BedConversionWorker", self.bed_conversion_queue, self.job_handlers_bed_conversion)
        )
        bed_conversion_worker.daemon = True
        self.bed_conversion_workers.append(bed_conversion_worker)
        bed_conversion_worker.start()
        
        # Create analysis worker (mgmt, cnv, target, fusion)
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
        
        # Wait for all workers to finish with periodic checks for shutdown
        all_workers = (self.preprocessing_workers + self.bed_conversion_workers + self.analysis_workers + 
                      self.classification_workers + self.slow_workers)
        
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
        all_workers = (self.preprocessing_workers + self.bed_conversion_workers + self.analysis_workers + 
                      self.classification_workers + self.slow_workers)
        
        # Wait for each worker with timeout
        for worker in all_workers:
            if worker.is_alive():
                worker.join(timeout=timeout)
                if worker.is_alive():
                    if self.verbose:
                        print(f"[WorkflowManager] Warning: Worker {worker.name} did not stop within timeout")
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
        analysis_queue_size = self.analysis_queue.qsize()
        classification_queue_size = self.classification_queue.qsize()
        slow_queue_size = self.slow_queue.qsize()
        
        # Get active jobs by worker type
        active_by_worker = {}
        for job_info in self.active_jobs.values():
            worker = job_info['worker']
            if worker not in active_by_worker:
                active_by_worker[worker] = []
            active_by_worker[worker].append({
                'job_type': job_info['job_type'],
                'filepath': job_info['filepath'],
                'duration': time.time() - job_info['start_time']
            })
        
        # Calculate total expected jobs based on completed + failed + active + queued
        # This gives us the current total of all jobs that have been processed or are in progress
        total_expected = (len(self.completed_jobs) + len(self.failed_jobs) + 
                         len(self.active_jobs) + 
                         preprocessing_queue_size + bed_conversion_queue_size + analysis_queue_size + 
                         classification_queue_size + slow_queue_size)
        
        # If we have no jobs in progress but have completed jobs, use the completed count as the total
        if total_expected == 0 and len(self.completed_jobs) > 0:
            total_expected = len(self.completed_jobs)
        
        # Ensure total_expected is at least as large as the number of completed jobs
        if total_expected < len(self.completed_jobs):
            total_expected = len(self.completed_jobs)
        
        # Calculate total jobs that were actually enqueued (excluding skipped jobs)
        total_actual_jobs = self.total_jobs_enqueued
        
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
            "queue_sizes": {
                "preprocessing": preprocessing_queue_size,
                "bed_conversion": bed_conversion_queue_size,
                "analysis": analysis_queue_size,
                "classification": classification_queue_size,
                "slow": slow_queue_size
            },
            "active_by_worker": active_by_worker
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
        verbose: bool = False,
        show_progress: bool = True
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
                disable=not self.show_progress
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
        process_existing: bool = True,
        show_progress: bool = True
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
            verbose=self.verbose,
            show_progress=show_progress
        )

        watcher.start(process_existing=process_existing)

        # Start progress monitoring if enabled
        progress_thread = None
        if show_progress:
            progress_thread = threading.Thread(
                target=self._monitor_progress,
                args=(watcher,),
                daemon=True
            )
            progress_thread.start()

        try:
            self.manager.run()
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
                    print("[WorkflowRunner] Warning: Some workers may not have stopped gracefully.")
                stats = self.manager.get_stats()
                print(f"[WorkflowRunner] Final stats: {stats}")
                print(f"[WorkflowRunner] Final completion by type: {stats.get('completed_by_type', {})}")
                print(f"[WorkflowRunner] Final failed by type: {stats.get('failed_by_type', {})}")

    def _monitor_progress(self, watcher) -> None:
        """Monitor and display worker progress in real-time."""
        import time
        
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
        analysis_pbar = tqdm(
            desc="Analysis", 
            unit="jobs",
            position=2,
            leave=True
        )
        classification_pbar = tqdm(
            desc="Classification",
            unit="jobs", 
            position=3,
            leave=True
        )
        slow_pbar = tqdm(
            desc="Slow",
            unit="jobs",
            position=4,
            leave=True
        )
        
        # Create overall progress bar
        overall_pbar = tqdm(
            desc="Overall Progress",
            unit="jobs",
            position=5,
            leave=True
        )
        
        progress_bars = {
            "preprocessing": preprocessing_pbar,
            "bed_conversion": bed_conversion_pbar,
            "analysis": analysis_pbar,
            "classification": classification_pbar,
            "slow": slow_pbar
        }
        
        # Track job type completion counts
        job_type_completion = {
            "preprocessing": 0,
            "bed_conversion": 0,
            "analysis": 0, 
            "classification": 0,
            "slow": 0
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
                            if queue_name == "preprocessing" and "PreprocessingWorker" in worker_name:
                                active_in_queue += 1
                                filename = job['filepath'].split('/')[-1]
                                duration = int(job['duration'])
                                active_jobs_info.append(f"{job['job_type']}:{filename}({duration}s)")
                            elif queue_name == "bed_conversion" and "BedConversionWorker" in worker_name:
                                active_in_queue += 1
                                filename = job['filepath'].split('/')[-1]
                                duration = int(job['duration'])
                                active_jobs_info.append(f"{job['job_type']}:{filename}({duration}s)")
                            elif queue_name == "analysis" and "AnalysisWorker" in worker_name:
                                active_in_queue += 1
                                filename = job['filepath'].split('/')[-1]
                                duration = int(job['duration'])
                                active_jobs_info.append(f"{job['job_type']}:{filename}({duration}s)")
                            elif queue_name == "classification" and "ClassificationWorker" in worker_name:
                                active_in_queue += 1
                                filename = job['filepath'].split('/')[-1]
                                duration = int(job['duration'])
                                active_jobs_info.append(f"{job['job_type']}:{filename}({duration}s)")
                            elif queue_name == "slow" and "SlowWorker" in worker_name:
                                active_in_queue += 1
                                filename = job['filepath'].split('/')[-1]
                                duration = int(job['duration'])
                                active_jobs_info.append(f"{job['job_type']}:{filename}({duration}s)")
                    
                    # Calculate completed jobs for this queue type (both successful and failed)
                    completed_in_queue = 0
                    for job_type, count in stats["completed_by_type"].items():
                        # Get queue type from the manager's job queue mapping
                        queue_type = self.manager.job_queue_mapping.get(job_type, "slow")
                        if queue_type == queue_name:
                            completed_in_queue += count
                    
                    # Add failed jobs to the completion count
                    for job_type, count in stats.get("failed_by_type", {}).items():
                        # Get queue type from the manager's job queue mapping
                        queue_type = self.manager.job_queue_mapping.get(job_type, "slow")
                        if queue_type == queue_name:
                            completed_in_queue += count
                    
                    # Update the progress bar description with queue info
                    pbar.set_description(f"{queue_name.title()} (Q:{queue_size} A:{active_in_queue} C:{completed_in_queue})")
                    
                    # Update active jobs info
                    if active_jobs_info:
                        pbar.set_postfix_str(" | ".join(active_jobs_info[:2]))  # Show first 2 active jobs
                
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