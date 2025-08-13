"""
Shared workflow state and event system for real-time GUI updates.

This module provides a centralized state management system that allows the GUI
to receive real-time updates about workflow progress, job status, and events.
"""

import threading
import time
import queue
from typing import Dict, List, Set, Optional, Any, Callable
from dataclasses import dataclass, field
from enum import Enum
import logging


class JobStatus(Enum):
    """Status of a job in the workflow."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class EventType(Enum):
    """Types of events that can occur in the workflow."""
    JOB_STARTED = "job_started"
    JOB_COMPLETED = "job_completed"
    JOB_FAILED = "job_failed"
    JOB_CANCELLED = "job_cancelled"
    WORKFLOW_STARTED = "workflow_started"
    WORKFLOW_COMPLETED = "workflow_completed"
    WORKFLOW_ERROR = "workflow_error"
    FILE_DETECTED = "file_detected"
    PROGRESS_UPDATE = "progress_update"
    LOG_MESSAGE = "log_message"


@dataclass
class WorkflowEvent:
    """Represents an event in the workflow."""
    event_type: EventType
    timestamp: float
    data: Dict[str, Any] = field(default_factory=dict)
    job_id: Optional[int] = None
    job_type: Optional[str] = None
    filepath: Optional[str] = None
    sample_id: Optional[str] = None
    error_message: Optional[str] = None


@dataclass
class JobInfo:
    """Information about a job in the workflow."""
    job_id: int
    job_type: str
    filepath: str
    sample_id: str
    status: JobStatus
    start_time: Optional[float] = None
    end_time: Optional[float] = None
    duration: Optional[float] = None
    worker_name: Optional[str] = None
    error_message: Optional[str] = None
    progress: float = 0.0
    # Add fields for storing actual job results and metadata
    job_metadata: Dict[str, Any] = field(default_factory=dict)
    job_results: Dict[str, Any] = field(default_factory=dict)


@dataclass
class QueueStatus:
    """Status of a specific queue."""
    queue_name: str
    pending_count: int = 0
    running_count: int = 0
    completed_count: int = 0
    failed_count: int = 0
    total_count: int = 0


class WorkflowState:
    """
    Centralized state management for workflow execution.
    
    This class maintains the current state of the workflow and provides
    methods for updating state and notifying listeners of changes.
    """
    
    def __init__(self):
        # Core state
        self.workflow_start_time: Optional[float] = None
        self.workflow_end_time: Optional[float] = None
        self.is_running: bool = False
        self.workflow_steps: List[str] = []
        self.monitored_directory: str = ""
        
        # Job tracking
        self.jobs: Dict[int, JobInfo] = {}
        self.next_job_id: int = 1
        
        # Sample tracking - NEW: Track all sample IDs seen during processing
        self.samples: Dict[str, Dict[str, Any]] = {}
        self.sample_jobs: Dict[str, List[int]] = {}  # sample_id -> list of job_ids
        
        # Queue tracking
        self.queues: Dict[str, QueueStatus] = {}
        
        # Event system
        self.event_listeners: List[Callable[[WorkflowEvent], None]] = []
        self.event_queue: queue.Queue = queue.Queue()
        
        # Logging
        self.logs: List[Dict[str, Any]] = []
        self.log_listeners: List[Callable[[Dict[str, Any]], None]] = []
        
        # Threading
        self._lock = threading.Lock()
        self._event_thread = None
        self._stop_event = threading.Event()
        
        # Start event processing
        self._start_event_processor()
    
    def _start_event_processor(self):
        """Start the event processing thread."""
        self.event_thread = threading.Thread(target=self._process_events, daemon=True)
        self.event_thread.start()
    
    def _process_events(self):
        """Process events from the event queue."""
        while True:
            try:
                event = self.event_queue.get(timeout=1.0)
                self._notify_listeners(event)
                self.event_queue.task_done()
            except queue.Empty:
                continue
            except Exception as e:
                logging.error(f"Error processing event: {e}")
    
    def add_event_listener(self, listener: Callable[[WorkflowEvent], None]):
        """Add an event listener."""
        with self._lock:
            self.event_listeners.append(listener)
    
    def remove_event_listener(self, listener: Callable[[WorkflowEvent], None]):
        """Remove an event listener."""
        with self._lock:
            if listener in self.event_listeners:
                self.event_listeners.remove(listener)
    
    def _notify_listeners(self, event: WorkflowEvent):
        """Notify all event listeners of an event."""
        with self._lock:
            for listener in self.event_listeners:
                try:
                    listener(event)
                except Exception as e:
                    logging.error(f"Error in event listener: {e}")
    
    def emit_event(self, event: WorkflowEvent):
        """Emit an event to be processed."""
        self.event_queue.put(event)
    
    def start_workflow(self, workflow_steps: List[str], monitored_directory: str):
        """Start a new workflow."""
        with self._lock:
            # Reset state for new workflow
            self.workflow_start_time = time.time()
            self.workflow_end_time = None
            self.is_running = True
            self.workflow_steps = workflow_steps.copy()
            self.monitored_directory = monitored_directory
            
            # Reset state
            self.jobs.clear()
            self.samples.clear()
            self.sample_jobs.clear()
            self.queues.clear()
            self.logs.clear()
            
            # Emit workflow started event
            self.emit_event(WorkflowEvent(
                event_type=EventType.WORKFLOW_STARTED,
                timestamp=time.time(),
                data={
                    "workflow_steps": workflow_steps,
                    "monitored_directory": monitored_directory
                }
            ))
    
    def stop_workflow(self):
        """Stop the current workflow."""
        with self._lock:
            self.is_running = False
            self.workflow_end_time = time.time()
            
            # Emit workflow completed event
            self.emit_event(WorkflowEvent(
                event_type=EventType.WORKFLOW_COMPLETED,
                timestamp=time.time(),
                data={
                    "duration": self.workflow_end_time - self.workflow_start_time if self.workflow_start_time else 0
                }
            ))
    
    def add_job(self, job_id: int, job_type: str, filepath: str, sample_id: str, queue_name: str):
        """Add a new job to the workflow state."""
        with self._lock:
            # Create job info
            job_info = JobInfo(
                job_id=job_id,
                job_type=job_type,
                filepath=filepath,
                sample_id=sample_id,
                status=JobStatus.PENDING,
                start_time=None,
                end_time=None,
                duration=None,
                worker_name=None,
                error_message=None,
                progress=0.0,
                job_metadata={},
                job_results={}
            )
            
            # Store job
            self.jobs[job_id] = job_info
            
            # Track sample - NEW: Register sample when job is added
            self._register_sample(sample_id, filepath, job_type)
            
            # Update queue status
            if queue_name not in self.queues:
                self.queues[queue_name] = QueueStatus(queue_name)
            
            queue_status = self.queues[queue_name]
            queue_status.pending_count += 1
            queue_status.total_count += 1
            
            # Track sample jobs
            if sample_id not in self.sample_jobs:
                self.sample_jobs[sample_id] = []
            self.sample_jobs[sample_id].append(job_id)
            
            # Emit event
            self.emit_event(WorkflowEvent(
                event_type=EventType.JOB_STARTED,
                timestamp=time.time(),
                job_id=job_id,
                job_type=job_type,
                filepath=filepath,
                sample_id=sample_id
            ))
            
            # Update next job ID
            if job_id >= self.next_job_id:
                self.next_job_id = job_id + 1
    
    def _register_sample(self, sample_id: str, filepath: str, job_type: str):
        """Register a new sample or update existing sample information."""
        with self._lock:
            if sample_id not in self.samples:
                # New sample
                self.samples[sample_id] = {
                    'sample_id': sample_id,
                    'first_seen': time.time(),
                    'last_seen': time.time(),
                    'filepaths': set([filepath]),
                    'job_types': set([job_type]),
                    'total_jobs': 0,
                    'completed_jobs': 0,
                    'failed_jobs': 0,
                    'current_status': 'pending',
                    'metadata': {},
                    'results': {}
                }
            else:
                # Update existing sample
                sample_info = self.samples[sample_id]
                sample_info['last_seen'] = time.time()
                sample_info['filepaths'].add(filepath)
                sample_info['job_types'].add(job_type)
            
            # Update total jobs count
            self.samples[sample_id]['total_jobs'] += 1
    
    def get_sample_info(self, sample_id: str) -> Optional[Dict[str, Any]]:
        """Get comprehensive information about a specific sample."""
        with self._lock:
            if sample_id not in self.samples:
                return None
            
            sample_info = self.samples[sample_id].copy()
            
            # Get all jobs for this sample
            sample_job_ids = self.sample_jobs.get(sample_id, [])
            sample_jobs = [self.jobs[job_id] for job_id in sample_job_ids if job_id in self.jobs]
            
            # Add job details
            sample_info['jobs'] = []
            for job in sample_jobs:
                job_detail = {
                    'job_id': job.job_id,
                    'job_type': job.job_type,
                    'filepath': job.filepath,
                    'status': job.status.value,
                    'start_time': job.start_time,
                    'end_time': job.end_time,
                    'duration': job.duration,
                    'progress': job.progress,
                    'error_message': job.error_message,
                    'metadata': job.job_metadata,
                    'results': job.job_results
                }
                sample_info['jobs'].append(job_detail)
            
            # Calculate sample statistics
            sample_info['completed_jobs'] = len([j for j in sample_jobs if j.status == JobStatus.COMPLETED])
            sample_info['failed_jobs'] = len([j for j in sample_jobs if j.status == JobStatus.FAILED])
            sample_info['running_jobs'] = len([j for j in sample_jobs if j.status == JobStatus.RUNNING])
            sample_info['pending_jobs'] = len([j for j in sample_jobs if j.status == JobStatus.PENDING])
            
            # Determine overall sample status
            if sample_info['failed_jobs'] > 0:
                sample_info['current_status'] = 'failed'
            elif sample_info['completed_jobs'] == sample_info['total_jobs']:
                sample_info['current_status'] = 'completed'
            elif sample_info['running_jobs'] > 0:
                sample_info['current_status'] = 'running'
            else:
                sample_info['current_status'] = 'pending'
            
            # Convert sets to lists for JSON serialization
            sample_info['filepaths'] = list(sample_info['filepaths'])
            sample_info['job_types'] = list(sample_info['job_types'])
            
            return sample_info
    
    def get_all_samples(self) -> List[Dict[str, Any]]:
        """Get a list of all samples with summary information."""
        with self._lock:
            samples_summary = []
            for sample_id in self.samples:
                sample_info = self.samples[sample_id]
                summary = {
                    'sample_id': sample_id,
                    'first_seen': sample_info['first_seen'],
                    'last_seen': sample_info['last_seen'],
                    'total_jobs': sample_info['total_jobs'],
                    'completed_jobs': sample_info['completed_jobs'],
                    'failed_jobs': sample_info['failed_jobs'],
                    'current_status': sample_info['current_status'],
                    'job_types': list(sample_info['job_types']),
                    'filepaths': list(sample_info['filepaths'])
                }
                samples_summary.append(summary)
            
            # Sort by last seen (most recent first)
            samples_summary.sort(key=lambda x: x['last_seen'], reverse=True)
            return samples_summary
    
    def get_sample_count(self) -> int:
        """Get the total number of samples tracked."""
        with self._lock:
            return len(self.samples)
    
    def start_job(self, job_id: int, worker_name: str):
        """Mark a job as started."""
        with self._lock:
            if job_id not in self.jobs:
                return
            
            job = self.jobs[job_id]
            job.status = JobStatus.RUNNING
            job.start_time = time.time()
            job.worker_name = worker_name
            
            # Update sample tracking - NEW: Update sample status
            if job.sample_id in self.samples:
                # Update sample status to running if it was pending
                if self.samples[job.sample_id]['current_status'] == 'pending':
                    self.samples[job.sample_id]['current_status'] = 'running'
            
            # Update queue status
            queue_name = self._get_queue_name_for_job_type(job.job_type)
            if queue_name in self.queues:
                queue_status = self.queues[queue_name]
                queue_status.pending_count -= 1
                queue_status.running_count += 1
            
            # Emit event
            self.emit_event(WorkflowEvent(
                event_type=EventType.JOB_STARTED,
                timestamp=time.time(),
                job_id=job_id,
                job_type=job.job_type,
                filepath=job.filepath,
                sample_id=job.sample_id
            ))
    
    def complete_job(self, job_id: int):
        """Mark a job as completed."""
        with self._lock:
            if job_id not in self.jobs:
                return
            
            job = self.jobs[job_id]
            job.status = JobStatus.COMPLETED
            job.end_time = time.time()
            if job.start_time:
                job.duration = job.end_time - job.start_time
            job.progress = 1.0
            
            # Update sample tracking - NEW: Update sample statistics
            if job.sample_id in self.samples:
                self.samples[job.sample_id]['completed_jobs'] += 1
            
            # Update queue status
            queue_name = self._get_queue_name_for_job_type(job.job_type)
            if queue_name in self.queues:
                queue_status = self.queues[queue_name]
                queue_status.running_count -= 1
                queue_status.completed_count += 1
            
            # Emit event
            self.emit_event(WorkflowEvent(
                event_type=EventType.JOB_COMPLETED,
                timestamp=time.time(),
                job_id=job_id,
                job_type=job.job_type,
                filepath=job.filepath,
                sample_id=job.sample_id
            ))
    
    def store_job_metadata(self, job_id: int, metadata: Dict[str, Any]):
        """Store metadata for a specific job."""
        with self._lock:
            if job_id in self.jobs:
                self.jobs[job_id].job_metadata.update(metadata)
    
    def store_job_results(self, job_id: int, results: Dict[str, Any]):
        """Store results for a specific job."""
        with self._lock:
            if job_id in self.jobs:
                self.jobs[job_id].job_results.update(results)
    
    def get_job_metadata(self, job_id: int) -> Optional[Dict[str, Any]]:
        """Get metadata for a specific job."""
        with self._lock:
            if job_id in self.jobs:
                return self.jobs[job_id].job_metadata
            return None
    
    def get_job_results(self, job_id: int) -> Optional[Dict[str, Any]]:
        """Get results for a specific job."""
        with self._lock:
            if job_id in self.jobs:
                return self.jobs[job_id].job_results
            return None
    
    def fail_job(self, job_id: int, error_message: str):
        """Mark a job as failed."""
        with self._lock:
            if job_id not in self.jobs:
                return
            
            job = self.jobs[job_id]
            job.status = JobStatus.FAILED
            job.end_time = time.time()
            if job.start_time:
                job.duration = job.end_time - job.start_time
            job.error_message = error_message
            
            # Update sample tracking - NEW: Update sample statistics
            if job.sample_id in self.samples:
                self.samples[job.sample_id]['failed_jobs'] += 1
            
            # Update queue status
            queue_name = self._get_queue_name_for_job_type(job.job_type)
            if queue_name in self.queues:
                queue_status = self.queues[queue_name]
                queue_status.running_count -= 1
                queue_status.failed_count += 1
            
            # Emit event
            self.emit_event(WorkflowEvent(
                event_type=EventType.JOB_FAILED,
                timestamp=time.time(),
                job_id=job_id,
                job_type=job.job_type,
                filepath=job.filepath,
                sample_id=job.sample_id,
                error_message=error_message
            ))
    
    def update_job_progress(self, job_id: int, progress: float):
        """Update the progress of a job (0.0 to 1.0)."""
        with self._lock:
            if job_id in self.jobs:
                self.jobs[job_id].progress = max(0.0, min(1.0, progress))
                
                # Emit progress update event
                self.emit_event(WorkflowEvent(
                    event_type=EventType.PROGRESS_UPDATE,
                    timestamp=time.time(),
                    job_id=job_id,
                    job_type=self.jobs[job_id].job_type,
                    filepath=self.jobs[job_id].filepath,
                    sample_id=self.jobs[job_id].sample_id,
                    data={"progress": progress}
                ))
    
    def add_file(self, filepath: str):
        """Add a detected file to the workflow."""
        with self._lock:
            # self.detected_files.add(filepath) # This line was removed from the original file
            
            # Emit file detected event
            self.emit_event(WorkflowEvent(
                event_type=EventType.FILE_DETECTED,
                timestamp=time.time(),
                filepath=filepath
            ))
    
    def add_log_message(self, level: str, message: str, job_id: Optional[int] = None, 
                       job_type: Optional[str] = None, filepath: Optional[str] = None):
        """Add a log message to the workflow state."""
        with self._lock:
            log_entry = {
                "timestamp": time.time(),
                "level": level.upper(),
                "message": message,
                "job_id": job_id,
                "job_type": job_type,
                "filepath": filepath
            }
            
            self.logs.append(log_entry)
            
            # Limit log size
            if len(self.logs) > 1000:  # Keep last 1000 log entries
                self.logs = self.logs[-1000:]
            
            # Notify log listeners
            for listener in self.log_listeners:
                try:
                    listener(log_entry)
                except Exception as e:
                    logging.error(f"Error in log listener: {e}")
    
    def _get_queue_name_for_job_type(self, job_type: str) -> str:
        """Get the queue name for a given job type."""
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
    
    def get_workflow_summary(self) -> Dict[str, Any]:
        """Get a comprehensive summary of the workflow state."""
        with self._lock:
            # Calculate job counts by status
            pending_jobs = len([j for j in self.jobs.values() if j.status == JobStatus.PENDING])
            running_jobs = len([j for j in self.jobs.values() if j.status == JobStatus.RUNNING])
            completed_jobs = len([j for j in self.jobs.values() if j.status == JobStatus.COMPLETED])
            failed_jobs = len([j for j in self.jobs.values() if j.status == JobStatus.FAILED])
            total_jobs = len(self.jobs)
            
            # Calculate queue summaries
            queue_summaries = {}
            for queue_name, queue_status in self.queues.items():
                queue_summaries[queue_name] = {
                    'pending': queue_status.pending_count,
                    'running': queue_status.running_count,
                    'completed': queue_status.completed_count,
                    'failed': queue_status.failed_count,
                    'total': queue_status.total_count
                }
            
            # Sample tracking summary - NEW: Include sample information
            sample_summary = {
                'total_samples': len(self.samples),
                'completed_samples': len([s for s in self.samples.values() if s['current_status'] == 'completed']),
                'running_samples': len([s for s in self.samples.values() if s['current_status'] == 'running']),
                'failed_samples': len([s for s in self.samples.values() if s['current_status'] == 'failed']),
                'pending_samples': len([s for s in self.samples.values() if s['current_status'] == 'pending'])
            }
            
            return {
                'is_running': self.is_running,
                'start_time': self.workflow_start_time,
                'end_time': self.workflow_end_time,
                'workflow_steps': self.workflow_steps,
                'monitored_directory': self.monitored_directory,
                'pending_jobs': pending_jobs,
                'running_jobs': running_jobs,
                'completed_jobs': completed_jobs,
                'failed_jobs': failed_jobs,
                'total_jobs': total_jobs,
                'queue_summaries': queue_summaries,
                'sample_summary': sample_summary,  # NEW: Sample tracking summary
                'recent_logs': self.get_recent_logs(10)
            }
    
    def get_jobs_by_status(self, status) -> List[JobInfo]:
        """Get all jobs with a specific status.
        
        Args:
            status: Can be either a JobStatus enum value or a string representation
                   (e.g., "completed", "running", "failed", "pending")
        """
        with self._lock:
            # Handle both string and enum inputs
            if isinstance(status, str):
                # Convert string to JobStatus enum
                try:
                    status = JobStatus(status.lower())
                except ValueError:
                    # If invalid status string, return empty list
                    return []
            
            return [job for job in self.jobs.values() if job.status == status]
    
    def get_jobs_by_type(self, job_type: str) -> List[JobInfo]:
        """Get all jobs of a specific type."""
        with self._lock:
            return [job for job in self.jobs.values() if job.job_type == job_type]
    
    def get_jobs_by_sample(self, sample_id: str) -> List[JobInfo]:
        """Get all jobs for a specific sample."""
        with self._lock:
            return [job for job in self.jobs.values() if job.sample_id == sample_id]
    
    def get_jobs_by_file(self, filepath: str) -> List[JobInfo]:
        """Get all jobs for a specific file."""
        with self._lock:
            return [job for job in self.jobs.values() if job.filepath == filepath]
    
    def get_recent_logs(self, count: int = 100) -> List[Dict[str, Any]]:
        """Get the most recent log messages."""
        with self._lock:
            return self.logs[-count:] if self.logs else []
    
    def get_logs_by_level(self, level: str) -> List[Dict[str, Any]]:
        """Get all log messages of a specific level."""
        with self._lock:
            return [log for log in self.logs if log["level"] == level]


# Global workflow state instance
workflow_state = WorkflowState()
