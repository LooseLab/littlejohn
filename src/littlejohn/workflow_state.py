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
    """Information about a specific job."""
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
    progress: float = 0.0  # 0.0 to 1.0


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
        self.jobs_by_type: Dict[str, List[int]] = {}
        self.jobs_by_sample: Dict[str, List[int]] = {}
        self.jobs_by_file: Dict[str, List[int]] = {}
        
        # Queue status
        self.queue_status: Dict[str, QueueStatus] = {
            "preprocessing": QueueStatus("preprocessing"),
            "bed_conversion": QueueStatus("bed_conversion"),
            "mgmt": QueueStatus("mgmt"),
            "cnv": QueueStatus("cnv"),
            "target": QueueStatus("target"),
            "fusion": QueueStatus("fusion"),
            "classification": QueueStatus("classification"),
            "slow": QueueStatus("slow"),
        }
        
        # File tracking
        self.detected_files: Set[str] = set()
        self.processed_files: Set[str] = set()
        self.failed_files: Set[str] = set()
        
        # Event system
        self.event_listeners: List[Callable[[WorkflowEvent], None]] = []
        self.event_queue: queue.Queue = queue.Queue()
        self.event_thread: Optional[threading.Thread] = None
        
        # Logging
        self.log_messages: List[Dict[str, Any]] = []
        self.max_log_messages: int = 1000
        
        # Thread safety
        self._lock = threading.RLock()
        
        # Start event processing thread
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
            self.workflow_start_time = time.time()
            self.workflow_end_time = None
            self.is_running = True
            self.workflow_steps = workflow_steps.copy()
            self.monitored_directory = monitored_directory
            
            # Reset state
            self.jobs.clear()
            self.jobs_by_type.clear()
            self.jobs_by_sample.clear()
            self.jobs_by_file.clear()
            self.detected_files.clear()
            self.processed_files.clear()
            self.failed_files.clear()
            self.log_messages.clear()
            
            # Reset queue status
            for queue_status in self.queue_status.values():
                queue_status.pending_count = 0
                queue_status.running_count = 0
                queue_status.completed_count = 0
                queue_status.failed_count = 0
                queue_status.total_count = 0
            
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
        """Add a new job to the workflow."""
        with self._lock:
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
                progress=0.0
            )
            
            self.jobs[job_id] = job_info
            
            # Update job type tracking
            if job_type not in self.jobs_by_type:
                self.jobs_by_type[job_type] = []
            self.jobs_by_type[job_type].append(job_id)
            
            # Update sample tracking
            if sample_id not in self.jobs_by_sample:
                self.jobs_by_sample[sample_id] = []
            self.jobs_by_sample[sample_id].append(job_id)
            
            # Update file tracking
            if filepath not in self.jobs_by_file:
                self.jobs_by_file[filepath] = []
            self.jobs_by_file[filepath].append(job_id)
            
            # Update queue status
            if queue_name in self.queue_status:
                self.queue_status[queue_name].pending_count += 1
                self.queue_status[queue_name].total_count += 1
            
            # Emit job started event
            self.emit_event(WorkflowEvent(
                event_type=EventType.JOB_STARTED,
                timestamp=time.time(),
                job_id=job_id,
                job_type=job_type,
                filepath=filepath,
                sample_id=sample_id,
                data={"queue_name": queue_name}
            ))
    
    def start_job(self, job_id: int, worker_name: str):
        """Mark a job as started."""
        with self._lock:
            if job_id in self.jobs:
                job = self.jobs[job_id]
                job.status = JobStatus.RUNNING
                job.start_time = time.time()
                job.worker_name = worker_name
                
                # Update queue status
                queue_name = self._get_queue_name_for_job_type(job.job_type)
                if queue_name in self.queue_status:
                    self.queue_status[queue_name].pending_count -= 1
                    self.queue_status[queue_name].running_count += 1
    
    def complete_job(self, job_id: int):
        """Mark a job as completed."""
        with self._lock:
            if job_id in self.jobs:
                job = self.jobs[job_id]
                job.status = JobStatus.COMPLETED
                job.end_time = time.time()
                job.duration = job.end_time - job.start_time if job.start_time else 0
                job.progress = 1.0
                
                # Update queue status
                queue_name = self._get_queue_name_for_job_type(job.job_type)
                if queue_name in self.queue_status:
                    self.queue_status[queue_name].running_count -= 1
                    self.queue_status[queue_name].completed_count += 1
                
                # Mark file as processed
                self.processed_files.add(job.filepath)
                
                # Emit job completed event
                self.emit_event(WorkflowEvent(
                    event_type=EventType.JOB_COMPLETED,
                    timestamp=time.time(),
                    job_id=job_id,
                    job_type=job.job_type,
                    filepath=job.filepath,
                    sample_id=job.sample_id,
                    data={"duration": job.duration}
                ))
    
    def fail_job(self, job_id: int, error_message: str):
        """Mark a job as failed."""
        with self._lock:
            if job_id in self.jobs:
                job = self.jobs[job_id]
                job.status = JobStatus.FAILED
                job.end_time = time.time()
                job.duration = job.end_time - job.start_time if job.start_time else 0
                job.error_message = error_message
                
                # Update queue status
                queue_name = self._get_queue_name_for_job_type(job.job_type)
                if queue_name in self.queue_status:
                    self.queue_status[queue_name].running_count -= 1
                    self.queue_status[queue_name].failed_count += 1
                
                # Mark file as failed
                self.failed_files.add(job.filepath)
                
                # Emit job failed event
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
            self.detected_files.add(filepath)
            
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
                "level": level,
                "message": message,
                "job_id": job_id,
                "job_type": job_type,
                "filepath": filepath
            }
            
            self.log_messages.append(log_entry)
            
            # Keep only the most recent messages
            if len(self.log_messages) > self.max_log_messages:
                self.log_messages = self.log_messages[-self.max_log_messages:]
            
            # Emit log message event
            self.emit_event(WorkflowEvent(
                event_type=EventType.LOG_MESSAGE,
                timestamp=time.time(),
                job_id=job_id,
                job_type=job_type,
                filepath=filepath,
                data={"level": level, "message": message}
            ))
    
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
        """Get a summary of the current workflow state."""
        with self._lock:
            total_jobs = len(self.jobs)
            completed_jobs = len([j for j in self.jobs.values() if j.status == JobStatus.COMPLETED])
            failed_jobs = len([j for j in self.jobs.values() if j.status == JobStatus.FAILED])
            running_jobs = len([j for j in self.jobs.values() if j.status == JobStatus.RUNNING])
            pending_jobs = len([j for j in self.jobs.values() if j.status == JobStatus.PENDING])
            
            duration = 0
            if self.workflow_start_time:
                if self.workflow_end_time:
                    duration = self.workflow_end_time - self.workflow_start_time
                else:
                    duration = time.time() - self.workflow_start_time
            
            return {
                "is_running": self.is_running,
                "start_time": self.workflow_start_time,
                "end_time": self.workflow_end_time,
                "duration": duration,
                "workflow_steps": self.workflow_steps,
                "monitored_directory": self.monitored_directory,
                "total_jobs": total_jobs,
                "completed_jobs": completed_jobs,
                "failed_jobs": failed_jobs,
                "running_jobs": running_jobs,
                "pending_jobs": pending_jobs,
                "detected_files": len(self.detected_files),
                "processed_files": len(self.processed_files),
                "failed_files": len(self.failed_files),
                "queue_status": {name: {
                    "pending": qs.pending_count,
                    "running": qs.running_count,
                    "completed": qs.completed_count,
                    "failed": qs.failed_count,
                    "total": qs.total_count
                } for name, qs in self.queue_status.items()}
            }
    
    def get_jobs_by_status(self, status: JobStatus) -> List[JobInfo]:
        """Get all jobs with a specific status."""
        with self._lock:
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
            return self.log_messages[-count:] if self.log_messages else []
    
    def get_logs_by_level(self, level: str) -> List[Dict[str, Any]]:
        """Get all log messages of a specific level."""
        with self._lock:
            return [log for log in self.log_messages if log["level"] == level]


# Global workflow state instance
workflow_state = WorkflowState()
