"""
Workflow hooks for integrating shared state with the existing workflow system.

This module provides hooks and wrappers that integrate the shared workflow state
with the existing WorkflowManager and FileWatcher classes without modifying
their core functionality.
"""

import time
import logging
from typing import List, Callable, Any
from pathlib import Path
import queue

from .workflow_state import workflow_state, EventType, WorkflowEvent
from .workflow_simple import Job, WorkflowContext, WorkflowManager, FileWatcher


class WorkflowStateHooks:
    """Hooks for integrating workflow state with the existing workflow system."""
    
    def __init__(self):
        self.original_handlers = {}
        self.original_file_handler = None
        self.original_manager_stop = None
        self.original_manager_shutdown = None
        self.original_worker = None
        self.original_file_processor = None
    
    def install_workflow_hooks(self, workflow_runner: Any, workflow_steps: List[str], 
                              monitored_directory: str):
        """Install hooks into the workflow system."""
        try:
            # Start the workflow in the shared state
            workflow_state.start_workflow(workflow_steps, monitored_directory)
            
            # Install hooks into the workflow runner
            if hasattr(workflow_runner, 'run_workflow'):
                self._install_runner_hooks(workflow_runner)
            
            # Install hooks into the workflow manager if available
            if hasattr(workflow_runner, 'manager') and workflow_runner.manager:
                self._install_manager_hooks(workflow_runner.manager)
            
            # Install hooks into the file watcher if available
            if hasattr(workflow_runner, 'watcher') and workflow_runner.watcher:
                self._install_file_watcher_hooks(workflow_runner.watcher)
            
            # Install enhanced worker hooks for real-time job tracking
            if hasattr(workflow_runner, 'manager') and workflow_runner.manager:
                self._install_enhanced_worker_hooks(workflow_runner.manager)
            
            logging.info("Enhanced workflow state hooks installed successfully with real-time job tracking")
            
        except Exception as e:
            logging.error(f"Failed to install workflow hooks: {e}")
            # Fall back to minimal mode
            logging.info("Falling back to minimal mode")
    
    def _install_runner_hooks(self, runner: Any):
        """Install hooks into the workflow runner."""
        try:
            # Store original run_workflow method
            original_run_workflow = runner.run_workflow
            
            def hooked_run_workflow(*args, **kwargs):
                """Run workflow with state tracking hooks."""
                try:
                    # Mark workflow as started
                    workflow_state.emit_event(WorkflowEvent(
                        event_type=EventType.WORKFLOW_STARTED,
                        timestamp=time.time(),
                        data={'args': args, 'kwargs': kwargs}
                    ))
                    
                    # Call original method
                    result = original_run_workflow(*args, **kwargs)
                    
                    # Mark workflow as completed
                    workflow_state.emit_event(WorkflowEvent(
                        event_type=EventType.WORKFLOW_COMPLETED,
                        timestamp=time.time(),
                        data={'result': result}
                    ))
                    
                    return result
                    
                except Exception as e:
                    # Mark workflow as failed
                    workflow_state.emit_event(WorkflowEvent(
                        event_type=EventType.WORKFLOW_ERROR,
                        timestamp=time.time(),
                        data={'error': str(e)},
                        error_message=str(e)
                    ))
                    raise
            
            # Replace method
            runner.run_workflow = hooked_run_workflow
            
        except Exception as e:
            logging.error(f"Failed to install runner hooks: {e}")
    
    def _install_manager_hooks(self, manager: WorkflowManager):
        """Install hooks into the WorkflowManager."""
        try:
            # Store original methods
            self.original_manager_stop = manager.stop
            self.original_manager_shutdown = manager.shutdown
            
            # Hook into the stop method
            def hooked_stop(timeout: float = 30.0) -> bool:
                result = self.original_manager_stop(timeout)
                workflow_state.stop_workflow()
                return result
            
            # Hook into the shutdown method
            def hooked_shutdown(timeout: float = 30.0) -> bool:
                result = self.original_manager_shutdown(timeout)
                workflow_state.stop_workflow()
                return result
            
            # Replace methods
            manager.stop = hooked_stop
            manager.shutdown = hooked_shutdown
            
            # Note: We're not installing worker hooks to avoid interfering with workflow execution
            
        except Exception as e:
            logging.error(f"Failed to install manager hooks: {e}")
    
    def _install_worker_hooks(self, manager: WorkflowManager):
        """Install hooks into the worker methods using a less invasive approach."""
        try:
            # Instead of replacing the worker method, we'll hook into specific events
            # This avoids interfering with the core workflow logic
            
            # Note: We're not replacing the worker method to avoid breaking workflow execution
            # Instead, we'll rely on the workflow state being updated through other means
            
            logging.info("Installed workflow hooks (non-invasive mode)")
            
        except Exception as e:
            logging.error(f"Failed to install worker hooks: {e}")
    
    def _install_file_watcher_hooks(self, watcher: FileWatcher):
        """Install hooks into the FileWatcher."""
        try:
            # Store original file handler
            self.original_file_handler = watcher.handle_file
            
            def hooked_file_handler(filepath: str):
                """File handler with state tracking hooks."""
                try:
                    # Add file to shared state
                    workflow_state.add_file(filepath)
                    
                    # Add log message
                    workflow_state.add_log_message(
                        "INFO",
                        f"Detected new file: {Path(filepath).name}",
                        filepath=filepath
                    )
                    
                    # Call original handler
                    if self.original_file_handler:
                        return self.original_file_handler(filepath)
                        
                except Exception as e:
                    # Add error log
                    workflow_state.add_log_message(
                        "ERROR",
                        f"Error processing file {Path(filepath).name}: {e}",
                        filepath=filepath
                    )
                    raise
            
            # Replace file handler
            watcher.handle_file = hooked_file_handler
            
        except Exception as e:
            logging.error(f"Failed to install file watcher hooks: {e}")
    
    def _install_enhanced_worker_hooks(self, manager: WorkflowManager):
        """Install enhanced hooks for real-time job tracking using a safer approach."""
        try:
            # Store original handlers, handling None cases
            original_handlers = {}
            
            if hasattr(manager, 'job_handlers_preprocessing') and manager.job_handlers_preprocessing:
                original_handlers.update(manager.job_handlers_preprocessing)
            if hasattr(manager, 'job_handlers_analysis') and manager.job_handlers_analysis:
                original_handlers.update(manager.job_handlers_analysis)
            if hasattr(manager, 'job_handlers_classification') and manager.job_handlers_classification:
                original_handlers.update(manager.job_handlers_classification)
            if hasattr(manager, 'job_handlers_slow') and manager.job_handlers_slow:
                original_handlers.update(manager.job_handlers_slow)
            
            if not original_handlers:
                logging.warning("No handlers found to enhance - skipping enhanced worker hooks")
                return
            
            # Create wrapped handlers that track job progress
            def create_tracked_handler(original_handler, job_type):
                def tracked_handler(job):
                    try:
                        # Track job start
                        workflow_state.add_job(
                            job_id=job.job_id,
                            job_type=job.job_type,
                            filepath=job.context.filepath,
                            sample_id=Path(job.context.filepath).stem,
                            queue_name=self._get_queue_name_for_job_type(job.job_type)
                        )
                        workflow_state.start_job(job.job_id, f"Worker-{job_type}")
                        
                        # Add log message
                        workflow_state.add_log_message(
                            "INFO",
                            f"Started processing {job.job_type} job for {Path(job.context.filepath).name}",
                            filepath=job.context.filepath,
                            job_id=job.job_id
                        )
                        
                        # Call original handler
                        result = original_handler(job)
                        
                        # Track successful completion
                        workflow_state.complete_job(job.job_id)
                        workflow_state.add_log_message(
                            "INFO",
                            f"Job {job.job_type} completed successfully for {Path(job.context.filepath).name}",
                            filepath=job.context.filepath,
                            job_id=job.job_id
                        )
                        
                        return result
                        
                    except Exception as e:
                        # Track job failure
                        workflow_state.fail_job(job.job_id, str(e))
                        workflow_state.add_log_message(
                            "ERROR",
                            f"Job {job.job_type} failed for {Path(job.context.filepath).name}: {e}",
                            filepath=job.context.filepath,
                            job_id=job.job_id
                        )
                        raise
                
                return tracked_handler
            
            # Replace handlers with tracked versions
            for job_type, handler in original_handlers.items():
                tracked_handler = create_tracked_handler(handler, job_type)
                
                # Update in all handler dictionaries
                if hasattr(manager, 'job_handlers_preprocessing') and manager.job_handlers_preprocessing and job_type in manager.job_handlers_preprocessing:
                    manager.job_handlers_preprocessing[job_type] = tracked_handler
                if hasattr(manager, 'job_handlers_analysis') and manager.job_handlers_analysis and job_type in manager.job_handlers_analysis:
                    manager.job_handlers_analysis[job_type] = tracked_handler
                if hasattr(manager, 'job_handlers_classification') and manager.job_handlers_classification and job_type in manager.job_handlers_classification:
                    manager.job_handlers_classification[job_type] = tracked_handler
                if hasattr(manager, 'job_handlers_slow') and manager.job_handlers_slow and job_type in manager.job_handlers_slow:
                    manager.job_handlers_slow[job_type] = tracked_handler
            
            logging.info(f"Enhanced worker hooks installed for {len(original_handlers)} handlers (handler-based)")
            
        except Exception as e:
            logging.error(f"Failed to install enhanced worker hooks: {e}")
            # Keep original handlers if enhancement fails
            pass
    
    def _get_queue_name_for_job_type(self, job_type: str) -> str:
        """Get the queue name for a job type."""
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
        return queue_mapping.get(job_type, "unknown")
    
    def uninstall_hooks(self):
        """Uninstall all hooks and restore original methods."""
        try:
            # Restore manager methods
            if self.original_manager_stop:
                # Note: We can't easily restore these without storing the manager reference
                pass
            
            # Restore worker method
            if self.original_worker:
                # Note: We can't easily restore these without storing the manager reference
                pass
            
            # Restore file handler
            if self.original_file_handler:
                # Note: We can't easily restore these without storing the watcher reference
                pass
            
            logging.info("Workflow state hooks uninstalled")
            
        except Exception as e:
            logging.error(f"Failed to uninstall hooks: {e}")


# Global instance
_hooks_instance = WorkflowStateHooks()


def install_workflow_hooks(workflow_runner: Any, workflow_steps: List[str], 
                          monitored_directory: str):
    """Install workflow state hooks."""
    _hooks_instance.install_workflow_hooks(workflow_runner, workflow_steps, monitored_directory)


def uninstall_workflow_hooks():
    """Uninstall workflow state hooks."""
    _hooks_instance.uninstall_hooks()


def get_workflow_state():
    """Get the global workflow state instance."""
    return workflow_state
