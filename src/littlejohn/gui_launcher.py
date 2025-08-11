"""
GUI launcher for LittleJohn workflow monitoring.

This module provides a clean interface for launching the workflow GUI
that doesn't interfere with CLI output and properly integrates with
the workflow state system.
"""

import threading
import time
import logging
from typing import Optional, Dict, Any
from pathlib import Path

try:
    from nicegui import ui, app
except ImportError:
    ui = None
    app = None


class GUILauncher:
    """Launcher for the LittleJohn workflow GUI."""
    
    def __init__(self, host: str = "localhost", port: int = 8081):
        self.host = host
        self.port = port
        self.gui_thread = None
        self.is_running = False
        self.workflow_runner = None
        self.workflow_steps = []
        self.monitored_directory = ""
        
    def launch_gui(self, workflow_runner: Any = None, workflow_steps: list = None, 
                   monitored_directory: str = "") -> bool:
        """Launch the GUI in a background thread."""
        if ui is None:
            logging.error("NiceGUI is not available")
            return False
            
        self.workflow_runner = workflow_runner
        self.workflow_steps = workflow_steps or []
        self.monitored_directory = monitored_directory
        
        try:
            # Start GUI in background thread
            self.gui_thread = threading.Thread(
                target=self._run_gui_worker,
                daemon=True
            )
            self.gui_thread.start()
            
            # Wait a moment for GUI to start
            time.sleep(1)
            self.is_running = True
            
            logging.info(f"GUI launched successfully on http://{self.host}:{self.port}")
            return True
            
        except Exception as e:
            logging.error(f"Failed to launch GUI: {e}")
            return False
    
    def _run_gui_worker(self):
        """Run the GUI worker in a separate thread."""
        try:
            # Create the main interface
            @ui.page('/')
            def workflow_monitor():
                self._create_workflow_monitor()
            
            # Start the NiceGUI server
            ui.run(
                host=self.host,
                port=self.port,
                show=False,  # Don't show browser automatically
                reload=False  # Disable auto-reload for production
            )
            
        except Exception as e:
            logging.error(f"GUI worker failed: {e}")
            self.is_running = False
    
    def _create_workflow_monitor(self):
        """Create the workflow monitoring interface."""
        ui.add_head_html('<title>LittleJohn Workflow Monitor</title>')
        
        # Initialize UI state variables
        status_label = ui.label('Initializing...').classes('text-lg font-semibold text-gray-600')
        started_time_label = ui.label('Started: --').classes('text-sm text-gray-600')
        active_jobs_label = ui.label('Active Jobs: 0').classes('text-sm text-gray-600')
        completed_jobs_label = ui.label('Completed Jobs: 0').classes('text-sm text-gray-600')
        failed_jobs_label = ui.label('Failed Jobs: 0').classes('text-sm text-gray-600')
        
        overall_progress = ui.linear_progress().classes('w-full mb-2')
        progress_label = ui.label('0% Complete').classes('text-center text-sm text-gray-600')
        
        preprocessing_progress = ui.label('Preprocessing: 0/0').classes('text-sm text-gray-600')
        analysis_progress = ui.label('Analysis: 0/0').classes('text-sm text-gray-600')
        classification_progress = ui.label('Classification: 0/0').classes('text-sm text-gray-600')
        
        log_area = ui.textarea('Workflow logs will appear here...').classes('w-full h-32').props('readonly')
        
        detected_files_label = ui.label('BAM Files: 0').classes('text-sm text-gray-600')
        processed_files_label = ui.label('Output Files: 0').classes('text-sm text-gray-600')
        failed_files_label = ui.label('Temp Files: 0').classes('text-sm text-gray-600')
        
        workflow_config_label = ui.label('Workflow: --').classes('text-sm text-gray-600')
        log_level_label = ui.label('Log Level: --').classes('text-sm text-gray-600')
        workers_label = ui.label('Analysis Workers: --').classes('text-sm text-gray-600')
        
        current_activity_label = ui.label('No active jobs').classes('text-sm text-gray-600 italic')
        
        with ui.column().classes('w-full h-screen'):
            # Header
            with ui.row().classes('w-full bg-blue-600 text-white p-4'):
                ui.html('<h1 class="text-2xl font-bold">LittleJohn Workflow Monitor</h1>')
            
            # Main content area with splitter
            with ui.row().classes('w-full flex-1'):
                # Left panel with vertical tabs
                with ui.column().classes('w-80 bg-gray-100 p-2'):
                    ui.html('<h2 class="text-lg font-semibold mb-3">Workflow Sections</h2>')
                    
                    # Vertical tabs container
                    with ui.tabs().classes('w-full').props('vertical') as tabs:
                        ui.tab('Status', icon='dashboard')
                        ui.tab('Progress', icon='trending_up')
                        ui.tab('Logs', icon='description')
                        ui.tab('Files', icon='folder')
                        ui.tab('Settings', icon='settings')
                    
                    # Tab panels
                    with ui.tab_panels(tabs, value='Status').classes('w-full'):
                        # Status tab
                        with ui.tab_panel('Status'):
                            status_label
                            started_time_label
                            active_jobs_label
                            completed_jobs_label
                            failed_jobs_label
                        
                        # Progress tab
                        with ui.tab_panel('Progress'):
                            ui.label('Overall Progress').classes('text-lg font-semibold mb-2')
                            overall_progress
                            progress_label
                            
                            ui.label('Queue Progress').classes('text-lg font-semibold mt-4 mb-2')
                            preprocessing_progress
                            analysis_progress
                            classification_progress
                        
                        # Logs tab
                        with ui.tab_panel('Logs'):
                            ui.label('Live Logs').classes('text-lg font-semibold mb-2')
                            log_area
                            
                            with ui.row().classes('w-full justify-center gap-2 mt-2'):
                                ui.button('Clear Logs', on_click=lambda: log_area.set_value('')).classes('bg-gray-500 hover:bg-gray-600 text-white')
                                ui.button('Export Logs', on_click=lambda: None).classes('bg-blue-500 hover:bg-blue-600 text-white')
                        
                        # Files tab
                        with ui.tab_panel('Files'):
                            ui.label('Monitored Files').classes('text-lg font-semibold mb-2')
                            detected_files_label
                            processed_files_label
                            failed_files_label
                        
                        # Settings tab
                        with ui.tab_panel('Settings'):
                            ui.label('GUI Settings').classes('text-lg font-semibold mb-2')
                            ui.switch('Auto-refresh', value=True).classes('mb-2')
                            ui.switch('Show notifications', value=True).classes('mb-2')
                            ui.switch('Dark mode', value=False).classes('mb-2')
                
                # Right panel with detailed information
                with ui.column().classes('flex-1 p-4'):
                    ui.html('<h2 class="text-xl font-semibold mb-4">Workflow Details</h2>')
                    
                    # Workflow configuration
                    with ui.card().classes('w-full mb-4'):
                        ui.html('<h3 class="text-lg font-semibold mb-2">Configuration</h3>')
                        workflow_config_label
                        log_level_label
                        workers_label
                    
                    # Current activity
                    with ui.card().classes('w-full mb-4'):
                        ui.html('<h3 class="text-lg font-semibold mb-2">Current Activity</h3>')
                        current_activity_label
                    
                    # Recent events
                    with ui.card().classes('w-full'):
                        ui.html('<h3 class="text-lg font-semibold mb-2">Recent Events</h3>')
                        ui.label('Workflow started').classes('text-sm text-gray-600')
                        ui.label('Monitoring directory for BAM files').classes('text-sm text-gray-600')
            
            # Footer
            with ui.row().classes('w-full bg-gray-200 p-2 justify-center'):
                ui.label('LittleJohn Workflow Monitor - Running').classes('text-sm text-gray-600')
            
            # Function to update UI with real-time data
            def update_ui():
                try:
                    # Try to get workflow state if available
                    try:
                        from littlejohn.workflow_state import workflow_state
                        summary = workflow_state.get_workflow_summary()
                        
                        # Update status
                        if summary['is_running']:
                            status_label.set_text('Workflow Status: Running')
                            status_label.classes('text-lg font-semibold text-green-600')
                        else:
                            status_label.set_text('Workflow Status: Stopped')
                            status_label.classes('text-lg font-semibold text-red-600')
                        
                        # Update timing
                        if summary['start_time']:
                            start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(summary['start_time']))
                            started_time_label.set_text(f'Started: {start_time}')
                        
                        # Update job counts
                        active_jobs_label.set_text(f'Active Jobs: {summary["running_jobs"]}')
                        completed_jobs_label.set_text(f'Completed Jobs: {summary["completed_jobs"]}')
                        failed_jobs_label.set_text(f'Failed Jobs: {summary["failed_jobs"]}')
                        
                        # Update progress
                        total_jobs = summary['total_jobs']
                        if total_jobs > 0:
                            progress = summary['completed_jobs'] / total_jobs
                            overall_progress.set_value(progress)
                            progress_label.set_text(f'{int(progress * 100)}% Complete')
                        
                        # Update queue progress
                        queue_status = summary['queue_status']
                        preprocessing_progress.set_text(f'Preprocessing: {queue_status["preprocessing"]["running"]}/{queue_status["preprocessing"]["total"]}')
                        analysis_progress.set_text(f'Analysis: {queue_status["mgmt"]["running"] + queue_status["cnv"]["running"] + queue_status["target"]["running"] + queue_status["fusion"]["running"]}/{queue_status["mgmt"]["total"] + queue_status["cnv"]["total"] + queue_status["target"]["total"] + queue_status["fusion"]["total"]}')
                        classification_progress.set_text(f'Classification: {queue_status["classification"]["running"]}/{queue_status["classification"]["total"]}')
                        
                        # Update file counts
                        detected_files_label.set_text(f'BAM Files: {summary["detected_files"]}')
                        processed_files_label.set_text(f'Output Files: {summary["processed_files"]}')
                        failed_files_label.set_text(f'Failed Files: {summary["failed_files"]}')
                        
                        # Update configuration
                        if summary['workflow_steps']:
                            workflow_config_label.set_text(f'Workflow: {", ".join(summary["workflow_steps"])}')
                        
                        # Update recent logs
                        recent_logs = workflow_state.get_recent_logs(20)
                        if recent_logs:
                            log_text = '\n'.join([f'[{time.strftime("%H:%M:%S", time.localtime(log["timestamp"]))}] {log["level"]}: {log["message"]}' for log in recent_logs])
                            log_area.set_value(log_text)
                        
                        # Update current activity
                        running_jobs = workflow_state.get_jobs_by_status(workflow_state.JobStatus.RUNNING)
                        if running_jobs:
                            activity_text = '\n'.join([f'{job.job_type}: {Path(job.filepath).name}' for job in running_jobs[:5]])
                            current_activity_label.set_text(activity_text)
                        else:
                            current_activity_label.set_text('No active jobs')
                            
                    except ImportError:
                        # Workflow state not available, show static info
                        status_label.set_text('Workflow Status: Unknown (hooks not installed)')
                        workflow_config_label.set_text(f'Workflow: {", ".join(self.workflow_steps) if self.workflow_steps else "Unknown"}')
                        detected_files_label.set_text(f'Monitored Directory: {self.monitored_directory}')
                        
                except Exception as e:
                    # If workflow state is not available, show static info
                    try:
                        status_label.set_text('Workflow Status: Running (Limited Info)')
                        workflow_config_label.set_text(f'Workflow: {", ".join(self.workflow_steps) if self.workflow_steps else "Unknown"}')
                        detected_files_label.set_text(f'Monitored Directory: {self.monitored_directory}')
                        
                        # Show that we're in limited mode
                        current_activity_label.set_text('Real-time updates not available - check CLI output for progress')
                        
                        # Add some basic log information
                        try:
                            # Use the correct NiceGUI API for textarea values
                            current_value = log_area.value if hasattr(log_area, 'value') else 'Workflow logs will appear here...'
                            if current_value == 'Workflow logs will appear here...':
                                log_area.set_value('GUI is running in limited mode.\nReal-time workflow updates not available.\nCheck the CLI output for detailed progress information.\n\nThis is normal when workflow hooks are not fully integrated.')
                        except Exception:
                            # If we can't access the textarea value, just set it directly
                            log_area.set_value('GUI is running in limited mode.\nReal-time workflow updates not available.\nCheck the CLI output for detailed progress information.\n\nThis is normal when workflow hooks are not fully integrated.')
                    except Exception as gui_error:
                        # If even the basic GUI updates fail, just log the error
                        print(f"GUI update error: {gui_error}")
                        pass
            
            # Set up auto-refresh timer
            ui.timer(2.0, update_ui)
            
            # Initial update
            update_ui()
    
    def stop_gui(self):
        """Stop the GUI."""
        self.is_running = False
        if self.gui_thread and self.gui_thread.is_alive():
            # Note: NiceGUI doesn't have a clean shutdown method
            # The thread will terminate when the main process ends
            pass
    
    def is_gui_running(self) -> bool:
        """Check if the GUI is running."""
        return self.is_running and self.gui_thread and self.gui_thread.is_alive()
    
    def get_gui_url(self) -> str:
        """Get the URL where the GUI is accessible."""
        return f"http://{self.host}:{self.port}"


def launch_gui(host: str = "localhost", port: int = 8081, show: bool = False, 
               workflow_runner: Any = None, workflow_steps: list = None, 
               monitored_directory: str = "") -> GUILauncher:
    """Launch the LittleJohn workflow GUI."""
    if ui is None:
        raise ImportError("NiceGUI is not available. Please install it with: pip install nicegui")
    
    launcher = GUILauncher(host, port)
    success = launcher.launch_gui(workflow_runner, workflow_steps, monitored_directory)
    
    if not success:
        raise RuntimeError("Failed to launch GUI")
    
    if show:
        import webbrowser
        webbrowser.open(launcher.get_gui_url())
    
    return launcher
