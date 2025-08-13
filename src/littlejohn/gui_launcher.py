"""
GUI launcher for LittleJohn workflow monitoring.

This module provides a clean interface for launching the workflow GUI
that runs in a separate thread but is completely isolated from
workflow execution to avoid any blocking.
"""

import threading
import time
import logging
import queue
from collections import deque
import csv
import pandas as pd
from typing import Optional, Dict, Any, List
from pathlib import Path
from dataclasses import dataclass
from enum import Enum

try:
    from nicegui import ui, app
except ImportError:
    ui = None
    app = None


class UpdateType(Enum):
    """Types of updates that can be sent to the GUI."""
    WORKFLOW_STATUS = "workflow_status"
    JOB_UPDATE = "job_update"
    QUEUE_UPDATE = "queue_update"
    LOG_MESSAGE = "log_message"
    PROGRESS_UPDATE = "progress_update"
    ERROR_UPDATE = "error_update"
    SAMPLES_UPDATE = "samples_update"


@dataclass
class GUIUpdate:
    """A single update message for the GUI."""
    update_type: UpdateType
    timestamp: float
    data: Dict[str, Any]
    priority: int = 0  # Higher priority updates are processed first


class GUILauncher:
    """Launcher for the LittleJohn workflow GUI using isolated threading with message queue."""
    
    def __init__(self, host: str = "localhost", port: int = 8081):
        self.host = host
        self.port = port
        self.gui_thread = None
        self.is_running = False
        self.workflow_runner = None
        self.workflow_steps = []
        self.monitored_directory = ""
        
        # Message queue for non-blocking communication
        self.update_queue = queue.PriorityQueue()
        self.gui_ready = threading.Event()
        self.shutdown_event = threading.Event()
        
        # GUI update thread
        self.update_thread = None  # not used anymore; updates processed on UI thread via timer

        # Debug counters
        self.total_updates_enqueued = 0
        self.total_updates_processed = 0

        # Runtime state
        self._start_time: Optional[float] = None
        self._is_running: bool = False

        # Log ring buffer (last 1000 entries)
        self._log_buffer = deque(maxlen=1000)
        
        # Cached samples data for persistence across navigation
        self._last_samples_rows: List[Dict[str, Any]] = []
        self._current_sample_id: Optional[str] = None
        self._selected_sample_id: Optional[str] = None
        self._known_sample_ids: set[str] = set()
    
    def launch_gui(self, workflow_runner: Any = None, workflow_steps: list = None, 
                   monitored_directory: str = "") -> bool:
        """Launch the GUI in a completely isolated background thread."""
        if ui is None:
            logging.error("NiceGUI is not available")
            return False
            
        self.workflow_runner = workflow_runner
        self.workflow_steps = workflow_steps or []
        self.monitored_directory = monitored_directory
        
        try:
            # Start GUI in completely isolated background thread
            self.gui_thread = threading.Thread(
                target=self._run_gui_worker,
                daemon=True,
                name="LittleJohn-GUI-Thread"
            )
            self.gui_thread.start()
            
            # NOTE: Update processing now happens inside the UI thread via ui.timer for thread-safety
            
            # Wait a moment for GUI to start
            time.sleep(2)
            
            # Check if thread is still running
            if self.gui_thread.is_alive():
                self.is_running = True
                logging.info(f"GUI launched successfully on http://{self.host}:{self.port}")
                return True
            else:
                logging.error("GUI thread failed to start")
                return False
            
        except Exception as e:
            logging.error(f"Failed to launch GUI: {e}")
            return False
    
    def send_update(self, update_type: UpdateType, data: Dict[str, Any], priority: int = 0):
        """Send an update to the GUI without blocking the workflow."""
        try:
            update = GUIUpdate(
                update_type=update_type,
                timestamp=time.time(),
                data=data,
                priority=priority
            )
            
            # Use negative priority so higher priority updates come first
            self.update_queue.put((-priority, update))
            self.total_updates_enqueued += 1
            logging.info(f"[GUI] Enqueued update #{self.total_updates_enqueued}: {update.update_type.value}")
            
        except Exception as e:
            # Don't let update failures affect the workflow
            logging.debug(f"Failed to send GUI update: {e}")
    
    def _drain_updates_on_ui(self):
        """Drain queued updates and apply them on the UI thread (called by ui.timer)."""
        if not self.gui_ready.is_set():
            return
        processed_any = False
        try:
            while True:
                try:
                    _, update = self.update_queue.get_nowait()
                except queue.Empty:
                    break
                self._handle_update(update)
                self.total_updates_processed += 1
                processed_any = True
                logging.info(f"[GUI] Processed update #{self.total_updates_processed}: {update.update_type.value}")
        except Exception as e:
            logging.debug(f"[GUI] Error draining updates on UI: {e}")
        finally:
            if processed_any:
                try:
                    ui.update()
                except Exception:
                    pass
    
    def _handle_update(self, update: GUIUpdate):
        """Handle a single update message."""
        try:
            if update.update_type == UpdateType.WORKFLOW_STATUS:
                self._update_workflow_status(update.data)
            elif update.update_type == UpdateType.JOB_UPDATE:
                self._update_job_status(update.data)
            elif update.update_type == UpdateType.QUEUE_UPDATE:
                self._update_queue_status(update.data)
            elif update.update_type == UpdateType.LOG_MESSAGE:
                self._update_logs(update.data)
            elif update.update_type == UpdateType.PROGRESS_UPDATE:
                self._update_progress(update.data)
            elif update.update_type == UpdateType.ERROR_UPDATE:
                self._update_errors(update.data)
            elif update.update_type == UpdateType.SAMPLES_UPDATE:
                self._update_samples_table(update.data)
                
        except Exception as e:
            logging.debug(f"Error handling GUI update: {e}")
    
    def _update_workflow_status(self, data: Dict[str, Any]):
        """Update workflow status in the GUI."""
        try:
            if hasattr(self, 'status_indicator') and hasattr(self, 'status_label'):
                if data.get('is_running', False):
                    self.status_indicator.classes('text-2xl text-green-400')
                    self.status_label.set_text('Workflow Status: Running')
                    self._is_running = True
                else:
                    self.status_indicator.classes('text-2xl text-red-400')
                    self.status_label.set_text('Workflow Status: Stopped')
                    self._is_running = False
                    
                # Update timing
                if data.get('start_time'):
                    self._start_time = float(data['start_time'])
                    if hasattr(self, 'workflow_start_time'):
                        start_str = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(self._start_time))
                        self.workflow_start_time.set_text(f'Started: {start_str}')
                    
                if hasattr(self, 'workflow_duration') and self._start_time:
                    elapsed_seconds = int(time.time() - self._start_time)
                    self.workflow_duration.set_text(f'Duration: {self._format_duration(elapsed_seconds)}')
                    
        except Exception as e:
            logging.debug(f"Error updating workflow status: {e}")
    
    def _update_job_status(self, data: Dict[str, Any]):
        """Update job status in the GUI."""
        try:
            if hasattr(self, 'active_jobs_table'):
                # Update active jobs table
                job_rows = []
                for job in data.get('active_jobs', []):
                    job_rows.append({
                        'job_id': str(job.get('job_id', ''))[:8],
                        'job_type': job.get('job_type', ''),
                        'filepath': Path(job.get('filepath', '')).name,
                        'worker': job.get('worker_name', ''),
                        'duration': self._format_duration(job.get('duration', 0)),
                        'progress': f"{int(job.get('progress', 0) * 100)}%"
                    })
                
                # Update table
                self.active_jobs_table.clear()
                if job_rows:
                    self.active_jobs_table.add_rows(job_rows)
                    
        except Exception as e:
            logging.debug(f"Error updating job status: {e}")
    
    def _update_queue_status(self, data: Dict[str, Any]):
        """Update queue status in the GUI."""
        try:
            # Update queue status displays
            if hasattr(self, 'preprocessing_status'):
                queue_data = data.get('preprocessing', {})
                self.preprocessing_status.set_text(f"{queue_data.get('running', 0)}/{queue_data.get('total', 0)}")
                
            if hasattr(self, 'analysis_status'):
                # Combine analysis queues
                analysis_running = sum(data.get(q, {}).get('running', 0) for q in ['mgmt', 'cnv', 'target', 'fusion'])
                analysis_total = sum(data.get(q, {}).get('total', 0) for q in ['mgmt', 'cnv', 'target', 'fusion'])
                self.analysis_status.set_text(f"{analysis_running}/{analysis_total}")
                
            if hasattr(self, 'classification_status'):
                classification_data = data.get('classification', {})
                self.classification_status.set_text(f"{classification_data.get('running', 0)}/{classification_data.get('total', 0)}")
                
            if hasattr(self, 'other_status'):
                other_data = data.get('other', {})
                self.other_status.set_text(f"{other_data.get('running', 0)}/{other_data.get('total', 0)}")
                
        except Exception as e:
            logging.debug(f"Error updating queue status: {e}")
    
    def _update_logs(self, data: Dict[str, Any]):
        """Update logs in the GUI."""
        try:
            if hasattr(self, 'log_area'):
                log_message = data.get('message', '')
                log_level = data.get('level', 'INFO')
                timestamp = time.strftime('%H:%M:%S')
                
                new_line = f'[{timestamp}] {log_level}: {log_message}\n'
                self._log_buffer.append(new_line)
                self.log_area.set_value(''.join(self._log_buffer))
                
        except Exception as e:
            logging.debug(f"Error updating logs: {e}")
    
    def _update_progress(self, data: Dict[str, Any]):
        """Update progress in the GUI."""
        try:
            if hasattr(self, 'progress_bar') and hasattr(self, 'progress_label'):
                progress = data.get('progress', 0.0)
                
                pct = max(0.0, min(100.0, round(progress * 100.0, 1)))
                
                self.progress_bar.set_value(pct)
                # Drop .0 for integers like 81.0 -> 81
                pct_str = f"{pct:.1f}" if pct % 1 else f"{int(pct)}"
                self.progress_label.set_text(f'{pct_str}% Complete')
                # Optional counts
                if hasattr(self, 'completed_count') and 'completed' in data:
                    self.completed_count.set_text(str(data['completed']))
                if hasattr(self, 'failed_count') and 'failed' in data:
                    self.failed_count.set_text(str(data['failed']))
                if hasattr(self, 'total_count') and 'total' in data:
                    self.total_count.set_text(str(data['total']))
                
        except Exception as e:
            logging.debug(f"Error updating progress: {e}")
    
    def _update_errors(self, data: Dict[str, Any]):
        """Update error information in the GUI."""
        try:
            # Update error counts if available
            if hasattr(self, 'preprocessing_errors'):
                self.preprocessing_errors.set_text(str(data.get('preprocessing_errors', 0)))
                
            if hasattr(self, 'analysis_errors'):
                self.analysis_errors.set_text(str(data.get('analysis_errors', 0)))
                
            if hasattr(self, 'classification_errors'):
                self.classification_errors.set_text(str(data.get('classification_errors', 0)))
                
        except Exception as e:
            logging.debug(f"Error updating error information: {e}")
    
    def _format_duration(self, seconds):
        """Format duration in seconds to human readable format."""
        if seconds < 60:
            return f"{seconds}s"
        elif seconds < 3600:
            minutes = seconds // 60
            return f"{minutes}m {seconds % 60}s"
        else:
            hours = seconds // 3600
            minutes = (seconds % 3600) // 60
            return f"{hours}h {minutes}m"
    
    def _run_gui_worker(self):
        """Run the GUI in a completely isolated thread."""
        try:
            # Set thread name for identification
            threading.current_thread().name = "LittleJohn-GUI-Thread"
            
            # Create the main workflow monitor page
            @ui.page('/')
            def welcome_page():
                """Welcome page at root route."""
                self._create_welcome_page()
            
            # Create the workflow monitoring page
            @ui.page('/littlejohn')
            def workflow_monitor():
                """Workflow monitoring page under /littlejohn route."""
                self._create_workflow_monitor()
            
            # Create the samples overview page
            @ui.page('/live_data')
            def samples_overview():
                """Samples overview page showing all tracked samples."""
                self._create_samples_overview()
            
            # Create individual sample detail pages
            @ui.page('/live_data/{sample_id}')
            def sample_detail(sample_id: str):
                """Individual sample detail page."""
                self._create_sample_detail_page(sample_id)
            
            # Enable global update processing regardless of which page is open
            try:
                self.gui_ready.set()
                ui.timer(0.3, self._drain_updates_on_ui, active=True)
            except Exception:
                pass
            
            # Start the GUI
            ui.run(
                host=self.host,
                port=self.port,
                show=False,
                reload=False
            )
        except Exception as e:
            print(f"❌ GUI worker error: {e}")
            import traceback
            traceback.print_exc()
    
    def _create_welcome_page(self):
        """Create the welcome page."""
        # Background and main container
        with ui.column().classes('w-full h-full items-center justify-center bg-gradient-to-br from-blue-50 to-indigo-100 p-8'):
            # Main title and description
            ui.label('🧬 LittleJohn').classes('text-6xl font-bold mb-4 text-blue-700')
            ui.label('Advanced BAM File Analysis & Workflow Management').classes('text-2xl text-gray-700 mb-8 text-center max-w-3xl')
            
            # Description card
            with ui.card().classes('w-full max-w-4xl mb-8 bg-white shadow-lg'):
                with ui.column().classes('p-6'):
                    ui.label('What is LittleJohn?').classes('text-xl font-semibold mb-4 text-blue-800')
                    ui.label('LittleJohn is a comprehensive bioinformatics workflow system designed for processing and analyzing BAM files. It provides automated preprocessing, multiple analysis pipelines, and real-time monitoring capabilities.').classes('text-gray-700 mb-4')
                    
                    with ui.row().classes('w-full justify-center gap-8 mt-6'):
                        with ui.column().classes('text-center'):
                            ui.label('🔬').classes('text-3xl mb-2')
                            ui.label('Preprocessing').classes('text-sm font-medium text-gray-600')
                        with ui.column().classes('text-center'):
                            ui.label('🧬').classes('text-3xl mb-2')
                            ui.label('MGMT Analysis').classes('text-sm font-medium text-gray-600')
                        with ui.column().classes('text-center'):
                            ui.label('📊').classes('text-3xl mb-2')
                            ui.label('CNV Detection').classes('text-sm font-medium text-gray-600')
                        with ui.column().classes('text-center'):
                            ui.label('🎯').classes('text-3xl mb-2')
                            ui.label('Target Analysis').classes('text-sm font-medium text-gray-600')
                        with ui.column().classes('text-center'):
                            ui.label('🔗').classes('text-3xl mb-2')
                            ui.label('Fusion Detection').classes('text-sm font-medium text-gray-600')
            
            # Action buttons
            with ui.row().classes('gap-6'):
                ui.link('📊 Open Workflow Monitor', '/littlejohn').classes('bg-blue-600 hover:bg-blue-700 text-white text-lg font-semibold px-8 py-4 rounded-lg shadow-lg transition-colors')
                
                # Placeholder buttons for future functionality
                ui.button('🚀 Launch New Workflow', on_click=lambda: self._launch_workflow_button_clicked()).classes('bg-green-600 hover:bg-green-700 text-white text-lg font-semibold px-8 py-4 rounded-lg shadow-lg transition-colors')
                ui.button('📋 View Documentation', on_click=lambda: self._view_docs_button_clicked()).classes('bg-purple-600 hover:bg-purple-700 text-white text-lg font-semibold px-8 py-4 rounded-lg shadow-lg transition-colors')
            
            # Sample Tracking Preview - SIMPLIFIED: Show basic info without workflow_state access
            with ui.card().classes('w-full max-w-4xl mt-8 bg-white shadow-lg'):
                with ui.column().classes('p-6'):
                    ui.label('🧬 Live Sample Tracking').classes('text-xl font-semibold mb-4 text-blue-800')
                    ui.label('Sample tracking will be available once the workflow is running.').classes('text-sm text-gray-600')
                    
                    with ui.row().classes('w-full justify-center gap-4 mt-4'):
                        ui.link('�� View All Samples', '/live_data').classes('bg-blue-600 hover:bg-blue-700 text-white text-sm px-4 py-2 rounded-lg')
                        ui.link('📊 Workflow Monitor', '/littlejohn').classes('bg-green-600 hover:bg-green-700 text-white text-sm px-4 py-2 rounded-lg')
            
            # Footer information
            with ui.row().classes('mt-12 text-center text-gray-500'):
                ui.label('LittleJohn Workflow Monitor - Advanced Bioinformatics Analysis Platform').classes('text-sm')
    
    def _create_samples_overview(self):
        """Create the samples overview page showing all tracked samples."""
        # Page title and navigation
        with ui.row().classes('w-full bg-blue-600 text-white p-4 items-center justify-between'):
            with ui.row().classes('items-center'):
                ui.label('🧬 Sample Tracking Overview').classes('text-2xl font-bold')
                ui.label('All samples processed by LittleJohn').classes('text-sm ml-4 opacity-80')
            
            # Navigation links
            with ui.row().classes('gap-4'):
                ui.link('🏠 Welcome', '/').classes('text-white hover:text-blue-200 text-sm')
                ui.link('📊 Workflow Monitor', '/littlejohn').classes('text-white hover:text-blue-200 text-sm')
                ui.label('📋 Sample Overview').classes('text-white text-sm font-semibold')
        
        # Main content area
        with ui.column().classes('w-full p-4 gap-4'):
            # Sample statistics
            with ui.card().classes('w-full bg-gradient-to-r from-blue-50 to-indigo-50'):
                ui.label('📊 Sample Statistics').classes('text-lg font-semibold mb-4 text-blue-800')
                
                # SIMPLIFIED: Show basic info without workflow_state access
                ui.label('Sample statistics will be available once the workflow is running.').classes('text-sm text-gray-600')
            
            # Samples table
            with ui.card().classes('w-full'):
                ui.label('📋 All Tracked Samples').classes('text-lg font-semibold mb-4')
                with ui.row().classes('items-center gap-2 mb-2'):
                    self.view_sample_button = ui.button(
                        'View',
                        on_click=lambda: ui.navigate.to(f"/live_data/{self._selected_sample_id}") if self._selected_sample_id else ui.notify('Select a sample first', type='warning')
                    ).props('color=primary')
                    self.view_sample_button.disable()

                # Create a placeholder table that will be updated later
                self.samples_table = ui.table(
                    columns=[
                        {'name': 'sample_id', 'label': 'Sample ID', 'field': 'sample_id'},
                        {'name': 'active_jobs', 'label': 'Active', 'field': 'active_jobs'},
                        {'name': 'total_jobs', 'label': 'Total', 'field': 'total_jobs'},
                        {'name': 'completed_jobs', 'label': 'Completed', 'field': 'completed_jobs'},
                        {'name': 'failed_jobs', 'label': 'Failed', 'field': 'failed_jobs'},
                        {'name': 'job_types', 'label': 'Job Types', 'field': 'job_types'},
                        {'name': 'last_seen', 'label': 'Last Activity', 'field': 'last_seen'}
                    ],
                    rows=[],
                    row_key='sample_id',
                    selection='single',
                    pagination=20
                ).classes('w-full')

                # Selection handler to enable the external View button
                try:
                    self.samples_table.on('selection', self._on_sample_selected)
                except Exception:
                    pass

                # If we have cached rows, show them immediately for instant UX
                if self._last_samples_rows:
                    try:
                        self.samples_table.rows = list(self._last_samples_rows)
                        self.samples_table.update()
                        # Auto-select if only one sample
                        if len(self._last_samples_rows) == 1:
                            self._selected_sample_id = self._last_samples_rows[0].get('sample_id')
                            self.view_sample_button.enable()
                    except Exception:
                        pass

    def _update_samples_table(self, data: Dict[str, Any]):
        """Update the samples overview table with new data."""
        try:
            if not hasattr(self, 'samples_table'):
                return
            samples = data.get('samples', [])
            # Deduplicate by sample_id taking the newest last_seen
            by_id: Dict[str, Dict[str, Any]] = {}
            for s in samples:
                sid = s.get('sample_id', '') or 'unknown'
                last_seen = float(s.get('last_seen', time.time()))
                existing = by_id.get(sid)
                if not existing or last_seen >= existing.get('_last_seen_raw', 0):
                    by_id[sid] = {
                        'sample_id': sid,
                        'active_jobs': s.get('active_jobs', 0),
                        'total_jobs': s.get('total_jobs', 0),
                        'completed_jobs': s.get('completed_jobs', 0),
                        'failed_jobs': s.get('failed_jobs', 0),
                        'job_types': ','.join(sorted(set(s.get('job_types', [])))) if isinstance(s.get('job_types', []), list) else str(s.get('job_types', '')),
                        'last_seen': time.strftime('%H:%M:%S', time.localtime(last_seen)),
                        'actions': 'View',
                        '_last_seen_raw': last_seen,
                    }

            rows = list(by_id.values())
            # Replace rows to avoid duplicates
            self.samples_table.rows = rows
            self.samples_table.update()
            # Cache for persistence
            self._last_samples_rows = rows
            self._known_sample_ids = {r.get('sample_id') for r in rows if r.get('sample_id')}
            # Track currently most active/recent sample
            if rows:
                rows_sorted = sorted(rows, key=lambda r: (r.get('active_jobs', 0), r.get('_last_seen_raw', 0)), reverse=True)
                self._current_sample_id = rows_sorted[0].get('sample_id')
            # Update external button state
            if self._selected_sample_id and any(r.get('sample_id') == self._selected_sample_id for r in rows):
                self.view_sample_button.enable()
            elif rows:
                if len(rows) == 1:
                    self._selected_sample_id = rows[0].get('sample_id')
                    self.view_sample_button.enable()
                else:
                    self.view_sample_button.disable()
            else:
                self._selected_sample_id = None
                self.view_sample_button.disable()
        except Exception as e:
            logging.debug(f"Error updating samples table: {e}")

    def _on_sample_selected(self, event) -> None:
        try:
            # NiceGUI passes {'rows': [selected_rows...]}
            rows = None
            if hasattr(event, 'args') and isinstance(event.args, dict):
                rows = event.args.get('rows')
            elif isinstance(event, dict):
                rows = event.get('rows')
            if rows and isinstance(rows, list) and len(rows) > 0 and isinstance(rows[0], dict):
                self._selected_sample_id = rows[0].get('sample_id')
                if self._selected_sample_id:
                    self.view_sample_button.enable()
                else:
                    self.view_sample_button.disable()
        except Exception:
            pass
    
    def _create_sample_detail_page(self, sample_id: str):
        """Create the individual sample detail page."""
        # Guard: unknown sample -> show message and back button; also redirect
        if self._known_sample_ids and sample_id not in self._known_sample_ids:
            with ui.column().classes('w-full items-center justify-center p-8'):
                ui.label(f'Unknown sample: {sample_id}').classes('text-xl font-semibold text-red-600')
                ui.label('This sample ID has not been seen yet in the current session.').classes('text-sm text-gray-600')
                ui.button('Back to Samples', on_click=lambda: ui.navigate.to('/live_data')).props('color=primary')
            # Soft redirect after short delay
            try:
                ui.timer(1.5, lambda: ui.navigate.to('/live_data'), once=True)
            except Exception:
                pass
            return

        sample_dir = Path(self.monitored_directory) / sample_id if self.monitored_directory else None

        # Page title and navigation
        with ui.row().classes('w-full bg-blue-600 text-white p-4 items-center justify-between'):
            with ui.row().classes('items-center'):
                ui.label(f'🧬 Sample: {sample_id}').classes('text-2xl font-bold')
                ui.label('Detailed sample information and job history').classes('text-sm ml-4 opacity-80')
            with ui.row().classes('gap-4'):
                ui.link('📋 Sample Overview', '/live_data').classes('text-white hover:text-blue-200 text-sm')
                ui.link('📊 Workflow Monitor', '/littlejohn').classes('text-white hover:text-blue-200 text-sm')

        with ui.column().classes('w-full p-4 gap-4'):
            # Files in output directory
            with ui.card().classes('w-full'):
                ui.label('📁 Output Files').classes('text-lg font-semibold mb-2')
                files_table = ui.table(
                    columns=[
                        {'name': 'name', 'label': 'File', 'field': 'name'},
                        {'name': 'size', 'label': 'Size (bytes)', 'field': 'size'},
                        {'name': 'mtime', 'label': 'Last Modified', 'field': 'mtime'},
                    ],
                    rows=[],
                    pagination=20
                ).classes('w-full')

            # master.csv summary
            with ui.card().classes('w-full'):
                ui.label('📊 master.csv Summary').classes('text-lg font-semibold mb-2')
                summary_table = ui.table(
                    columns=[
                        {'name': 'key', 'label': 'Field', 'field': 'key'},
                        {'name': 'value', 'label': 'Value', 'field': 'value'},
                    ],
                    rows=[],
                    pagination=0
                ).classes('w-full')

            # Classification tab: per-tool charts (only shown if files exist)
            with ui.card().classes('w-full'):
                ui.label('🧪 Classification').classes('text-lg font-semibold mb-2')
                # Configure per-tool CSV and scaling mode
                # mode: 'fraction' -> values in 0-1; 'percent' -> values in 0-100
                tool_to_file = {
                    'Sturgeon': {'file': 'sturgeon_scores.csv', 'mode': 'fraction'},
                    'NanoDX': {'file': 'NanoDX_scores.csv', 'mode': 'fraction'},
                    'PanNanoDX': {'file': 'PanNanoDX_scores.csv', 'mode': 'fraction'},
                    'Random Forest': {'file': 'random_forest_scores.csv', 'mode': 'percent'},
                }
                charts: Dict[str, Dict[str, Any]] = {}
                for tool_name, cfg in tool_to_file.items():
                    exp = ui.expansion(tool_name, icon='analytics').classes('w-full')
                    with exp:
                        summary_labels = None
                        if tool_name == 'Sturgeon':
                            with ui.card().classes('w-full bg-gradient-to-r from-blue-50 to-indigo-50 mb-2 p-2'):
                                with ui.row().classes('gap-6 items-center'):
                                    st_class = ui.label('Sturgeon classification: Unknown').classes('text-sm font-semibold text-blue-800')
                                    st_conf = ui.label('Confidence: --%').classes('text-sm text-gray-700')
                                    st_probes = ui.label('Probes: --').classes('text-sm text-gray-700')
                                summary_labels = {'class': st_class, 'conf': st_conf, 'probes': st_probes}
                        elif tool_name in ('NanoDX', 'PanNanoDX'):
                            with ui.card().classes('w-full bg-gradient-to-r from-blue-50 to-indigo-50 mb-2 p-2'):
                                with ui.row().classes('gap-6 items-center'):
                                    ndx_class = ui.label(f'{tool_name} classification: Unknown').classes('text-sm font-semibold text-blue-800')
                                    ndx_conf = ui.label('Confidence: --%').classes('text-sm text-gray-700')
                                    ndx_feats = ui.label('Probes: --').classes('text-sm text-gray-700')
                                summary_labels = {'class': ndx_class, 'conf': ndx_conf, 'probes': ndx_feats}
                        elif tool_name == 'Random Forest':
                            with ui.card().classes('w-full bg-gradient-to-r from-blue-50 to-indigo-50 mb-2 p-2'):
                                with ui.row().classes('gap-6 items-center'):
                                    rf_class = ui.label('Forest classification: Unknown').classes('text-sm font-semibold text-blue-800')
                                    rf_conf = ui.label('Confidence: --%').classes('text-sm text-gray-700')
                                    rf_feats = ui.label('Features: --').classes('text-sm text-gray-700')
                                summary_labels = {'class': rf_class, 'conf': rf_conf, 'probes': rf_feats}
                        ui.label(f'{tool_name} current classification').classes('text-sm text-gray-700')
                        # Bar chart with richer style
                        bar = ui.echart({
                            'backgroundColor': 'transparent',
                            'title': {'text': f'{tool_name} (Top classes)', 'left': 'center', 'top': 10, 'textStyle': {'fontSize': 16, 'color': '#000'}},
                            'tooltip': {'trigger': 'axis', 'axisPointer': {'type': 'shadow'}, 'formatter': '{b}: {c}%'},
                            'grid': {'left': '5%', 'right': '5%', 'bottom': '5%', 'top': '25%', 'containLabel': True},
                            'xAxis': {'type': 'value', 'min': 0, 'max': 100, 'interval': 20, 'axisLabel': {'formatter': '{value}%'}},
                            'yAxis': {'type': 'category', 'inverse': True, 'data': []},
                            'series': [{
                                'type': 'bar', 'data': [], 'barMaxWidth': '60%',
                                'itemStyle': {'color': '#007AFF', 'borderRadius': [0, 4, 4, 0]},
                                'label': {'show': True, 'position': 'right', 'formatter': '{c}%'}
                            }],
                        }).classes('w-full h-60')
                        ui.label(f'{tool_name} confidence over time').classes('text-sm text-gray-700 mt-2')
                        ts = ui.echart({
                            'backgroundColor': 'transparent',
                            'title': {'text': f'{tool_name} (time series)', 'left': 'center', 'top': 5, 'textStyle': {'fontSize': 16, 'color': '#000'}},
                            'tooltip': {'trigger': 'axis'},
                            'legend': {'type': 'scroll', 'top': 45},
                            'grid': {'left': '5%', 'right': '5%', 'bottom': '5%', 'top': '30%', 'containLabel': True},
                            'xAxis': {'type': 'time'},
                            'yAxis': {
                                'type': 'value', 'min': 0, 'max': 100,
                                'axisLabel': {'formatter': '{value}%'},
                                'splitLine': {'show': True, 'lineStyle': {'type': 'dashed', 'color': '#E0E0E0'}},
                            },
                            'series': [],
                        }).classes('w-full h-64')
                        charts[tool_name] = {
                            'bar': bar,
                            'ts': ts,
                            'file': cfg['file'],
                            'mode': cfg['mode'],
                            'summary': summary_labels,
                            'expansion': exp,
                            'last_mtime': None,
                        }

            # Periodic refresher
            def _read_scores_csv(csv_path: Path, mode: str):
                try:
                    with csv_path.open('r', newline='') as fh:
                        reader = csv.DictReader(fh)
                        rows = list(reader)
                        if not rows:
                            return None
                        # Determine which columns are numeric scores
                        numeric_keys = []
                        for k in reader.fieldnames or []:
                            if k.lower() in {'timestamp', 'number_probes'}:
                                continue
                            try:
                                float(rows[-1].get(k, ''))
                                numeric_keys.append(k)
                            except Exception:
                                continue
                        if not numeric_keys:
                            return None
                        # Last row scores
                        last_scores = {}
                        for k in numeric_keys:
                            try:
                                v = float(rows[-1].get(k, 0))
                                # Scale based on mode
                                if mode == 'fraction':
                                    v *= 100.0
                                last_scores[k] = round(v, 2)
                            except Exception:
                                pass
                        # Time series (use index or timestamp if present)
                        time_key = 'timestamp' if 'timestamp' in (reader.fieldnames or []) else None
                        x_labels = []
                        series_map: Dict[str, List[List[Any]]] = {k: [] for k in numeric_keys}
                        for idx, r in enumerate(rows):
                            if time_key:
                                raw_x = r.get(time_key)
                                try:
                                    x_val = float(raw_x)
                                except Exception:
                                    x_val = raw_x
                            else:
                                x_val = idx + 1
                            x_labels.append(x_val)
                            for k in numeric_keys:
                                try:
                                    vv = float(r.get(k, 0))
                                    if mode == 'fraction':
                                        vv *= 100.0
                                    series_map[k].append([x_val, round(vv, 2)])
                                except Exception:
                                    series_map[k].append([x_val, 0])
                        # Ensure chronological order when timestamp column is present
                        if time_key:
                            def _sort_key(p):
                                try:
                                    return float(p[0])
                                except Exception:
                                    return str(p[0])
                            for k in series_map:
                                series_map[k] = sorted(series_map[k], key=_sort_key)
                        return {'last': last_scores, 'x': x_labels, 'series': series_map, 'has_time': bool(time_key)}
                except Exception:
                    return None

            def _update_charts_from_file(tool_name: str, file_name: str):
                try:
                    file_path = sample_dir / file_name if sample_dir else None
                    if not file_path or not file_path.exists():
                        return
                    # Special handling for Random Forest CSVs to match previous implementation
                    if tool_name == 'Random Forest':
                        try:
                            df = pd.read_csv(file_path, index_col=0)
                        except Exception:
                            return
                        if df.empty:
                            return
                        # Determine number_probes and drop for plotting
                        lastrow = df.iloc[-1]
                        n_features = lastrow.get('number_probes', 0)
                        if 'number_probes' in df.columns:
                            df_plot = df.drop(columns=['number_probes'])
                        else:
                            df_plot = df.copy()
                        # Filter columns where any value > 0.5 (percentage units) over time
                        cols_keep = [c for c in df_plot.columns if (df_plot[c] > 0.5).any()]
                        if cols_keep:
                            df_plot = df_plot[cols_keep]
                        # Bar: top 10 from latest row
                        latest = df_plot.iloc[-1]
                        top_series = latest.sort_values(ascending=False).head(10)
                        bar = charts[tool_name]['bar']
                        bar.options['yAxis']['data'] = list(top_series.index[::-1])
                        bar.options['series'][0]['data'] = [round(float(v), 2) for v in list(top_series.values[::-1])]
                        bar.update()
                        # Time series: use top 10 diagnoses
                        ts = charts[tool_name]['ts']
                        top10 = latest.nlargest(10).index.tolist()
                        filtered_df = df_plot[top10]
                        # Build time-axis series: [timestamp_ms, value]
                        raw_index = list(filtered_df.index.tolist())
                        x_vals = list(raw_index)

                        ts.options['series'] = []
                        for col in filtered_df.columns:
                            y_vals = [round(float(v), 2) for v in filtered_df[col].tolist()]
                            series_data = [[x_vals[i], y_vals[i]] for i in range(len(y_vals))]
                            ts.options['series'].append({
                                'name': col, 'type': 'line', 'smooth': True, 'animation': False,
                                'data': series_data,
                            })
                        # Threshold lines 65/85
                        ts.options['series'].append({'name': 'thresholds', 'type': 'line', 'data': [], 'markLine': {
                            'silent': True, 'lineStyle': {'type': 'dashed', 'color': '#999'}, 'data': [{'yAxis': 65}, {'yAxis': 85}],
                        }})
                        ts.update()
                        # Summary
                        if charts[tool_name].get('summary'):
                            labels_map = charts[tool_name]['summary']
                            if labels_map and len(top_series) > 0:
                                best_label = top_series.index[0]
                                best_value = float(top_series.values[0])
                                labels_map['class'].set_text(f'Forest classification: {best_label}')
                                labels_map['conf'].set_text(f'Confidence: {best_value:.2f}%')
                                if n_features is not None:
                                    try:
                                        labels_map['probes'].set_text(f'Features: {int(float(n_features))}')
                                    except Exception:
                                        pass
                        return

                    mode = charts[tool_name].get('mode', 'percent')
                    data = _read_scores_csv(file_path, mode)
                    if not data:
                        return
                    bar = charts[tool_name]['bar']
                    ts = charts[tool_name]['ts']
                    # Update bar chart with top 10
                    last = data['last']
                    top = sorted(last.items(), key=lambda kv: kv[1], reverse=True)[:10]
                    labels = [k for k, _ in top][::-1]
                    values = [round(v, 2) for _, v in top][::-1]
                    bar.options['yAxis']['data'] = labels
                    bar.options['series'][0]['data'] = values
                    bar.update()
                    # Update time series (limit series to top 5 for readability)
                    ts.options['series'] = []
                    if data.get('has_time'):
                        ts.options['xAxis'] = {'type': 'time'}
                        for k, _ in top[:5]:
                            ts.options['series'].append({
                                'name': k, 'type': 'line', 'smooth': True, 'animation': False,
                                'data': data['series'][k],
                            })
                    else:
                        ts.options['xAxis'] = {'type': 'category', 'data': data['x']}
                        for k, _ in top[:5]:
                            ts.options['series'].append({
                                'name': k, 'type': 'line', 'smooth': True, 'animation': False,
                                'data': [y for _, y in data['series'][k]],
                            })
                    # Add threshold lines for Sturgeon and Random Forest
                    if tool_name == 'Sturgeon':
                        ts.options.setdefault('series', [])
                        ts.options['yAxis'].setdefault('axisLabel', {'formatter': '{value}%'} )
                        ts.options.setdefault('grid', {'containLabel': True})
                        ts.options.setdefault('visualMap', None)
                        ts.options.setdefault('markLine', None)
                        # Use a markLine via a fake series overlay
                        ts.options['series'].append({
                            'name': 'thresholds', 'type': 'line', 'data': [], 'markLine': {
                                'silent': True,
                                'lineStyle': {'type': 'dashed', 'color': '#999'},
                                'data': [ {'yAxis': 60}, {'yAxis': 80} ],
                            }
                        })
                    if tool_name == 'Random Forest':
                        ts.options.setdefault('series', [])
                        ts.options['series'].append({
                            'name': 'thresholds', 'type': 'line', 'data': [], 'markLine': {
                                'silent': True,
                                'lineStyle': {'type': 'dashed', 'color': '#999'},
                                'data': [ {'yAxis': 65}, {'yAxis': 85} ],
                            }
                        })
                    ts.update()
                    # Update summary card for Sturgeon
                    if tool_name == 'Sturgeon' and charts[tool_name].get('summary'):
                        labels_map = charts[tool_name]['summary']
                        if labels_map:
                            # Highest class and value
                            if top:
                                best_label, best_value = top[0]
                                labels_map['class'].set_text(f'Sturgeon classification: {best_label}')
                                labels_map['conf'].set_text(f'Confidence: {best_value:.2f}%')
                            # number_probes if present
                            try:
                                # 'number_probes' column (raw, not multiplied)
                                file_path = sample_dir / file_name if sample_dir else None
                                with file_path.open('r', newline='') as fh:
                                    reader = csv.DictReader(fh)
                                    rows = list(reader)
                                if rows:
                                    npv = rows[-1].get('number_probes') or rows[-1].get('number_probes'.lower())
                                    if npv is not None:
                                        labels_map['probes'].set_text(f'Probes: {int(float(npv))}')
                            except Exception:
                                pass
                        # Also update expansion header to show summary inline
                        try:
                            exp = charts[tool_name].get('expansion')
                            if exp and top:
                                best_label, best_value = top[0]
                                header = f"{tool_name} — {best_label} ({round(best_value, 2)}%)"
                                exp.props(f'label="{header}"')
                        except Exception:
                            pass

                    # Update summary card for NanoDX / PanNanoDX
                    if tool_name in ('NanoDX', 'PanNanoDX') and charts[tool_name].get('summary'):
                        labels_map = charts[tool_name]['summary']
                        if labels_map:
                            if top:
                                best_label, best_value = top[0]
                                labels_map['class'].set_text(f'{tool_name} classification: {best_label}')
                                labels_map['conf'].set_text(f'Confidence: {best_value:.2f}%')
                            # number_probes if present
                            try:
                                file_path = sample_dir / file_name if sample_dir else None
                                with file_path.open('r', newline='') as fh:
                                    reader = csv.DictReader(fh)
                                    rows = list(reader)
                                if rows:
                                    npv = rows[-1].get('number_probes') or rows[-1].get('number_probes'.lower())
                                    if npv is not None:
                                        labels_map['probes'].set_text(f'Probes: {int(float(npv))}')
                            except Exception:
                                pass
                        # Inline header
                        try:
                            exp = charts[tool_name].get('expansion')
                            if exp and top:
                                best_label, best_value = top[0]
                                header = f"{tool_name} — {best_label} ({round(best_value, 2)}%)"
                                exp.props(f'label="{header}"')
                        except Exception:
                            pass
                    # Update summary for Random Forest
                    if tool_name == 'Random Forest' and charts[tool_name].get('summary'):
                        labels_map = charts[tool_name]['summary']
                        if labels_map:
                            if top:
                                best_label, best_value = top[0]
                                labels_map['class'].set_text(f'Forest classification: {best_label}')
                                labels_map['conf'].set_text(f'Confidence: {best_value:.2f}%')
                            try:
                                with file_path.open('r', newline='') as fh:
                                    reader = csv.DictReader(fh)
                                    rows = list(reader)
                                if rows:
                                    npv = rows[-1].get('number_probes') or rows[-1].get('number_probes'.lower())
                                    if npv is not None:
                                        labels_map['probes'].set_text(f'Features: {int(float(npv))}')
                            except Exception:
                                pass
                        # Update expansion header summary for RF as well
                        try:
                            exp = charts[tool_name].get('expansion')
                            if exp and top:
                                best_label, best_value = top[0]
                                header = f"{tool_name} — {best_label} ({round(best_value, 2)}%)"
                                exp.props(f'label="{header}"')
                        except Exception:
                            pass
                except Exception:
                    pass

            def _refresh_sample_detail() -> None:
                # Refresh files list
                try:
                    rows = []
                    if sample_dir and sample_dir.exists():
                        for f in sorted(sample_dir.iterdir()):
                            if f.is_file():
                                try:
                                    stat = f.stat()
                                    rows.append({
                                        'name': f.name,
                                        'size': stat.st_size,
                                        'mtime': time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(stat.st_mtime)),
                                    })
                                except Exception:
                                    continue
                    files_table.rows = rows
                    files_table.update()
                except Exception:
                    pass

                # Refresh master.csv summary
                try:
                    if sample_dir:
                        csv_path = sample_dir / 'master.csv'
                        if csv_path.exists():
                            # Read first row of CSV
                            with csv_path.open('r', newline='') as fh:
                                reader = csv.DictReader(fh)
                                first_row = next(reader, None)
                            if first_row:
                                # Show a compact selection first; fall back to all
                                preferred_keys = [
                                    'counter_bam_passed', 'counter_bam_failed', 'counter_bases_count',
                                    'counter_mapped_count', 'counter_unmapped_count',
                                    'run_info_run_time', 'run_info_device', 'run_info_model', 'run_info_flow_cell',
                                    'bam_tracking_counter', 'bam_tracking_total_files',
                                ]
                                rows = []
                                for k in preferred_keys:
                                    if k in first_row:
                                        rows.append({'key': k, 'value': first_row.get(k, '')})
                                # Add a few more dynamic fields if present
                                additional = [
                                    ('devices', 'devices'),
                                    ('basecall_models', 'basecall_models'),
                                    ('flowcell_ids', 'flowcell_ids'),
                                ]
                                for k, _ in additional:
                                    if k in first_row:
                                        rows.append({'key': k, 'value': first_row.get(k, '')})
                                if not rows:
                                    rows = [{'key': k, 'value': v} for k, v in first_row.items()]
                                summary_table.rows = rows
                                summary_table.update()
                        else:
                            summary_table.rows = [{'key': 'Status', 'value': 'master.csv not found'}]
                            summary_table.update()
                except Exception:
                    # Avoid breaking the UI; skip on CSV parse errors
                    pass

                # Refresh classification charts only when source file changed (or first load)
                for tool_name, info in charts.items():
                    file_path = sample_dir / info['file'] if sample_dir else None
                    try:
                        if file_path and file_path.exists():
                            mtime = file_path.stat().st_mtime
                            if info.get('last_mtime') is None or mtime > info.get('last_mtime', 0):
                                _update_charts_from_file(tool_name, info['file'])
                                info['last_mtime'] = mtime
                    except Exception:
                        pass

            ui.timer(2.0, _refresh_sample_detail, active=True)

    
    def _create_workflow_monitor(self):
        """Create the main workflow monitoring page."""
        # Page title and navigation
        with ui.row().classes('w-full bg-blue-600 text-white p-4 items-center justify-between'):
            with ui.row().classes('items-center'):
                ui.label('📊 LittleJohn Workflow Monitor').classes('text-2xl font-bold')
                ui.label('Real-time workflow monitoring and control').classes('text-sm ml-4 opacity-80')
            
            # Navigation links
            with ui.row().classes('gap-4'):
                ui.link('🏠 Welcome', '/').classes('text-white hover:text-blue-200 text-sm')
                ui.link('📋 Sample Overview', '/live_data').classes('text-white hover:text-blue-200 text-sm')
                ui.label('📊 Workflow Monitor').classes('text-white text-sm font-semibold')
        
        # Main content area
        with ui.column().classes('w-full p-4 gap-4'):
            # Workflow Status Overview
            with ui.card().classes('w-full bg-gradient-to-r from-blue-50 to-indigo-50'):
                ui.label('🚀 Workflow Status Overview').classes('text-lg font-semibold mb-4 text-blue-800')
                
                # Status indicator
                with ui.row().classes('w-full items-center gap-4'):
                    self.status_indicator = ui.label('🟢').classes('text-2xl')
                    self.status_label = ui.label('Workflow Status: Running').classes('text-sm font-medium text-green-600')
                
                # Timing information
                with ui.row().classes('w-full gap-8 mt-4'):
                    self.workflow_start_time = ui.label('Started: --').classes('text-sm text-gray-600')
                    self.workflow_duration = ui.label('Duration: --').classes('text-sm text-gray-600')
                
                # Progress bar (hide internal float value text; use external formatted label below)
                self.progress_bar = ui.linear_progress(0.0).classes('w-full mt-4').style('color: transparent')
                self.progress_label = ui.label('0% Complete').classes('text-sm text-center text-gray-600')

                # Counts summary
                with ui.row().classes('w-full gap-6 mt-2'):
                    with ui.row().classes('items-center gap-2'):
                        ui.label('Completed:').classes('text-xs text-gray-600')
                        self.completed_count = ui.label('0').classes('text-xs font-semibold')
                    with ui.row().classes('items-center gap-2'):
                        ui.label('Failed:').classes('text-xs text-gray-600')
                        self.failed_count = ui.label('0').classes('text-xs font-semibold')
                    with ui.row().classes('items-center gap-2'):
                        ui.label('Total:').classes('text-xs text-gray-600')
                        self.total_count = ui.label('0').classes('text-xs font-semibold')
            
            # Queue Status
            with ui.card().classes('w-full'):
                ui.label('📋 Queue Status').classes('text-lg font-semibold mb-4')
                
                # Queue status grid
                with ui.grid(columns=4).classes('w-full gap-4'):
                    # Preprocessing
                    with ui.card().classes('bg-green-50 p-4'):
                        ui.label('🔬 Preprocessing').classes('text-sm font-medium text-green-800')
                        self.preprocessing_status = ui.label('0/0').classes('text-2xl font-bold text-green-600')
                    
                    # Analysis
                    with ui.card().classes('bg-blue-50 p-4'):
                        ui.label('🧬 Analysis').classes('text-sm font-medium text-blue-800')
                        self.analysis_status = ui.label('0/0').classes('text-2xl font-bold text-blue-600')
                    
                    # Classification
                    with ui.card().classes('bg-purple-50 p-4'):
                        ui.label('🎯 Classification').classes('text-sm font-medium text-purple-800')
                        self.classification_status = ui.label('0/0').classes('text-2xl font-bold text-purple-800')
                    
                    # Other
                    with ui.card().classes('bg-gray-50 p-4'):
                        ui.label('⚙️ Other').classes('text-sm font-medium text-gray-800')
                        self.other_status = ui.label('0/0').classes('text-2xl font-bold text-gray-600')
            
            # Active Jobs
            with ui.card().classes('w-full'):
                ui.label('⚡ Active Jobs').classes('text-lg font-semibold mb-4')
                
                # Active jobs table
                self.active_jobs_table = ui.table(
                    columns=[
                        {'name': 'job_id', 'label': 'Job ID', 'field': 'job_id'},
                        {'name': 'job_type', 'label': 'Type', 'field': 'job_type'},
                        {'name': 'filepath', 'label': 'File', 'field': 'filepath'},
                        {'name': 'worker', 'label': 'Worker', 'field': 'worker'},
                        {'name': 'duration', 'label': 'Duration', 'field': 'duration'},
                        {'name': 'progress', 'label': 'Progress', 'field': 'progress'}
                    ],
                    rows=[],
                    pagination=10
                ).classes('w-full')
                
                # Placeholder for when no jobs are active
                ui.label('No active jobs at the moment.').classes('text-sm text-gray-500 mt-2')
            
            # Live Logs
            with ui.card().classes('w-full'):
                ui.label('📝 Live Logs').classes('text-lg font-semibold mb-4')
                
                # Log controls
                with ui.row().classes('w-full justify-between items-center mb-2'):
                    with ui.row().classes('gap-2'):
                        ui.button('Clear', on_click=self._clear_logs).classes('bg-gray-500 hover:bg-gray-600 text-white text-xs')
                        ui.button('Export', on_click=lambda: self._export_logs()).classes('bg-blue-500 hover:bg-blue-600 text-white text-xs')
                
                # Log area
                self.log_area = ui.textarea('Workflow logs will appear here...').classes('w-full h-40').props('readonly')
            
            # Configuration
            with ui.card().classes('w-full'):
                ui.label('⚙️ Workflow Configuration').classes('text-lg font-semibold mb-4')
                
                # Configuration details
                with ui.grid(columns=2).classes('w-full gap-4'):
                    with ui.column():
                        ui.label('Monitored Directory:').classes('text-sm font-medium')
                        ui.label(self.monitored_directory or 'Not specified').classes('text-sm text-gray-600')
                        
                        ui.label('Workflow Steps:').classes('text-sm font-medium mt-2')
                        ui.label(', '.join(self.workflow_steps) if self.workflow_steps else 'Not specified').classes('text-sm text-gray-600')
                    
                    with ui.column():
                        ui.label('Log Level:').classes('text-sm font-medium')
                        ui.label('--').classes('text-sm text-gray-600')
                        
                        ui.label('Analysis Workers:').classes('text-sm font-medium mt-2')
                        ui.label('--').classes('text-sm text-gray-600')
            
            # Error Summary & Troubleshooting
            with ui.card().classes('w-full'):
                ui.label('⚠️ Error Summary & Troubleshooting').classes('text-lg font-semibold mb-2')
                
                # Error counts by type
                with ui.row().classes('w-full justify-between'):
                    with ui.column().classes('text-center'):
                        self.preprocessing_errors = ui.label('0').classes('text-xl font-bold text-red-600')
                        ui.label('Preprocessing').classes('text-xs text-gray-600')
                    with ui.column().classes('text-center'):
                        self.analysis_errors = ui.label('0').classes('text-xl font-bold text-red-600')
                        ui.label('Analysis').classes('text-xs text-gray-600')
                    with ui.column().classes('text-center'):
                        self.classification_errors = ui.label('0').classes('text-xl font-bold text-red-600')
                        ui.label('Classification').classes('text-xs text-gray-600')
                
                # Common error messages
                ui.separator()
                ui.label('Recent Errors').classes('text-sm font-medium mt-2')
                self.error_summary_label = ui.label('No errors detected').classes('text-xs text-gray-600')
            
            # Footer
            with ui.row().classes('w-full bg-gray-200 p-2 justify-center'):
                ui.label('LittleJohn Workflow Monitor - Running').classes('text-sm text-gray-600')
            
            # Signal that GUI is ready to receive updates
            self.gui_ready.set()
            logging.info("[GUI] UI created and ready to receive updates")
            # Start a periodic UI-thread drain of the update queue
            ui.timer(0.3, self._drain_updates_on_ui, active=True)
            # Start duration refresher
            ui.timer(1.0, lambda: self._refresh_duration(), active=True)

    def _refresh_duration(self):
        if self._is_running and self._start_time and hasattr(self, 'workflow_duration'):
            try:
                elapsed_seconds = int(time.time() - self._start_time)
                self.workflow_duration.set_text(f'Duration: {self._format_duration(elapsed_seconds)}')
            except Exception:
                pass
    
    def _export_logs(self):
        """Export logs to a file."""
        try:
            # Simple log export functionality
            log_content = ''.join(self._log_buffer)
            if log_content:
                # Create a simple download
                ui.download(log_content, 'workflow_logs.txt')
            else:
                ui.notify('No logs to export', type='warning')
        except Exception as e:
            ui.notify(f'Export failed: {e}', type='error')

    def _clear_logs(self):
        try:
            self._log_buffer.clear()
            self.log_area.set_value('')
        except Exception:
            pass
    
    def _launch_workflow_button_clicked(self):
        """Handle launch workflow button click."""
        ui.notify('Launch workflow functionality not implemented yet', type='info')
    
    def _view_docs_button_clicked(self):
        """Handle view documentation button click."""
        ui.notify('Documentation not available yet', type='info')
    
    def stop_gui(self):
        """Stop the GUI thread."""
        self.is_running = False
        # Note: NiceGUI doesn't have a clean shutdown method
        # The thread will terminate when the main process ends
        logging.info("GUI shutdown requested")
    
    def is_gui_running(self) -> bool:
        """Check if the GUI thread is running."""
        return self.is_running and self.gui_thread and self.gui_thread.is_alive()
    
    def get_gui_url(self) -> str:
        """Get the URL where the GUI is running."""
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
    
    # Store global reference for workflow updates
    global _current_gui_launcher
    _current_gui_launcher = launcher
    logging.info(f"[GUI] Global launcher set (id={id(_current_gui_launcher)})")
    
    if show:
        import webbrowser
        webbrowser.open(launcher.get_gui_url())
    
    return launcher


def get_gui_launcher() -> Optional[GUILauncher]:
    """Get the current GUI launcher instance for sending updates."""
    try:
        launcher = _current_gui_launcher
        logging.info(f"[GUI] get_gui_launcher -> {id(launcher) if launcher else 'None'}")
        return launcher
    except NameError:
        logging.info("[GUI] get_gui_launcher -> NameError (no global set)")
        return None


def send_gui_update(update_type: UpdateType, data: Dict[str, Any], priority: int = 0):
    """Send an update to the GUI without blocking the workflow."""
    launcher = get_gui_launcher()
    if launcher:
        launcher.send_update(update_type, data, priority)
    else:
        logging.info("[GUI] No GUI launcher available for update (update dropped)")


# Global reference to current GUI launcher
_current_gui_launcher = None
