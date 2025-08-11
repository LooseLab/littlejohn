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

# Import the global workflow state instance
try:
    from .workflow_state import workflow_state
except ImportError:
    try:
        from littlejohn.workflow_state import workflow_state
    except ImportError:
        workflow_state = None


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
        """Create the comprehensive workflow monitoring interface."""
        # Main header with status indicator
        with ui.header().classes('bg-blue-600 text-white p-4'):
            with ui.row().classes('w-full justify-between items-center'):
                ui.label('LittleJohn Workflow Monitor').classes('text-2xl font-bold')
                self.status_indicator = ui.label('●').classes('text-2xl text-yellow-400')
        
        # Main content area with comprehensive monitoring
        with ui.row().classes('w-full h-full'):
            # Left column - Workflow Overview & Status
            with ui.column().classes('w-1/3 p-4 border-r'):
                ui.label('Workflow Overview').classes('text-xl font-bold mb-4 text-blue-700')
                
                # Real-time Status Card
                with ui.card().classes('w-full mb-4 bg-gradient-to-r from-blue-50 to-indigo-50'):
                    ui.label('🔄 Real-time Status').classes('text-lg font-semibold mb-2 text-blue-800')
                    self.status_label = ui.label('Initializing...').classes('text-sm font-medium')
                    
                    # Status indicator
                    with ui.row().classes('w-full justify-center mb-2'):
                        self.status_indicator
                    
                    # Timing information
                    self.started_time_label = ui.label('Started: --').classes('text-xs text-gray-600')
                    self.elapsed_time_label = ui.label('Elapsed: --').classes('text-xs text-gray-600')
                
                # Progress Overview Card
                with ui.card().classes('w-full mb-4'):
                    ui.label('📊 Progress Overview').classes('text-lg font-semibold mb-2')
                    
                    # Overall progress
                    ui.label('Overall Progress').classes('text-sm font-medium mt-2')
                    self.overall_progress = ui.linear_progress().classes('w-full mb-1')
                    self.progress_label = ui.label('0% Complete').classes('text-xs text-gray-600 text-center')
                    
                    # Job summary
                    with ui.row().classes('w-full justify-between mt-3'):
                        with ui.column().classes('text-center'):
                            self.active_jobs_label = ui.label('0').classes('text-2xl font-bold text-blue-600')
                            ui.label('Active').classes('text-xs text-gray-600')
                        with ui.column().classes('text-center'):
                            self.completed_jobs_label = ui.label('0').classes('text-2xl font-bold text-green-600')
                            ui.label('Completed').classes('text-xs text-gray-600')
                        with ui.column().classes('text-center'):
                            self.failed_jobs_label = ui.label('0').classes('text-2xl font-bold text-red-600')
                            ui.label('Failed').classes('text-xs text-gray-600')
                
                # Queue Performance Card
                with ui.card().classes('w-full mb-4'):
                    ui.label('⚡ Queue Performance').classes('text-lg font-semibold mb-2')
                    
                    # Queue progress bars
                    self.preprocessing_progress = ui.label('Preprocessing: 0/0').classes('text-sm')
                    self.analysis_progress = ui.label('Analysis: 0/0').classes('text-sm')
                    self.classification_progress = ui.label('Classification: 0/0').classes('text-sm')
                    
                    # Performance metrics
                    ui.separator()
                    ui.label('Performance Metrics').classes('text-sm font-medium mt-2')
                    self.jobs_per_second_label = ui.label('Jobs/sec: --').classes('text-xs text-gray-600')
                    self.avg_job_time_label = ui.label('Avg Job Time: --').classes('text-xs text-gray-600')
                
                # File Monitoring Card
                with ui.card().classes('w-full mb-4'):
                    ui.label('📁 File Monitoring').classes('text-lg font-semibold mb-2')
                    
                    self.detected_files_label = ui.label('BAM Files: 0').classes('text-sm')
                    self.processed_files_label = ui.label('Output Files: 0').classes('text-sm')
                    self.failed_files_label = ui.label('Failed Files: 0').classes('text-sm')
                    
                    # File processing rate
                    ui.separator()
                    ui.label('Processing Rate').classes('text-sm font-medium mt-2')
                    self.files_per_minute_label = ui.label('Files/min: --').classes('text-xs text-gray-600')
                
                # Configuration Card
                with ui.card().classes('w-full mb-4'):
                    ui.label('⚙️ Configuration').classes('text-lg font-semibold mb-2')
                    self.workflow_config_label = ui.label('Workflow: Loading...').classes('text-xs text-gray-600')
                    
                    # Queue priorities
                    ui.separator()
                    ui.label('Queue Priorities').classes('text-sm font-medium mt-2')
                    self.priority_info_label = ui.label('Loading...').classes('text-xs text-gray-600')
            
            # Right column - Detailed Monitoring & Control
            with ui.column().classes('w-2/3 p-4'):
                ui.label('Detailed Monitoring & Control').classes('text-2xl font-bold mb-6 text-blue-700')
                
                # Current Activity & Job Details
                with ui.card().classes('w-full mb-6'):
                    ui.label('🎯 Current Activity').classes('text-lg font-semibold mb-2')
                    self.current_activity_label = ui.label('No active jobs').classes('text-sm text-gray-600')
                    
                    # Active jobs table
                    ui.separator()
                    ui.label('Active Jobs').classes('text-sm font-medium mt-2')
                    self.active_jobs_table = ui.table(
                        columns=[
                            {'name': 'job_id', 'label': 'Job ID', 'field': 'job_id'},
                            {'name': 'type', 'label': 'Type', 'field': 'type'},
                            {'name': 'file', 'label': 'File', 'field': 'file'},
                            {'name': 'queue', 'label': 'Queue', 'field': 'queue'},
                            {'name': 'start_time', 'label': 'Started', 'field': 'start_time'},
                            {'name': 'duration', 'label': 'Duration', 'field': 'duration'}
                        ],
                        rows=[],
                        pagination=5
                    ).classes('w-full')
                
                # Live Logs with Filtering
                with ui.card().classes('w-full mb-6'):
                    with ui.row().classes('w-full justify-between items-center mb-2'):
                        ui.label('📝 Live Logs').classes('text-lg font-semibold')
                        with ui.row().classes('gap-2'):
                            self.log_level_filter = ui.select(
                                ['ALL', 'INFO', 'WARNING', 'ERROR', 'DEBUG'],
                                value='ALL',
                                label='Log Level'
                            ).classes('w-24')
                            ui.button('Clear', on_click=lambda: self.log_area.set_value('')).classes('bg-gray-500 hover:bg-gray-600 text-white text-xs')
                            ui.button('Export', on_click=lambda: self._export_logs()).classes('bg-blue-500 hover:bg-blue-600 text-white text-xs')
                    
                    self.log_area = ui.textarea('Workflow logs will appear here...').classes('w-full h-40').props('readonly')
                
                # Recent Events & Notifications
                with ui.card().classes('w-full mb-6'):
                    ui.label('🔔 Recent Events & Notifications').classes('text-lg font-semibold mb-2')
                    
                    # Events list
                    self.events_list = ui.column().classes('w-full')
                    ui.label('Workflow started').classes('text-sm text-gray-600')
                    ui.label('Monitoring directory for BAM files').classes('text-sm text-gray-600')
                
                # Metadata Monitoring Section
                with ui.card().classes('w-full mb-6'):
                    ui.label('📋 Preprocessing & Workflow Metadata').classes('text-lg font-semibold mb-2')
                    
                    # Workflow context metadata display
                    ui.label('Workflow Context Metadata').classes('text-sm font-medium mt-2')
                    self.workflow_context_display = ui.textarea('No workflow context metadata available...').classes('w-full h-32').props('readonly')
                    
                    # Metadata summary
                    ui.separator()
                    ui.label('Preprocessing Summary').classes('text-sm font-medium mt-2')
                    self.metadata_summary = ui.label('No preprocessing metadata available').classes('text-xs text-gray-600')
                
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
            
            # Function to update UI with comprehensive real-time data
            def update_ui():
                try:
                    # Try to get workflow state if available
                    try:
                        # from littlejohn.workflow_state import workflow_state # This line is removed
                        summary = workflow_state.get_workflow_summary()
                        
                        # Update status indicator
                        if summary['is_running']:
                            self.status_indicator.classes('text-2xl text-green-400')
                            self.status_label.set_text('Workflow Status: Running')
                            self.status_label.classes('text-sm font-medium text-green-600')
                        else:
                            self.status_indicator.classes('text-2xl text-red-400')
                            self.status_label.set_text('Workflow Status: Stopped')
                            self.status_label.classes('text-sm font-medium text-red-600')
                        
                        # Update timing information
                        if summary['start_time']:
                            start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(summary['start_time']))
                            self.started_time_label.set_text(f'Started: {start_time}')
                            
                            # Calculate elapsed time
                            elapsed_seconds = int(time.time() - summary['start_time'])
                            elapsed_time = self._format_duration(elapsed_seconds)
                            self.elapsed_time_label.set_text(f'Elapsed: {elapsed_time}')
                        
                        # Update job counts with enhanced styling
                        self.active_jobs_label.set_text(str(summary['running_jobs']))
                        self.completed_jobs_label.set_text(str(summary['completed_jobs']))
                        self.failed_jobs_label.set_text(str(summary['failed_jobs']))
                        
                        # Update progress
                        total_jobs = summary['total_jobs']
                        if total_jobs > 0:
                            progress = summary['completed_jobs'] / total_jobs
                            self.overall_progress.set_value(progress)
                            self.progress_label.set_text(f'{int(progress * 100)}% Complete')
                        
                        # Update queue progress with detailed counts
                        queue_status = summary['queue_status']
                        self.preprocessing_progress.set_text(f'Preprocessing: {queue_status["preprocessing"]["running"]}/{queue_status["preprocessing"]["total"]}')
                        
                        analysis_running = queue_status["mgmt"]["running"] + queue_status["cnv"]["running"] + queue_status["target"]["running"] + queue_status["fusion"]["running"]
                        analysis_total = queue_status["mgmt"]["total"] + queue_status["cnv"]["total"] + queue_status["target"]["total"] + queue_status["fusion"]["total"]
                        self.analysis_progress.set_text(f'Analysis: {analysis_running}/{analysis_total}')
                        
                        self.classification_progress.set_text(f'Classification: {queue_status["classification"]["running"]}/{queue_status["classification"]["total"]}')
                        
                        # Calculate and update performance metrics
                        if summary['start_time'] and summary['completed_jobs'] > 0:
                            elapsed = time.time() - summary['start_time']
                            jobs_per_second = summary['completed_jobs'] / elapsed if elapsed > 0 else 0
                            avg_job_time = elapsed / summary['completed_jobs'] if summary['completed_jobs'] > 0 else 0
                            
                            self.jobs_per_second_label.set_text(f'Jobs/sec: {jobs_per_second:.2f}')
                            self.avg_job_time_label.set_text(f'Avg Job Time: {self._format_duration(int(avg_job_time))}')
                        
                        # Update file counts
                        self.detected_files_label.set_text(f'BAM Files: {summary["detected_files"]}')
                        self.processed_files_label.set_text(f'Output Files: {summary["processed_files"]}')
                        self.failed_files_label.set_text(f'Failed Files: {summary["failed_files"]}')
                        
                        # Calculate file processing rate
                        if summary['start_time'] and summary['processed_files'] > 0:
                            elapsed_minutes = (time.time() - summary['start_time']) / 60
                            files_per_minute = summary['processed_files'] / elapsed_minutes if elapsed_minutes > 0 else 0
                            self.files_per_minute_label.set_text(f'Files/min: {files_per_minute:.1f}')
                        
                        # Update configuration
                        if summary['workflow_steps']:
                            self.workflow_config_label.set_text(f'Workflow: {", ".join(summary["workflow_steps"])}')
                        
                        # Update queue priorities if available
                        if hasattr(workflow_state, 'queue_priorities'):
                            priority_text = ', '.join([f'{q}:{p}' for q, p in workflow_state.queue_priorities.items()])
                            self.priority_info_label.set_text(priority_text[:100] + '...' if len(priority_text) > 100 else priority_text)
                        
                        # Update active jobs table
                        running_jobs = workflow_state.get_jobs_by_status("running")
                        if running_jobs:
                            table_rows = []
                            for job in running_jobs[:10]:  # Show up to 10 active jobs
                                job_start = getattr(job, 'start_time', time.time())
                                duration = int(time.time() - job_start) if job_start else 0
                                
                                table_rows.append({
                                    'job_id': str(job.job_id)[:8],
                                    'type': job.job_type,
                                    'file': Path(job.filepath).name if hasattr(job, 'filepath') else 'Unknown',
                                    'queue': getattr(job, 'queue_name', 'Unknown'),
                                    'start_time': time.strftime('%H:%M:%S', time.localtime(job_start)) if job_start else '--',
                                    'duration': self._format_duration(duration)
                                })
                            
                            # Clear existing rows and add new ones
                            self.active_jobs_table.clear()
                            if table_rows:
                                self.active_jobs_table.add_rows(table_rows)
                        
                        # Update recent logs with filtering
                        recent_logs = workflow_state.get_recent_logs(50)
                        if recent_logs:
                            # Apply log level filter
                            selected_level = self.log_level_filter.value
                            if selected_level != 'ALL':
                                filtered_logs = [log for log in recent_logs if log.get('level', 'INFO') == selected_level]
                            else:
                                filtered_logs = recent_logs
                            
                            log_text = '\n'.join([f'[{time.strftime("%H:%M:%S", time.localtime(log["timestamp"]))}] {log["level"]}: {log["message"]}' for log in filtered_logs[-30:]])
                            self.log_area.set_value(log_text)
                        
                        # Update current activity
                        if running_jobs:
                            activity_text = '\n'.join([f'{job.job_type}: {Path(job.filepath).name if hasattr(job, "filepath") else "Unknown"}' for job in running_jobs[:5]])
                            self.current_activity_label.set_text(activity_text)
                        else:
                            self.current_activity_label.set_text('No active jobs')
                        
                        # Update error summary
                        failed_jobs = workflow_state.get_jobs_by_status("failed")
                        if failed_jobs:
                            # Count errors by type
                            error_counts = {'preprocessing': 0, 'analysis': 0, 'classification': 0}
                            for job in failed_jobs:
                                if job.job_type in ['preprocessing', 'bed_conversion']:
                                    error_counts['preprocessing'] += 1
                                elif job.job_type in ['mgmt', 'cnv', 'target', 'fusion']:
                                    error_counts['analysis'] += 1
                                elif job.job_type in ['sturgeon', 'nanodx', 'pannanodx', 'random_forest']:
                                    error_counts['classification'] += 1
                            
                            self.preprocessing_errors.set_text(str(error_counts['preprocessing']))
                            self.analysis_errors.set_text(str(error_counts['analysis']))
                            self.classification_errors.set_text(str(error_counts['classification']))
                            
                            # Show recent error messages
                            recent_errors = [job for job in failed_jobs[-5:]]  # Last 5 failed jobs
                            error_text = '\n'.join([f'{job.job_type}: {getattr(job, "error_message", "Unknown error")}' for job in recent_errors])
                            self.error_summary_label.set_text(error_text[:200] + '...' if len(error_text) > 200 else error_text)
                        else:
                            self.preprocessing_errors.set_text('0')
                            self.analysis_errors.set_text('0')
                            self.classification_errors.set_text('0')
                            self.error_summary_label.set_text('No errors detected')
                        
                        # Update workflow context metadata
                        # Monitor master.csv file for preprocessing metadata
                        if self.monitored_directory:
                            master_csv_path = Path(self.monitored_directory) / "master.csv"
                            if master_csv_path.exists():
                                try:
                                    # Read and display master.csv content
                                    import pandas as pd
                                    df = pd.read_csv(master_csv_path)
                                    
                                    if not df.empty:
                                        # Display preprocessing metadata from master.csv
                                        metadata_text = f"Preprocessing Metadata from master.csv\n"
                                        metadata_text += f"File: {master_csv_path.name}\n"
                                        metadata_text += f"Total samples: {len(df)}\n"
                                        metadata_text += f"Last updated: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(master_csv_path.stat().st_mtime))}\n\n"
                                        
                                        # Show column information
                                        metadata_text += f"Available columns:\n"
                                        for col in df.columns:
                                            metadata_text += f"  - {col}\n"
                                        
                                        metadata_text += f"\nSample Data (first 5 samples):\n"
                                        metadata_text += "-" * 50 + "\n"
                                        
                                        # Display first 5 samples with key preprocessing data
                                        for idx, row in df.head().iterrows():
                                            metadata_text += f"Sample {idx + 1}:\n"
                                            
                                            # Show key preprocessing fields if available
                                            if 'sample_id' in df.columns:
                                                metadata_text += f"  Sample ID: {row.get('sample_id', 'N/A')}\n"
                                            if 'state' in df.columns:
                                                metadata_text += f"  State: {row.get('state', 'N/A')}\n"
                                            if 'mapped_reads' in df.columns:
                                                mapped_reads = row.get('mapped_reads', 0)
                                                if pd.notna(mapped_reads):
                                                    metadata_text += f"  Mapped reads: {int(mapped_reads):,}\n"
                                            if 'unmapped_reads' in df.columns:
                                                unmapped_reads = row.get('unmapped_reads', 0)
                                                if pd.notna(unmapped_reads):
                                                    metadata_text += f"  Unmapped reads: {int(unmapped_reads):,}\n"
                                            if 'yield_tracking' in df.columns:
                                                yield_tracking = row.get('yield_tracking', 0)
                                                if pd.notna(yield_tracking):
                                                    metadata_text += f"  Total yield: {int(yield_tracking):,} bases\n"
                                            if 'file_size' in df.columns:
                                                file_size = row.get('file_size', 0)
                                                if pd.notna(file_size):
                                                    metadata_text += f"  File size: {int(file_size):,} bytes\n"
                                            if 'has_supplementary_reads' in df.columns:
                                                has_supp = row.get('has_supplementary_reads', False)
                                                if pd.notna(has_supp):
                                                    metadata_text += f"  Has supplementary reads: {has_supp}\n"
                                            if 'has_mgmt_reads' in df.columns:
                                                has_mgmt = row.get('has_mgmt_reads', False)
                                                if pd.notna(has_mgmt):
                                                    metadata_text += f"  Has MGMT reads: {has_mgmt}\n"
                                            
                                            metadata_text += "-" * 30 + "\n"
                                        
                                        # Show summary statistics
                                        metadata_text += f"\nSummary Statistics:\n"
                                        if 'mapped_reads' in df.columns:
                                            total_mapped = df['mapped_reads'].sum()
                                            if pd.notna(total_mapped):
                                                metadata_text += f"Total mapped reads: {int(total_mapped):,}\n"
                                        if 'yield_tracking' in df.columns:
                                            total_yield = df['yield_tracking'].sum()
                                            if pd.notna(total_yield):
                                                metadata_text += f"Total yield: {int(total_yield):,} bases\n"
                                        if 'state' in df.columns:
                                            pass_count = len(df[df['state'] == 'pass'])
                                            fail_count = len(df[df['state'] == 'fail'])
                                            metadata_text += f"Pass samples: {pass_count}, Fail samples: {fail_count}\n"
                                        
                                        self.workflow_context_display.set_value(metadata_text)
                                        self.metadata_summary.set_text(f"Preprocessing metadata loaded from master.csv - {len(df)} samples processed")
                                        
                                    else:
                                        self.workflow_context_display.set_value('master.csv exists but is empty. No preprocessing data available yet.')
                                        self.metadata_summary.set_text('master.csv empty - no preprocessing data')
                                        
                                except Exception as e:
                                    error_msg = f"Error reading master.csv: {e}\n\n"
                                    error_msg += f"File: {master_csv_path}\n"
                                    error_msg += f"Size: {master_csv_path.stat().st_size:,} bytes\n"
                                    error_msg += f"Last modified: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(master_csv_path.stat().st_mtime))}"
                                    
                                    self.workflow_context_display.set_value(error_msg)
                                    self.metadata_summary.set_text('Error reading master.csv')
                            else:
                                # Fallback to workflow state if master.csv doesn't exist
                                preprocessing_jobs = workflow_state.get_jobs_by_status("completed")
                                preprocessing_metadata = []
                                
                                for job_info in preprocessing_jobs:
                                    if job_info.job_type == "preprocessing":
                                        # Get the actual preprocessing metadata and results
                                        job_metadata = workflow_state.get_job_metadata(job_info.job_id) or {}
                                        job_results = workflow_state.get_job_results(job_info.job_id) or {}
                                        
                                        # Extract preprocessing metadata from job context
                                        bam_metadata = job_metadata.get('bam_metadata', {})
                                        preprocessing_result = job_results.get('preprocessing', {})
                                        
                                        # Create comprehensive metadata display
                                        job_display = {
                                            'job_id': job_info.job_id,
                                            'sample_id': bam_metadata.get('sample_id', job_info.sample_id),
                                            'filepath': job_info.filepath,
                                            'status': 'completed',
                                            'timestamp': job_info.start_time,
                                            'bam_metadata': bam_metadata,
                                            'preprocessing_result': preprocessing_result
                                        }
                                        preprocessing_metadata.append(job_display)
                                
                                if preprocessing_metadata:
                                    # Display comprehensive preprocessing metadata
                                    metadata_text = "Preprocessing Metadata:\n\n"
                                    for meta in preprocessing_metadata:
                                        metadata_text += f"Job ID: {meta['job_id']}\n"
                                        metadata_text += f"Sample ID: {meta['sample_id']}\n"
                                        metadata_text += f"File: {Path(meta['filepath']).name}\n"
                                        metadata_text += f"Status: {meta['status']}\n"
                                        if meta['timestamp']:
                                            metadata_text += f"Completed: {time.strftime('%H:%M:%S', time.localtime(meta['timestamp']))}\n"
                                        
                                        # Display BAM metadata if available
                                        bam_metadata = meta.get('bam_metadata', {})
                                        if bam_metadata:
                                            metadata_text += f"\nBAM Metadata:\n"
                                            metadata_text += f"  State: {bam_metadata.get('state', 'unknown')}\n"
                                            metadata_text += f"  Mapped reads: {bam_metadata.get('mapped_reads', 0):,}\n"
                                            metadata_text += f"  Unmapped reads: {bam_metadata.get('unmapped_reads', 0):,}\n"
                                            metadata_text += f"  Total yield: {bam_metadata.get('yield_tracking', 0):,} bases\n"
                                            metadata_text += f"  File size: {bam_metadata.get('file_size', 0):,} bytes\n"
                                            
                                            # Show supplementary read info
                                            if bam_metadata.get('has_supplementary_reads', False):
                                                metadata_text += f"  Supplementary reads: {bam_metadata.get('supplementary_reads', 0):,}\n"
                                                metadata_text += f"  Reads with supplementary: {bam_metadata.get('reads_with_supplementary', 0):,}\n"
                                            
                                            # Show MGMT read info
                                            if bam_metadata.get('has_mgmt_reads', False):
                                                metadata_text += f"  MGMT reads: {bam_metadata.get('mgmt_read_count', 0):,}\n"
                                        
                                        # Display preprocessing results if available
                                        preprocessing_result = meta.get('preprocessing_result', {})
                                        if preprocessing_result:
                                            metadata_text += f"\nPreprocessing Results:\n"
                                            metadata_text += f"  Status: {preprocessing_result.get('status', 'unknown')}\n"
                                            metadata_text += f"  Sample ID: {preprocessing_result.get('sample_id', 'unknown')}\n"
                                            metadata_text += f"  State: {preprocessing_result.get('state', 'unknown')}\n"
                                            metadata_text += f"  Mapped reads: {preprocessing_result.get('mapped_reads', 0):,}\n"
                                            metadata_text += f"  Unmapped reads: {preprocessing_result.get('unmapped_reads', 0):,}\n"
                                            metadata_text += f"  Total yield: {preprocessing_result.get('total_yield', 0):,} bases\n"
                                        
                                        metadata_text += "-" * 50 + "\n"
                                    
                                    self.workflow_context_display.set_value(metadata_text)
                                    self.metadata_summary.set_text(f"Preprocessing metadata available - {len(preprocessing_metadata)} completed jobs with detailed BAM analysis")
                                else:
                                    self.workflow_context_display.set_value('No preprocessing metadata available yet.\n\nPreprocessing jobs will appear here as they complete with detailed BAM analysis results.')
                                    self.metadata_summary.set_text('No preprocessing metadata available')
                        else:
                            self.workflow_context_display.set_value('No monitored directory specified. Cannot access preprocessing metadata.')
                            self.metadata_summary.set_text('No monitored directory')
                            
                    except ImportError:
                        # Workflow state not available, show static info
                        self._show_limited_mode_ui()
                        
                except Exception as e:
                    # If workflow state is not available, show static info
                    self._show_limited_mode_ui()
                    print(f"GUI update error: {e}")
            
            # Set up auto-refresh timer
            ui.timer(2.0, update_ui)
            
            # Initial update
            update_ui()
    
    def _show_limited_mode_ui(self):
        """Show limited mode UI when workflow state is not available."""
        try:
            # Update basic status
            self.status_indicator.classes('text-2xl text-yellow-400')
            self.status_label.set_text('Workflow Status: Limited Mode')
            self.status_label.classes('text-sm font-medium text-yellow-600')
            
            # Show static configuration
            self.workflow_config_label.set_text(f'Workflow: {", ".join(self.workflow_steps) if self.workflow_steps else "Unknown"}')
            self.detected_files_label.set_text(f'Monitored Directory: {self.monitored_directory}')
            
            # Show limited mode message
            self.current_activity_label.set_text('Real-time updates not available - check CLI output for progress')
            
            # Add basic log information
            self.log_area.set_value('GUI is running in limited mode.\nReal-time workflow updates not available.\nCheck the CLI output for detailed progress information.\n\nThis is normal when workflow hooks are not fully integrated.')
            
        except Exception as e:
            print(f"Limited mode UI error: {e}")
    
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
    
    def _export_logs(self):
        """Export current logs to a file."""
        try:
            from pathlib import Path
            import tempfile
            
            # Create a temporary file
            temp_dir = Path(tempfile.gettempdir())
            log_file = temp_dir / f"littlejohn_workflow_logs_{int(time.time())}.txt"
            
            # Get current log content
            log_content = self.log_area.value if hasattr(self.log_area, 'value') else 'No logs available'
            
            # Write to file
            with open(log_file, 'w') as f:
                f.write(f"LittleJohn Workflow Logs\n")
                f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Workflow: {', '.join(self.workflow_steps) if self.workflow_steps else 'Unknown'}\n")
                f.write(f"Monitored Directory: {self.monitored_directory}\n")
                f.write("-" * 50 + "\n\n")
                f.write(log_content)
            
            # Show success message
            ui.notify(f'Logs exported to: {log_file}', type='positive')
            
        except Exception as e:
            ui.notify(f'Failed to export logs: {e}', type='negative')
    
    def _update_metadata_summary(self, metadata_files):
        """Update the preprocessing metadata summary information."""
        try:
            if not metadata_files:
                self.metadata_summary.set_text('No preprocessing metadata available')
                return
            
            # For workflow context, show a simple summary
            if metadata_files and metadata_files[0].get('type') == 'workflow_context':
                self.metadata_summary.set_text('Workflow context metadata available - showing preprocessing and job results')
                return
            
            # Count files by workflow step
            step_counts = {}
            total_size = 0
            
            for file_info in metadata_files:
                file_type = file_info.get('type', 'general')
                step_counts[file_type] = step_counts.get(file_type, 0) + 1
                total_size += file_info['size']
            
            # Create preprocessing-focused summary
            summary_parts = []
            summary_parts.append(f"{len(metadata_files)} metadata files")
            
            # Show preprocessing files first
            if 'preprocessing' in step_counts:
                summary_parts.append(f"Preprocessing: {step_counts['preprocessing']}")
            
            # Show analysis files
            analysis_files = sum(step_counts.get(step, 0) for step in ['fusion_analysis', 'cnv_analysis', 'target_analysis', 'sturgeon_analysis'])
            if analysis_files > 0:
                summary_parts.append(f"Analysis: {analysis_files}")
            
            # Show workflow files
            if 'workflow' in step_counts or 'execution' in step_counts:
                workflow_files = step_counts.get('workflow', 0) + step_counts.get('execution', 0)
                summary_parts.append(f"Workflow: {workflow_files}")
            
            # Show total size
            if total_size > 0:
                size_mb = total_size / (1024 * 1024)
                summary_parts.append(f"Total size: {size_mb:.1f} MB")
            
            # Add most recent preprocessing file info
            preprocessing_files = [f for f in metadata_files if f.get('type') == 'preprocessing']
            if preprocessing_files:
                most_recent = preprocessing_files[0]
                recent_time = time.strftime('%H:%M:%S', time.localtime(most_recent['modified']))
                summary_parts.append(f"Latest preprocessing: {most_recent['name']} ({recent_time})")
            
            summary_text = ' | '.join(summary_parts)
            self.metadata_summary.set_text(summary_text)
            
        except Exception as e:
            print(f"Metadata summary error: {e}")
            self.metadata_summary.set_text('Error updating preprocessing metadata summary')
    
    def _on_metadata_file_change(self, event):
        """Handle metadata file selection change."""
        try:
            if event.value and event.value != 'No metadata files detected':
                self._display_metadata_content()
        except Exception as e:
            print(f"Metadata file change error: {e}")
    
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
