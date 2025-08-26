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
from datetime import datetime
from littlejohn.analysis.master_csv_manager import MasterCSVManager
from typing import Optional, Dict, Any, List
from pathlib import Path
from dataclasses import dataclass
from enum import Enum
import os

from littlejohn.gui import theme, images

from littlejohn.gui.components.news_feed import NewsFeed

from littlejohn.reporting.report import create_pdf
from littlejohn.reporting.sections.disclaimer_text import EXTENDED_DISCLAIMER_TEXT


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

    def __init__(self, host: str = "0.0.0.0", port: int = 8081):
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

        self.news_feed = None
        # GUI update thread
        self.update_thread = (
            None  # not used anymore; updates processed on UI thread via timer
        )

        # Debug counters
        self.total_updates_enqueued = 0
        self.total_updates_processed = 0
        # Monotonic tiebreaker for priority queue to avoid tuple comparison of GUIUpdate
        self._update_seq: int = 0

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
        self._preexisting_sample_ids: set[str] = set()
        self._preexisting_scanned: bool = False
        # MGMT per-sample cache (latest count seen)
        self._mgmt_state: Dict[str, Dict[str, Any]] = {}
        # Coverage per-sample cache (file mtimes, computed metrics)
        self._coverage_state: Dict[str, Dict[str, Any]] = {}
        # CNV per-sample cache
        self._cnv_state: Dict[str, Dict[str, Any]] = {}
        # Cache last seen queue status so we can populate immediately on page creation
        self._last_queue_status: Dict[str, Any] = {}

    def launch_gui(
        self,
        workflow_runner: Any = None,
        workflow_steps: list = None,
        monitored_directory: str = "",
    ) -> bool:
        """Launch the GUI in a completely isolated background thread."""
        if ui is None:
            logging.error("NiceGUI is not available")
            return False

        self.workflow_runner = workflow_runner
        self.workflow_steps = workflow_steps or []
        
        # Store absolute monitored directory to avoid relative path issues
        try:
            self.monitored_directory = (
                str(Path(monitored_directory).resolve()) if monitored_directory else ""
            )
        except Exception:
            self.monitored_directory = monitored_directory

        try:
            # Start GUI in completely isolated background thread
            self.gui_thread = threading.Thread(
                target=self._run_gui_worker, daemon=True, name="LittleJohn-GUI-Thread"
            )
            self.gui_thread.start()

            # NOTE: Update processing now happens inside the UI thread via ui.timer for thread-safety

            # Wait a moment for GUI to start
            time.sleep(1)

            # Check if thread is still running
            if self.gui_thread.is_alive():
                self.is_running = True
                logging.info(
                    f"GUI launched successfully on http://{self.host}:{self.port}"
                )
                return True
            else:
                logging.error("GUI thread failed to start")
                return False

        except Exception as e:
            logging.error(f"Failed to launch GUI: {e}")
            return False

    def send_update(
        self, update_type: UpdateType, data: Dict[str, Any], priority: int = 0
    ):
        """Send an update to the GUI without blocking the workflow."""
        try:
            update = GUIUpdate(
                update_type=update_type,
                timestamp=time.time(),
                data=data,
                priority=priority,
            )

            # Use negative priority so higher priority updates come first.
            # Include a monotonically increasing sequence as a tiebreaker so
            # heap comparisons never fallback to comparing GUIUpdate objects.
            self._update_seq += 1
            self.update_queue.put((-priority, self._update_seq, update))
            self.total_updates_enqueued += 1
            logging.info(
                f"[GUI] Enqueued update #{self.total_updates_enqueued}: {update.update_type.value}"
            )

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
                    # Entries are (-priority, seq, GUIUpdate)
                    _, _, update = self.update_queue.get_nowait()
                except queue.Empty:
                    break
                self._handle_update(update)
                self.total_updates_processed += 1
                processed_any = True
                logging.info(
                    f"[GUI] Processed update #{self.total_updates_processed}: {update.update_type.value}"
                )
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
            if hasattr(self, "status_indicator") and hasattr(self, "status_label"):
                if data.get("is_running", False):
                    self.status_indicator.classes("text-2xl text-green-400")
                    self.status_label.set_text("Workflow Status: Running")
                    self._is_running = True
                else:
                    self.status_indicator.classes("text-2xl text-red-400")
                    self.status_label.set_text("Workflow Status: Stopped")
                    self._is_running = False

                # Update timing
                if data.get("start_time"):
                    self._start_time = float(data["start_time"])
                    if hasattr(self, "workflow_start_time"):
                        start_str = time.strftime(
                            "%Y-%m-%d %H:%M:%S", time.localtime(self._start_time)
                        )
                        self.workflow_start_time.set_text(f"Started: {start_str}")

                if hasattr(self, "workflow_duration") and self._start_time:
                    elapsed_seconds = int(time.time() - self._start_time)
                    self.workflow_duration.set_text(
                        f"Duration: {self._format_duration(elapsed_seconds)}"
                    )

        except Exception as e:
            logging.debug(f"Error updating workflow status: {e}")

    def _update_job_status(self, data: Dict[str, Any]):
        """Update job status in the GUI with search/sort/filter support."""
        try:
            if hasattr(self, "active_jobs_table"):
                # Build all rows
                job_rows: List[Dict[str, Any]] = []
                for job in data.get("active_jobs", []):
                    job_rows.append(
                        {
                            "job_id": str(job.get("job_id", ""))[:8],
                            "job_type": job.get("job_type", ""),
                            "filepath": Path(job.get("filepath", "")).name,
                            "worker": job.get("worker_name", ""),
                            "duration": self._format_duration(job.get("duration", 0)),
                            "progress": f"{int(job.get('progress', 0) * 100)}%",
                        }
                    )

                # Cache for filtering
                self._active_jobs_all_rows = job_rows

                # Update filter options dynamically
                try:
                    type_options = ["All"] + sorted(
                        {r.get("job_type", "") for r in job_rows if r.get("job_type")}
                    )
                    worker_options = ["All"] + sorted(
                        {r.get("worker", "") for r in job_rows if r.get("worker")}
                    )
                    if hasattr(self, "active_jobs_type_filter"):
                        self.active_jobs_type_filter.set_options(type_options)
                    if hasattr(self, "active_jobs_worker_filter"):
                        self.active_jobs_worker_filter.set_options(worker_options)
                except Exception:
                    pass

                # Apply filters and refresh table
                self._apply_active_jobs_filters_and_update()

        except Exception as e:
            logging.debug(f"Error updating job status: {e}")

    def _update_queue_status(self, data: Dict[str, Any]):
        """Update queue status in the GUI."""
        try:
            # Cache the latest payload so newly created pages can show current values
            self._last_queue_status = data or {}
            # Emit a concise log line to the Live Logs pane for visibility
            try:
                pr = data.get("preprocessing", {})
                ar = {
                    "mgmt": data.get("mgmt", {}),
                    "cnv": data.get("cnv", {}),
                    "target": data.get("target", {}),
                    "fusion": data.get("fusion", {}),
                }
                an_running = sum(int(v.get("running", 0) or 0) for v in ar.values())
                an_total = sum(int(v.get("total", 0) or 0) for v in ar.values())
                cl = data.get("classification", {})
                ot = data.get("other", {})
                msg = (
                    f"Queues | Pre:{pr.get('running',0)}/{pr.get('total',0)} "
                    f"| An:{an_running}/{an_total} "
                    f"| Cl:{cl.get('running',0)}/{cl.get('total',0)} "
                    f"| Ot:{ot.get('running',0)}/{ot.get('total',0)}"
                )
                self._log_buffer.append(f"[queue] {msg}\n")
                if hasattr(self, "log_area"):
                    self.log_area.set_value("".join(self._log_buffer))
            except Exception:
                pass
            # Update queue status displays
            if hasattr(self, "preprocessing_status"):
                queue_data = data.get("preprocessing", {})
                self.preprocessing_status.set_text(
                    f"{queue_data.get('running', 0)}/{queue_data.get('total', 0)}"
                )

            if hasattr(self, "analysis_status"):
                # Combine analysis queues
                analysis_running = sum(
                    data.get(q, {}).get("running", 0)
                    for q in ["mgmt", "cnv", "target", "fusion"]
                )
                analysis_total = sum(
                    data.get(q, {}).get("total", 0)
                    for q in ["mgmt", "cnv", "target", "fusion"]
                )
                self.analysis_status.set_text(f"{analysis_running}/{analysis_total}")

            if hasattr(self, "classification_status"):
                classification_data = data.get("classification", {})
                self.classification_status.set_text(
                    f"{classification_data.get('running', 0)}/{classification_data.get('total', 0)}"
                )

            if hasattr(self, "other_status"):
                other_data = data.get("other", {})
                self.other_status.set_text(
                    f"{other_data.get('running', 0)}/{other_data.get('total', 0)}"
                )

        except Exception as e:
            logging.debug(f"Error updating queue status: {e}")

    def _apply_active_jobs_filters_and_update(self) -> None:
        """Apply dropdown filters and current search to the Active Jobs table."""
        try:
            all_rows = getattr(self, "_active_jobs_all_rows", [])
            if not hasattr(self, "active_jobs_table"):
                return
            selected_type = "All"
            selected_worker = "All"
            try:
                if hasattr(self, "active_jobs_type_filter"):
                    selected_type = self.active_jobs_type_filter.value or "All"
                if hasattr(self, "active_jobs_worker_filter"):
                    selected_worker = self.active_jobs_worker_filter.value or "All"
            except Exception:
                pass

            def _keep(row: Dict[str, Any]) -> bool:
                if selected_type != "All" and row.get("job_type") != selected_type:
                    return False
                if selected_worker != "All" and row.get("worker") != selected_worker:
                    return False
                return True

            filtered = [r for r in all_rows if _keep(r)]
            self.active_jobs_table.rows = filtered
            self.active_jobs_table.update()
        except Exception:
            pass

    def _update_logs(self, data: Dict[str, Any]):
        """Update logs in the GUI."""
        try:
            if hasattr(self, "log_area"):
                log_message = data.get("message", "")
                log_level = data.get("level", "INFO")
                timestamp = time.strftime("%H:%M:%S")

                new_line = f"[{timestamp}] {log_level}: {log_message}\n"
                self._log_buffer.append(new_line)
                self.log_area.set_value("".join(self._log_buffer))

        except Exception as e:
            logging.debug(f"Error updating logs: {e}")

    def _update_progress(self, data: Dict[str, Any]):
        """Update progress in the GUI."""
        try:
            if hasattr(self, "progress_bar") and hasattr(self, "progress_label"):
                progress = data.get("progress", 0.0)

                pct = max(0.0, min(100.0, round(progress * 100.0, 1)))

                self.progress_bar.set_value(round(progress, 2))
                # Drop .0 for integers like 81.0 -> 81
                pct_str = f"{pct:.1f}" if pct % 1 else f"{int(pct)}"
                self.progress_label.set_text(f"{pct_str}% Complete")
                # Optional counts
                if hasattr(self, "completed_count") and "completed" in data:
                    self.completed_count.set_text(str(data["completed"]))
                if hasattr(self, "failed_count") and "failed" in data:
                    self.failed_count.set_text(str(data["failed"]))
                if hasattr(self, "total_count") and "total" in data:
                    self.total_count.set_text(str(data["total"]))

        except Exception as e:
            logging.debug(f"Error updating progress: {e}")

    def _update_errors(self, data: Dict[str, Any]):
        """Update error information in the GUI."""
        try:
            # Update error counts if available
            if hasattr(self, "preprocessing_errors"):
                self.preprocessing_errors.set_text(
                    str(data.get("preprocessing_errors", 0))
                )

            if hasattr(self, "analysis_errors"):
                self.analysis_errors.set_text(str(data.get("analysis_errors", 0)))

            if hasattr(self, "classification_errors"):
                self.classification_errors.set_text(
                    str(data.get("classification_errors", 0))
                )

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
            @ui.page("/")
            def welcome_page():
                """Welcome page at root route."""
                self._create_welcome_page()

            # Create the workflow monitoring page
            @ui.page("/littlejohn")
            def workflow_monitor():
                """Workflow monitoring page under /littlejohn route."""
                self._create_workflow_monitor()

            # Create the samples overview page
            @ui.page("/live_data")
            def samples_overview():
                """Samples overview page showing all tracked samples."""
                self._create_samples_overview()

            # Create individual sample detail pages
            @ui.page("/live_data/{sample_id}")
            def sample_detail(sample_id: str):
                """Individual sample detail page."""
                self._create_sample_detail_page(sample_id)

            # Enable global update processing regardless of which page is open
            try:
                self.gui_ready.set()
                ui.timer(0.3, self._drain_updates_on_ui, active=True)
                # Seed pre-existing samples shortly after startup, then poll for new ones
                ui.timer(
                    0.5,
                    lambda: self._scan_and_seed_samples(preexisting=True),
                    once=True,
                )
                ui.timer(3.0, self._scan_for_new_samples, active=True)
            except Exception:
                pass
            ui.add_css(
                """
                .shadows-into light-regular {
                    font-family: "Shadows Into Light", cursive;
                    font-weight: 800;
                    font-style: normal;
                }
            """
            )
            # Register fonts from the GUI package if available
            try:
                fonts_dir = Path(__file__).parent / "gui" / "fonts"
                if fonts_dir.exists():
                    app.add_static_files("/fonts", str(fonts_dir))
                else:
                    logging.debug(f"Fonts directory not found: {fonts_dir}")
            except Exception as e:
                logging.debug(f"Could not register fonts static dir: {e}")


            try:
                iconfile = os.path.join(
                    os.path.dirname(os.path.abspath(images.__file__)), "favicon.ico"
                )
                if not os.path.exists(iconfile):
                    logging.warning(f"Favicon file not found: {iconfile}")
            except Exception as e:
                logging.error(f"Error locating favicon: {str(e)}")
                iconfile = None
            # Start the GUI
            ui.run(
                host=self.host,
                port=self.port,
                show=False,
                reload=False,
                title="ROBIN",
                storage_secret="robin",
                favicon=iconfile,
            )
        except Exception as e:
            print(f"❌ GUI worker error: {e}")
            import traceback

            traceback.print_exc()

    def _create_welcome_page(self):
        """Create the welcome page."""
        # Background and main container
        with theme.frame(
            "Welcome",
            smalltitle="Welcome",
            batphone=False,
        ):
            with ui.card().classes("w-full shadow-lg rounded-xl"):
                with ui.row().classes(
                    "w-full items-center justify-between p-4 border-b border-gray-200"
                ):
                    with ui.row().classes("items-center gap-2"):
                        # ui.icon('feed', color='primary').classes('text-xl')
                        ui.label("Welcome to R.O.B.I.N").classes(
                            "text-sky-600 dark:text-white"
                        ).style("font-size: 150%; font-weight: 300").tailwind(
                            "drop-shadow", "font-bold"
                        )
                with ui.row().classes(
                    "w-full items-center justify-between p-4 border-b border-gray-200"
                ):
                    ui.label(
                        "This tool enables classification of brain tumours in real time from Oxford Nanopore Data."
                    ).classes("text-black-600 dark:text-white").style(
                        "font-size: 100%; font-weight: 300"
                    ).tailwind(
                        "drop-shadow", "font-bold"
                    )

                    # Description card
                    with ui.card().classes("w-full bg-white shadow-lg"):
                        with ui.column().classes("p-6"):
                            ui.label("What is R.O.B.I.N?").classes(
                                "text-xl font-semibold mb-4"
                            )
                            ui.label(
                                "ROBIN is a comprehensive bioinformatics workflow system designed for processing and analyzing BAM files. It provides automated preprocessing, multiple analysis pipelines, and real-time monitoring capabilities. It now encorporates Little John to help with the heavy lifting."
                            ).classes("text-gray-700 mb-4")

                            with ui.row().classes("w-full justify-center gap-8 mt-6"):
                                ui.link("�� View All Samples", "/live_data").classes(
                                    "bg-blue-600 hover:bg-blue-700 text-white text-lg font-semibold px-8 py-4 rounded-lg shadow-lg transition-colors"
                                )
                                ui.link(
                                    "📊 Open Workflow Monitor", "/littlejohn"
                                ).classes(
                                    "bg-blue-600 hover:bg-blue-700 text-white text-lg font-semibold px-8 py-4 rounded-lg shadow-lg transition-colors"
                                )

                                # Placeholder buttons for future functionality
                                ui.button(
                                    "🚀 Launch New Workflow",
                                    on_click=lambda: self._launch_workflow_button_clicked(),
                                ).classes(
                                    "bg-green-600 hover:bg-green-700 text-white text-lg font-semibold px-8 py-4 rounded-lg shadow-lg transition-colors"
                                )
                                ui.link("📋 View Documentation", "https://looselab.github.io/ROBIN/").classes(
                                    "bg-blue-600 hover:bg-blue-700 text-white text-lg font-semibold px-8 py-4 rounded-lg shadow-lg transition-colors"
                                )
                                
                            with ui.column().classes("w-full justify-center gap-8 mt-6"):
                                # Initialize news feed only if it hasn't been initialized yet
                                if self.news_feed is None:
                                    self.news_feed = NewsFeed()
                                    self.news_feed.start_update_timer()
                                # Create the news element
                                self.news_feed.create_news_element()

                            with ui.row().classes("w-full justify-center gap-8 mt-6"):
                                with ui.column().classes("text-center"):
                                    ui.label("🔬").classes("text-3xl mb-2")
                                    ui.label("Preprocessing").classes(
                                        "text-sm font-medium text-gray-600"
                                    )
                                with ui.column().classes("text-center"):
                                    ui.label("🧬").classes("text-3xl mb-2")
                                    ui.label("Methylation Classification").classes(
                                        "text-sm font-medium text-gray-600"
                                    )
                                with ui.column().classes("text-center"):
                                    ui.label("🧬").classes("text-3xl mb-2")
                                    ui.label("MGMT Analysis").classes(
                                        "text-sm font-medium text-gray-600"
                                    )
                                with ui.column().classes("text-center"):
                                    ui.label("📊").classes("text-3xl mb-2")
                                    ui.label("CNV Detection").classes(
                                        "text-sm font-medium text-gray-600"
                                    )
                                with ui.column().classes("text-center"):
                                    ui.label("🎯").classes("text-3xl mb-2")
                                    ui.label("Target Analysis").classes(
                                        "text-sm font-medium text-gray-600"
                                    )
                                with ui.column().classes("text-center"):
                                    ui.label("🔗").classes("text-3xl mb-2")
                                    ui.label("Fusion Detection").classes(
                                        "text-sm font-medium text-gray-600"
                                    )

                    # Action buttons

            # ToDo: Reimplement this.
            """
            with ui.row().classes("w-full no-wrap"):
                with ui.column().classes("w-1/4"):
                    with ui.card().classes("w-full shadow-lg rounded-xl"):
                        with ui.row().classes(
                            "w-full items-center justify-between p-4 border-b border-gray-200"
                        ):
                            with ui.row().classes("items-center gap-2"):
                                # ui.icon('feed', color='primary').classes('text-xl')
                                ui.label("View Data").classes(
                                    "text-sky-600 dark:text-white"
                                ).style("font-size: 150%; font-weight: 300").tailwind(
                                    "drop-shadow", "font-bold"
                                )
                        with ui.button(on_click=lambda: ui.navigate.to("/live")).props(
                            "color=green"
                        ):
                            ui.label("View Live Data")
                            ui.image(
                                os.path.join(
                                    os.path.dirname(os.path.abspath(images.__file__)),
                                    "ROBIN_logo_small.png",
                                )
                            ).classes("rounded-full w-16 h-16 ml-4")
                        with ui.button(
                            on_click=lambda: ui.navigate.to("/browse")
                        ).props("color=green"):
                            ui.label("Browse Historic Data")
                            ui.image(
                                os.path.join(
                                    os.path.dirname(os.path.abspath(images.__file__)),
                                    "ROBIN_logo_small.png",
                                )
                            ).classes("rounded-full w-16 h-16 ml-4")
                
                with ui.column().classes("w-2/4"):
                    # Initialize news feed only if it hasn't been initialized yet
                    if self.news_feed is None:
                        self.news_feed = NewsFeed()
                        self.news_feed.start_update_timer()
                    # Create the news element
                    self.news_feed.create_news_element()
                with ui.column().classes("w-1/4"):
                    # Initialize telemetry if enabled
                    if self.telemetry:
                        logging.info("Adding telemetry map to splash screen")
                        self.telemetry.create_map_element()
                    else:
                        logging.warning(
                            "No telemetry instance available for map display"
                        )
            """

    def _create_samples_overview(self):
        """Create the samples overview page showing all tracked samples."""
        # Page title and navigation
        with theme.frame(
            "R.O.B.I.N - Sample Tracking Overview",
            smalltitle="Samples",
            batphone=False,
        ):
            # Main content area
            with ui.column().classes("w-full p-4 gap-4"):
                # Sample statistics
                with ui.card().classes(
                    "w-full bg-gradient-to-r from-blue-50 to-indigo-50"
                ):
                    ui.label("📊 Sample Statistics").classes(
                        "text-lg font-semibold mb-4 text-blue-800"
                    )

                    # SIMPLIFIED: Show basic info without workflow_state access
                    ui.label(
                        "Sample statistics will be available once the workflow is running."
                    ).classes("text-sm text-gray-600")

                # Samples table
                with ui.card().classes("w-full"):
                    ui.label("📋 All Tracked Samples").classes(
                        "text-lg font-semibold mb-4"
                    )
                    
                    # Filters row
                    with ui.row().classes("items-center gap-2 mb-2"):

                        # Initialize filters model
                        self._samples_filters = getattr(
                            self, "_samples_filters", None
                        ) or {"query": "", "origin": "All"}

                        # Global search box
                        self.samples_search = (
                            ui.input(placeholder="Search…")
                            .props("clearable dense")
                            .on(
                                "update:model-value",
                                lambda e: self._set_samples_query(
                                    str(getattr(e, "value", "") or "")
                                ),
                            )
                            .classes("ml-auto")
                        )

                        # Origin filter
                        self.origin_filter = (
                            ui.select(
                                options=["All", "Live", "Pre-existing", "Complete"],
                                value=self._samples_filters.get("origin", "All"),
                                label="Origin",
                            )
                            .props("dense clearable")
                            .on(
                                "update:model-value",
                                lambda e: self._set_samples_origin_filter(
                                    e.args['label'] if e.args and isinstance(e.args, dict) and 'label' in e.args else "All"
                                ),
                            )
                        )

                    # Create a placeholder table that will be updated later
                    from littlejohn.gui.theme import styled_table
                    _samples_container, self.samples_table = styled_table(
                            columns=[
                                {
                                    "name": "actions",
                                    "label": "Actions",
                                    "field": "actions",
                                },
                                {
                                    "name": "sample_id",
                                    "label": "Sample ID",
                                    "field": "sample_id",
                                    "sortable": True,
                                },
                                {
                                    "name": "origin",
                                    "label": "Origin",
                                    "field": "origin",
                                    "sortable": True,
                                },
                                {
                                    "name": "run_start",
                                    "label": "Run Start",
                                    "field": "run_start",
                                    "sortable": True,
                                },
                                {
                                    "name": "device",
                                    "label": "Device",
                                    "field": "device",
                                    "sortable": True,
                                },
                                {
                                    "name": "flowcell",
                                    "label": "Flowcell",
                                    "field": "flowcell",
                                    "sortable": True,
                                },
                                {
                                    "name": "active_jobs",
                                    "label": "Active",
                                    "field": "active_jobs",
                                    "sortable": True,
                                },
                                {
                                    "name": "total_jobs",
                                    "label": "Total",
                                    "field": "total_jobs",
                                    "sortable": True,
                                },
                                {
                                    "name": "completed_jobs",
                                    "label": "Completed",
                                    "field": "completed_jobs",
                                    "sortable": True,
                                },
                                {
                                    "name": "failed_jobs",
                                    "label": "Failed",
                                    "field": "failed_jobs",
                                    "sortable": True,
                                },
                                {
                                    "name": "job_types",
                                    "label": "Job Types",
                                    "field": "job_types",
                                    "sortable": True,
                                },
                                {
                                    "name": "last_seen",
                                    "label": "Last Activity",
                                    "field": "last_seen",
                                    "sortable": True,
                                },
                                {
                                    "name": "export",
                                    "label": "Export",
                                    "field": "export",
                                },
                            ],
                            rows=[],
                            pagination=20,
                            class_size="table-xs",
                            row_key="sample_id",
                        )
                    try:
                        self.samples_table.props("rows-per-page-options=[10,20,50,0]")
                    except Exception:
                        pass

                    # Per-row action button to view sample
                    try:
                        self.samples_table.add_slot(
                            "body-cell-actions",
                            """
<q-td key=\"actions\" :props=\"props\">
  <q-btn color=\"primary\" size=\"sm\" label=\"View\"
         @click=\"$parent.$emit('action', props.row.sample_id)\" />
</q-td>
""",
                        )
                    except Exception:
                        pass

                    # Add export checkbox as rightmost column
                    try:
                        self.samples_table.add_slot(
                            "body-cell-export",
                            """
<q-td key=\"export\" :props=\"props\">
  <q-checkbox size=\"sm\"
              :model-value=\"props.row.export === true\"
              @update:model-value=\"$parent.$emit('export-toggled', { id: props.row.sample_id, value: $event })\" />
</q-td>
""",
                        )
                    except Exception:
                        pass

                    # Handle action emitted from table slot
                    try:
                        self.samples_table.on(
                            "action",
                            lambda e: ui.navigate.to(f"/live_data/{e.args}")
                            if isinstance(getattr(e, "args", None), str)
                            else ui.notify("Invalid button payload", type="warning"),
                        )
                    except Exception:
                        pass

                    # Track multi-selection for batch export via custom checkbox column
                    try:
                        self._selected_sample_ids = set()

                        def _on_export_toggled(event):
                            try:
                                payload = None
                                if hasattr(event, "args"):
                                    payload = getattr(event, "args", None)
                                elif isinstance(event, dict):
                                    payload = event
                                if isinstance(payload, dict):
                                    sid = payload.get("id")
                                    val = bool(payload.get("value"))
                                    if sid:
                                        if val:
                                            self._selected_sample_ids.add(sid)
                                        else:
                                            self._selected_sample_ids.discard(sid)
                                        # reflect state back into rows
                                        try:
                                            for r in (self.samples_table.rows or []):
                                                if r.get("sample_id") == sid:
                                                    r["export"] = sid in self._selected_sample_ids
                                            self.samples_table.update()
                                        except Exception:
                                            pass
                                        if self._selected_sample_ids:
                                            self.export_reports_button.enable()
                                        else:
                                            self.export_reports_button.disable()
                            except Exception:
                                pass

                        self.samples_table.on("export-toggled", _on_export_toggled)
                    except Exception:
                        pass

                    # Export selected reports button moved to the right of the table
                    with ui.row().classes("w-full justify-end mt-2"):
                        self.export_reports_button = ui.button(
                            "Export Reports",
                            on_click=lambda: None,
                        ).props("color=primary")
                        self.export_reports_button.disable()

                    # Batch export selected reports
                    try:
                        async def _export_selected_reports(state: Dict[str, Any]):
                            try:
                                selected = list(getattr(self, "_selected_sample_ids", set()) or [])
                                if not selected:
                                    ui.notify("No samples selected", type="warning")
                                    return
                                try:
                                    from nicegui import run as ng_run  # type: ignore
                                except Exception:
                                    ng_run = None  # type: ignore
                                for sid in selected:
                                    try:
                                        sample_dir = (
                                            Path(self.monitored_directory) / sid
                                            if self.monitored_directory
                                            else None
                                        )
                                        if not sample_dir or not sample_dir.exists():
                                            ui.notify(f"Missing output for {sid}", type="warning")
                                            continue
                                        filename = f"{sid}_run_report.pdf"
                                        pdf_path = os.path.join(str(sample_dir), filename)
                                        os.makedirs(str(sample_dir), exist_ok=True)
                                        export_csv_dir = None
                                        if bool(state.get("export_csv", False)):
                                            export_csv_dir = os.path.join(str(sample_dir), "report_csv")
                                        if ng_run is not None:
                                            pdf_file = await ng_run.io_bound(
                                                create_pdf,
                                                pdf_path,
                                                str(sample_dir),
                                                state.get("type", "detailed"),
                                                export_csv_dir=export_csv_dir,
                                                export_xlsx=False,
                                                export_zip=bool(state.get("export_csv", False)),
                                            )
                                        else:
                                            pdf_file = create_pdf(
                                                pdf_path,
                                                str(sample_dir),
                                                state.get("type", "detailed"),
                                                export_csv_dir=export_csv_dir,
                                                export_xlsx=False,
                                                export_zip=bool(state.get("export_csv", False)),
                                            )
                                        ui.download(pdf_file)
                                        # Also offer CSV ZIP if requested
                                        if bool(state.get("export_csv", False)) and export_csv_dir:
                                            zip_path = os.path.join(export_csv_dir, f"{sid}_report_data.zip")
                                            if os.path.exists(zip_path):
                                                ui.download(zip_path)
                                    except Exception as e:
                                        ui.notify(f"Export failed for {sid}: {e}", type="error")
                                ui.notify("Export complete")
                            except Exception:
                                pass

                        async def _confirm_bulk_export():
                            # Ensure there is at least one selection before opening
                            if not getattr(self, "_selected_sample_ids", None):
                                ui.notify("No samples selected", type="warning")
                                return

                            report_types = {
                                "summary": "Summary Only",
                                "detailed": "Detailed",
                            }
                            state: Dict[str, Any] = {"type": "detailed", "export_csv": False}

                            with ui.dialog() as dialog, ui.card().classes("w-96 p-4"):
                                ui.label("Generate Report").classes("text-h6 font-bold mb-4")

                                # Report type selector
                                with ui.column().classes("mb-4"):
                                    ui.label("Report Type").classes("font-bold mb-2")
                                    ui.toggle(
                                        report_types,
                                        value="detailed",
                                        on_change=lambda e: state.update({"type": e.value}),
                                    )

                                # Data export options
                                with ui.column().classes("mb-4"):
                                    ui.label("Include Data").classes("font-bold mb-2")
                                    ui.checkbox(
                                        "CSV data (ZIP)",
                                        value=False,
                                        on_change=lambda e: state.update({"export_csv": bool(e.value)}),
                                    )

                                # Disclaimer section (match single-report dialog)
                                with ui.column().classes("mb-4"):
                                    ui.label("Disclaimer").classes("font-bold mb-2")
                                    formatted_text = EXTENDED_DISCLAIMER_TEXT.replace("\n\n", "<br><br>").replace("\n", " ")
                                    ui.label(formatted_text).classes("text-sm text-gray-600 mb-4")

                                ui.label("Are you sure you want to generate a report?").classes("mb-4")

                                # Buttons
                                with ui.row().classes("justify-end gap-2"):
                                    ui.button("Cancel", on_click=lambda: dialog.submit(("No", None))).props("flat")
                                    ui.button(
                                        "Export",
                                        on_click=lambda: dialog.submit(("Yes", state)),
                                    ).props("color=primary")

                            dialog_result = await dialog
                            if dialog_result is None:
                                return
                            result, state_val = dialog_result
                            if result == "Yes" and isinstance(state_val, dict):
                                await _export_selected_reports(state_val)

                        # Wire the button now that handlers exist
                        self.export_reports_button.on_click(_confirm_bulk_export)
                    except Exception:
                        pass

                    # If we have cached rows, show them immediately for instant UX
                    if self._last_samples_rows:
                        try:
                            self._apply_samples_table_filters()
                            # No selection behavior needed; per-row buttons handle navigation
                        except Exception:
                            pass
                    


    def _update_samples_table(self, data: Dict[str, Any]):
        """Update the samples overview table with new data."""
        try:
            if not hasattr(self, "samples_table"):
                return
            samples = data.get("samples", [])
            
            # Deduplicate by sample_id taking the newest last_seen
            by_id: Dict[str, Dict[str, Any]] = {}
            for s in samples:
                sid = s.get("sample_id", "") or "unknown"
                last_seen = float(s.get("last_seen", time.time()))
                existing = by_id.get(sid)
                if not existing or last_seen >= existing.get("_last_seen_raw", 0):
                    origin_value = (
                        "Pre-existing"
                        if sid in self._preexisting_sample_ids
                        else "Live"
                    )
                    # Flip Live samples to Complete if inactive for 60 minutes
                    try:
                        if origin_value == "Live" and (time.time() - last_seen) >= 3600:
                            origin_value = "Complete"
                    except Exception:
                        pass
                    by_id[sid] = {
                        "sample_id": sid,
                        "origin": origin_value,
                        # Persisted run info will be patched in below from master.csv
                        "run_start": "",
                        "device": "",
                        "flowcell": "",
                        "active_jobs": s.get("active_jobs", 0),
                        "total_jobs": s.get("total_jobs", 0),
                        "completed_jobs": s.get("completed_jobs", 0),
                        "failed_jobs": s.get("failed_jobs", 0),
                        "job_types": (
                            ",".join(sorted(set(s.get("job_types", []))))
                            if isinstance(s.get("job_types", []), list)
                            else str(s.get("job_types", ""))
                        ),
                        "last_seen": time.strftime(
                            "%H:%M:%S", time.localtime(last_seen)
                        ),
                        "actions": "View",
                        "_last_seen_raw": last_seen,
                    }

            # Patch from master.csv and persist the new overview values for later reload
            try:
                base = Path(self.monitored_directory) if self.monitored_directory else None
                manager = (
                    MasterCSVManager(str(base)) if base and base.exists() else None
                )
            except Exception:
                base, manager = None, None

            for sid, row in by_id.items():
                try:
                    # Persist overview numbers to master.csv so we can restore later
                    if manager is not None:
                        persist_payload = {
                            "active_jobs": int(row.get("active_jobs", 0)),
                            "total_jobs": int(row.get("total_jobs", 0)),
                            "completed_jobs": int(row.get("completed_jobs", 0)),
                            "failed_jobs": int(row.get("failed_jobs", 0)),
                            "job_types": row.get("job_types", ""),
                            "last_seen": float(row.get("_last_seen_raw", time.time())),
                        }
                        manager.update_sample_overview(sid, persist_payload)

                    # Read run info from master.csv to display in table
                    if base is not None:
                        csv_path = base / sid / "master.csv"
                        if csv_path.exists():
                            with csv_path.open("r", newline="") as fh:
                                reader = csv.DictReader(fh)
                                first_row = next(reader, None)
                            if first_row:
                                row["run_start"] = self._format_timestamp_for_display(
                                    first_row.get("run_info_run_time", "")
                                )
                                row["device"] = first_row.get("run_info_device", "")
                                row["flowcell"] = first_row.get("run_info_flow_cell", "")
                except Exception:
                    pass

            # Merge with preexisting scans if any
            existing_rows_by_id = {
                r["sample_id"]: r for r in (self._last_samples_rows or [])
            }
            for sid, row in by_id.items():
                existing_rows_by_id[sid] = row
            rows = list(existing_rows_by_id.values())
            # Replace rows to avoid duplicates then apply filters
            self._last_samples_rows = rows
            self._apply_samples_table_filters()
            
            # Update known IDs from unfiltered cache
            self._known_sample_ids = {
                r.get("sample_id")
                for r in (self._last_samples_rows or [])
                if r.get("sample_id") and r.get("sample_id") != "unknown"
            }
            # Track currently most active/recent sample (for convenience elsewhere)
            visible_rows = self.samples_table.rows or []
            if visible_rows:
                rows_sorted = sorted(
                    visible_rows,
                    key=lambda r: (r.get("active_jobs", 0), r.get("_last_seen_raw", 0)),
                    reverse=True,
                )
                self._current_sample_id = rows_sorted[0].get("sample_id")
        except Exception as e:
            logging.debug(f"Error updating samples table: {e}")

    # -------- Samples table helpers: search, filter, sort --------
    def _format_timestamp_for_display(self, value: Any) -> str:
        """Convert timestamps (ISO 8601 or epoch) to 'YYYY-MM-DD HH:MM'."""
        try:
            if value is None:
                return ""
            if isinstance(value, (int, float)):
                return datetime.fromtimestamp(float(value)).strftime("%Y-%m-%d %H:%M")
            s = str(value).strip()
            if not s:
                return ""
            if s.endswith("Z"):
                s = s[:-1] + "+00:00"
            try:
                return datetime.fromisoformat(s).strftime("%Y-%m-%d %H:%M")
            except Exception:
                try:
                    return datetime.fromtimestamp(float(s)).strftime("%Y-%m-%d %H:%M")
                except Exception:
                    return s
        except Exception:
            try:
                return str(value)
            except Exception:
                return ""
                
    def _normalize_rows_for_display(
        self, rows: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        normalized: List[Dict[str, Any]] = []
        for r in rows or []:
            sid = r.get("sample_id")
            if not sid or sid == "unknown":
                continue
            jt = r.get("job_types")
            if isinstance(jt, set):
                r = dict(r)
                r["job_types"] = ",".join(sorted(jt))
            normalized.append(r)
        return normalized

    def _set_samples_query(self, query: str) -> None:
        try:
            self._samples_filters["query"] = (query or "").strip().lower()
            self._apply_samples_table_filters()
        except Exception:
            pass

    def _set_samples_origin_filter(self, origin_value: str) -> None:
        try:
            self._samples_filters["origin"] = origin_value or "All"
            self._apply_samples_table_filters()
        except Exception:
            pass

    def _apply_samples_table_filters(self) -> None:
        try:
            base_rows = getattr(self, "_last_samples_rows", []) or []
            rows = self._normalize_rows_for_display(base_rows)

            # Origin filter, compute dynamic 'Complete' for display if needed
            now_ts = time.time()
            for r in rows:
                try:
                    if r.get("origin") == "Live":
                        last_raw = float(r.get("_last_seen_raw", 0))
                        if last_raw and (now_ts - last_raw) >= 3600:
                            r["origin"] = "Complete"
                except Exception:
                    pass
            
            origin = (self._samples_filters or {}).get("origin", "All")
            if origin and origin != "All":
                rows = [r for r in rows if (r.get("origin") == origin)]

            # Global query filter
            q = (self._samples_filters or {}).get("query", "")
            if q:
                ql = q.lower()

                def match_any(r: Dict[str, Any]) -> bool:
                    return any(
                        (str(r.get(k, "")).lower().find(ql) >= 0)
                        for k in ["sample_id", "origin", "job_types", "last_seen"]
                    ) or any(
                        (str(r.get(k, 0)).lower().find(ql) >= 0)
                        for k in [
                            "active_jobs",
                            "total_jobs",
                            "completed_jobs",
                            "failed_jobs",
                        ]
                    )

                rows = [r for r in rows if match_any(r)]

            # annotate export selection state per row for rightmost checkbox column
            try:
                selected = getattr(self, "_selected_sample_ids", set()) or set()
            except Exception:
                selected = set()
            for r in rows:
                try:
                    r["export"] = r.get("sample_id") in selected
                except Exception:
                    r["export"] = False

            self.samples_table.rows = rows
            self.samples_table.update()
        except Exception:
            pass

    def _create_sample_detail_page(self, sample_id: str):
        """Create the individual sample detail page."""
        with theme.frame(
            f"R.O.B.I.N - Sample {sample_id}",
            smalltitle="Samples",
            batphone=False,
        ):
            # Guard: unknown sample -> show message and back button; also redirect
            if self._known_sample_ids and sample_id not in self._known_sample_ids:
                with ui.column().classes("w-full items-center justify-center p-8"):
                    ui.label(f"Unknown sample: {sample_id}").classes(
                        "text-xl font-semibold text-red-600"
                    )
                    ui.label(
                        "This sample ID has not been seen yet in the current session."
                    ).classes("text-sm text-gray-600")
                    ui.button(
                        "Back to Samples", on_click=lambda: ui.navigate.to("/live_data")
                    ).props("color=primary")
                # Soft redirect after short delay
                try:
                    ui.timer(1.5, lambda: ui.navigate.to("/live_data"), once=True)
                except Exception:
                    pass
                return

            sample_dir = (
                Path(self.monitored_directory) / sample_id
                if self.monitored_directory
                else None
            )

            # Debug visibility: entering sample page and directory availability
            try:
                ui.notify(f"Opening sample {sample_id}", type="info")
            except Exception:
                pass
            if not sample_dir or not sample_dir.exists():
                try:
                    ui.notify(
                        f"Sample output directory not found for {sample_id}",
                        type="warning",
                    )
                except Exception:
                    pass

            # Page title and navigation
            with ui.row().classes("w-full"):  #'w-full bg-blue-600 text-white p-4 items-center justify-between'):
                with ui.row().classes("items-center"):
                    ui.label(f"{sample_id}").classes("text-2xl font-bold")
                    ui.label("Detailed sample information.").classes(
                        "text-sm ml-4 opacity-80"
                    )
                with ui.row().classes("w-full"):
                    # ui.link('📋 Sample Overview', '/live_data').classes('text-white hover:text-blue-200 text-sm')
                    # ui.link('📊 Workflow Monitor', '/littlejohn').classes('text-white hover:text-blue-200 text-sm')
                    # Report generation (confirmation dialog + download)
                    async def confirm_report_generation():
                        """Show a confirmation dialog before generating the report."""
                        report_types = {
                            "summary": "Summary Only",
                            "detailed": "Detailed",
                        }
                        state: Dict[str, Any] = {"type": "detailed", "export_csv": False}

                        with ui.dialog() as dialog, ui.card().classes("w-96 p-4"):
                            ui.label("Generate Report").classes(
                                "text-h6 font-bold mb-4"
                            )

                            # Report type selector
                            with ui.column().classes("mb-4"):
                                ui.label("Report Type").classes("font-bold mb-2")
                                ui.toggle(
                                    report_types,
                                    value="detailed",
                                    on_change=lambda e: state.update({"type": e.value}),
                                )

                            # Disclaimer section
                            with ui.column().classes("mb-4"):
                                ui.label("Disclaimer").classes("font-bold mb-2")
                                formatted_text = EXTENDED_DISCLAIMER_TEXT.replace(
                                    "\n\n", "<br><br>"
                                ).replace("\n", " ")
                                ui.label(formatted_text).classes(
                                    "text-sm text-gray-600 mb-4"
                                )

                            # Data export options
                            with ui.column().classes("mb-4"):
                                ui.label("Include Data").classes("font-bold mb-2")
                                ui.checkbox(
                                    "CSV data (ZIP)",
                                    value=False,
                                    on_change=lambda e: state.update({"export_csv": bool(e.value)}),
                                )

                            ui.label(
                                "Are you sure you want to generate a report?"
                            ).classes("mb-4")

                            # Buttons
                            with ui.row().classes("justify-end gap-2"):
                                ui.button(
                                    "No", on_click=lambda: dialog.submit(("No", None))
                                ).props("flat")
                                ui.button(
                                    "Yes",
                                    on_click=lambda: dialog.submit(("Yes", state)),
                                ).props("color=primary")

                        dialog_result = await dialog
                        if dialog_result is None:
                            return
                        result, state_val = dialog_result
                        if result == "Yes" and isinstance(state_val, dict):
                            await download_report(state_val)

                    async def download_report(state: Dict[str, Any]):
                        """Generate and download the report for this sample."""
                        try:
                            # Import here to avoid global dependency if GUI isn't used
                            from nicegui import run as ng_run  # type: ignore

                            if not sample_dir or not sample_dir.exists():
                                ui.notify(
                                    "Output directory not available for this sample",
                                    type="warning",
                                )
                                return
                            ui.notify("Generating report…")
                            filename = f"{sample_id}_run_report.pdf"
                            pdf_path = os.path.join(str(sample_dir), filename)
                            os.makedirs(str(sample_dir), exist_ok=True)
                            export_csv_dir = None
                            if bool(state.get("export_csv", False)):
                                export_csv_dir = os.path.join(str(sample_dir), "report_csv")
                            pdf_file = await ng_run.io_bound(
                                create_pdf,
                                pdf_path,
                                str(sample_dir),
                                state.get("type", "detailed"),
                                export_csv_dir=export_csv_dir,
                                export_xlsx=False,
                                export_zip=bool(state.get("export_csv", False)),
                            )
                            ui.download(pdf_file)
                            # Also offer CSV ZIP if requested
                            if bool(state.get("export_csv", False)) and export_csv_dir:
                                zip_path = os.path.join(export_csv_dir, f"{sample_id}_report_data.zip")
                                if os.path.exists(zip_path):
                                    ui.download(zip_path)
                            ui.notify("Report downloaded")
                        except Exception as e:
                            ui.notify(f"Report generation failed: {e}", type="error")

                    ui.button(
                        "Generate Report", on_click=confirm_report_generation
                    ).classes("object-right ml-auto text-sm font-semibold px-3 py-1 rounded")

            # Main content area with loading state
            with ui.column().classes("w-full p-4 gap-4"):
                # Loading container that will be hidden when data is ready
                loading_container = ui.column().classes("w-full items-center justify-center p-8")
                with loading_container:
                    ui.spinner('bars', size='4em').classes("mb-4")
                    ui.label("Loading sample data...").classes("text-lg text-gray-600")
                    ui.label("This may take a moment for large datasets").classes("text-sm text-gray-500")
                
                # Content container that will be shown when data is ready
                content_container = ui.column().classes("w-full gap-4").style("display: none")
                
                with content_container:
                    # Summary section (new component)
                    try:
                        from .gui.components.summary import add_summary_section  # type: ignore

                        add_summary_section(sample_dir, sample_id)
                    except Exception as e:
                        logging.exception(f"[GUI] Summary section failed: {e}")
                        try:
                            ui.notify(f"Summary section failed: {e}", type="warning")
                        except Exception:
                            pass

                    # Classification section (refactored component)
                    try:
                        from .gui.components.classification import add_classification_section  # type: ignore

                        add_classification_section(sample_dir)
                    except Exception as e:
                        logging.exception(f"[GUI] Classification section failed: {e}")
                        try:
                            ui.notify(f"Classification section failed: {e}", type="warning")
                        except Exception:
                            pass

                    # Coverage section (refactored component)
                    try:
                        from .gui.components.coverage import add_coverage_section  # type: ignore

                        add_coverage_section(self, sample_dir)
                    except Exception as e:
                        try:
                            ui.notify(f"Coverage section failed: {e}", type="warning")
                        except Exception:
                            pass

                    # MGMT section (refactored component)
                    try:
                        from .gui.components.mgmt import add_mgmt_section  # type: ignore

                        add_mgmt_section(self, sample_dir)
                    except Exception as e:
                        logging.exception(f"[GUI] MGMT section failed: {e}")
                        try:
                            ui.notify(f"MGMT section failed: {e}", type="warning")
                        except Exception:
                            pass

                    # CNV section (refactored component)
                    try:
                        from .gui.components.cnv import add_cnv_section  # type: ignore

                        # Pass launcher for shared state access (launcher._cnv_state)
                        add_cnv_section(self, sample_dir)
                    except Exception as e:
                        logging.exception(f"[GUI] CNV section failed: {e}")
                        try:
                            ui.notify(f"CNV section failed: {e}", type="warning")
                        except Exception:
                            pass

                    # Fusion section (target and genome-wide; excludes full SV UI)
                    try:
                        from .gui.components.fusion import add_fusion_section  # type: ignore

                        add_fusion_section(self, sample_dir)
                    except Exception as e:
                        logging.exception(f"[GUI] Fusion section failed: {e}")
                        try:
                            ui.notify(f"Fusion section failed: {e}", type="warning")
                        except Exception:
                            pass

                    # Files in output directory
                    with ui.card().classes("w-full"):
                        ui.label("📁 Output Files").classes("text-lg font-semibold mb-2")
                        with ui.row().classes("items-center gap-3 mb-2"):
                            files_search = ui.input("Search files…").props(
                                "borderless dense clearable"
                            )
                        from littlejohn.gui.theme import styled_table
                        _files_container, files_table = styled_table(
                            columns=[
                                {
                                    "name": "name",
                                    "label": "File",
                                    "field": "name",
                                    "sortable": True,
                                },
                                {
                                    "name": "size",
                                    "label": "Size (bytes)",
                                    "field": "size",
                                    "sortable": True,
                                },
                                {
                                    "name": "mtime",
                                    "label": "Last Modified",
                                    "field": "mtime",
                                    "sortable": True,
                                },
                            ],
                            rows=[],
                            pagination=20,
                            class_size="table-xs",
                        )
                        try:
                            files_table.props(
                                'multi-sort rows-per-page-options="[10,20,50,0]"'
                            )
                            files_search.bind_value(files_table, "filter")
                        except Exception:
                            pass

                    # master.csv summary
                    with ui.card().classes("w-full"):
                        ui.label("📊 master.csv Summary").classes(
                            "text-lg font-semibold mb-2"
                        )
                        with ui.row().classes("items-center gap-3 mb-2"):
                            summary_search = ui.input("Search fields…").props(
                                "borderless dense clearable"
                            )
                        from littlejohn.gui.theme import styled_table
                        _summary_container, summary_table = styled_table(
                            columns=[
                                {
                                    "name": "key",
                                    "label": "Field",
                                    "field": "key",
                                    "sortable": True,
                                },
                                {
                                    "name": "value",
                                    "label": "Value",
                                    "field": "value",
                                    "sortable": True,
                                },
                            ],
                            rows=[],
                            pagination=0,
                            class_size="table-xs",
                        )
                        try:
                            summary_table.props("multi-sort")
                            summary_search.bind_value(summary_table, "filter")
                        except Exception:
                            pass

                # Periodic refresher for files table and master.csv summary
                _notify_state = {"files_error": False, "csv_error": False}

                def _refresh_sample_detail() -> None:
                    # Refresh files list
                    try:
                        rows = []
                        if sample_dir and sample_dir.exists():
                            for f in sorted(sample_dir.iterdir()):
                                if f.is_file():
                                    try:
                                        stat = f.stat()
                                        rows.append(
                                            {
                                                "name": f.name,
                                                "size": stat.st_size,
                                                "mtime": time.strftime(
                                                    "%Y-%m-%d %H:%M:%S",
                                                    time.localtime(stat.st_mtime),
                                                ),
                                            }
                                        )
                                    except Exception:
                                        continue
                        files_table.rows = rows
                        files_table.update()
                    except Exception as e:
                        if not _notify_state["files_error"]:
                            try:
                                ui.notify(
                                    f"Failed to list output files for {sample_id}: {e}",
                                    type="warning",
                                )
                            except Exception:
                                pass
                            _notify_state["files_error"] = True

                    # Refresh master.csv summary
                    try:
                        if sample_dir:
                            csv_path = sample_dir / "master.csv"
                            if csv_path.exists():
                                with csv_path.open("r", newline="") as fh:
                                    reader = csv.DictReader(fh)
                                    first_row = next(reader, None)
                                if first_row:
                                    preferred_keys = [
                                        "counter_bam_passed",
                                        "counter_bam_failed",
                                        "counter_bases_count",
                                        "counter_mapped_count",
                                        "counter_unmapped_count",
                                        "run_info_run_time",
                                        "run_info_device",
                                        "run_info_model",
                                        "run_info_flow_cell",
                                        "bam_tracking_counter",
                                        "bam_tracking_total_files",
                                    ]
                                    rows2 = []
                                    for k in preferred_keys:
                                        if k in first_row:
                                            rows2.append(
                                                {"key": k, "value": first_row.get(k, "")}
                                            )
                                    if not rows2:
                                        rows2 = [
                                            {"key": k, "value": v}
                                            for k, v in first_row.items()
                                        ]
                                    summary_table.rows = rows2
                                    summary_table.update()
                            else:
                                summary_table.rows = [
                                    {"key": "Status", "value": "master.csv not found"}
                                ]
                                summary_table.update()
                    except Exception as e:
                        if not _notify_state["csv_error"]:
                            try:
                                ui.notify(
                                    f"Failed to read master.csv for {sample_id}: {e}",
                                    type="warning",
                                )
                            except Exception:
                                pass
                            _notify_state["csv_error"] = True

                # Show content and hide loading after initial data load
                def _show_content():
                    loading_container.style("display: none")
                    content_container.style("display: flex")
                
                # Initial data load with loading state
                try:
                    # Use a timer to simulate async loading and then show content
                    ui.timer(0.1, _show_content, once=True)
                except Exception:
                    # Fallback: show content immediately if timer fails
                    _show_content()

                # Start periodic refresh
                try:
                    ui.timer(5.0, _refresh_sample_detail)
                except Exception:
                    pass

    def _create_workflow_monitor(self):
        """Create the main workflow monitoring page."""
        with theme.frame(
            "R.O.B.I.N - LittleJohn Workflow Monitor",
            smalltitle="Samples",
            batphone=False,
        ):
            # Page title and navigation
            with ui.row().classes(
                "w-full p-4 items-center justify-between"
            ):
                with ui.row().classes("items-center"):
                    ui.label("📊 LittleJohn Workflow Monitor").classes(
                        "text-2xl font-bold"
                    )
                    ui.label("Real-time workflow monitoring and control").classes(
                        "text-sm ml-4 opacity-80"
                    )


            # Main content area
            with ui.column().classes("w-full p-4 gap-4"):
                # Workflow Status Overview
                with ui.card().classes(
                    "w-full bg-gradient-to-r from-blue-50 to-indigo-50"
                ):
                    ui.label("🚀 Workflow Status Overview").classes(
                        "text-lg font-semibold mb-4 text-blue-800"
                    )

                    # Status indicator
                    with ui.row().classes("w-full items-center gap-4"):
                        self.status_indicator = ui.label("🟢").classes("text-2xl")
                        self.status_label = ui.label(
                            "Workflow Status: Running"
                        ).classes("text-sm font-medium text-green-600")

                    # Timing information
                    with ui.row().classes("w-full gap-8 mt-4"):
                        self.workflow_start_time = ui.label("Started: --").classes(
                            "text-sm text-gray-600"
                        )
                        self.workflow_duration = ui.label("Duration: --").classes(
                            "text-sm text-gray-600"
                        )

                    # Progress bar (hide internal float value text; use external formatted label below)
                    self.progress_bar = (
                        ui.linear_progress(0.0)
                        .classes("w-full mt-4")
                        .style("color: transparent")
                    )
                    self.progress_label = ui.label("0% Complete").classes(
                        "text-sm text-center text-gray-600"
                    )

                    # Counts summary
                    with ui.row().classes("w-full gap-6 mt-2"):
                        with ui.row().classes("items-center gap-2"):
                            ui.label("Completed:").classes("text-xs text-gray-600")
                            self.completed_count = ui.label("0").classes(
                                "text-xs font-semibold"
                            )
                        with ui.row().classes("items-center gap-2"):
                            ui.label("Failed:").classes("text-xs text-gray-600")
                            self.failed_count = ui.label("0").classes(
                                "text-xs font-semibold"
                            )
                        with ui.row().classes("items-center gap-2"):
                            ui.label("Total:").classes("text-xs text-gray-600")
                            self.total_count = ui.label("0").classes(
                                "text-xs font-semibold"
                            )

                # Queue Status
                with ui.card().classes("w-full"):
                    ui.label("📋 Queue Status").classes("text-lg font-semibold mb-4")

                    # Queue status grid
                    with ui.grid(columns=4).classes("w-full gap-4"):
                        # Preprocessing
                        with ui.card().classes("bg-green-50 p-4"):
                            ui.label("🔬 Preprocessing").classes(
                                "text-sm font-medium text-green-800"
                            )
                            self.preprocessing_status = ui.label("0/0").classes(
                                "text-2xl font-bold text-green-600"
                            )

                        # Analysis
                        with ui.card().classes("bg-blue-50 p-4"):
                            ui.label("🧬 Analysis").classes(
                                "text-sm font-medium text-blue-800"
                            )
                            self.analysis_status = ui.label("0/0").classes(
                                "text-2xl font-bold text-blue-600"
                            )

                        # Classification
                        with ui.card().classes("bg-purple-50 p-4"):
                            ui.label("🎯 Classification").classes(
                                "text-sm font-medium text-purple-800"
                            )
                            self.classification_status = ui.label("0/0").classes(
                                "text-2xl font-bold text-purple-800"
                            )

                        # Other
                        with ui.card().classes("bg-gray-50 p-4"):
                            ui.label("⚙️ Other").classes(
                                "text-sm font-medium text-gray-800"
                            )
                            self.other_status = ui.label("0/0").classes(
                                "text-2xl font-bold text-gray-600"
                            )

                    # If we have a cached queue status from before the page was created, apply it now
                    try:
                        if self._last_queue_status:
                            self._update_queue_status(self._last_queue_status)
                    except Exception:
                        pass

                # Active Jobs
                with ui.card().classes("w-full"):
                    ui.label("⚡ Active Jobs").classes("text-lg font-semibold mb-4")

                    with ui.row().classes("items-center gap-3 mb-2"):
                        self.active_jobs_search = ui.input("Search…").props(
                            "borderless dense clearable"
                        )
                        self.active_jobs_type_filter = (
                            ui.select(options=["All"], value="All", label="Type")
                            .props("dense clearable")
                            .classes("w-40")
                        )
                        self.active_jobs_worker_filter = (
                            ui.select(options=["All"], value="All", label="Worker")
                            .props("dense clearable")
                            .classes("w-40")
                        )

                    # Active jobs table
                    from littlejohn.gui.theme import styled_table
                    _jobs_container, self.active_jobs_table = styled_table(
                        columns=[
                            {
                                "name": "job_id",
                                "label": "Job ID",
                                "field": "job_id",
                                "sortable": True,
                            },
                            {
                                "name": "job_type",
                                "label": "Type",
                                "field": "job_type",
                                "sortable": True,
                            },
                            {
                                "name": "filepath",
                                "label": "File",
                                "field": "filepath",
                                "sortable": True,
                            },
                            {
                                "name": "worker",
                                "label": "Worker",
                                "field": "worker",
                                "sortable": True,
                            },
                            {
                                "name": "duration",
                                "label": "Duration",
                                "field": "duration",
                                "sortable": True,
                            },
                            {
                                "name": "progress",
                                "label": "Progress",
                                "field": "progress",
                                "sortable": True,
                            },
                        ],
                        rows=[],
                        pagination=10,
                        class_size="table-xs",
                    )
                    try:
                        self.active_jobs_table.props(
                            'multi-sort rows-per-page-options="[10,20,50,0]"'
                        )
                        self.active_jobs_search.bind_value(
                            self.active_jobs_table, "filter"
                        )
                    except Exception:
                        pass
                    try:

                        def _on_active_jobs_filter_change(_=None):
                            self._apply_active_jobs_filters_and_update()

                        self.active_jobs_type_filter.on(
                            "update:model-value", _on_active_jobs_filter_change
                        )
                        self.active_jobs_worker_filter.on(
                            "update:model-value", _on_active_jobs_filter_change
                        )
                    except Exception:
                        pass

                    # Placeholder for when no jobs are active
                    ui.label("No active jobs at the moment.").classes(
                        "text-sm text-gray-500 mt-2"
                    )

                # Live Logs
                with ui.card().classes("w-full"):
                    ui.label("📝 Live Logs").classes("text-lg font-semibold mb-4")

                    # Log controls
                    with ui.row().classes("w-full justify-between items-center mb-2"):
                        with ui.row().classes("gap-2"):
                            ui.button("Clear", on_click=self._clear_logs).classes(
                                "bg-gray-500 hover:bg-gray-600 text-white text-xs"
                            )
                            ui.button(
                                "Export", on_click=lambda: self._export_logs()
                            ).classes(
                                "bg-blue-500 hover:bg-blue-600 text-white text-xs"
                            )

                    # Log area
                    self.log_area = (
                        ui.textarea("Workflow logs will appear here...")
                        .classes("w-full h-40")
                        .props("readonly")
                    )

                # Configuration
                with ui.card().classes("w-full"):
                    ui.label("⚙️ Workflow Configuration").classes(
                        "text-lg font-semibold mb-4"
                    )

                    # Configuration details
                    with ui.grid(columns=2).classes("w-full gap-4"):
                        with ui.column():
                            ui.label("Monitored Directory:").classes(
                                "text-sm font-medium"
                            )
                            ui.label(
                                self.monitored_directory or "Not specified"
                            ).classes("text-sm text-gray-600")

                            ui.label("Workflow Steps:").classes(
                                "text-sm font-medium mt-2"
                            )
                            ui.label(
                                ", ".join(self.workflow_steps)
                                if self.workflow_steps
                                else "Not specified"
                            ).classes("text-sm text-gray-600")

                        with ui.column():
                            ui.label("Log Level:").classes("text-sm font-medium")
                            ui.label("--").classes("text-sm text-gray-600")

                            ui.label("Analysis Workers:").classes(
                                "text-sm font-medium mt-2"
                            )
                            ui.label("--").classes("text-sm text-gray-600")

                # Error Summary & Troubleshooting
                with ui.card().classes("w-full"):
                    ui.label("⚠️ Error Summary & Troubleshooting").classes(
                        "text-lg font-semibold mb-2"
                    )

                    # Error counts by type
                    with ui.row().classes("w-full justify-between"):
                        with ui.column().classes("text-center"):
                            self.preprocessing_errors = ui.label("0").classes(
                                "text-xl font-bold text-red-600"
                            )
                            ui.label("Preprocessing").classes("text-xs text-gray-600")
                        with ui.column().classes("text-center"):
                            self.analysis_errors = ui.label("0").classes(
                                "text-xl font-bold text-red-600"
                            )
                            ui.label("Analysis").classes("text-xs text-gray-600")
                        with ui.column().classes("text-center"):
                            self.classification_errors = ui.label("0").classes(
                                "text-xl font-bold text-red-600"
                            )
                            ui.label("Classification").classes("text-xs text-gray-600")

                    # Common error messages
                    ui.separator()
                    ui.label("Recent Errors").classes("text-sm font-medium mt-2")
                    self.error_summary_label = ui.label("No errors detected").classes(
                        "text-xs text-gray-600"
                    )

                # Footer
                with ui.row().classes("w-full bg-gray-200 p-2 justify-center"):
                    ui.label("LittleJohn Workflow Monitor - Running").classes(
                        "text-sm text-gray-600"
                    )

                # Signal that GUI is ready to receive updates
                self.gui_ready.set()
                logging.info("[GUI] UI created and ready to receive updates")
                # Start a periodic UI-thread drain of the update queue
                ui.timer(0.3, self._drain_updates_on_ui, active=True)
                # Start duration refresher
                ui.timer(1.0, lambda: self._refresh_duration(), active=True)

    def _refresh_duration(self):
        if self._is_running and self._start_time and hasattr(self, "workflow_duration"):
            try:
                elapsed_seconds = int(time.time() - self._start_time)
                self.workflow_duration.set_text(
                    f"Duration: {self._format_duration(elapsed_seconds)}"
                )
            except Exception:
                pass

    def _export_logs(self):
        """Export logs to a file."""
        try:
            # Simple log export functionality
            log_content = "".join(self._log_buffer)
            if log_content:
                # Create a simple download
                ui.download(log_content, "workflow_logs.txt")
            else:
                ui.notify("No logs to export", type="warning")
        except Exception as e:
            ui.notify(f"Export failed: {e}", type="error")

    def _clear_logs(self):
        try:
            self._log_buffer.clear()
            self.log_area.set_value("")
        except Exception:
            pass

    def _scan_and_seed_samples(self, preexisting: bool = False) -> None:
        try:
            base = Path(self.monitored_directory) if self.monitored_directory else None
            if not base or not base.exists():
                return
            rows: List[Dict[str, Any]] = []
            for sample_dir in base.iterdir():
                if not sample_dir.is_dir():
                    continue
                master = sample_dir / "master.csv"
                if master.exists():
                    sid = sample_dir.name
                    if preexisting:
                        self._preexisting_sample_ids.add(sid)
                    # Determine last_seen from file mtime
                    last_seen = master.stat().st_mtime
                    # Try to load persisted overview and run info
                    run_start = ""
                    device = ""
                    flowcell = ""
                    ov_active = 0
                    ov_total = 0
                    ov_completed = 0
                    ov_failed = 0
                    ov_job_types = ""
                    try:
                        with master.open("r", newline="") as fh:
                            reader = csv.DictReader(fh)
                            first_row = next(reader, None)
                        if first_row:
                            run_start = first_row.get("run_info_run_time", "")
                            device = first_row.get("run_info_device", "")
                            flowcell = first_row.get("run_info_flow_cell", "")
                            # Use saved last_seen if present
                            try:
                                saved_last = float(first_row.get("samples_overview_last_seen", 0.0))
                                if saved_last:
                                    last_seen = saved_last
                            except Exception:
                                pass
                            ov_active = int(first_row.get("samples_overview_active_jobs", 0) or 0)
                            ov_total = int(first_row.get("samples_overview_total_jobs", 0) or 0)
                            ov_completed = int(first_row.get("samples_overview_completed_jobs", 0) or 0)
                            ov_failed = int(first_row.get("samples_overview_failed_jobs", 0) or 0)
                            ov_job_types = str(first_row.get("samples_overview_job_types", "") or "")
                    except Exception:
                        pass
                    origin_value = (
                        "Pre-existing"
                        if sid in self._preexisting_sample_ids
                        else "Live"
                    )
                    try:
                        if origin_value == "Live" and (time.time() - last_seen) >= 3600:
                            origin_value = "Complete"
                    except Exception:
                        pass
                    rows.append(
                        {
                            "sample_id": sid,
                            "origin": origin_value,
                            "run_start": self._format_timestamp_for_display(run_start),
                            "device": device,
                            "flowcell": flowcell,
                            "active_jobs": ov_active,
                            "total_jobs": ov_total,
                            "completed_jobs": ov_completed,
                            "failed_jobs": ov_failed,
                            "job_types": ov_job_types,
                            "last_seen": time.strftime(
                                "%H:%M:%S", time.localtime(last_seen)
                            ),
                            "_last_seen_raw": last_seen,
                        }
                    )
            if rows:
                # Merge with any current rows and update table
                existing = {r["sample_id"]: r for r in (self._last_samples_rows or [])}
                for r in rows:
                    existing[r["sample_id"]] = r
                merged = list(existing.values())
                if hasattr(self, "samples_table"):
                    try:
                        self.samples_table.rows = merged
                        self.samples_table.update()
                    except Exception:
                        pass
                self._last_samples_rows = merged
                self._known_sample_ids = {r["sample_id"] for r in merged}
            if preexisting:
                self._preexisting_scanned = True
        except Exception:
            pass

    def _scan_for_new_samples(self) -> None:
        try:
            base = Path(self.monitored_directory) if self.monitored_directory else None
            if not base or not base.exists():
                return
            # Prepare mappings for efficient lookups and updates
            new_rows: List[Dict[str, Any]] = []
            updated_rows: List[Dict[str, Any]] = []
            existing_by_id: Dict[str, Dict[str, Any]] = {
                r.get("sample_id"): r
                for r in (self._last_samples_rows or [])
                if r.get("sample_id")
            }
            for sample_dir in base.iterdir():
                if not sample_dir.is_dir():
                    continue
                sid = sample_dir.name
                master = sample_dir / "master.csv"
                if master.exists():
                    try:
                        last_seen = master.stat().st_mtime
                        # Load persisted overview and run info
                        run_start = ""
                        device = ""
                        flowcell = ""
                        ov_active = 0
                        ov_total = 0
                        ov_completed = 0
                        ov_failed = 0
                        ov_job_types = ""
                        with master.open("r", newline="") as fh:
                            reader = csv.DictReader(fh)
                            first_row = next(reader, None)
                        if first_row:
                            run_start = first_row.get("run_info_run_time", "")
                            device = first_row.get("run_info_device", "")
                            flowcell = first_row.get("run_info_flow_cell", "")
                            try:
                                saved_last = float(first_row.get("samples_overview_last_seen", 0.0))
                                if saved_last:
                                    last_seen = saved_last
                            except Exception:
                                pass
                            ov_active = int(first_row.get("samples_overview_active_jobs", 0) or 0)
                            ov_total = int(first_row.get("samples_overview_total_jobs", 0) or 0)
                            ov_completed = int(first_row.get("samples_overview_completed_jobs", 0) or 0)
                            ov_failed = int(first_row.get("samples_overview_failed_jobs", 0) or 0)
                            ov_job_types = str(first_row.get("samples_overview_job_types", "") or "")
                    except Exception:
                        last_seen = None
                    if last_seen is None:
                        continue

                    existing_row = existing_by_id.get(sid)
                    if existing_row is None:
                        # New sample discovered → mark as Live
                        new_rows.append(
                            {
                                "sample_id": sid,
                                "origin": "Live",
                                "run_start": run_start,
                                "device": device,
                                "flowcell": flowcell,
                                "active_jobs": ov_active,
                                "total_jobs": ov_total,
                                "completed_jobs": ov_completed,
                                "failed_jobs": ov_failed,
                                "job_types": ov_job_types,
                                "last_seen": time.strftime(
                                    "%H:%M:%S", time.localtime(last_seen)
                                ),
                                "_last_seen_raw": last_seen,
                            }
                        )
                    else:
                        # Existing sample: if master.csv has a newer mtime, update and flip to Live
                        prev_seen = existing_row.get("_last_seen_raw") or 0
                        if last_seen > prev_seen:
                            updated = dict(existing_row)
                            # Preserve or refresh run info/overview from persisted values
                            formatted_run_start = (
                                self._format_timestamp_for_display(run_start) if run_start else ""
                            )
                            updated["run_start"] = formatted_run_start or existing_row.get("run_start", "")
                            updated["device"] = device or existing_row.get("device", "")
                            updated["flowcell"] = flowcell or existing_row.get("flowcell", "")
                            if ov_job_types:
                                updated["job_types"] = ov_job_types
                            updated["active_jobs"] = ov_active
                            updated["total_jobs"] = ov_total
                            updated["completed_jobs"] = ov_completed
                            updated["failed_jobs"] = ov_failed
                            updated["last_seen"] = time.strftime(
                                "%H:%M:%S", time.localtime(last_seen)
                            )
                            updated["_last_seen_raw"] = last_seen
                            # Determine origin based on inactivity threshold
                            try:
                                if (time.time() - last_seen) >= 3600:
                                    updated["origin"] = "Complete"
                                else:
                                    updated["origin"] = "Live"
                            except Exception:
                                updated["origin"] = "Live"
                            updated_rows.append(updated)

            if new_rows or updated_rows:
                merged_map: Dict[str, Dict[str, Any]] = {
                    r["sample_id"]: r for r in (self._last_samples_rows or [])
                }
                for r in new_rows:
                    merged_map[r["sample_id"]] = r
                for r in updated_rows:
                    merged_map[r["sample_id"]] = r
                merged = list(merged_map.values())
                if hasattr(self, "samples_table"):
                    try:
                        self.samples_table.rows = merged
                        self.samples_table.update()
                    except Exception:
                        pass
                self._last_samples_rows = merged
                self._known_sample_ids = {r["sample_id"] for r in merged}
        except Exception:
            pass

    def _launch_workflow_button_clicked(self):
        """Handle launch workflow button click."""
        ui.notify("Launch workflow functionality not implemented yet", type="info")

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

    def submit_sample_job(
        self, sample_dir: str, job_type: str, sample_id: str = None
    ) -> bool:
        """
        Submit a job for an existing sample directory.
        
        This allows users to manually trigger specific job types for samples
        that have already been processed or need reprocessing.
        
        Args:
            sample_dir: Path to the sample directory
            job_type: Type of job to run (e.g., 'igv_bam')
            sample_id: Optional sample ID (defaults to directory name)
            
        Returns:
            True if job was successfully submitted, False otherwise
        """
        try:
            if not hasattr(self, 'workflow_runner') or self.workflow_runner is None:
                print(f"[GUI] No workflow runner available for {job_type} job")
                return False
            
            if sample_id is None:
                sample_id = Path(sample_dir).name
            
            # Check if it's a Simple workflow or Ray workflow
            if hasattr(self.workflow_runner, 'submit_sample_job'):
                # Simple workflow
                return self.workflow_runner.submit_sample_job(sample_dir, job_type, sample_id)
            elif hasattr(self.workflow_runner, 'manager') and hasattr(self.workflow_runner.manager, 'submit_sample_job'):
                # Ray workflow - this is more complex, so we'll provide a fallback
                print(f"[GUI] Ray workflow detected for {job_type} job - manual submission not yet supported")
                return False
            else:
                print(f"[GUI] Unknown workflow type for {job_type} job")
                return False
                
        except Exception as e:
            print(f"[GUI] Failed to submit {job_type} job for sample {sample_id}: {e}")
            return False


def launch_gui(
    host="0.0.0.0",
    port: int = 8081,
    show: bool = False,
    workflow_runner: Any = None,
    workflow_steps: list = None,
    monitored_directory: str = "",
) -> GUILauncher:
    """Legacy launch function (kept for backward compatibility).

    Internally delegates to the refactored launcher in `littlejohn.gui.app`.
    """
    from .gui.app import launch_gui as _launch  # type: ignore

    return _launch(
        host=host,
        port=port,
        show=show,
        workflow_runner=workflow_runner,
        workflow_steps=workflow_steps,
        monitored_directory=monitored_directory,
    )


def get_gui_launcher() -> Optional[GUILauncher]:
    """Compatibility shim that delegates to littlejohn.gui.app.get_gui_launcher."""
    try:
        from .gui.app import get_gui_launcher as _get  # type: ignore

        return _get()
    except Exception:
        return None


def send_gui_update(update_type: UpdateType, data: Dict[str, Any], priority: int = 0):
    """Compatibility shim that delegates to littlejohn.gui.app.send_gui_update."""
    try:
        from .gui.app import send_gui_update as _send  # type: ignore

        _send(update_type, data, priority)
    except Exception:
        logging.info("[GUI] No GUI launcher available for update (update dropped)")


# Global reference retained for backward compatibility (now managed in gui.app)
_current_gui_launcher = None
