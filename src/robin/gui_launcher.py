"""
GUI launcher for robin workflow monitoring.

This module provides a clean interface for launching the workflow GUI
that runs in a separate thread but is completely isolated from
workflow execution to avoid any blocking.
"""

# Suppress pkg_resources deprecation warnings from sorted_nearest
import warnings
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
# Suppress matplotlib tight_layout warnings
warnings.filterwarnings(
    "ignore", message="The figure layout has changed to tight", category=UserWarning
)

import asyncio
import threading
import time
import logging
import queue
from collections import deque
import csv
from datetime import datetime
from robin.analysis.master_csv_manager import MasterCSVManager

from typing import Optional, Dict, Any, List, Set
from pathlib import Path
from dataclasses import dataclass, field, asdict
from enum import Enum
import os
import json
import pickle

from robin.gui import theme, images

from robin.gui.components.news_feed import NewsFeed

from robin.reporting.report import create_pdf
from robin.reporting.sections.disclaimer_text import EXTENDED_DISCLAIMER_TEXT


# Files that indicate an analysis step is complete.
COMPLETION_JOB_PATTERNS: Dict[str, List[str]] = {
    "fusion": [
        "fusion_candidates_master_processed.pkl",
        "fusion_candidates_all_processed.pkl",
        "fusion_summary.csv",
        "fusion_results.csv",
        "sv_count.txt",
    ],
    "cnv": [
        "CNV_dict.npy",
        "XYestimate.pkl",
        "cnv_results.csv",
        "cnv_summary.txt",
        "cnv_analysis_counter.txt",
        "cnv_analysis_results.pkl",
    ],
    "mgmt": ["final_mgmt.csv", "*_mgmt.csv"],
    "target": ["coverage_main.csv", "bed_coverage_main.csv"],
    "sturgeon": ["sturgeon_scores.csv", "sturgeon_results.csv", "sturgeon_summary.csv"],
    "nanodx": ["NanoDX_scores.csv", "nanodx_results.csv", "nanodx_summary.csv"],
    "pannanodx": ["PanNanoDX_scores.csv", "pannanodx_results.csv", "pannanodx_summary.csv"],
    "random_forest": [
        "random_forest_scores.csv",
        "random_forest_results.csv",
        "random_forest_summary.csv",
    ],
}

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
    WARNING_NOTIFICATION = "warning_notification"


@dataclass
class GUIUpdate:
    """A single update message for the GUI."""

    update_type: UpdateType
    timestamp: float
    data: Dict[str, Any]
    priority: int = 0  # Higher priority updates are processed first


@dataclass
class SampleRecord:
    """Master record for a single sample with all tracked information."""
    
    sample_id: str
    origin: str = "Live"  # Live, Pre-existing, or Complete
    run_start: str = ""
    device: str = ""
    flowcell: str = ""
    active_jobs: int = 0
    pending_jobs: int = 0
    total_jobs: int = 0
    completed_jobs: int = 0
    failed_jobs: int = 0
    job_types: str = ""
    last_seen: str = ""  # Formatted timestamp string
    _last_seen_raw: float = 0.0  # Raw timestamp for sorting
    _dirty: bool = False  # Flag to indicate record needs UI update
    _file_mtime: float = 0.0  # master.csv modification time
    files_seen: int = 0  # Total files seen/processed
    files_processed: int = 0  # Files completely processed through all analysis steps
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for table rows."""
        return {
            "sample_id": self.sample_id,
            "origin": self.origin,
            "run_start": self.run_start,
            "device": self.device,
            "flowcell": self.flowcell,
            "active_jobs": self.active_jobs,
            "pending_jobs": self.pending_jobs,
            "total_jobs": self.total_jobs,
            "completed_jobs": self.completed_jobs,
            "failed_jobs": self.failed_jobs,
            "job_types": self.job_types,
            "last_seen": self.last_seen,
            "files_seen": self.files_seen,
            "files_processed": self.files_processed,
            "file_progress": self.files_processed / self.files_seen if self.files_seen > 0 else 0.0,
            "actions": "View",
            "_last_seen_raw": self._last_seen_raw,
        }


class GUILauncher:
    """Launcher for the robin workflow GUI using isolated threading with message queue."""

    def __init__(self, host: str = "0.0.0.0", port: int = 8081, reload: bool = False):
        self.host = host
        self.port = port
        self.gui_thread = None
        self.is_running = False
        self.workflow_runner = None
        self.workflow_steps = []
        self.monitored_directory = ""
        self.reload = reload
        self.center = None  # Center ID for the analysis
        
        # Sample status transition timeout (in seconds)
        # Change to 60 for testing (1 minute), 3600 for production (60 minutes)
        self.completion_timeout_seconds = 60  # Set to 3600 for production
        
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
        
        # Adaptive throttling and update coalescing
        self._last_update_times: Dict[UpdateType, float] = {}
        self._last_queue_size_logged: int = 0
        self._last_log_time: float = 0.0

        # Runtime state
        self._start_time: Optional[float] = None
        self._is_running: bool = False
        
        # Track which samples have had target.bam finalization triggered
        self._finalized_samples: set = set()

        # SNP analysis concurrency control
        self._snp_analysis_lock = threading.Lock()
        self._snp_analysis_running_sample: Optional[str] = None

        # Log ring buffer (last 1000 entries)
        self._log_buffer = deque(maxlen=1000)

        # Cached samples data for persistence across navigation
        self._last_samples_rows: List[Dict[str, Any]] = []
        self._current_sample_id: Optional[str] = None
        self._selected_sample_id: Optional[str] = None
        self._known_sample_ids: set[str] = set()
        self._preexisting_sample_ids: set[str] = set()
        self._preexisting_scanned: bool = False
        self._pending_samples_data: Optional[Dict[str, Any]] = None
        
        # Caching for samples data
        self._last_cache_time: float = 0.0
        self._cache_duration: float = 30.0  # Cache for 30 seconds
        
        # Master record system - centralized source of truth for samples
        self._samples_master_record: Dict[str, SampleRecord] = {}
        self._samples_record_lock = threading.Lock()
        self._master_record_cache_file: Optional[Path] = None  # Set when monitored_directory is known
        self._background_scan_interval: float = 10.0  # Scan every 10 seconds
        self._background_scan_in_progress: bool = False
        # MGMT per-sample cache (latest count seen)
        self._mgmt_state: Dict[str, Dict[str, Any]] = {}
        # Coverage per-sample cache (file mtimes, computed metrics)
        self._coverage_state: Dict[str, Dict[str, Any]] = {}
        
        # Progress notification event
        from nicegui import Event
        self.progress_notification_event = Event[Dict[str, Any]]()
        # CNV per-sample cache
        self._cnv_state: Dict[str, Dict[str, Any]] = {}
        # Cache last seen queue status so we can populate immediately on page creation
        self._last_queue_status: Dict[str, Any] = {}

    # ========== Master Record System Methods ==========
    
    def _setup_master_record_cache(self) -> None:
        """Initialize cache file path when monitored_directory is known."""
        try:
            if self.monitored_directory:
                base = Path(self.monitored_directory)
                if base.exists():
                    # Store cache in the monitored directory
                    self._master_record_cache_file = base / ".samples_master_record.pkl"
                    logging.info(f"Master record cache file: {self._master_record_cache_file}")
        except Exception as e:
            logging.debug(f"Error setting up cache file: {e}")

    def _load_master_record_cache(self) -> bool:
        """Load master record from cache file. Returns True if successful."""
        try:
            if not self._master_record_cache_file or not self._master_record_cache_file.exists():
                return False
            
            with open(self._master_record_cache_file, 'rb') as f:
                cached_data = pickle.load(f)
            
            if isinstance(cached_data, dict):
                # Convert dicts back to SampleRecord objects
                with self._samples_record_lock:
                    self._samples_master_record = {
                        sid: SampleRecord(**data) if isinstance(data, dict) else data
                        for sid, data in cached_data.items()
                    }
                logging.info(f"Loaded {len(self._samples_master_record)} samples from cache")
                # Mark all as dirty to trigger UI refresh
                for record in self._samples_master_record.values():
                    record._dirty = True
                return True
        except Exception as e:
            logging.debug(f"Error loading master record cache: {e}")
        return False

    def _save_master_record_cache(self) -> None:
        """Save master record to cache file."""
        try:
            if not self._master_record_cache_file:
                return
            
            with self._samples_record_lock:
                # Convert SampleRecord objects to dicts for pickle
                cache_data = {
                    sid: asdict(record)
                    for sid, record in self._samples_master_record.items()
                }
            
            # Save to temporary file first, then rename (atomic operation)
            temp_file = self._master_record_cache_file.with_suffix('.pkl.tmp')
            with open(temp_file, 'wb') as f:
                pickle.dump(cache_data, f)
            temp_file.replace(self._master_record_cache_file)
            logging.debug(f"Saved master record cache ({len(cache_data)} samples)")
        except Exception as e:
            logging.debug(f"Error saving master record cache: {e}")

    def _background_scan_samples(self) -> None:
        """Schedule a background scan without blocking the UI thread."""
        try:
            import asyncio
            asyncio.create_task(self._background_scan_samples_async())
        except RuntimeError:
            # Fallback if no event loop is available
            self._background_scan_samples_sync()
        except Exception as e:
            logging.debug(f"Failed to schedule background scan: {e}")

    async def _background_scan_samples_async(self) -> None:
        """Run the background scan off the UI thread and refresh UI afterwards."""
        if getattr(self, "_background_scan_in_progress", False):
            return
        self._background_scan_in_progress = True
        try:
            import asyncio
            updated_any = await asyncio.to_thread(self._background_scan_samples_sync)
            if updated_any and hasattr(self, "samples_table"):
                self._refresh_table_from_master()
        except Exception as e:
            logging.debug(f"Background scan task failed: {e}")
        finally:
            self._background_scan_in_progress = False

    def _background_scan_samples_sync(self) -> bool:
        """Scan folder and update master record (file I/O). Returns True if updated."""
        try:
            if not self.monitored_directory:
                return False
            
            base = Path(self.monitored_directory)
            if not base.exists():
                return False
            
            # Scan all sample directories
            now_ts = time.time()
            updated_any = False
            
            with self._samples_record_lock:
                # Track which samples we've seen in this scan
                seen_sample_ids = set()
                
                for sample_dir in base.iterdir():
                    if not sample_dir.is_dir():
                        continue
                    
                    sid = sample_dir.name
                    seen_sample_ids.add(sid)
                    master_csv = sample_dir / "master.csv"
                    
                    if not master_csv.exists():
                        continue
                    
                    try:
                        # Get file modification time
                        file_mtime = master_csv.stat().st_mtime
                        
                        # Get or create record
                        record = self._samples_master_record.get(sid)
                        if record is None:
                            # New sample - create record
                            record = SampleRecord(sample_id=sid)
                            record._dirty = True
                            updated_any = True
                            # Check if it's pre-existing
                            if sid in self._preexisting_sample_ids:
                                record.origin = "Pre-existing"
                            else:
                                record.origin = "Live"
                        elif file_mtime > record._file_mtime:
                            # File has changed - update record
                            record._dirty = True
                            updated_any = True
                        
                        # Update file mtime
                        record._file_mtime = file_mtime
                        
                        # Read master.csv for metadata
                        try:
                            with master_csv.open("r", newline="") as fh:
                                reader = csv.DictReader(fh)
                                first_row = next(reader, None)
                            if first_row:
                                # Update run info
                                record.run_start = self._format_timestamp_for_display(
                                    first_row.get("run_info_run_time", "")
                                )
                                record.device = first_row.get("run_info_device", "") or ""
                                record.flowcell = first_row.get("run_info_flow_cell", "") or ""
                                
                                # Update last_seen from saved value or use file mtime
                                try:
                                    saved_last = float(
                                        first_row.get("samples_overview_last_seen", 0.0) or 0.0
                                    )
                                    if saved_last > 0:
                                        record._last_seen_raw = saved_last
                                    else:
                                        record._last_seen_raw = file_mtime
                                except Exception:
                                    record._last_seen_raw = file_mtime
                                
                                # Update job counts from persisted overview
                                record.active_jobs = int(
                                    first_row.get("samples_overview_active_jobs", 0) or 0
                                )
                                record.pending_jobs = int(
                                    first_row.get("samples_overview_pending_jobs", 0) or 0
                                )
                                record.total_jobs = int(
                                    first_row.get("samples_overview_total_jobs", 0) or 0
                                )
                                record.completed_jobs = int(
                                    first_row.get("samples_overview_completed_jobs", 0) or 0
                                )
                                record.failed_jobs = int(
                                    first_row.get("samples_overview_failed_jobs", 0) or 0
                                )
                                record.job_types = str(
                                    first_row.get("samples_overview_job_types", "") or ""
                                )
                        except Exception as e:
                            logging.debug(f"Error reading master.csv for {sid}: {e}")
                        
                        # Update last_seen formatted string
                        record.last_seen = time.strftime(
                            "%Y-%m-%d %H:%M:%S", time.localtime(record._last_seen_raw)
                        )
                        
                        # Update origin based on inactivity timeout AND active jobs status
                        prev_origin = record.origin
                        # Only mark as Complete if timeout passed AND no active jobs
                        if record.origin == "Live" and (now_ts - record._last_seen_raw) >= self.completion_timeout_seconds:
                            if record.active_jobs == 0 and record.pending_jobs == 0:
                                record.origin = "Complete"
                                record._dirty = True
                                # Trigger finalization if transitioning from Live to Complete
                                if prev_origin == "Live" and sid not in self._finalized_samples:
                                    self._trigger_target_bam_finalization(sid)
                                    self._finalized_samples.add(sid)
                            # If there are active jobs, keep as Live even if timeout passed
                        elif record.origin == "Pre-existing" and (now_ts - record._last_seen_raw) >= self.completion_timeout_seconds:
                            # Keep as Pre-existing if it was pre-existing and still inactive
                            pass
                        elif record.origin == "Complete":
                            # Reactivate if file was modified recently OR if there are active jobs
                            if (now_ts - record._last_seen_raw) < self.completion_timeout_seconds or record.active_jobs > 0 or record.pending_jobs > 0:
                                record.origin = "Live"
                                record._dirty = True
                        
                        # Save to master record
                        self._samples_master_record[sid] = record
                        
                        # Persist overview to master.csv via MasterCSVManager
                        try:
                            manager = MasterCSVManager(str(base))
                            persist_payload = {
                                "active_jobs": int(record.active_jobs),
                                "pending_jobs": int(record.pending_jobs),
                                "total_jobs": int(record.total_jobs),
                                "completed_jobs": int(record.completed_jobs),
                                "failed_jobs": int(record.failed_jobs),
                                "job_types": record.job_types,
                                "last_seen": float(record._last_seen_raw),
                            }
                            manager.update_sample_overview(sid, persist_payload)
                        except Exception as e:
                            logging.debug(f"Error persisting overview for {sid}: {e}")
                    
                    except Exception as e:
                        logging.debug(f"Error processing sample {sid}: {e}")
                
                # Remove samples that no longer exist
                to_remove = set(self._samples_master_record.keys()) - seen_sample_ids
                if to_remove:
                    for sid in to_remove:
                        del self._samples_master_record[sid]
                    updated_any = True
                    logging.debug(f"Removed {len(to_remove)} deleted samples from master record")
            
            # Save cache if anything changed
            if updated_any:
                self._save_master_record_cache()
            
            return updated_any
        
        except Exception as e:
            logging.error(f"Error in background sample scan: {e}")
            return False

    def _refresh_table_from_master(self) -> None:
        """Fast UI refresh - reads from master record and applies filters.
        This is called from the background scanner and does NOT do file I/O."""
        try:
            if not hasattr(self, "samples_table"):
                return
            
            # Get all records from master
            with self._samples_record_lock:
                records = list(self._samples_master_record.values())
            
            # Convert to row dicts
            rows = [record.to_dict() for record in records]
            
            # Update cached rows
            self._last_samples_rows = rows
            self._last_cache_time = time.time()
            
            # Apply filters and update table
            self._apply_samples_table_filters()
            
        except Exception as e:
            logging.error(f"Error refreshing table from master: {e}")

    def _update_master_record_from_workflow(self, samples_data: List[Dict[str, Any]]) -> None:
        """Update master record from workflow polling data.
        This merges workflow stats into the master record without doing file I/O."""
        try:
            with self._samples_record_lock:
                for s in samples_data:
                    sid = s.get("sample_id", "") or "unknown"
                    if sid == "unknown":
                        continue
                    
                    last_seen = float(s.get("last_seen", time.time()))
                    
                    # Get or create record
                    record = self._samples_master_record.get(sid)
                    if record is None:
                        # Create new record (will be populated by background scan)
                        record = SampleRecord(
                            sample_id=sid,
                            _last_seen_raw=last_seen,
                            last_seen=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(last_seen)),
                        )
                        if sid in self._preexisting_sample_ids:
                            record.origin = "Pre-existing"
                        else:
                            record.origin = "Live"
                    
                    # Update job counts from workflow data (these are more current than file scans)
                    record.active_jobs = s.get("active_jobs", 0)
                    record.pending_jobs = s.get("pending_jobs", 0)
                    record.total_jobs = s.get("total_jobs", 0)
                    record.completed_jobs = s.get("completed_jobs", 0)
                    record.failed_jobs = s.get("failed_jobs", 0)
                    
                    # Update file progress from job counts (same data source as other columns)
                    record.files_seen = record.total_jobs
                    record.files_processed = record.completed_jobs
                    
                    # Update job types
                    job_types = s.get("job_types", [])
                    if isinstance(job_types, list):
                        record.job_types = ",".join(sorted(set(job_types)))
                    else:
                        record.job_types = str(job_types or "")
                    
                    # Update last_seen if this is newer
                    if last_seen > record._last_seen_raw:
                        record._last_seen_raw = last_seen
                        record.last_seen = time.strftime(
                            "%Y-%m-%d %H:%M:%S", time.localtime(last_seen)
                        )
                        record._dirty = True
                    
                    # Mark as dirty to trigger UI update
                    record._dirty = True
                    self._samples_master_record[sid] = record
                
                # Persist updates to master.csv
                try:
                    base = Path(self.monitored_directory) if self.monitored_directory else None
                    if base and base.exists():
                        manager = MasterCSVManager(str(base))
                        for record in self._samples_master_record.values():
                            if record._dirty:
                                persist_payload = {
                                    "active_jobs": int(record.active_jobs),
                                    "pending_jobs": int(record.pending_jobs),
                                    "total_jobs": int(record.total_jobs),
                                    "completed_jobs": int(record.completed_jobs),
                                    "failed_jobs": int(record.failed_jobs),
                                    "job_types": record.job_types,
                                    "last_seen": float(record._last_seen_raw),
                                }
                                manager.update_sample_overview(record.sample_id, persist_payload)
                                record._dirty = False
                except Exception as e:
                    logging.debug(f"Error persisting workflow updates: {e}")
            
            # Trigger UI refresh
            if hasattr(self, "samples_table"):
                self._refresh_table_from_master()
        
        except Exception as e:
            logging.error(f"Error updating master record from workflow: {e}")

    def launch_gui(
        self,
        workflow_runner: Any = None,
        workflow_steps: list = None,
        monitored_directory: str = "",
        reload: bool = False,
        center: str = None,
    ) -> bool:
        """Launch the GUI in a completely isolated background thread.
        
        Note: When reload=True, the GUI must run in the main thread because
        signal handlers can only be set in the main thread. In this case,
        this method will block until the GUI is stopped.
        """
        if ui is None:
            logging.error("NiceGUI is not available")
            return False

        self.workflow_runner = workflow_runner
        self.workflow_steps = workflow_steps or []
        self.center = center

        # Store absolute monitored directory to avoid relative path issues
        try:
            self.monitored_directory = (
                str(Path(monitored_directory).resolve()) if monitored_directory else ""
            )
        except Exception:
            self.monitored_directory = monitored_directory
        
        # Setup master record cache when directory is known
        if self.monitored_directory:
            self._setup_master_record_cache()
            # Try to load cache immediately for fast startup
            if self._load_master_record_cache():
                logging.info("Loaded samples from cache - table will populate immediately")

        try:
            # When reload is requested, we must run in the main thread
            # because signal handlers can only be set in the main thread
            if self.reload:
                # Check if we're in the main thread
                if threading.current_thread() is threading.main_thread():
                    # Run directly in main thread (blocking)
                    logging.info("Running GUI with reload in main thread")
                    self.is_running = True
                    try:
                        self._run_gui_worker()
                    except KeyboardInterrupt:
                        logging.info("GUI stopped by user")
                    finally:
                        self.is_running = False
                    return True
                else:
                    # Not in main thread - can't use reload
                    logging.warning(
                        "Reload requested but not in main thread. "
                        "Disabling reload functionality. "
                        "Run from main thread to enable reload."
                    )
                    self.reload = False
            
            # Start GUI in completely isolated background thread (when reload is False)
            self.gui_thread = threading.Thread(
                target=self._run_gui_worker, daemon=True, name="robin-GUI-Thread"
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
            current_time = time.time()
            queue_size = self.update_queue.qsize()
            # Adaptive threshold: higher threshold for large workloads
            # Threshold increases based on queue size to prevent dropping updates during bursts
            threshold = 500 if queue_size < 200 else 1000
            
            # Rate limiting: Skip low-priority updates if queue is getting too large
            if queue_size > threshold and priority < 5:
                # Only log when queue size changes significantly to reduce log spam
                if abs(queue_size - self._last_queue_size_logged) > 50 or (current_time - self._last_log_time) > 10.0:
                    logging.warning(
                        f"[GUI] Skipping low-priority update due to queue size: "
                        f"{queue_size} (threshold: {threshold})"
                    )
                    self._last_queue_size_logged = queue_size
                    self._last_log_time = current_time
                return
            
            # Update coalescing: Skip duplicate low-priority updates of the same type within 0.5 seconds
            # This prevents queue buildup from rapid duplicate updates
            last_time = self._last_update_times.get(update_type, 0)
            time_since_last = current_time - last_time
            
            if priority < 5 and time_since_last < 0.5:
                # Skip duplicate low-priority updates within coalescing window
                return
            
            # Create and enqueue the current update
            update = GUIUpdate(
                update_type=update_type,
                timestamp=current_time,
                data=data,
                priority=priority,
            )
            
            self._last_update_times[update_type] = current_time

            # Use negative priority so higher priority updates come first.
            # Include a monotonically increasing sequence as a tiebreaker so
            # heap comparisons never fallback to comparing GUIUpdate objects.
            self._update_seq += 1
            self.update_queue.put((-priority, self._update_seq, update))
            self.total_updates_enqueued += 1
            
            # Reduced logging: only log every 100th update or important updates
            if self.total_updates_enqueued % 100 == 0 or priority >= 8:
                logging.debug(
                    f"[GUI] Enqueued update #{self.total_updates_enqueued}: {update_type.value} (queue: {queue_size})"
                )

        except Exception as e:
            # Don't let update failures affect the workflow
            logging.debug(f"Failed to send GUI update: {e}")

    def _process_progress_queue(self):
        """Process queued progress updates for report generation."""
        try:
            from robin.gui.report_progress import progress_manager
            
            # Process updates in the UI context
            while not progress_manager.progress_queue.empty():
                update = progress_manager.progress_queue.get_nowait()
                self._handle_progress_update(update)
                
        except queue.Empty:
            pass
        except Exception as e:
            logging.debug(f"Error processing progress queue: {e}")
    
    def _setup_notification_system(self, container):
        """Set up the notification system with the provided container."""
        try:
            # Set up event subscription in the proper UI context
            self.progress_notification_event.subscribe(
                lambda data: self._show_notification_in_container(container, data)
            )
            logging.info("Notification system set up successfully")
        except Exception as e:
            logging.error(f"Error setting up notification system: {e}")

    def _show_notification_in_container(self, container, data: Dict[str, Any]):
        """Show notification in the dedicated container."""
        try:
            message = data.get('message', '')
            notification_type = data.get('type', 'info')
            timeout = data.get('timeout', 5000)
            
            # Create notification in the container
            with container:
                ui.notify(
                    message,
                    type=notification_type,
                    timeout=timeout,
                    position="top-right"
                )
                
        except Exception as e:
            logging.error(f"Error showing notification in container: {e}")

    def _handle_progress_notification_event(self, event_data: Dict[str, Any]):
        """Handle progress notification events."""
        try:
            from nicegui import ui

            message = event_data.get('message', '')
            notification_type = event_data.get('type', 'info')
            timeout = event_data.get('timeout', 5000)

            # Show notification in the proper UI context
            ui.notify(
                message,
                type=notification_type,
                timeout=timeout,
                position="top-right"
            )

        except Exception as e:
            logging.error(f"Error handling progress notification event: {e}")

    def _handle_progress_update(self, update: Dict[str, Any]):
        """Handle a single progress update in the UI context."""
        try:
            from nicegui import ui
            
            sample_id = update['sample_id']
            update_type = update['type']
            
            if update_type == 'start':
                # Emit event for initial notification
                self.progress_notification_event.emit({
                    'sample_id': sample_id,
                    'message': f"[{sample_id}] Starting report generation...",
                    'type': 'info',
                    'timeout': 3000  # 3 seconds
                })
                
                logging.debug(f"Started report generation for {sample_id}")
                
            elif update_type == 'update':
                stage = update['stage']
                message = update['message']
                progress = update.get('progress')
                
                # Calculate progress percentage
                progress_percent = int((progress or 0.0) * 100) if progress is not None else ""
                progress_text = f" ({progress_percent}%)" if progress_percent else ""
                
                # Emit event for progress notification
                self.progress_notification_event.emit({
                    'sample_id': sample_id,
                    'message': f"[{sample_id}] {message}{progress_text}",
                    'type': 'info',
                    'timeout': 2000  # 2 seconds for progress updates
                })
                
                logging.debug(f"Updated progress for {sample_id}: {stage} - {message}")
                
            elif update_type == 'complete':
                filename = update.get('filename')
                
                # Show completion notification
                completion_message = f"[{sample_id}] Report generation completed"
                if filename:
                    completion_message += f": {filename}"
                
                # Emit event for completion notification
                self.progress_notification_event.emit({
                    'sample_id': sample_id,
                    'message': completion_message,
                    'type': 'positive',
                    'timeout': 5000
                })
                
                logging.info(f"Completed report generation for {sample_id}")
                
            elif update_type == 'error':
                error_message = update['error_message']
                
                # Emit event for error notification
                self.progress_notification_event.emit({
                    'sample_id': sample_id,
                    'message': f"[{sample_id}] Report generation failed: {error_message}",
                    'type': 'negative',
                    'timeout': 10000
                })
                
                logging.error(f"Report generation failed for {sample_id}: {error_message}")
                
        except Exception as e:
            logging.error(f"Error handling progress update: {e}")

    def _drain_updates_on_ui(self):
        """Drain queued updates and apply them on the UI thread (called by ui.timer).
        
        Uses adaptive throttling: processes more updates when queue is large,
        fewer when queue is small. Batches UI updates for better performance.
        """
        if not self.gui_ready.is_set():
            return
        
        processed_count = 0
        max_updates_per_cycle = 50  # Limit updates per cycle to prevent UI blocking
        
        try:
            queue_size = self.update_queue.qsize()
            
            # Adaptive processing: process more updates when queue is large
            if queue_size > 200:
                max_updates_per_cycle = 100  # Process more aggressively when backlogged
            elif queue_size > 100:
                max_updates_per_cycle = 75
            elif queue_size < 20:
                max_updates_per_cycle = 25  # Process fewer when queue is small
            
            # Process updates in batch
            while processed_count < max_updates_per_cycle:
                try:
                    # Entries are (-priority, seq, GUIUpdate)
                    _, _, update = self.update_queue.get_nowait()
                except queue.Empty:
                    break
                
                self._handle_update(update)
                self.total_updates_processed += 1
                processed_count += 1
                
                # Reduced logging: only log every 50th update or when queue size changes significantly
                if self.total_updates_processed % 50 == 0:
                    logging.debug(
                        f"[GUI] Processed update #{self.total_updates_processed}: {update.update_type.value} (queue: {queue_size})"
                    )
            
            # Note: Timer interval is fixed at 0.5s for consistent performance
            # Adaptive processing (max_updates_per_cycle) handles queue size variations
                
        except Exception as e:
            logging.debug(f"[GUI] Error draining updates on UI: {e}")
        finally:
            # Batch UI update: only call ui.update() once per drain cycle
            if processed_count > 0:
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
                # Update master record from workflow polling data
                # This merges workflow stats into master record without file I/O
                self._pending_samples_data = update.data
                try:
                    samples_list = update.data.get("samples", [])
                    if samples_list:
                        self._update_master_record_from_workflow(samples_list)
                except Exception as e:
                    logging.error(f"[GUI] Samples master record update failed: {e}")
            elif update.update_type == UpdateType.WARNING_NOTIFICATION:
                self._show_warning_notification(update.data)

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

                # Show/hide the "no active jobs" placeholder
                if hasattr(self, "no_active_jobs_label"):
                    if job_rows:
                        self.no_active_jobs_label.set_visibility(False)
                    else:
                        self.no_active_jobs_label.set_visibility(True)

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

    def _update_file_progress(self):
        """Update file progress display in workflow monitor from current samples data."""
        try:
            if not hasattr(self, "sample_files_progress_container") or not hasattr(self, "_last_samples_rows"):
                return
            
            rows = self._last_samples_rows or []
            total_files_seen = 0
            total_files_processed = 0
            sample_progress_data = []
            
            # Calculate totals and collect per-sample data
            for row in rows:
                files_seen = row.get("files_seen", 0) or 0
                files_processed = row.get("files_processed", 0) or 0
                sample_id = row.get("sample_id", "")
                
                if files_seen > 0:
                    total_files_seen += files_seen
                    total_files_processed += files_processed
                    sample_progress_data.append({
                        "sample_id": sample_id,
                        "files_seen": files_seen,
                        "files_processed": files_processed,
                        "progress": files_processed / files_seen if files_seen > 0 else 0.0
                    })
            
            # Update overall progress
            if hasattr(self, "overall_files_progress") and hasattr(self, "overall_files_label"):
                overall_progress = total_files_processed / total_files_seen if total_files_seen > 0 else 0.0
                self.overall_files_progress.set_value(round(overall_progress, 2))
                self.overall_files_label.set_text(f"{total_files_processed}/{total_files_seen} files processed")
            
            # Update per-sample progress bars (limit to first 20 to avoid UI overload)
            if hasattr(self, "sample_files_progress_container"):
                try:
                    # Clear existing content
                    self.sample_files_progress_container.clear()
                    
                    if sample_progress_data:
                        # Sort by progress (lowest first) to show samples needing attention
                        sample_progress_data.sort(key=lambda x: x["progress"])
                        
                        # Show up to 20 samples with most activity
                        for sample_info in sample_progress_data[:20]:
                            sample_id = sample_info["sample_id"]
                            progress = sample_info["progress"]
                            files_seen = sample_info["files_seen"]
                            files_processed = sample_info["files_processed"]
                            
                            with self.sample_files_progress_container:
                                with ui.row().classes("w-full items-center gap-2 mb-1"):
                                    ui.label(sample_id).classes("text-xs font-medium min-w-[120px]")
                                    progress_bar = ui.linear_progress(progress).classes("flex-1")
                                    if progress >= 1.0:
                                        progress_bar.props("color=positive")
                                    ui.label(f"{files_processed}/{files_seen}").classes("text-xs text-gray-600 min-w-[60px]")
                    else:
                        # Show placeholder when no data
                        with self.sample_files_progress_container:
                            ui.label("No file progress data available yet").classes("text-xs text-gray-500 italic")
                except Exception as e:
                    logging.debug(f"Error updating per-sample file progress: {e}")
                            
        except Exception as e:
            logging.debug(f"Error updating file progress: {e}")

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

    def _show_warning_notification(self, data: Dict[str, Any]):
        """Show a warning notification in the GUI with dismiss button."""
        try:
            if ui is None:
                return
            
            message = data.get("message", "")
            title = data.get("title", "Warning")
            sample_id = data.get("sample_id", "")
            filename = data.get("filename", "")
            level = data.get("level", "warning")
            
            # Build notification message
            if sample_id:
                notification_msg = f"[{sample_id}] {message}"
            elif filename:
                notification_msg = f"[{filename}] {message}"
            else:
                notification_msg = message
            
            # Show dismissible notification
            # NiceGUI notifications are dismissible by default with a close button
            # timeout=0 makes it persistent until manually dismissed
            ui.notify(
                notification_msg,
                type=level,
                timeout=0,  # Persistent until manually dismissed (close button available)
                position="top-right",
            )
            
        except Exception as e:
            logging.debug(f"Error showing warning notification: {e}")

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
            threading.current_thread().name = "robin-GUI-Thread"

            # Create the main workflow monitor page
            @ui.page("/")
            def welcome_page():
                """Welcome page at root route."""
                _setup_global_resources()
                self._create_welcome_page()

            # Create the workflow monitoring page
            @ui.page("/robin")
            def workflow_monitor():
                """Workflow monitoring page under /robin route."""
                _setup_global_resources()
                self._create_workflow_monitor()

            # Create the samples overview page
            @ui.page("/live_data")
            def samples_overview():
                """Samples overview page showing all tracked samples."""
                _setup_global_resources()
                self._create_samples_overview()

            # Create individual sample detail pages
            @ui.page("/live_data/{sample_id}")
            def sample_detail(sample_id: str):
                """Individual sample detail page."""
                _setup_global_resources()

                # Clear cached state BEFORE creating the page so plots refresh on load
                # This ensures all graphs load properly on each page visit
                self._refresh_sample_plots(sample_id)
                
                # Add a page visit handler to refresh plots on reconnect
                def on_page_visit():
                    """Handle page visit - refresh plots if needed."""
                    try:
                        logging.info(
                            f"Page visit detected for sample {sample_id} - triggering plot refresh"
                        )
                        # Small delay to ensure components are loaded
                        ui.timer(
                            1.0,
                            lambda: self._refresh_sample_plots(sample_id),
                            once=True,
                        )
                    except Exception as e:
                        logging.debug(f"Page visit handler failed for {sample_id}: {e}")

                # Register the page visit handler
                ui.context.client.on_connect(on_page_visit)

                self._create_sample_detail_page(sample_id)

            # Create sample details page
            @ui.page("/live_data/{sample_id}/details")
            def sample_details(sample_id: str):
                """Sample details page with comprehensive information."""
                _setup_global_resources()
                self._create_sample_details_page(sample_id)

            # Download API endpoint
            @ui.page("/api/download/{sample_id}/{filename}")
            def download_file(sample_id: str, filename: str):
                """Download a file from a sample directory."""
                try:
                    # Security: Only allow alphanumeric characters and common file extensions
                    import re
                    if not re.match(r'^[a-zA-Z0-9._-]+$', filename):
                        ui.notify("Invalid filename", type="error")
                        return
                    
                    # Find the sample directory
                    base_dir = Path(self.monitored_directory) if self.monitored_directory else None
                    if not base_dir or not base_dir.exists():
                        ui.notify("Sample directory not found", type="error")
                        return
                    
                    sample_dir = base_dir / sample_id
                    if not sample_dir.exists():
                        ui.notify(f"Sample {sample_id} not found", type="error")
                        return
                    
                    file_path = sample_dir / filename
                    if not file_path.exists() or not file_path.is_file():
                        ui.notify(f"File {filename} not found", type="error")
                        return
                    
                    # Read file content
                    with open(file_path, 'rb') as f:
                        content = f.read()
                    
                    # Use NiceGUI's download functionality
                    ui.download(
                        content,
                        filename=filename,
                        media_type='application/octet-stream'
                    )
                    
                except Exception as e:
                    ui.notify(f"Download failed: {e}", type="error")

            # Setup global CSS and static files - moved to a helper function
            def _setup_global_resources():
                """Setup global CSS and static file resources."""
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

            # Setup global timers and processing - moved to a helper function
            def _setup_global_timers():
                """Setup global timers and update processing."""
                try:
                    self.gui_ready.set()
                    # Faster drain rate: 0.5s to keep up with high-volume updates
                    # The _drain_updates_on_ui method adaptively processes more/fewer updates per cycle
                    app.timer(0.5, self._drain_updates_on_ui, active=True)
                    
                    # Background master record scanner - scans folder periodically
                    # Runs every 10 seconds to keep master record up to date
                    # First scan happens immediately if cache was loaded, otherwise after 1 second
                    initial_delay = 1.0 if not self._samples_master_record else 0.1
                    app.timer(
                        initial_delay,
                        lambda: self._background_scan_samples(),
                        once=True,
                    )
                    # Subsequent scans every 10 seconds
                    app.timer(self._background_scan_interval, self._background_scan_samples, active=True)
                    
                    # If cache was not loaded, do initial preexisting scan after GUI is ready
                    if not self._samples_master_record:
                        app.timer(
                            2.0,
                            lambda: self._background_scan_samples(),  # Initial scan
                            once=True,
                        )
                    
                    # Process progress queue for report generation
                    app.timer(0.1, self._process_progress_queue, active=True)
                    
                except Exception:
                    pass

            try:
                iconfile = os.path.join(
                    os.path.dirname(os.path.abspath(images.__file__)), "favicon.ico"
                )
                if not os.path.exists(iconfile):
                    logging.warning(f"Favicon file not found: {iconfile}")
            except Exception as e:
                logging.error(f"Error locating favicon: {str(e)}")
                iconfile = None
            
            # Setup global timers once when the app starts
            _setup_global_timers()
            
            # Start the GUI
            ui.run(
                host=self.host,
                port=self.port,
                show=False,
                reload=self.reload,
                title="ROBIN",
                storage_secret="robin",
                favicon=iconfile,
                reconnect_timeout=60,
            )
        except Exception as e:
            print(f"GUI worker error: {e}")
            import traceback

            traceback.print_exc()

    def _create_welcome_page(self):
        """Create the welcome page."""
        # Background and main container
        with theme.frame(
            "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>",
            smalltitle="<strong>R.O.B.I.N</strong>",
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            # Main content container with full width and proper centering
            with ui.column().classes("w-full min-h-screen p-4 md:p-8"):
                # Header section
                with ui.card().classes(
                    "w-full max-w-6xl mx-auto shadow-lg rounded-xl mb-8"
                ):
                    with ui.column().classes("p-8 text-center"):
                        ui.label("Welcome to R.O.B.I.N").classes(
                            "text-sky-600 dark:text-white text-4xl md:text-5xl font-bold mb-4"
                        ).style("font-weight: 600")

                        ui.label(
                            "This tool enables real time analysis of data from Oxford Nanopore Technologies sequencers."
                        ).classes(
                            "text-gray-700 dark:text-gray-300 text-lg md:text-xl mb-6"
                        )

                # Main content section
                with ui.card().classes(
                    "w-full max-w-6xl mx-auto shadow-lg rounded-xl mb-8"
                ):
                    with ui.column().classes("p-8"):
                        ui.label("What is R.O.B.I.N?").classes(
                            "text-2xl font-semibold mb-4 text-center"
                        )
                        ui.label(
                            "ROBIN (Rapid nanopOre Brain intraoperatIve classificatioN) is a comprehensive bioinformatics workflow system designed for processing and analyzing BAM files. It provides automated preprocessing, multiple analysis pipelines, and real-time monitoring capabilities. It now incorporates LITTLE JOHN (Lightweight Infrastructure for Task Tracking and Logging with Extensible Job Orchestration for High-throughput aNalysis), which handles the heavy lifting behind the scenes."
                        ).classes("text-gray-700 mb-8 text-center text-lg")

                        # Action buttons - responsive grid layout
                        with ui.grid(columns=2).classes("w-full gap-4 mb-8"):
                            ui.link("View All Samples", "/live_data").classes(
                                "bg-blue-600 hover:bg-blue-700 text-white text-lg font-semibold px-6 py-4 rounded-lg shadow-lg transition-colors text-center"
                            )
                            ui.link("Open Workflow Monitor", "/robin").classes(
                                "bg-blue-600 hover:bg-blue-700 text-white text-lg font-semibold px-6 py-4 rounded-lg shadow-lg transition-colors text-center"
                            )
                            ui.button(
                                "Launch New Workflow",
                                on_click=lambda: self._launch_workflow_button_clicked(),
                            ).classes(
                                "bg-green-600 hover:bg-green-700 text-white text-lg font-semibold px-6 py-4 rounded-lg shadow-lg transition-colors"
                            )
                            ui.link(
                                "View Documentation",
                                "https://looselab.github.io/ROBIN/",
                            ).classes(
                                "bg-blue-600 hover:bg-blue-700 text-white text-lg font-semibold px-6 py-4 rounded-lg shadow-lg transition-colors text-center"
                            )

                # News feed section
                with ui.card().classes(
                    "w-full max-w-6xl mx-auto shadow-lg rounded-xl mb-8"
                ):
                    with ui.column().classes("w-full"):
                        # Initialize news feed only if it hasn't been initialized yet
                        if self.news_feed is None:
                            self.news_feed = NewsFeed()
                            self.news_feed.start_update_timer()
                        # Create the news element
                        self.news_feed.create_news_element()

    def _create_samples_overview(self):
        """Create the samples overview page showing all tracked samples."""
        # Page title and navigation
        with theme.frame(
            "R.O.B.I.N - Sample Tracking Overview",
            smalltitle="Samples",
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            # Main content area
            with ui.column().classes("w-full p-4 gap-4"):
                # Sample statistics
                with ui.card().classes(
                    "w-full bg-gradient-to-r from-blue-50 to-indigo-50"
                ):
                    ui.label("Sample Statistics").classes(
                        "text-lg font-semibold mb-4 text-blue-800"
                    )

                    # SIMPLIFIED: Show basic info without workflow_state access
                    ui.label(
                        "Sample statistics will be available once the workflow is running."
                    ).classes("text-sm text-gray-600")

                # Samples table section - removed card wrapper to enable touch scrolling on mobile
                with ui.column().classes("w-full gap-2"):
                    with ui.row().classes("w-full items-center justify-between mb-4"):
                        ui.label("All Tracked Samples").classes(
                            "text-lg font-semibold"
                        )
                        
                        # Loading indicator
                        self.samples_loading_indicator = ui.spinner(size="sm", color="primary")
                        self.samples_loading_indicator.set_visibility(False)
                        
                        #ui.button(
                        #    "Refresh All Plots",
                        #    on_click=lambda: self._refresh_all_sample_plots(),
                        #).classes("q-btn--secondary")

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
                            .on("change", lambda e: self._set_samples_query(self._extract_event_value(e)))
                            .on("update:model-value", lambda e: self._set_samples_query(self._extract_event_value(e)))
                            .classes("ml-auto")
                        )
                        
                        # Initialize search input with current filter value
                        try:
                            current_query = self._samples_filters.get("query", "")
                            if current_query:
                                self.samples_search.value = current_query
                        except Exception:
                            pass

                        # Origin filter
                        self.origin_filter = (
                            ui.select(
                                options=["All", "Live", "Pre-existing", "Complete"],
                                value=self._samples_filters.get("origin", "All"),
                                label="Origin",
                            )
                            .props("dense clearable")
                            .on("change", lambda e: self._set_samples_origin_filter(self._extract_event_value(e)))
                            .on("update:model-value", lambda e: self._set_samples_origin_filter(self._extract_event_value(e)))
                        )

                    # Loading state container
                    self.samples_loading_container = ui.card().classes("w-full p-8 text-center")
                    with self.samples_loading_container:
                        ui.spinner(size="lg", color="primary")
                        ui.label("Loading samples...").classes("ml-2 text-lg")
                        ui.label("This may take a moment for large directories").classes("text-sm text-gray-500 mt-2")
                    
                    # Create a placeholder table that will be updated later
                    from robin.gui.theme import styled_table

                    # Create samples table
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
                                "name": "file_progress",
                                "label": "File Progress",
                                "field": "file_progress",
                                "sortable": True,
                            },
                            {
                                "name": "active_jobs",
                                "label": "Active",
                                "field": "active_jobs",
                                "sortable": True,
                            },
                            {
                                "name": "pending_jobs",
                                "label": "Pending",
                                "field": "pending_jobs",
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
                    
                    # Set default sorting by last activity in reverse chronological order
                    try:
                        # Try multiple approaches to set default sorting
                        self.samples_table.props("default-sort=last_seen desc")
                        # Alternative approach using sort property
                        self.samples_table.props("sort=last_seen desc")
                        # Set initial sort on the column
                        for col in self.samples_table.columns:
                            if col.get("field") == "last_seen":
                                col["sort"] = "desc"
                                break
                    except Exception:
                        pass
                    
                    # Populate table from master record if available (fast initial load)
                    try:
                        if self._samples_master_record:
                            logging.info("Populating table from master record cache")
                            self._refresh_table_from_master()
                    except Exception as e:
                        logging.debug(f"Error populating table from master record: {e}")

                    # Per-row action buttons (Finalize only shows for Complete/Pre-existing samples)
                    try:
                        self.samples_table.add_slot(
                            "body-cell-actions",
                            """
<q-td key=\"actions\" :props=\"props\">
  <q-btn-group>
    <q-btn color=\"primary\" size=\"sm\" label=\"View\"
           @click=\"$parent.$emit('action', props.row.sample_id)\" />
    <q-btn v-if=\"props.row.origin === 'Complete' || props.row.origin === 'Pre-existing'\"
           color=\"secondary\" size=\"sm\" label=\"Finalize\" icon=\"merge_type\"
           @click=\"$parent.$emit('finalize-target', props.row.sample_id)\" />
  </q-btn-group>
</q-td>
""",
                        )
                    except Exception:
                        pass
                    
                    # Add file progress column with linear progress bar
                    try:
                        self.samples_table.add_slot(
                            "body-cell-file_progress",
                            """
<q-td key=\"file_progress\" :props=\"props\">
  <div style=\"min-width: 120px;\">
    <q-linear-progress 
      :value=\"props.row.file_progress || 0\" 
      :color=\"props.row.file_progress >= 1 ? 'positive' : 'primary'\" 
      size=\"12px\" 
      rounded
      class=\"q-mb-xs\" />
    <div style=\"font-size: 10px; color: #666; text-align: center;\">
      {{ props.row.files_processed || 0 }}/{{ props.row.files_seen || 0 }} files
    </div>
  </div>
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
                            lambda e: (
                                ui.navigate.to(f"/live_data/{e.args}")
                                if isinstance(getattr(e, "args", None), str)
                                else ui.notify("Invalid button payload", type="warning")
                            ),
                        )
                    except Exception:
                        pass
                    
                    # Handle finalize-target action
                    def _on_finalize_target(event):
                        try:
                            sample_id = getattr(event, "args", None) if hasattr(event, "args") else None
                            if isinstance(sample_id, str):
                                self._trigger_target_bam_finalization(sample_id, trigger_snp=True)
                                ui.notify(
                                    f"Triggered run finalization and SNP analysis for {sample_id}",
                                    type="positive",
                                )
                            else:
                                ui.notify("Invalid sample ID", type="warning")
                        except Exception as e:
                            ui.notify(f"Error triggering finalization: {e}", type="negative")
                    
                    try:
                        self.samples_table.on("finalize-target", _on_finalize_target)
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
                                            for r in self.samples_table.rows or []:
                                                if r.get("sample_id") == sid:
                                                    r["export"] = (
                                                        sid in self._selected_sample_ids
                                                    )
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

                        async def _export_selected_reports(state: Dict[str, Any], progress_dialog, files_to_download, download_complete, progress_callback, progress_updates):
                            try:
                                selected = list(
                                    getattr(self, "_selected_sample_ids", set()) or []
                                )
                                if not selected:
                                    ui.notify("No samples selected", type="warning")
                                    return
                                
                                total_samples = len(selected)
                                
                                try:
                                    from nicegui import run as ng_run  # type: ignore
                                except Exception:
                                    ng_run = None  # type: ignore
                                
                                for idx, sid in enumerate(selected):
                                    try:
                                        # Update overall progress
                                        overall_progress = (idx / total_samples) * 0.9
                                        
                                        # Emit progress update showing current sample and mark as starting
                                        current_sample_msg = f"Generating {idx + 1}/{total_samples} - {sid}"
                                        progress_updates.put({
                                            'stage': 'processing_sections',
                                            'message': 'Starting...',
                                            'progress': 0.0,
                                            'sample_id': sid
                                        })
                                        
                                        sample_dir = (
                                            Path(self.monitored_directory) / sid
                                            if self.monitored_directory
                                            else None
                                        )
                                        if not sample_dir or not sample_dir.exists():
                                            logging.warning(f"Missing output for {sid}")
                                            continue
                                        
                                        # Don't use notification system - only update dialog
                                        
                                        filename = f"{sid}_run_report.pdf"
                                        pdf_path = os.path.join(
                                            str(sample_dir), filename
                                        )
                                        os.makedirs(str(sample_dir), exist_ok=True)
                                        export_csv_dir = None
                                        if bool(state.get("export_csv", False)):
                                            export_csv_dir = os.path.join(
                                                str(sample_dir), "report_csv"
                                            )
                                        
                                        # Don't use the notification system - use only our dialog callback
                                        if ng_run is not None:
                                            # Use custom callback that updates dialog only
                                            def sample_progress_callback(data: Dict[str, Any]):
                                                data['sample_id'] = sid  # Add sample ID to track which bar to update
                                                progress_callback(data)  # This updates the dialog via progress_updates queue
                                            
                                            pdf_file = await ng_run.io_bound(
                                                create_pdf,
                                                pdf_path,
                                                str(sample_dir),
                                                self.center or "Unknown",
                                                report_type=state.get("type", "detailed"),
                                                export_csv_dir=export_csv_dir,
                                                export_xlsx=False,
                                                export_zip=bool(
                                                    state.get("export_csv", False)
                                                ),
                                                progress_callback=sample_progress_callback,  # Pass the dialog-only callback
                                                workflow_steps=self.workflow_steps if hasattr(self, 'workflow_steps') else None,
                                            )
                                        else:
                                            # Use custom callback that updates dialog only
                                            def sample_progress_callback(data: Dict[str, Any]):
                                                data['sample_id'] = sid  # Add sample ID to track which bar to update
                                                progress_callback(data)  # This updates the dialog via progress_updates queue
                                            
                                            pdf_file = create_pdf(
                                                pdf_path,
                                                str(sample_dir),
                                                self.center or "Unknown",
                                                report_type=state.get("type", "detailed"),
                                                export_csv_dir=export_csv_dir,
                                                export_xlsx=False,
                                                export_zip=bool(
                                                    state.get("export_csv", False)
                                                ),
                                                progress_callback=sample_progress_callback,  # Pass the dialog-only callback
                                                workflow_steps=self.workflow_steps if hasattr(self, 'workflow_steps') else None,
                                            )
                                        
                                        # Queue file for download instead of downloading immediately
                                        files_to_download.append(pdf_file)
                                        
                                        # Also offer CSV ZIP if requested
                                        if (
                                            bool(state.get("export_csv", False))
                                            and export_csv_dir
                                        ):
                                            zip_path = os.path.join(
                                                export_csv_dir, f"{sid}_report_data.zip"
                                            )
                                            if os.path.exists(zip_path):
                                                files_to_download.append(zip_path)
                                        
                                        # Mark sample as complete
                                        progress_updates.put({
                                            'stage': 'completed',
                                            'message': 'Completed',
                                            'progress': 1.0,
                                            'sample_id': sid
                                        })
                                        
                                    except Exception as e:
                                        # Report generation failed
                                        logging.error(f"Export failed for {sid}: {e}")
                                        # Mark sample as failed
                                        progress_updates.put({
                                            'stage': 'error',
                                            'message': f'Failed: {str(e)[:50]}',
                                            'progress': 1.0,
                                            'sample_id': sid
                                        })
                                
                                # Mark as complete
                                download_complete["done"] = True
                                logging.info(f"Bulk export complete. {len(files_to_download)} file(s) ready for download.")
                            except Exception as e:
                                logging.error(f"Error in bulk export: {e}")
                                download_complete["done"] = True

                        async def _confirm_bulk_export():
                            # Ensure there is at least one selection before opening
                            if not getattr(self, "_selected_sample_ids", None):
                                ui.notify("No samples selected", type="warning")
                                return

                            selected_ids = list(getattr(self, "_selected_sample_ids", set()) or [])
                            num_selected = len(selected_ids)

                            report_types = {
                                "summary": "Summary Only",
                                "detailed": "Detailed",
                            }
                            state: Dict[str, Any] = {
                                "type": "detailed",
                                "export_csv": False,
                            }

                            with ui.dialog().props("persistent") as dialog:
                                with ui.card().classes("w-96 p-4"):
                                    ui.label("Export Reports").classes(
                                        "text-h6 font-bold mb-4"
                                    )

                                    with ui.column():
                                        # Report type selector
                                        with ui.column().classes("mb-4"):
                                            ui.label("Report Type").classes("font-bold mb-2")
                                            ui.toggle(
                                                report_types,
                                                value="detailed",
                                                on_change=lambda e: state.update(
                                                    {"type": e.value}
                                                ),
                                            )

                                        # Data export options
                                        with ui.column().classes("mb-4"):
                                            ui.label("Include Data").classes("font-bold mb-2")
                                            ui.checkbox(
                                                "CSV data (ZIP)",
                                                value=False,
                                                on_change=lambda e: state.update(
                                                    {"export_csv": bool(e.value)}
                                                ),
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

                                        ui.label(
                                            f"Are you sure you want to export reports for {num_selected} sample(s)?"
                                        ).classes("mb-4")

                                        # Buttons
                                        with ui.row().classes("justify-end gap-2"):
                                            ui.button(
                                                "Cancel",
                                                on_click=lambda: dialog.submit("Cancel"),
                                            ).props("flat")
                                            ui.button(
                                                "Export",
                                                on_click=lambda: dialog.submit("Export"),
                                            ).props("color=primary")

                                    # Display selected samples below buttons
                                    ui.label(
                                        f"Exporting reports for {num_selected} sample(s): {', '.join(selected_ids[:3])}"
                                        + (f" and {num_selected - 3} more" if num_selected > 3 else "")
                                    ).classes("text-sm font-medium text-gray-700 mt-4")

                            dialog_result = await dialog
                            if dialog_result != "Export":
                                return

                            # Now show the progress dialog
                            with ui.dialog().props("persistent") as progress_dialog:
                                with ui.card().classes("w-96 p-4"):
                                    # Title
                                    ui.label("Exporting Reports").classes(
                                        "text-h6 font-bold mb-4"
                                    )

                                    # Sample count
                                    ui.label(f"Exporting {num_selected} report(s)").classes(
                                        "text-sm font-medium text-gray-700 mb-4"
                                    )

                                    # Report type and output displays
                                    ui.label(
                                        f"Report Type: {state.get('type', 'detailed').title()}"
                                    ).classes("text-sm mb-2")
                                    
                                    ui.label(
                                        f"Output: PDF{' + CSV (ZIP)' if state.get('export_csv', False) else ''}"
                                    ).classes("text-sm mb-4")
                                    
                                    # Create individual progress bars for each sample
                                    sample_progress_bars = {}
                                    sample_progress_labels = {}
                                    
                                    with ui.column().classes("w-full"):
                                        for sid in selected_ids:
                                            with ui.column().classes("mb-3 w-full"):
                                                ui.label(sid).classes("text-xs font-medium text-gray-700 mb-1")
                                                progress_bar = ui.linear_progress(0.0).classes("mb-1")
                                                progress_label = ui.label("Waiting...").classes("text-xs text-gray-500")
                                                sample_progress_bars[sid] = progress_bar
                                                sample_progress_labels[sid] = progress_label
                                    
                                    # Messages container - use label with newlines for multiple messages
                                    messages_label = ui.label("").classes("text-xs text-gray-500")

                                    # Track messages
                                    progress_updates = queue.Queue()
                                    messages_list = []
                                    
                                    # Track current sample being processed
                                    current_sample = {"id": None}
                                    
                                    # Timer to process progress updates on UI thread
                                    def process_progress_updates():
                                        """Process queued progress updates."""
                                        try:
                                            while not progress_updates.empty():
                                                update = progress_updates.get_nowait()
                                                stage = update.get("stage", "unknown")
                                                message = update.get("message", "")
                                                progress = update.get("progress", 0.0)
                                                sample_id = update.get("sample_id")
                                                
                                                # Update the current sample's progress bar
                                                if sample_id and sample_id in sample_progress_bars:
                                                    if progress is not None:
                                                        sample_progress_bars[sample_id].value = progress
                                                        sample_progress_labels[sample_id].text = f"{int(progress * 100)}% - {message}"
                                                    else:
                                                        sample_progress_labels[sample_id].text = message
                                                    current_sample["id"] = sample_id
                                                
                                                # Add message to messages list
                                                messages_list.append(message)
                                                # Keep only last 10 messages
                                                if len(messages_list) > 10:
                                                    messages_list.pop(0)
                                                
                                                # Update messages label
                                                messages_label.text = "\n".join(messages_list[-5:])
                                                
                                        except queue.Empty:
                                            pass
                                        except Exception as e:
                                            logging.debug(f"Error processing progress updates: {e}")
                                    
                                    # Set up timer to process updates
                                    update_timer = ui.timer(0.1, process_progress_updates)
                                    
                                    def progress_callback(progress_data: Dict[str, Any]):
                                        """Custom progress callback to update dialog (called from background thread)."""
                                        try:
                                            # Queue the update instead of directly updating UI
                                            progress_updates.put(progress_data)
                                        except Exception as e:
                                            logging.debug(f"Error in progress callback: {e}")
                                    
                                    # Track if still generating
                                    is_generating = {"active": True}

                                    # Storage for the files to download
                                    files_to_download = []
                                    download_complete = {"done": False}
                                    
                                    # Timer to handle downloads once background task is done
                                    def handle_downloads():
                                        """Handle downloads in UI context once generation is complete."""
                                        if download_complete["done"] and files_to_download:
                                            # Log that export is complete
                                            file_count = len([f for f in files_to_download if f is not None])
                                            logging.info(f"Bulk export complete. {file_count} file(s) ready for download.")
                                            
                                            for file_path in files_to_download:
                                                if file_path is not None:
                                                    logging.debug(f"Initiating download: {file_path}")
                                                    ui.download(file_path)
                                            # Close dialog after 3 seconds
                                            ui.timer(3.0, lambda: progress_dialog.submit(None), once=True)
                                            download_timer.deactivate()
                                    
                                    download_timer = ui.timer(0.1, handle_downloads)
                                    
                                    # Start report generation
                                    async def complete_export():
                                        """Complete the export and close dialog."""
                                        try:
                                            await _export_selected_reports(state, progress_dialog, files_to_download, download_complete, progress_callback, progress_updates)
                                        finally:
                                            # Clean up timer
                                            update_timer.deactivate()
                                            is_generating["active"] = False
                                    
                                    asyncio.create_task(complete_export())

                            await progress_dialog

                        # Wire the button now that handlers exist
                        self.export_reports_button.on_click(_confirm_bulk_export)
                    except Exception:
                        pass

                    # If we have cached rows, show them immediately for instant UX
                    if self._last_samples_rows:
                        try:
                            self._apply_samples_table_filters()
                            # Hide loading container if we have data
                            if hasattr(self, "samples_loading_container"):
                                self.samples_loading_container.set_visibility(False)
                            # No selection behavior needed; per-row buttons handle navigation
                        except Exception:
                            pass
                    else:
                        # Show loading container if no cached data
                        if hasattr(self, "samples_loading_container"):
                            self.samples_loading_container.set_visibility(True)


    def _update_samples_table_sync(self, data: Dict[str, Any]) -> None:
        """Synchronous version of samples table update"""
        try:
            if not hasattr(self, "samples_table"):
                return
            samples = data.get("samples", [])
            expected_job_types = self._get_expected_completion_job_types()
            base = (
                Path(self.monitored_directory) if self.monitored_directory else None
            )

            # Deduplicate by sample_id taking the newest last_seen
            by_id: Dict[str, Dict[str, Any]] = {}
            for s in samples:
                sid = s.get("sample_id", "") or "unknown"
                last_seen = float(s.get("last_seen", time.time()))
                existing = by_id.get(sid)
                if not existing or last_seen >= existing.get("_last_seen_raw", 0):
                    origin_value = (
                        "Pre-existing"
                        if sid in self._preexisting_sample_ids and (time.time() - last_seen) >= self.completion_timeout_seconds
                        else "Live"
                    )
                    # Flip Live samples to Complete if inactive for configured timeout
                    prev_origin = None
                    if existing:
                        prev_origin = existing.get("origin")
                    else:
                        # Check previous rows for this sample
                        existing_rows = self._last_samples_rows or []
                        for prev_row in existing_rows:
                            if prev_row.get("sample_id") == sid:
                                prev_origin = prev_row.get("origin")
                                break
                    try:
                        active_jobs_count = s.get("active_jobs", 0)
                        pending_jobs_count = s.get("pending_jobs", 0)
                        # Only mark as Complete if timeout passed AND no active jobs
                        if origin_value == "Live" and (time.time() - last_seen) >= self.completion_timeout_seconds:
                            if active_jobs_count == 0 and pending_jobs_count == 0:
                                should_complete = True
                                if expected_job_types and base is not None:
                                    sample_dir = base / sid
                                    should_complete = self._expected_jobs_completed(
                                        sample_dir, expected_job_types
                                    )
                                if should_complete:
                                    origin_value = "Complete"
                                    # Trigger finalization if transitioning from Live to Complete
                                    if prev_origin == "Live" or prev_origin is None:
                                        self._trigger_target_bam_finalization(sid)
                            # If there are active jobs, keep as Live even if timeout passed
                    except Exception:
                        pass
                    # Get job counts
                    total_jobs = s.get("total_jobs", 0)
                    completed_jobs = s.get("completed_jobs", 0)
                    
                    # Set file progress directly from job counts (same data source as other columns)
                    files_seen = total_jobs
                    files_processed = completed_jobs
                    file_progress = completed_jobs / total_jobs if total_jobs > 0 else 0.0
                    
                    by_id[sid] = {
                        "sample_id": sid,
                        "origin": origin_value,
                        # Persisted run info will be patched in below from master.csv
                        "run_start": "",
                        "device": "",
                        "flowcell": "",
                        "active_jobs": s.get("active_jobs", 0),
                        "pending_jobs": s.get("pending_jobs", 0),
                        "total_jobs": total_jobs,
                        "completed_jobs": completed_jobs,
                        "failed_jobs": s.get("failed_jobs", 0),
                        "job_types": (
                            ",".join(sorted(set(s.get("job_types", []))))
                            if isinstance(s.get("job_types", []), list)
                            else str(s.get("job_types", ""))
                        ),
                        "last_seen": time.strftime(
                            "%Y-%m-%d %H:%M:%S", time.localtime(last_seen)
                        ),
                        "files_seen": files_seen,
                        "files_processed": files_processed,
                        "file_progress": file_progress,
                        "actions": "View",
                        "_last_seen_raw": last_seen,
                    }

            # Patch from master.csv and persist the new overview values for later reload
            try:
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
                            "pending_jobs": int(row.get("pending_jobs", 0)),
                            "total_jobs": int(row.get("total_jobs", 0)),
                            "completed_jobs": int(row.get("completed_jobs", 0)),
                            "failed_jobs": int(row.get("failed_jobs", 0)),
                            "job_types": row.get("job_types", ""),
                            "last_seen": float(row.get("_last_seen_raw", time.time())),
                        }
                        manager.update_sample_overview(sid, persist_payload)

                    # Read run info and file progress from master.csv to display in table
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
                                row["flowcell"] = first_row.get(
                                    "run_info_flow_cell", ""
                                )
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
            
            # Update UI on main thread
            if rows and hasattr(self, "samples_table"):
                self._apply_samples_table_filters()
                # Also update file progress display in workflow monitor
                self._update_file_progress()

        except Exception as e:
            logging.error(f"Error in samples table update: {e}")

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

    def _extract_event_value(self, event, default=""):
        """Extract value from NiceGUI event arguments."""
        try:
            # Handle different event types
            if hasattr(event, "value"):
                return str(event.value or default)
            elif hasattr(event, "args"):
                if isinstance(event.args, dict):
                    return str(event.args.get("value", event.args.get("label", default)))
                elif isinstance(event.args, list) and len(event.args) > 0:
                    return str(event.args[0] or default)
                else:
                    return str(event.args or default)
            else:
                return str(default)
        except Exception:
            return str(default)

    def _set_samples_query(self, query: str) -> None:
        try:
            # Ensure we have a filters dict
            if not hasattr(self, "_samples_filters"):
                self._samples_filters = {"query": "", "origin": "All"}
            
            # Normalize the query
            normalized_query = (query or "").strip().lower()
            self._samples_filters["query"] = normalized_query
            
            # Apply filters immediately
            self._apply_samples_table_filters()
            
            # Debug logging
            logging.debug(f"Search query updated: '{normalized_query}'")
        except Exception as e:
            logging.error(f"Error setting samples query: {e}")
            pass

    def _set_samples_origin_filter(self, origin_value: str) -> None:
        try:
            # Ensure we have a filters dict
            if not hasattr(self, "_samples_filters"):
                self._samples_filters = {"query": "", "origin": "All"}
            
            # Normalize the origin value
            normalized_origin = origin_value or "All"
            self._samples_filters["origin"] = normalized_origin
            
            # Apply filters immediately
            self._apply_samples_table_filters()
            
            # Debug logging
            logging.debug(f"Origin filter updated: '{normalized_origin}'")
        except Exception as e:
            logging.error(f"Error setting samples origin filter: {e}")
            pass

    def _apply_samples_table_filters(self) -> None:
        try:
            base_rows = getattr(self, "_last_samples_rows", []) or []
            rows = self._normalize_rows_for_display(base_rows)
            
            logging.debug(f"Applying filters to {len(rows)} base rows")

            # Origin filter, compute dynamic 'Complete' for display if needed
            now_ts = time.time()
            for r in rows:
                    try:
                        if r.get("origin") == "Live":
                            last_raw = float(r.get("_last_seen_raw", 0))
                            active_jobs_count = r.get("active_jobs", 0)
                            pending_jobs_count = r.get("pending_jobs", 0)
                            # Only mark as Complete if timeout passed AND no active jobs
                            if last_raw and (now_ts - last_raw) >= self.completion_timeout_seconds:
                                if active_jobs_count == 0 and pending_jobs_count == 0:
                                    r["origin"] = "Complete"
                                # If there are active jobs, keep as Live even if timeout passed
                    except Exception:
                        pass

            origin = (self._samples_filters or {}).get("origin", "All")
            if origin and origin != "All":
                rows = [r for r in rows if (r.get("origin") == origin)]
                logging.debug(f"After origin filter ({origin}): {len(rows)} rows")

            # Global query filter
            q = (self._samples_filters or {}).get("query", "")
            if q:
                ql = q.lower()
                logging.debug(f"Applying search filter: '{ql}'")

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
                logging.debug(f"After search filter: {len(rows)} rows")

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

            # Sort rows by last activity in reverse chronological order (newest first)
            try:
                rows.sort(key=lambda r: float(r.get("_last_seen_raw", 0)), reverse=True)
            except Exception:
                # Fallback to sorting by last_seen string if _last_seen_raw is not available
                try:
                    rows.sort(key=lambda r: r.get("last_seen", ""), reverse=True)
                except Exception:
                    pass

            # Update the table
            if hasattr(self, "samples_table"):
                self.samples_table.rows = rows
                self.samples_table.update()
                logging.debug(f"Updated samples table with {len(rows)} filtered rows")
            else:
                logging.warning("samples_table not found, cannot update")
                
        except Exception as e:
            logging.error(f"Error applying samples table filters: {e}")
            pass

    def _refresh_sample_plots(self, sample_id: str):
        """Refresh all plots for a specific sample by clearing all cached state."""
        try:
            if not self.monitored_directory:
                return

            sample_dir = Path(self.monitored_directory) / sample_id
            if not sample_dir.exists():
                return

            # Clear ALL cached state to force complete regeneration
            key = str(sample_dir)

            # Completely clear coverage state for this sample
            if hasattr(self, "_coverage_state") and key in self._coverage_state:
                old_state = self._coverage_state[key]
                logging.info(
                    f"Clearing coverage state for {sample_id}: {list(old_state.keys())}"
                )
                self._coverage_state[key] = {}
                logging.info(f"Cleared all coverage state for sample {sample_id}")

            # Clear other component states if they exist
            for attr in ["_mgmt_state", "_cnv_state", "_fusion_state"]:
                if hasattr(self, attr):
                    state_dict = getattr(self, attr)
                    if key in state_dict:
                        old_state = state_dict[key]
                        logging.info(
                            f"Clearing {attr} state for {sample_id}: {list(old_state.keys())}"
                        )
                        # Clear all state for these components
                        state_dict[key] = {}
                        logging.info(f"Cleared all {attr} state for sample {sample_id}")

        except Exception as e:
            logging.debug(f"Failed to refresh sample plots for {sample_id}: {e}")

    '''
    def _refresh_all_sample_plots(self):
        """Refresh all sample plots by resetting modification time caches."""
        try:
            if not self.monitored_directory:
                return

            # Get all sample directories
            monitored_path = Path(self.monitored_directory)
            if not monitored_path.exists():
                return

            sample_dirs = [d for d in monitored_path.iterdir() if d.is_dir()]

            for sample_dir in sample_dirs:
                sample_id = sample_dir.name
                self._refresh_sample_plots(sample_id)

            ui.notify(
                f"Refreshed plots for {len(sample_dirs)} samples", type="positive"
            )
            logging.info(f"Refreshed plots for {len(sample_dirs)} samples")

        except Exception as e:
            ui.notify(f"Failed to refresh all plots: {e}", type="negative")
            logging.exception("Failed to refresh all sample plots")
    '''
    
    def _create_sample_detail_page(self, sample_id: str):
        """Create the individual sample detail page."""
        
        # Check if this is a page refresh vs navigation
        # Use browser storage to detect if this is a refresh
        is_page_refresh = False
        try:
            # Check if we have a flag indicating this page was recently loaded
            if hasattr(ui, 'storage') and hasattr(ui.storage, 'browser'):
                last_load_time = ui.storage.browser.get(f'sample_{sample_id}_last_load', 0)
                current_time = time.time()
                # If last load was very recent (< 2 seconds), it's likely a page refresh
                is_page_refresh = (current_time - last_load_time) < 2.0
                
                # Update the last load time
                ui.storage.browser[f'sample_{sample_id}_last_load'] = current_time
        except Exception:
            # If storage is not available, assume it's navigation (show loading)
            is_page_refresh = False
        
        # Show loading spinner unless it's a page refresh
        show_loading = not is_page_refresh
        
        async def confirm_report_generation():
            """Show a confirmation dialog before generating the report."""
            report_types = {
                "summary": "Summary Only",
                "detailed": "Detailed",
            }
            state: Dict[str, Any] = {
                "type": "detailed",
                "export_csv": False,
            }

            with ui.dialog().props("persistent") as dialog:
                with ui.card().classes("w-96 p-4"):
                    # Title
                    title_label = ui.label("Generate Report").classes(
                        "text-h6 font-bold mb-4"
                    )

                    # Container for initial form content
                    with ui.column():
                        # Report type selector
                        with ui.column().classes("mb-4"):
                            ui.label("Report Type").classes("font-bold mb-2")
                            type_toggle = ui.toggle(
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
                            csv_checkbox = ui.checkbox(
                                "CSV data (ZIP)",
                                value=False,
                                on_change=lambda e: state.update(
                                    {"export_csv": bool(e.value)}
                                ),
                            )

                        ui.label(
                            "Are you sure you want to generate a report?"
                        ).classes("mb-4")

                        # Buttons
                        with ui.row().classes("justify-end gap-2"):
                            ui.button(
                                "No", on_click=lambda: dialog.submit("No")
                            ).props("flat")
                            ui.button(
                                "Yes",
                                on_click=lambda: dialog.submit("Yes"),
                            ).props("color=primary")
                    
                    # Run name below buttons
                    ui.label(f"Exporting data for run: {sample_id}").classes(
                        "text-sm font-medium text-gray-700 mt-4"
                    )

            dialog_result = await dialog
            if dialog_result != "Yes":
                return

            # Now show the progress dialog
            with ui.dialog().props("persistent") as progress_dialog:
                with ui.card().classes("w-96 p-4"):
                    # Title
                    ui.label("Generating Report").classes(
                        "text-h6 font-bold mb-2"
                    )
                    
                    # Divider
                    ui.separator().classes("mb-4")
                    
                    # Run name
                    ui.label(f"Exporting data for run: {sample_id}").classes(
                        "text-sm font-medium text-gray-700 mb-4"
                    )

                    # Report type and output displays
                    report_type_display = ui.label(
                        f"Report Type: {state.get('type', 'detailed').title()}"
                    ).classes("text-sm mb-2")
                    
                    output_display = ui.label(
                        f"Output: PDF{' + CSV (ZIP)' if state.get('export_csv', False) else ''}"
                    ).classes("text-sm mb-4")
                    
                    # Progress indicator
                    progress_bar = ui.linear_progress(0.0).classes("mb-2")
                    
                    # Progress text
                    progress_text = ui.label("Preparing...").classes(
                        "text-xs text-gray-600 text-center mb-2"
                    )
                    
                    # Messages container - use label with newlines for multiple messages
                    messages_label = ui.label("").classes("text-xs text-gray-500")

                    # Track messages
                    progress_updates = queue.Queue()
                    messages_list = []
                    
                    # Timer to process progress updates on UI thread
                    def process_progress_updates():
                        """Process queued progress updates."""
                        try:
                            while not progress_updates.empty():
                                update = progress_updates.get_nowait()
                                stage = update.get("stage", "unknown")
                                message = update.get("message", "")
                                progress = update.get("progress", 0.0)
                                
                                if progress is not None:
                                    progress_bar.value = progress
                                    progress_text.text = f"{int(progress * 100)}% - {message}"
                                else:
                                    progress_text.text = message
                                
                                # Add message to messages list
                                messages_list.append(message)
                                # Keep only last 10 messages
                                if len(messages_list) > 10:
                                    messages_list.pop(0)
                                
                                # Update messages label
                                messages_label.text = "\n".join(messages_list[-5:])
                                
                        except queue.Empty:
                            pass
                        except Exception as e:
                            logging.debug(f"Error processing progress updates: {e}")
                    
                    # Set up timer to process updates
                    update_timer = ui.timer(0.1, process_progress_updates)
                    
                    def progress_callback(progress_data: Dict[str, Any]):
                        """Custom progress callback to update dialog (called from background thread)."""
                        try:
                            # Queue the update instead of directly updating UI
                            progress_updates.put(progress_data)
                        except Exception as e:
                            logging.debug(f"Error in progress callback: {e}")
                    
                    # Track if still generating
                    is_generating = {"active": True}

                    # Storage for the files to download (will be populated by background task)
                    files_to_download = []
                    download_complete = {"done": False}
                    
                    # Timer to handle downloads once background task is done
                    def handle_downloads():
                        """Handle downloads in UI context once generation is complete."""
                        if download_complete["done"] and files_to_download:
                            # Log that report generation is complete and downloads are available
                            file_count = len([f for f in files_to_download if f is not None])
                            logging.info(f"Report generation complete for {sample_id}. {file_count} file(s) ready for download.")
                            
                            for file_path in files_to_download:
                                if file_path is not None:
                                    logging.debug(f"Initiating download: {file_path}")
                                    ui.download(file_path)
                            # Close dialog after 3 seconds
                            ui.timer(3.0, lambda: progress_dialog.submit(None), once=True)
                            download_timer.deactivate()
                    
                    download_timer = ui.timer(0.1, handle_downloads)
                    
                    # Start report generation
                    async def complete_generation():
                        """Complete the report generation and close dialog."""
                        try:
                            await generate_and_download_report(state, progress_callback, progress_dialog, is_generating, files_to_download)
                        finally:
                            # Clean up timer
                            update_timer.deactivate()
                            is_generating["active"] = False
                            download_complete["done"] = True
                    
                    asyncio.create_task(complete_generation())

            await progress_dialog

        async def generate_and_download_report(state: Dict[str, Any], progress_callback, progress_dialog, is_generating, files_to_download):
            """Generate report and update progress in dialog."""
            try:
                from nicegui import run as ng_run  # type: ignore
                
                if not sample_dir or not sample_dir.exists():
                    # Queue error notification
                    files_to_download.append(None)  # Signal error
                    ui.timer(0.1, lambda: ui.notify(
                        "Output directory not available for this sample",
                        type="warning",
                    ), once=True)
                    return
                
                filename = f"{sample_id}_run_report.pdf"
                pdf_path = os.path.join(str(sample_dir), filename)
                os.makedirs(str(sample_dir), exist_ok=True)
                export_csv_dir = None
                if bool(state.get("export_csv", False)):
                    export_csv_dir = os.path.join(
                        str(sample_dir), "report_csv"
                    )
                
                # Use only our custom callback, not the notification system
                def combined_callback(progress_data: Dict[str, Any]):
                    """Custom callback for dialog updates only (no notifications)."""
                    progress_callback(progress_data)
                
                pdf_file = await ng_run.io_bound(
                    create_pdf,
                    pdf_path,
                    str(sample_dir),
                    self.center or "Unknown",
                    report_type=state.get("type", "detailed"),
                    export_csv_dir=export_csv_dir,
                    export_xlsx=False,
                    export_zip=bool(state.get("export_csv", False)),
                    progress_callback=combined_callback,
                    workflow_steps=self.workflow_steps if hasattr(self, 'workflow_steps') else None,
                )
                
                # Mark report as completed
                from robin.gui.report_progress import progress_manager
                progress_manager.complete_report(sample_id, filename)
                
                # Queue files for download in UI context
                files_to_download.append(pdf_file)
                
                # Also offer CSV ZIP if requested
                if bool(state.get("export_csv", False)) and export_csv_dir:
                    zip_path = os.path.join(
                        export_csv_dir, f"{sample_id}_report_data.zip"
                    )
                    if os.path.exists(zip_path):
                        files_to_download.append(zip_path)
                
            except Exception as e:
                # Mark report as failed
                from robin.gui.report_progress import progress_manager
                progress_manager.error_report(sample_id, str(e))
                
                # Queue error notification in UI context
                ui.timer(0.1, lambda: ui.notify(
                    f"Error generating report: {str(e)}",
                    type="negative",
                ), once=True)
                files_to_download.append(None)  # Signal error

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
                
                # Show initial notification with more explicit styling
                ui.notify(
                    f"[{sample_id}] Starting report generation...",
                    type="info",
                    timeout=0,  # Persistent notification
                    position="top-right"
                )
                
                filename = f"{sample_id}_run_report.pdf"
                pdf_path = os.path.join(str(sample_dir), filename)
                os.makedirs(str(sample_dir), exist_ok=True)
                export_csv_dir = None
                if bool(state.get("export_csv", False)):
                    export_csv_dir = os.path.join(
                        str(sample_dir), "report_csv"
                    )
                
                # Create progress callback
                from robin.gui.report_progress import create_progress_callback
                progress_callback = create_progress_callback(sample_id)
                
                pdf_file = await ng_run.io_bound(
                    create_pdf,
                    pdf_path,
                    str(sample_dir),
                    self.center or "Unknown",
                    report_type=state.get("type", "detailed"),
                    export_csv_dir=export_csv_dir,
                    export_xlsx=False,
                    export_zip=bool(state.get("export_csv", False)),
                    progress_callback=progress_callback,
                    workflow_steps=self.workflow_steps if hasattr(self, 'workflow_steps') else None,
                )
                
                # Mark report as completed
                from robin.gui.report_progress import progress_manager
                progress_manager.complete_report(sample_id, filename)
                
                ui.download(pdf_file)
                # Also offer CSV ZIP if requested
                if bool(state.get("export_csv", False)) and export_csv_dir:
                    zip_path = os.path.join(
                        export_csv_dir, f"{sample_id}_report_data.zip"
                    )
                    if os.path.exists(zip_path):
                        ui.download(zip_path)
                
            except Exception as e:
                # Mark report as failed
                from robin.gui.report_progress import progress_manager
                progress_manager.error_report(sample_id, str(e))
                
        with theme.frame(
            f"R.O.B.I.N - Sample {sample_id}",
            smalltitle=sample_id,
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
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
                    ui.label(
                        "Redirecting to sample list..."
                    ).classes("text-sm text-gray-500 mt-2")
                    ui.button(
                        "Back to Samples", on_click=lambda: ui.navigate.to("/live_data")
                    ).props("color=primary")
                
                # Soft redirect after short delay
                def redirect_to_samples():
                    """Redirect to the sample list page."""
                    try:
                        ui.navigate.to("/live_data")
                    except Exception as e:
                        logging.warning(f"Failed to redirect to /live_data: {e}")
                        # Fallback: try using JavaScript redirect
                        try:
                            ui.run_javascript('window.location.href = "/live_data"')
                        except Exception as e2:
                            logging.error(f"Failed JavaScript redirect: {e2}")
                
                # Create and activate the timer
                redirect_timer = ui.timer(2.0, redirect_to_samples, once=True)
                redirect_timer.activate()
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
            
            # Check directory existence asynchronously to avoid blocking
            async def check_directory_and_notify():
                try:
                    if not sample_dir or not sample_dir.exists():
                        ui.notify(
                            f"Sample output directory not found for {sample_id}",
                            type="warning",
                        )
                except Exception:
                    pass
            
            # Start directory check in background
            ui.timer(0.1, check_directory_and_notify, once=True)

            with ui.card().classes("w-full").style("border: 2px solid var(--md-primary)"):
                with ui.row().classes("w-full flex justify-between items-center"):
                    with ui.column():
                        ui.label(f"{sample_id}").classes("text-2xl font-bold")
                        ui.label("Detailed sample information.").classes(
                            "text-sm ml-4 opacity-80"
                        )
                    with ui.column().classes("flex gap-2 items-center"):
                        # Check if target.bam exists to show "More Details" button
                        target_bam_exists = False
                        if sample_dir and sample_dir.exists():
                            target_bam = sample_dir / "target.bam"
                            target_bam_exists = target_bam.exists() and target_bam.is_file()
                        
                        if target_bam_exists:
                            ui.button(
                                "More Details",
                                on_click=lambda: ui.navigate.to(f"/live_data/{sample_id}/details")
                            ).classes(
                                "text-sm font-semibold px-3 py-1 rounded bg-secondary text-white"
                            )
                        ui.button(
                            "Generate Report", on_click=confirm_report_generation
                        ).classes(
                            "object-right ml-auto text-sm font-semibold px-3 py-1 rounded"
                        )
                ui.separator().classes().style("border: 1px solid var(--md-primary)")
                # Main content area with conditional loading state
                with ui.column().classes("w-full p-4 gap-4"):
                    if show_loading:
                        # Loading container that will be hidden when data is ready
                        loading_container = ui.column().classes(
                            "w-full items-center justify-center p-8"
                        )
                        with loading_container:
                            ui.spinner("bars", size="4em").classes("mb-4")
                            ui.label("Loading sample data...").classes("text-lg text-gray-600")
                            ui.label("This may take a moment for large datasets").classes(
                                "text-sm text-gray-500"
                            )

                        # Content container that will be shown when data is ready
                        content_container = (
                            ui.column().classes("w-full gap-4").style("display: none")
                        )
                    else:
                        # For page refreshes, show content immediately
                        content_container = ui.column().classes("w-full gap-4")
                        loading_container = None

                    with content_container:
                        # Summary section (new component) - create UI immediately, load data async
                        try:
                            try:
                                from .gui.components.summary import add_summary_section  # type: ignore
                            except ImportError:
                                # Try absolute import if relative fails
                                from robin.gui.components.summary import add_summary_section

                            # Create the UI components immediately on the main thread
                            add_summary_section(sample_dir, sample_id, self)
                        except Exception as e:
                            logging.exception(f"[GUI] Summary section failed: {e}")
                            try:
                                ui.notify(f"Summary section failed: {e}", type="warning")
                            except Exception:
                                pass

                        # Import section visibility helper
                        try:
                            from robin.gui.config import is_section_enabled, get_enabled_classification_steps
                        except ImportError:
                            is_section_enabled = lambda name, steps: True
                            get_enabled_classification_steps = lambda steps: {"sturgeon", "nanodx", "random_forest", "pannanodx"}
                        
                        workflow_steps = self.workflow_steps if hasattr(self, 'workflow_steps') else None

                        # Defer heavy analysis sections until after initial render
                        analysis_container = ui.column().classes("w-full gap-4")

                        def _build_analysis_sections():
                            try:
                                analysis_container.clear()
                                with analysis_container:
                                    # MNP-Flex section
                                    try:
                                        try:
                                            from .gui.components.mnpflex import add_mnpflex_section  # type: ignore
                                        except ImportError:
                                            from robin.gui.components.mnpflex import add_mnpflex_section

                                        add_mnpflex_section(self, sample_dir, sample_id)
                                    except Exception as e:
                                        logging.exception(f"[GUI] MNP-Flex section failed: {e}")
                                        try:
                                            ui.notify(
                                                f"MNP-Flex section failed: {e}", type="warning"
                                            )
                                        except Exception:
                                            pass

                                    # Classification section
                                    enabled_classification_steps = get_enabled_classification_steps(workflow_steps)
                                    if not workflow_steps or enabled_classification_steps:
                                        try:
                                            try:
                                                from .gui.components.classification import add_classification_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.classification import (
                                                    add_classification_section,
                                                )

                                            add_classification_section(sample_dir, self)
                                        except Exception as e:
                                            logging.exception(f"[GUI] Classification section failed: {e}")
                                            try:
                                                ui.notify(
                                                    f"Classification section failed: {e}", type="warning"
                                                )
                                            except Exception:
                                                pass

                                    # Coverage section (target)
                                    if not workflow_steps or is_section_enabled("target", workflow_steps):
                                        try:
                                            try:
                                                from .gui.components.coverage import add_coverage_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.coverage import (
                                                    add_coverage_section,
                                                    add_igv_viewer,
                                                )

                                            add_coverage_section(self, sample_dir)
                                        except Exception as e:
                                            try:
                                                ui.notify(f"Coverage section failed: {e}", type="warning")
                                            except Exception:
                                                pass
                                    
                                    # MGMT section
                                    if not workflow_steps or is_section_enabled("mgmt", workflow_steps):
                                        try:
                                            try:
                                                from .gui.components.mgmt import add_mgmt_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.mgmt import add_mgmt_section

                                            add_mgmt_section(self, sample_dir)
                                        except Exception as e:
                                            logging.exception(f"[GUI] MGMT section failed: {e}")
                                            try:
                                                ui.notify(f"MGMT section failed: {e}", type="warning")
                                            except Exception:
                                                pass
                                    
                                    # CNV section
                                    if not workflow_steps or is_section_enabled("cnv", workflow_steps):
                                        try:
                                            try:
                                                from .gui.components.cnv import add_cnv_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.cnv import add_cnv_section

                                            add_cnv_section(self, sample_dir)
                                        except Exception as e:
                                            logging.exception(f"[GUI] CNV section failed: {e}")
                                            try:
                                                ui.notify(f"CNV section failed: {e}", type="warning")
                                            except Exception:
                                                pass

                                    # Fusion section + BED coverage
                                    if not workflow_steps or is_section_enabled("fusion", workflow_steps):
                                        try:
                                            try:
                                                from .gui.components.fusion import add_fusion_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.fusion import add_fusion_section

                                            add_fusion_section(self, sample_dir)
                                        except Exception as e:
                                            logging.exception(f"[GUI] Fusion section failed: {e}")
                                            try:
                                                ui.notify(f"Fusion section failed: {e}", type="warning")
                                            except Exception:
                                                pass
                                        
                                        try:
                                            try:
                                                from .gui.components.bed_coverage import add_bed_coverage_section  # type: ignore
                                            except ImportError:
                                                from robin.gui.components.bed_coverage import add_bed_coverage_section
                                            
                                            add_bed_coverage_section(self, sample_dir)
                                        except Exception as e:
                                            logging.exception(f"[GUI] BED Coverage section failed: {e}")
                                            try:
                                                ui.notify(f"BED Coverage section failed: {e}", type="warning")
                                            except Exception:
                                                pass
                            except Exception as e:
                                logging.exception(f"[GUI] Failed to build analysis sections: {e}")

                        # Delay heavy UI creation to allow websocket handshake to complete
                        ui.timer(0.5, _build_analysis_sections, once=True)

                        # Files in output directory
                        with ui.card().classes("w-full"):
                            ui.label("Output Files").classes(
                                "text-lg font-semibold mb-2"
                            )
                            with ui.row().classes("items-center gap-3 mb-2"):
                                files_search = ui.input("Search files…").props(
                                    "borderless dense clearable"
                                )
                            from robin.gui.theme import styled_table

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
                        {
                            "name": "actions",
                            "label": "Download",
                            "field": "actions",
                            "sortable": False,
                            "align": "center",
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
                            
                            # Add download button slot that emits an event to trigger Python download
                            try:
                                files_table.add_slot(
                                    "body-cell-actions",
                                    """
<q-td key="actions" :props="props">
  <q-btn color="primary" size="sm" icon="download" 
         @click="() => $parent.$emit('download-file', props.row.name)" 
         title="Download file" />
</q-td>
""",
                                )
                            except Exception:
                                pass
                            
                            # Add event handler for download button clicks
                            try:
                                files_table.on('download-file', lambda event: _download_file(event.args))
                            except Exception:
                                pass

                        # master.csv summary
                        with ui.card().classes("w-full"):
                            ui.label("master.csv Summary").classes(
                                "text-lg font-semibold mb-2"
                            )
                            with ui.row().classes("items-center gap-3 mb-2"):
                                summary_search = ui.input("Search fields…").props(
                                    "borderless dense clearable"
                                )
                            from robin.gui.theme import styled_table

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

                    # Download method for individual files
                    def _download_file(filename: str):
                        """Download a single file from the sample directory."""
                        try:
                            if not sample_dir or not sample_dir.exists():
                                ui.notify("Sample directory not found", type="error")
                                return

                            file_path = sample_dir / filename
                            if not file_path.exists() or not file_path.is_file():
                                ui.notify(f"File {filename} not found", type="error")
                                return

                            with open(file_path, 'rb') as f:
                                content = f.read()

                            ui.download(
                                content,
                                filename=filename,
                                media_type='application/octet-stream'
                            )

                        except Exception as e:
                            ui.notify(f"Download failed: {e}", type="error")

                    # Periodic refresher for files table and master.csv summary
                    _notify_state = {"files_error": False, "csv_error": False}

                    def _refresh_files_list_sync() -> List[Dict[str, Any]]:
                        """Synchronous file list refresh - runs in background thread"""
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
                                                "actions": f.name,  # Store filename for actions
                                            }
                                        )
                                    except Exception:
                                        continue
                        return rows

                    def _refresh_csv_summary_sync() -> List[Dict[str, Any]]:
                        """Synchronous CSV summary refresh - runs in background thread"""
                        if not sample_dir:
                            return [{"key": "Status", "value": "No sample directory"}]
                            
                        csv_path = sample_dir / "master.csv"
                        if not csv_path.exists():
                            return [{"key": "Status", "value": "master.csv not found"}]
                            
                        try:
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
                                            {
                                                "key": k,
                                                "value": first_row.get(k, ""),
                                            }
                                        )
                                if not rows2:
                                    rows2 = [
                                        {"key": k, "value": v}
                                        for k, v in first_row.items()
                                    ]
                                return rows2
                        except Exception:
                            return [{"key": "Error", "value": "Failed to read CSV"}]

                    async def _refresh_sample_detail_async() -> None:
                        """Asynchronous version of sample detail refresh"""
                        try:
                            # Run file operations in background threads
                            #import concurrent.futures
                            # Run file I/O operations in background threads to avoid blocking GUI
                            import asyncio
                            files_result = await asyncio.to_thread(_refresh_files_list_sync)
                            csv_result = await asyncio.to_thread(_refresh_csv_summary_sync)
                            
                            
                            # Update UI with results
                            files_table.rows = files_result
                            files_table.update()
                            
                            summary_table.rows = csv_result
                            summary_table.update()
                            
                            # Reset error states on success
                            _notify_state["files_error"] = False
                            _notify_state["csv_error"] = False
                            
                        except Exception as e:
                            logging.error(f"Error in async sample detail refresh: {e}")
                            # Only show error notification once per error type
                            if not _notify_state["files_error"]:
                                try:
                                    ui.notify(
                                        f"Failed to refresh sample data for {sample_id}: {e}",
                                        type="warning",
                                    )
                                except Exception:
                                    pass
                                _notify_state["files_error"] = True

                    def _refresh_sample_detail() -> None:
                        """Synchronous version - kept for backward compatibility"""
                        # Refresh files list
                        try:
                            rows = _refresh_files_list_sync()
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

                        # Refresh master.csv summary (non-blocking)
                        def update_csv_summary():
                            try:
                                rows2 = _refresh_csv_summary_sync()
                                summary_table.rows = rows2
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
                        
                        # Run in background thread to avoid blocking GUI
                        import threading
                        thread = threading.Thread(target=update_csv_summary, daemon=True)
                        thread.start()

                    # Show content and hide loading after initial data load (only when showing loading)
                    if show_loading and loading_container:
                        async def _load_initial_data_and_show():
                            """Load initial data asynchronously then show content"""
                            try:
                                # Load initial data
                                await _refresh_sample_detail_async()
                                
                                # Show content and hide loading
                                loading_container.style("display: none")
                                content_container.style("display: flex")
                                
                            except Exception as e:
                                logging.error(f"Error loading initial data: {e}")
                                # Show content anyway to avoid infinite loading
                                loading_container.style("display: none")
                                content_container.style("display: flex")

                        # Start initial data loading
                        try:
                            ui.timer(0.1, _load_initial_data_and_show, once=True)
                        except Exception:
                            # Fallback: show content immediately if timer fails
                            loading_container.style("display: none")
                            content_container.style("display: flex")
                    else:
                        # For page refreshes, load data immediately
                        try:
                            ui.timer(0.1, _refresh_sample_detail_async, once=True)
                        except Exception:
                            pass

                    # Start periodic refresh with async version
                    try:
                        refresh_timer = ui.timer(30.0, _refresh_sample_detail_async)
                        try:
                            ui.context.client.on_disconnect(
                                lambda: refresh_timer.deactivate()
                            )
                        except Exception:
                            pass
                    except Exception:
                        pass

    def _create_sample_details_page(self, sample_id: str):
        """Create the sample details page with comprehensive information."""
        
        # Get sample directory
        sample_dir = (
            Path(self.monitored_directory) / sample_id
            if self.monitored_directory
            else None
        )
        
        # Check if sample is known
        if self._known_sample_ids and sample_id not in self._known_sample_ids:
            with theme.frame(
                f"R.O.B.I.N - Sample Details",
                smalltitle="Details",
                batphone=False,
                center=self.center,
                setup_notifications=self._setup_notification_system,
            ):
                with ui.column().classes("w-full items-center justify-center p-8"):
                    ui.label(f"Unknown sample: {sample_id}").classes(
                        "text-xl font-semibold text-red-600"
                    )
                    ui.label(
                        "This sample ID has not been seen yet in the current session."
                    ).classes("text-sm text-gray-600")
                    ui.button(
                        "Back to Sample", 
                        on_click=lambda: ui.navigate.to(f"/live_data/{sample_id}")
                    ).props("color=primary").classes("mt-4")
                    ui.button(
                        "Back to Samples", 
                        on_click=lambda: ui.navigate.to("/live_data")
                    ).classes("mt-2")
            return
        
        # Create the page with theme frame
        with theme.frame(
            f"R.O.B.I.N - Sample Details: {sample_id}",
            smalltitle=f"{sample_id} Details",
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            # Main content - full width like sample page
            with ui.column().classes("w-full p-4 gap-4"):
                # Header section
                with ui.row().classes("w-full flex justify-between items-center flex-wrap gap-2"):
                    with ui.column():
                        ui.label(f"Sample Details: {sample_id}").classes("text-2xl font-bold")
                        ui.label("IGV Viewer, sample information and analysis results.").classes(
                            "text-sm ml-4 opacity-80"
                        )
                    with ui.column():
                        ui.button(
                            "Back to Sample", 
                            on_click=lambda: ui.navigate.to(f"/live_data/{sample_id}")
                        ).classes("bg-primary text-white rounded-md mobile-button")
                
                ui.separator().classes().style("border: 1px solid var(--md-primary)")
                
                # Placeholder content
                with ui.column().classes("w-full gap-4"):
                    ui.label("Sample Details Page").classes("text-headline-medium text-center")
                    # Sample information section
                    with ui.card().classes("w-full elevation-2 rounded-lg p-4"):
                        ui.label("Sample Information").classes("text-headline-small font-bold mb-4")
                        
                        if sample_dir and sample_dir.exists():
                            ui.label(f"Sample Directory: {sample_dir}").classes("text-body-medium mb-2")
                            ui.label("Status: Directory found").classes("text-body-medium text-green-600")
                        else:
                            ui.label("Status: Directory not found").classes("text-body-medium text-red-600")
                            if sample_dir:
                                ui.label(f"Expected path: {sample_dir}").classes("text-body-small text-gray-500")
                    
                    # Center information
                    if self.center:
                        with ui.card().classes("w-full elevation-2 rounded-lg p-4"):
                            ui.label("Analysis Center").classes("text-headline-small font-bold mb-4")
                            ui.label(f"Center: {self.center}").classes("text-body-medium")
                    
                    # IGV Viewer section - moved to top, before tables
                    if sample_dir and sample_dir.exists():
                        from robin.gui.components.coverage import add_igv_viewer
                        add_igv_viewer(self, sample_dir)
                    
                    # SNP Analysis section
                    if sample_dir and sample_dir.exists():
                        from robin.gui.components.snp import add_snp_section
                        add_snp_section(self, sample_dir)
                    
                    # Fusion Pairs Table section
                    if sample_dir and sample_dir.exists():
                        from robin.gui.components.fusion import (
                            _load_processed_pickle,
                            _cluster_fusion_reads
                        )
                        import pandas as pd
                        
                        # Load fusion data
                        fusion_data_loaded = False
                        fusion_data = None
                        try:
                            target_file = sample_dir / "fusion_candidates_master_processed.pkl"
                            genome_file = sample_dir / "fusion_candidates_all_processed.pkl"
                            
                            # Try to load target panel data first, then genome-wide
                            if target_file.exists():
                                fusion_data = _load_processed_pickle(target_file)
                                if fusion_data and fusion_data.get("annotated_data") is not None:
                                    fusion_data_loaded = True
                            elif genome_file.exists():
                                fusion_data = _load_processed_pickle(genome_file)
                                if fusion_data and fusion_data.get("annotated_data") is not None:
                                    fusion_data_loaded = True
                        except Exception as e:
                            logging.warning(f"Failed to load fusion data: {e}")
                        
                        if fusion_data_loaded and fusion_data:
                            annotated_data = fusion_data.get("annotated_data", pd.DataFrame())
                            goodpairs = fusion_data.get("goodpairs", pd.Series())
                            
                            if not annotated_data.empty:
                                # Filter to good pairs if available
                                if not goodpairs.empty and goodpairs.sum() > 0:
                                    aligned_goodpairs = goodpairs.reindex(annotated_data.index, fill_value=False)
                                    filtered_data = annotated_data[aligned_goodpairs]
                                else:
                                    filtered_data = annotated_data
                                
                                # Cluster fusion reads
                                clustered_data = _cluster_fusion_reads(
                                    filtered_data, 
                                    max_distance=10000, 
                                    use_breakpoint_validation=True
                                )
                                
                                if not clustered_data.empty:
                                    with ui.card().classes("w-full elevation-2 rounded-lg p-4"):
                                        ui.label("Fusion Pairs").classes("text-headline-small font-bold mb-4")
                                        ui.label("Click on a row to view the fusion pair in IGV.").classes(
                                            "text-body-small text-gray-600 mb-2"
                                        )
                                        
                                        # Create columns for the table
                                        from robin.gui.theme import styled_table
                                        
                                        columns = [
                                            {"name": "fusion_pair", "label": "Fusion Pair", "field": "fusion_pair", "sortable": True},
                                            {"name": "chr1", "label": "Chr 1", "field": "chr1", "sortable": True},
                                            {"name": "pos1", "label": "Breakpoint 1", "field": "pos1", "sortable": True},
                                            {"name": "chr2", "label": "Chr 2", "field": "chr2", "sortable": True},
                                            {"name": "pos2", "label": "Breakpoint 2", "field": "pos2", "sortable": True},
                                            {"name": "reads", "label": "Supporting Reads", "field": "reads", "sortable": True},
                                            {"name": "action", "label": "View in IGV", "field": "action", "sortable": False}
                                        ]
                                        
                                        # Format rows for display
                                        rows = []
                                        fusion_region_map = {}  # Store region info for each row
                                        
                                        # Import re for position parsing
                                        import re
                                        
                                        for idx, row in clustered_data.iterrows():
                                            # Try to get start/end coordinates from the row if available
                                            # (e.g., from breakpoint validation)
                                            if all(col in row for col in ["gene1_start", "gene1_end", "gene2_start", "gene2_end"]):
                                                # Use actual start/end coordinates if available
                                                start1_raw = int(row["gene1_start"])
                                                end1_raw = int(row["gene1_end"])
                                                start2_raw = int(row["gene2_start"])
                                                end2_raw = int(row["gene2_end"])
                                            else:
                                                # Parse from position strings (format: "start-end")
                                                pos1_str = str(row.get("gene1_position", ""))
                                                pos2_str = str(row.get("gene2_position", ""))
                                                
                                                # Parse range format "start-end" or just a single number
                                                pos1_match = re.match(r'(\d+)[-–—](\d+)', pos1_str.replace(',', ''))
                                                pos2_match = re.match(r'(\d+)[-–—](\d+)', pos2_str.replace(',', ''))
                                                
                                                if pos1_match and pos2_match:
                                                    # Extract both start and end from the range
                                                    start1_raw = int(pos1_match.group(1))
                                                    end1_raw = int(pos1_match.group(2))
                                                    start2_raw = int(pos2_match.group(1))
                                                    end2_raw = int(pos2_match.group(2))
                                                else:
                                                    # Fallback: try to extract single coordinates
                                                    pos1_single = re.search(r'(\d+)', pos1_str.replace(',', ''))
                                                    pos2_single = re.search(r'(\d+)', pos2_str.replace(',', ''))
                                                    if pos1_single and pos2_single:
                                                        # Single coordinate - use it as both start and end
                                                        start1_raw = end1_raw = int(pos1_single.group(1))
                                                        start2_raw = end2_raw = int(pos2_single.group(1))
                                                    else:
                                                        # Skip this row if we can't parse coordinates
                                                        continue
                                            
                                            # Use the actual range (start to end), then add/subtract 10kb padding
                                            padding = 10000
                                            min1 = min(start1_raw, end1_raw)
                                            max1 = max(start1_raw, end1_raw)
                                            min2 = min(start2_raw, end2_raw)
                                            max2 = max(start2_raw, end2_raw)
                                            
                                            # Subtract 10kb from min and add 10kb to max
                                            start1 = max(1, min1 - padding)
                                            end1 = max1 + padding
                                            start2 = max(1, min2 - padding)
                                            end2 = max2 + padding
                                            
                                            chr1 = str(row.get("chr1", "Unknown"))
                                            chr2 = str(row.get("chr2", "Unknown"))
                                            
                                            # Format region as "chr1:start-end chr2:start-end"
                                            region = f"{chr1}:{start1}-{end1} {chr2}:{start2}-{end2}"
                                            
                                            # For display, show the breakpoint range (not the padded version)
                                            # Use a single representative coordinate or the range midpoint
                                            display_pos1 = f"{min1:,}-{max1:,}" if min1 != max1 else f"{min1:,}"
                                            display_pos2 = f"{min2:,}-{max2:,}" if min2 != max2 else f"{min2:,}"
                                            
                                            formatted_row = {
                                                "fusion_pair": row.get("fusion_pair", ""),
                                                "chr1": chr1,
                                                "pos1": display_pos1,
                                                "chr2": chr2,
                                                "pos2": display_pos2,
                                                "reads": int(row.get("reads", 0)),
                                                "action": "🔍"  # Icon for viewing in IGV
                                            }
                                            
                                            rows.append(formatted_row)
                                            fusion_region_map[len(rows) - 1] = region
                                        
                                        if rows:
                                            # Store fusion regions mapped by fusion pair for easy lookup
                                            fusion_regions_by_pair = {}
                                            for idx, row_data in enumerate(rows):
                                                fusion_regions_by_pair[row_data["fusion_pair"]] = fusion_region_map[idx]
                                            
                                            # Create JavaScript map of regions for IGV navigation
                                            import json
                                            js_regions_json = json.dumps(fusion_regions_by_pair)
                                            
                                            # Function to navigate IGV to a fusion region
                                            def navigate_to_fusion_region(fusion_pair: str):
                                                """Navigate IGV browser to the specified fusion pair region."""
                                                if fusion_pair in fusion_regions_by_pair:
                                                    region = fusion_regions_by_pair[fusion_pair]
                                                    # Escape region string for JavaScript
                                                    escaped_region = region.replace('"', '\\"').replace("'", "\\'")
                                                    js_navigate = f"""
                                                        (function() {{
                                                            try {{
                                                                if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                    console.log('[IGV] Navigating to fusion region: {escaped_region}');
                                                                    window.lj_igv.search('{escaped_region}');
                                                                }} else {{
                                                                    console.warn('[IGV] Browser not ready yet, will navigate when ready');
                                                                    setTimeout(function() {{
                                                                        if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                            window.lj_igv.search('{escaped_region}');
                                                                        }}
                                                                    }}, 500);
                                                                }}
                                                            }} catch (error) {{
                                                                console.error('[IGV] Navigation error:', error);
                                                            }}
                                                        }})();
                                                    """
                                                    ui.run_javascript(js_navigate, timeout=5.0)
                                            
                                            # Initialize fusion regions map in window BEFORE creating table
                                            # This ensures it's available when the slot template renders
                                            js_init_regions = f"""
                                                (function() {{
                                                    window.fusionRegionsMap = window.fusionRegionsMap || {{}};
                                                    Object.assign(window.fusionRegionsMap, {js_regions_json});
                                                    console.log('[Fusion] Initialized fusion regions map with', Object.keys(window.fusionRegionsMap).length, 'regions');
                                                }})();
                                            """
                                            ui.run_javascript(js_init_regions, timeout=5.0)
                                            
                                            # Create styled table
                                            table_container, fusion_table = styled_table(
                                                columns=columns,
                                                rows=rows,
                                                pagination=20,
                                                class_size="table-xs"
                                            )
                                            
                                            # Add clickable action button column using slot that emits events to Python
                                            try:
                                                # Use Vue event emission which is more reliable than window functions
                                                fusion_table.add_slot(
                                                    "body-cell-action",
                                                    """
<q-td key="action" :props="props">
  <q-btn 
    icon="visibility" 
    size="sm" 
    dense 
    flat 
    color="primary"
    @click="$parent.$emit('fusion-view-igv', props.row.fusion_pair)"
    title="View in IGV"
  />
</q-td>
"""
                                                )
                                                
                                                # Handle the event from the slot
                                                def on_fusion_view_igv(e):
                                                    """Handle fusion view IGV event from table button."""
                                                    try:
                                                        fusion_pair = e.args if isinstance(e.args, str) else getattr(e, 'args', None)
                                                        if fusion_pair:
                                                            logging.debug(f"[Fusion] Button clicked for: {fusion_pair}")
                                                            navigate_to_fusion_region(fusion_pair)
                                                    except Exception as ex:
                                                        logging.warning(f"Error handling fusion view IGV event: {ex}")
                                                
                                                fusion_table.on("fusion-view-igv", on_fusion_view_igv)
                                                logging.debug("Added action button column slot with event handler")
                                            except Exception as slot_ex:
                                                logging.warning(f"Could not add action column slot: {slot_ex}")
                                                import traceback
                                                logging.warning(traceback.format_exc())
                                                
                                                # Fallback: Use JavaScript with inline region lookup
                                                try:
                                                    # Create a safer inline handler that doesn't rely on window functions
                                                    js_inline_handler = f"""
                                                        (function() {{
                                                            // Store regions as a constant in the closure
                                                            const fusionRegions = {js_regions_json};
                                                            
                                                            // Create handler function immediately
                                                            window.handleFusionNav = function(fusionPair) {{
                                                                console.log('[Fusion] Navigation handler called for:', fusionPair);
                                                                const region = fusionRegions[fusionPair];
                                                                if (region) {{
                                                                    console.log('[Fusion] Navigating to:', region);
                                                                    function nav() {{
                                                                        if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                            window.lj_igv.search(region);
                                                                            return true;
                                                                        }}
                                                                        return false;
                                                                    }}
                                                                    if (!nav()) {{
                                                                        setTimeout(nav, 500);
                                                                        setTimeout(nav, 1500);
                                                                    }}
                                                                }}
                                                            }};
                                                            console.log('[Fusion] Created navigation handler');
                                                        }})();
                                                    """
                                                    ui.run_javascript(js_inline_handler, timeout=5.0)
                                                    
                                                    # Wait a moment for JS to execute, then add slot
                                                    ui.timer(0.1, lambda: None, once=True)
                                                    
                                                    fusion_table.add_slot(
                                                        "body-cell-action",
                                                        """
<q-td key="action" :props="props">
  <q-btn 
    icon="visibility" 
    size="sm" 
    dense 
    flat 
    color="primary"
    @click="window.handleFusionNav && window.handleFusionNav(props.row.fusion_pair)"
    title="View in IGV"
  />
</q-td>
"""
                                                    )
                                                except Exception as fallback_ex:
                                                    logging.warning(f"Fallback approach also failed: {fallback_ex}")
                                            
                                            # Use JavaScript to attach click handlers to rows - improved with event delegation
                                            js_attach_handlers = f"""
                                                (function() {{
                                                    let attached = false;
                                                    
                                                    function attachFusionHandlers() {{
                                                        if (attached) return;
                                                        
                                                        console.log('[Fusion] Attaching click handlers...');
                                                        
                                                        // Use event delegation on the document body for better reliability
                                                        function handleFusionRowClick(e) {{
                                                            // Check if click is on a fusion table row
                                                            let target = e.target;
                                                            let row = target.closest('tbody tr');
                                                            
                                                            if (!row) return;
                                                            
                                                            // Check if this row is in a fusion table
                                                            const table = row.closest('table');
                                                            if (!table) return;
                                                            
                                                            const headers = table.querySelectorAll('thead th');
                                                            let isFusionTable = false;
                                                            for (let i = 0; i < headers.length; i++) {{
                                                                if ((headers[i].textContent || '').toLowerCase().includes('fusion pair')) {{
                                                                    isFusionTable = true;
                                                                    break;
                                                                }}
                                                            }}
                                                            
                                                            if (!isFusionTable) return;
                                                            
                                                            // Skip if clicking on the action button (let button handle it)
                                                            if (target.closest('button') || target.closest('.q-btn')) {{
                                                                return;
                                                            }}
                                                            
                                                            e.preventDefault();
                                                            e.stopPropagation();
                                                            
                                                            console.log('[Fusion] Row clicked');
                                                            
                                                            // Get fusion pair from first cell
                                                            const cells = row.querySelectorAll('td');
                                                            if (cells.length > 0) {{
                                                                const fusionPair = cells[0].textContent.trim();
                                                                console.log('[Fusion] Fusion pair:', fusionPair);
                                                                
                                                                if (fusionPair && window.fusionRegionsMap && window.fusionRegionsMap[fusionPair]) {{
                                                                    const region = window.fusionRegionsMap[fusionPair];
                                                                    console.log('[Fusion] Navigating to region:', region);
                                                                    
                                                                    // Visual feedback
                                                                    row.style.backgroundColor = '#e3f2fd';
                                                                    setTimeout(function() {{
                                                                        row.style.backgroundColor = '';
                                                                    }}, 400);
                                                                    
                                                                    // Navigate IGV
                                                                    function navIGV() {{
                                                                        if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                            try {{
                                                                                window.lj_igv.search(region);
                                                                                console.log('[Fusion] IGV navigation successful');
                                                                                return true;
                                                                            }} catch(err) {{
                                                                                console.error('[Fusion] IGV error:', err);
                                                                                return false;
                                                                            }}
                                                                        }}
                                                                        return false;
                                                                    }}
                                                                    
                                                                    if (!navIGV()) {{
                                                                        setTimeout(function() {{ navIGV(); }}, 500);
                                                                        setTimeout(function() {{ navIGV(); }}, 1500);
                                                                    }}
                                                                }}
                                                            }}
                                                        }}
                                                        
                                                        // Attach event listener to document (event delegation)
                                                        document.addEventListener('click', handleFusionRowClick, true);
                                                        
                                                        // Also make rows visually clickable
                                                        const tables = document.querySelectorAll('table');
                                                        tables.forEach(function(table) {{
                                                            const headers = table.querySelectorAll('thead th');
                                                            let isFusionTable = false;
                                                            for (let i = 0; i < headers.length; i++) {{
                                                                if ((headers[i].textContent || '').toLowerCase().includes('fusion pair')) {{
                                                                    isFusionTable = true;
                                                                    break;
                                                                }}
                                                            }}
                                                            
                                                            if (isFusionTable) {{
                                                                const tbody = table.querySelector('tbody');
                                                                if (tbody) {{
                                                                    const rows = tbody.querySelectorAll('tr');
                                                                    rows.forEach(function(row) {{
                                                                        row.style.cursor = 'pointer';
                                                                    }});
                                                                    console.log('[Fusion] Made', rows.length, 'rows clickable');
                                                                }}
                                                            }}
                                                        }});
                                                        
                                                        attached = true;
                                                        console.log('[Fusion] Click handlers attached successfully');
                                                    }}
                                                    
                                                    // Try multiple times to ensure table is rendered
                                                    attachFusionHandlers();
                                                    setTimeout(attachFusionHandlers, 300);
                                                    setTimeout(attachFusionHandlers, 800);
                                                    setTimeout(attachFusionHandlers, 1500);
                                                    setTimeout(attachFusionHandlers, 2500);
                                                }})();
                                            """
                                            
                                            # Execute JavaScript after table is created
                                            ui.timer(0.5, lambda: ui.run_javascript(js_attach_handlers, timeout=10.0), once=True)
                                            ui.timer(2.0, lambda: ui.run_javascript(js_attach_handlers, timeout=10.0), once=True)
                                            
                                            # Add summary
                                            total_fusions = len(rows)
                                            total_reads = sum(r["reads"] for r in rows)
                                            ui.label(
                                                f"Total fusions: {total_fusions} | Total supporting reads: {total_reads}"
                                            ).classes("text-xs text-gray-500 mt-2")
                                        else:
                                            ui.label("No fusion pairs found").classes("text-gray-600")
                                else:
                                    with ui.card().classes("w-full elevation-2 rounded-lg p-4"):
                                        ui.label("Fusion Pairs").classes("text-headline-small font-bold mb-4")
                                        ui.label("No validated fusion pairs found").classes("text-gray-600")
                            else:
                                with ui.card().classes("w-full elevation-2 rounded-lg p-4"):
                                    ui.label("Fusion Pairs").classes("text-headline-small font-bold mb-4")
                                    ui.label("Fusion data is empty").classes("text-gray-600")
                        else:
                            with ui.card().classes("w-full elevation-2 rounded-lg p-4"):
                                ui.label("Fusion Pairs").classes("text-headline-small font-bold mb-4")
                                ui.label("Fusion analysis data not available. Please run fusion analysis first.").classes(
                                    "text-gray-600"
                                )
                    
                    # Target Genes Table section
                    if sample_dir and sample_dir.exists():
                        target_coverage_file = sample_dir / "target_coverage.csv"
                        bed_coverage_file = sample_dir / "bed_coverage_main.csv"
                        
                        # Try target_coverage.csv first, fallback to bed_coverage_main.csv
                        coverage_file = None
                        if target_coverage_file.exists():
                            coverage_file = target_coverage_file
                        elif bed_coverage_file.exists():
                            coverage_file = bed_coverage_file
                        
                        if coverage_file:
                            # Load coverage data asynchronously to avoid blocking GUI
                            def load_coverage_data():
                                try:
                                    import pandas as pd
                                    return pd.read_csv(coverage_file)
                                except Exception as e:
                                    logging.error(f"Failed to load coverage file {coverage_file}: {e}")
                                    return None
                            
                            # Create placeholder table that will be updated when data loads
                            coverage_table_placeholder = ui.table(
                                columns=[
                                    {"name": "gene", "label": "Gene", "field": "gene"},
                                    {"name": "chrom", "label": "Chrom", "field": "chrom"},
                                    {"name": "startpos", "label": "Start", "field": "startpos"},
                                    {"name": "endpos", "label": "End", "field": "endpos"},
                                    {"name": "coverage", "label": "Coverage", "field": "coverage"},
                                ],
                                rows=[],
                            ).classes("w-full")
                            
                            # Load coverage data synchronously (fast enough for typical file sizes)
                            # Note: This runs during page creation, not during user interaction
                            try:
                                import pandas as pd
                                df = pd.read_csv(coverage_file)
                                
                                # Ensure we have the required columns
                                required_cols = ["chrom", "startpos", "endpos", "name"]
                                if all(col in df.columns for col in required_cols):
                                    # Calculate coverage if not present
                                    if "coverage" not in df.columns:
                                        if "length" in df.columns and "bases" in df.columns:
                                            df["coverage"] = df["bases"] / df["length"]
                                        elif "startpos" in df.columns and "endpos" in df.columns:
                                            df["length"] = df["endpos"] - df["startpos"] + 1
                                            if "bases" in df.columns:
                                                df["coverage"] = df["bases"] / df["length"]
                                            else:
                                                df["coverage"] = 0
                                    
                                    # Prepare table data
                                    table_data = []
                                    gene_regions_by_name = {}  # Store regions for navigation
                                    for _, row in df.iterrows():
                                        gene_name = str(row["name"])
                                        chrom = str(row["chrom"])
                                        startpos_raw = int(row["startpos"])
                                        endpos_raw = int(row["endpos"])
                                        
                                        # Store region info for navigation (with 10kb padding)
                                        padding = 10000
                                        startpos_nav = max(1, startpos_raw - padding)
                                        endpos_nav = endpos_raw + padding
                                        region = f"{chrom}:{startpos_nav}-{endpos_nav}"
                                        gene_regions_by_name[gene_name] = region
                                        
                                        table_data.append({
                                            "chrom": chrom,
                                            "startpos": f"{startpos_raw:,}",  # Format with commas for display
                                            "endpos": f"{endpos_raw:,}",  # Format with commas for display
                                            "name": gene_name,
                                            "coverage": float(row.get("coverage", 0)),
                                            "action": "🔍"  # Icon for viewing in IGV
                                        })
                                    
                                    if table_data:
                                        with ui.card().classes("w-full elevation-2 rounded-lg p-4"):
                                            ui.label("Target Genes").classes("text-headline-small font-bold mb-4")
                                            ui.label("Click on a gene row to view it in IGV.").classes(
                                                "text-body-small text-gray-600 mb-2"
                                            )
                                            
                                            from robin.gui.theme import styled_table
                                            
                                            columns = [
                                                {"name": "name", "label": "Gene Name", "field": "name", "sortable": True},
                                                {"name": "chrom", "label": "Chromosome", "field": "chrom", "sortable": True},
                                                {"name": "startpos", "label": "Start", "field": "startpos", "sortable": True},
                                                {"name": "endpos", "label": "End", "field": "endpos", "sortable": True},
                                                {"name": "coverage", "label": "Coverage (x)", "field": "coverage", "sortable": True},
                                                {"name": "action", "label": "View in IGV", "field": "action", "sortable": False}
                                            ]
                                            
                                            table_container, gene_table = styled_table(
                                                columns=columns,
                                                rows=table_data,
                                                pagination=25,
                                                class_size="table-xs"
                                            )
                                            
                                            # Add search functionality
                                            try:
                                                with gene_table.add_slot("top-right"):
                                                    with ui.input(placeholder="Search genes...").props("type=search").bind_value(
                                                        gene_table, "filter"
                                                    ).add_slot("append"):
                                                        ui.icon("search")
                                            except Exception:
                                                pass
                                            
                                            # Add colored coverage badges
                                            try:
                                                gene_table.add_slot(
                                                    "body-cell-coverage",
                                                    """
    <q-td key="coverage" :props="props">
    <q-badge :color="props.value >= 30 ? 'green' : props.value >= 20 ? 'blue' : props.value >= 10 ? 'orange' : 'red'">
        {{ Number(props.value).toFixed(2) }}
    </q-badge>
    </q-td>
    """,
                                                )
                                            except Exception:
                                                pass
                                            
                                            # Function to navigate IGV to a gene region
                                            import json
                                            js_gene_regions_json = json.dumps(gene_regions_by_name)
                                            
                                            def navigate_to_gene_region(gene_name: str):
                                                """Navigate IGV browser to the specified gene region."""
                                                if gene_name in gene_regions_by_name:
                                                    region = gene_regions_by_name[gene_name]
                                                    # Escape region string for JavaScript
                                                    escaped_region = region.replace('"', '\\"').replace("'", "\\'")
                                                    js_navigate = f"""
                                                        (function() {{
                                                            try {{
                                                                if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                    console.log('[IGV] Navigating to gene region: {escaped_region}');
                                                                    window.lj_igv.search('{escaped_region}');
                                                                }} else {{
                                                                    console.warn('[IGV] Browser not ready yet, will navigate when ready');
                                                                    setTimeout(function() {{
                                                                        if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                            window.lj_igv.search('{escaped_region}');
                                                                        }}
                                                                    }}, 500);
                                                                }}
                                                            }} catch (error) {{
                                                                console.error('[IGV] Navigation error:', error);
                                                            }}
                                                        }})();
                                                    """
                                                    ui.run_javascript(js_navigate, timeout=5.0)
                                            
                                            # Store regions in window object for JavaScript access
                                            js_init_gene_regions = f"""
                                                window.geneRegionsMap = window.geneRegionsMap || {{}};
                                                Object.assign(window.geneRegionsMap, {js_gene_regions_json});
                                                console.log('[Gene] Loaded', Object.keys(window.geneRegionsMap).length, 'gene regions');
                                            """
                                            ui.run_javascript(js_init_gene_regions, timeout=5.0)
                                            
                                            # Add clickable action button column using slot that emits events to Python
                                            try:
                                                # Use Vue event emission which is more reliable than window functions
                                                gene_table.add_slot(
                                                    "body-cell-action",
                                                    """
<q-td key="action" :props="props">
  <q-btn 
    icon="visibility" 
    size="sm" 
    dense 
    flat 
    color="primary"
    @click="$parent.$emit('gene-view-igv', props.row.name)"
    title="View in IGV"
  />
</q-td>
"""
                                                )
                                                
                                                # Handle the event from the slot
                                                def on_gene_view_igv(e):
                                                    """Handle gene view IGV event from table button."""
                                                    try:
                                                        gene_name = e.args if isinstance(e.args, str) else getattr(e, 'args', None)
                                                        if gene_name:
                                                            logging.debug(f"[Gene] Button clicked for: {gene_name}")
                                                            navigate_to_gene_region(gene_name)
                                                    except Exception as ex:
                                                        logging.warning(f"Error handling gene view IGV event: {ex}")
                                                
                                                gene_table.on("gene-view-igv", on_gene_view_igv)
                                                logging.debug("Added action button column slot with event handler")
                                            except Exception as slot_ex:
                                                logging.warning(f"Could not add action column slot: {slot_ex}")
                                                import traceback
                                                logging.warning(traceback.format_exc())
                                            
                                            # Also add row click handlers as a fallback
                                            try:
                                                # Create JavaScript to handle row clicks
                                                js_gene_table_handlers = f"""
                                                    (function() {{
                                                        // Store regions in window for JavaScript access
                                                        window.geneRegionsMap = window.geneRegionsMap || {{}};
                                                        Object.assign(window.geneRegionsMap, {js_gene_regions_json});
                                                        
                                                        function attachGeneTableHandlers() {{
                                                            // Find tables with "Gene Name" header
                                                            const tables = document.querySelectorAll('table');
                                                            
                                                            tables.forEach(function(table) {{
                                                                const headers = table.querySelectorAll('thead th');
                                                                let isGeneTable = false;
                                                                
                                                                for (let i = 0; i < headers.length; i++) {{
                                                                    const headerText = (headers[i].textContent || '').toLowerCase();
                                                                    if (headerText.includes('gene name') && headerText.includes('chromosome')) {{
                                                                        isGeneTable = true;
                                                                        break;
                                                                    }}
                                                                }}
                                                                
                                                                if (isGeneTable) {{
                                                                    const tbody = table.querySelector('tbody');
                                                                    if (tbody) {{
                                                                        const rows = tbody.querySelectorAll('tr');
                                                                        rows.forEach(function(row) {{
                                                                            // Skip if clicking on the action button (let button handle it)
                                                                            if (row.hasAttribute('data-gene-handled')) return;
                                                                            row.setAttribute('data-gene-handled', 'true');
                                                                            
                                                                            // Make row clickable (but button takes precedence)
                                                                            row.style.cursor = 'pointer';
                                                                            
                                                                            row.onclick = function(e) {{
                                                                                // Don't handle if clicking on button
                                                                                if (e.target.closest('button') || e.target.closest('.q-btn')) {{
                                                                                    return;
                                                                                }}
                                                                                
                                                                                e.preventDefault();
                                                                                e.stopPropagation();
                                                                                
                                                                                const cells = row.querySelectorAll('td');
                                                                                if (cells.length >= 5) {{
                                                                                    const geneName = cells[0].textContent.trim();
                                                                                    
                                                                                    if (geneName && window.geneRegionsMap && window.geneRegionsMap[geneName]) {{
                                                                                        const region = window.geneRegionsMap[geneName];
                                                                                        console.log('[Gene] Row clicked, navigating to region:', region);
                                                                                        
                                                                                        // Visual feedback
                                                                                        const origBg = row.style.backgroundColor;
                                                                                        row.style.backgroundColor = '#e3f2fd';
                                                                                        setTimeout(function() {{
                                                                                            row.style.backgroundColor = origBg;
                                                                                        }}, 400);
                                                                                        
                                                                                        // Navigate IGV
                                                                                        function navIGV() {{
                                                                                            if (window.lj_igv && window.lj_igv_browser_ready) {{
                                                                                                try {{
                                                                                                    window.lj_igv.search(region);
                                                                                                    console.log('[Gene] IGV navigation successful');
                                                                                                    return true;
                                                                                                }} catch(err) {{
                                                                                                    console.error('[Gene] IGV error:', err);
                                                                                                    return false;
                                                                                                }}
                                                                                            }}
                                                                                            return false;
                                                                                        }}
                                                                                        
                                                                                        if (!navIGV()) {{
                                                                                            setTimeout(function() {{ navIGV(); }}, 500);
                                                                                            setTimeout(function() {{ navIGV(); }}, 1500);
                                                                                        }}
                                                                                    }}
                                                                                }}
                                                                            }};
                                                                        }});
                                                                    }}
                                                                }}
                                                            }});
                                                        }}
                                                        
                                                        attachGeneTableHandlers();
                                                        setTimeout(attachGeneTableHandlers, 300);
                                                        setTimeout(attachGeneTableHandlers, 800);
                                                        setTimeout(attachGeneTableHandlers, 1500);
                                                    }})();
                                                """
                                                
                                                ui.timer(0.5, lambda: ui.run_javascript(js_gene_table_handlers, timeout=10.0), once=True)
                                                ui.timer(2.0, lambda: ui.run_javascript(js_gene_table_handlers, timeout=10.0), once=True)
                                            except Exception as e:
                                                logging.warning(f"Could not add gene table click handlers: {e}")
                                            
                                            # Add summary
                                            total_genes = len(table_data)
                                            ui.label(f"Total target genes: {total_genes}").classes("text-xs text-gray-500 mt-2")
                            except Exception as e:
                                logging.warning(f"Could not load target gene table: {e}")

    def _create_workflow_monitor(self):
        """Create the main workflow monitoring page."""
        with theme.frame(
            "R.O.B.I.N - robin Workflow Monitor",
            smalltitle="Samples",
            batphone=False,
            center=self.center,
            setup_notifications=self._setup_notification_system,
        ):
            # Page title and navigation
            with ui.row().classes("w-full p-4 items-center justify-between"):
                with ui.row().classes("items-center"):
                    ui.label("robin Workflow Monitor").classes("text-2xl font-bold")
                    ui.label("Real-time workflow monitoring and control").classes(
                        "text-sm ml-4 opacity-80"
                    )

            # Main content area
            with ui.column().classes("w-full p-4 gap-4"):
                # Workflow Status Overview
                with ui.card().classes(
                    "w-full bg-gradient-to-r from-blue-50 to-indigo-50"
                ):
                    ui.label("Workflow Status Overview").classes(
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

                # File Processing Progress (Per-Run)
                with ui.card().classes("w-full bg-gradient-to-r from-green-50 to-emerald-50"):
                    ui.label("File Processing Progress (Per Run)").classes(
                        "text-lg font-semibold mb-4 text-green-800"
                    )
                    with ui.column().classes("w-full gap-3"):
                        # Overall file progress summary
                        with ui.row().classes("w-full items-center gap-4"):
                            ui.label("Overall Files:").classes("text-sm font-medium text-gray-700")
                            self.overall_files_progress = ui.linear_progress(0.0).classes("flex-1")
                            self.overall_files_label = ui.label("0/0 files processed").classes(
                                "text-sm text-gray-600 min-w-[120px]"
                            )
                        
                        # Per-sample file progress will be shown in a scrollable container
                        ui.label("Per-Sample Progress:").classes("text-sm font-medium text-gray-700 mt-2")
                        with ui.scroll_area().classes("w-full").style("max-height: 300px;"):
                            self.sample_files_progress_container = ui.column().classes("w-full gap-2 p-2")
                            # Individual sample progress bars will be added here dynamically

                # Queue Status
                with ui.card().classes("w-full"):
                    ui.label("Queue Status").classes("text-lg font-semibold mb-4")

                    # Queue status grid
                    with ui.grid(columns=4).classes("w-full gap-4"):
                        # Preprocessing
                        with ui.card().classes("bg-green-50 p-4"):
                            ui.label("Preprocessing").classes(
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
                            ui.label("Classification").classes(
                                "text-sm font-medium text-purple-800"
                            )
                            self.classification_status = ui.label("0/0").classes(
                                "text-2xl font-bold text-purple-800"
                            )

                        # Other
                        with ui.card().classes("bg-gray-50 p-4"):
                            ui.label("Other").classes(
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
                    ui.label("Active Jobs").classes("text-lg font-semibold mb-4")

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
                    from robin.gui.theme import styled_table

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
                    self.no_active_jobs_label = ui.label(
                        "No active jobs at the moment."
                    ).classes("text-sm text-gray-500 mt-2")

                # Live Logs
                with ui.card().classes("w-full"):
                    ui.label("Live Logs").classes("text-lg font-semibold mb-4")

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
                    ui.label("Workflow Configuration").classes(
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
                    ui.label("Error Summary & Troubleshooting").classes(
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
                    ui.label("robin Workflow Monitor - Running").classes(
                        "text-sm text-gray-600"
                    )

                # Signal that GUI is ready to receive updates
                self.gui_ready.set()
                logging.info("[GUI] UI created and ready to receive updates")
                # Start a periodic UI-thread drain of the update queue
                app.timer(0.3, self._drain_updates_on_ui, active=True)
                # Start duration refresher
                app.timer(1.0, lambda: self._refresh_duration(), active=True)

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

    def _locate_reference_for_sample(self, sample_dir: Path) -> Optional[str]:
        """Locate a reference genome for SNP analysis, preferring CLI-specified value."""

        # 1. Prefer reference supplied on the command line / workflow runner
        try:
            runner = getattr(self, "workflow_runner", None)
            runner_reference = None

            if runner is not None:
                runner_reference = getattr(runner, "reference", None)

                # Some wrappers may store the reference on an inner object
                if not runner_reference and hasattr(runner, "manager"):
                    runner_reference = getattr(runner.manager, "reference", None)

            if runner_reference:
                ref_path = Path(str(runner_reference))
                if ref_path.exists():
                    logging.info(
                        "Using reference genome specified via command line: %s",
                        ref_path,
                    )
                    return str(ref_path)
                logging.warning(
                    "Reference genome provided via command line does not exist: %s",
                    ref_path,
                )
        except Exception as exc:
            logging.debug(f"Error checking workflow runner for reference genome: {exc}")

        # 2. Fallback to legacy heuristics (sample directory, monitored directory, env var)
        candidates: List[Path] = []

        try:
            candidates.extend([
                sample_dir / "reference.fasta",
                sample_dir / "reference.fa",
            ])

            base_dir = Path(self.monitored_directory) if self.monitored_directory else sample_dir.parent
            if base_dir:
                candidates.extend([
                    base_dir / "reference.fasta",
                    base_dir / "reference.fa",
                ])

            env_reference = os.environ.get("robin_REFERENCE")
            if env_reference:
                candidates.append(Path(env_reference))
        except Exception as exc:
            logging.debug(f"Error assembling fallback reference candidates for {sample_dir}: {exc}")

        for candidate in candidates:
            try:
                if candidate and candidate.exists():
                    return str(candidate)
            except Exception:
                continue

        return None

    def _trigger_target_bam_finalization(self, sample_id: str, *, trigger_snp: bool = False) -> None:
        """Trigger target.bam finalization for a sample. Optionally start SNP analysis afterwards."""
        if not self.monitored_directory:
            return
        
        try:
            from robin.analysis.target_analysis import finalize_accumulation_for_sample
            import csv
            
            # Read target_panel from master.csv
            sample_dir = Path(self.monitored_directory) / sample_id
            master_csv = sample_dir / "master.csv"
            target_panel = "rCNS2"  # Default
            
            if master_csv.exists():
                try:
                    with master_csv.open("r", newline="") as fh:
                        reader = csv.DictReader(fh)
                        first_row = next(reader, None)
                        if first_row:
                            panel = first_row.get("analysis_panel", "").strip()
                            if panel:
                                target_panel = panel
                except Exception:
                    pass  # Use default
            
            already_finalized = sample_id in self._finalized_samples

            # Trigger finalization (and optional SNP analysis) in background thread to avoid blocking GUI
            import threading
            def finalize_in_background():
                try:
                    finalization_succeeded = True

                    if already_finalized:
                        logging.info(f"Sample {sample_id} already finalized; skipping target.bam merge")
                    else:
                        logging.info(
                            f"Triggering target.bam finalization for sample {sample_id} (status: Complete)"
                        )
                        result = finalize_accumulation_for_sample(
                            sample_id=sample_id,
                            work_dir=str(self.monitored_directory),
                            target_panel=target_panel
                        )
                        if result.get("status") != "error":
                            self._finalized_samples.add(sample_id)
                            logging.info(
                                f"Successfully finalized target.bam for sample {sample_id}"
                            )
                        else:
                            finalization_succeeded = False
                            logging.warning(
                                f"Target.bam finalization returned error for {sample_id}: {result.get('error', 'Unknown')}"
                            )

                    if trigger_snp and finalization_succeeded:
                        # Ensure only one SNP analysis submission runs at a time
                        with self._snp_analysis_lock:
                            running_sample = self._snp_analysis_running_sample
                            if running_sample:
                                if running_sample == sample_id:
                                    message = "SNP analysis is already running for this sample."
                                else:
                                    message = (
                                        f"SNP analysis already running for sample {running_sample}. "
                                        "Please wait for it to finish before starting another."
                                    )
                                logging.info(message)
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis in progress",
                                        "message": message,
                                        "sample_id": sample_id,
                                        "level": "warning",
                                    },
                                    priority=6,
                                )
                                return

                            self._snp_analysis_running_sample = sample_id

                        try:
                            reference_path = self._locate_reference_for_sample(sample_dir)
                            if not reference_path:
                                message = (
                                    "Reference genome not available; skipping SNP analysis. "
                                    "Please provide a reference via the workflow CLI."
                                )
                                logging.warning(message)
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis skipped",
                                        "message": message,
                                        "sample_id": sample_id,
                                        "level": "warning",
                                    },
                                    priority=6,
                                )
                                return

                            target_bam = sample_dir / "target.bam"
                            if not target_bam.exists():
                                message = "target.bam not found; cannot run SNP analysis."
                                logging.warning(message)
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis skipped",
                                        "message": message,
                                        "sample_id": sample_id,
                                        "level": "warning",
                                    },
                                    priority=6,
                                )
                                return

                            workflow_runner = getattr(self, "workflow_runner", None)
                            if workflow_runner and hasattr(
                                workflow_runner, "submit_snp_analysis_job"
                            ):
                                submitted = workflow_runner.submit_snp_analysis_job(
                                    sample_dir=str(sample_dir),
                                    sample_id=sample_id,
                                    reference=reference_path,
                                    threads=4,
                                    force_regenerate=False,
                                )
                            elif workflow_runner and hasattr(
                                workflow_runner, "submit_sample_job"
                            ):
                                submitted = workflow_runner.submit_sample_job(
                                    sample_dir=str(sample_dir),
                                    job_type="snp_analysis",
                                    sample_id=sample_id,
                                )
                            else:
                                submitted = False

                            if submitted:
                                logging.info(
                                    f"SNP analysis job submitted for {sample_id} using reference {reference_path}"
                                )
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis queued",
                                        "message": "SNP analysis job submitted to the workflow queue.",
                                        "sample_id": sample_id,
                                        "level": "info",
                                    },
                                    priority=6,
                                )
                            else:
                                logging.warning(
                                    "No workflow runner available; SNP analysis not submitted."
                                )
                                self.send_update(
                                    UpdateType.WARNING_NOTIFICATION,
                                    {
                                        "title": "SNP analysis skipped",
                                        "message": "Workflow runner unavailable; SNP analysis not submitted.",
                                        "sample_id": sample_id,
                                        "level": "warning",
                                    },
                                    priority=6,
                                )
                        except Exception as snp_exc:
                            logging.error(
                                f"Error submitting SNP analysis for {sample_id}: {snp_exc}",
                                exc_info=True,
                            )
                            self.send_update(
                                UpdateType.WARNING_NOTIFICATION,
                                {
                                    "title": "SNP analysis failed",
                                    "message": f"Error submitting SNP analysis: {snp_exc}",
                                    "sample_id": sample_id,
                                    "level": "negative",
                                },
                                priority=6,
                            )
                        finally:
                            with self._snp_analysis_lock:
                                self._snp_analysis_running_sample = None
                except Exception as e:
                    logging.error(f"Error finalizing target.bam for {sample_id}: {e}")
            
            thread = threading.Thread(target=finalize_in_background, daemon=True)
            thread.start()
            
        except Exception as e:
            logging.error(f"Failed to trigger target.bam finalization for {sample_id}: {e}")
    
    def _scan_and_seed_samples(self, preexisting: bool = False) -> None:
        """Synchronous version - kept for backward compatibility and io_bound calls"""
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
                    ov_pending = 0
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
                                saved_last = float(
                                    first_row.get("samples_overview_last_seen", 0.0)
                                )
                                if saved_last:
                                    last_seen = saved_last
                            except Exception:
                                pass
                            ov_active = int(
                                first_row.get("samples_overview_active_jobs", 0) or 0
                            )
                            ov_pending = int(
                                first_row.get("samples_overview_pending_jobs", 0) or 0
                            )
                            ov_total = int(
                                first_row.get("samples_overview_total_jobs", 0) or 0
                            )
                            ov_completed = int(
                                first_row.get("samples_overview_completed_jobs", 0) or 0
                            )
                            ov_failed = int(
                                first_row.get("samples_overview_failed_jobs", 0) or 0
                            )
                            ov_job_types = str(
                                first_row.get("samples_overview_job_types", "") or ""
                            )
                    except Exception:
                        pass
                    origin_value = (
                        "Pre-existing"
                        if sid in self._preexisting_sample_ids and (time.time() - last_seen) >= self.completion_timeout_seconds
                        else "Live"
                    )
                    try:
                        # Only mark as Complete if timeout passed AND no active jobs
                        if origin_value == "Live" and (time.time() - last_seen) >= self.completion_timeout_seconds:
                            if ov_active == 0 and ov_pending == 0:
                                origin_value = "Complete"
                                # Check if this is a transition from Live to Complete
                                existing = self._last_samples_rows or []
                                prev_origin = None
                                for prev_row in existing:
                                    if prev_row.get("sample_id") == sid:
                                        prev_origin = prev_row.get("origin")
                                        break
                                # Trigger finalization if transitioning from Live to Complete
                                if prev_origin == "Live" or prev_origin is None:
                                    self._trigger_target_bam_finalization(sid)
                            # If there are active jobs, keep as Live even if timeout passed
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
                            "pending_jobs": ov_pending,
                            "total_jobs": ov_total,
                            "completed_jobs": ov_completed,
                            "failed_jobs": ov_failed,
                            "job_types": ov_job_types,
                            "last_seen": time.strftime(
                                "%Y-%m-%d %H:%M:%S", time.localtime(last_seen)
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
                # Sort rows by last activity in reverse chronological order (newest first)
                try:
                    merged.sort(key=lambda r: float(r.get("_last_seen_raw", 0)), reverse=True)
                except Exception:
                    # Fallback to sorting by last_seen string if _last_seen_raw is not available
                    try:
                        merged.sort(key=lambda r: r.get("last_seen", ""), reverse=True)
                    except Exception:
                        pass
                
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

    def _get_cached_samples(self) -> Optional[List[Dict[str, Any]]]:
        """Return cached sample data if available and recent"""
        if self._last_samples_rows and time.time() - self._last_cache_time < self._cache_duration:
            return self._last_samples_rows
        return None

    def _get_expected_completion_job_types(self) -> Set[str]:
        """Get the set of job types expected to complete for this workflow."""
        if not self.workflow_steps:
            return set()
        expected: Set[str] = set()
        for step in self.workflow_steps:
            step_name = step.split(":")[-1] if ":" in step else step
            if step_name in COMPLETION_JOB_PATTERNS:
                expected.add(step_name)
        return expected

    def _expected_jobs_completed(self, sample_dir: Path, expected_job_types: Set[str]) -> bool:
        """Return True if all expected job types have completion outputs."""
        if not expected_job_types:
            return True
        if not sample_dir.exists():
            return False
        counts = self._calculate_job_counts_from_files(
            sample_dir, expected_job_types=expected_job_types
        )
        total_jobs = int(counts.get("total_jobs", 0))
        completed_jobs = int(counts.get("completed_jobs", 0))
        return total_jobs > 0 and completed_jobs >= total_jobs

    def _calculate_job_counts_from_files(
        self, sample_dir: Path, expected_job_types: Optional[Set[str]] = None
    ) -> Dict[str, Any]:
        """Calculate job counts and types from actual analysis result files in the sample directory."""
        try:
            total_jobs = 0
            completed_jobs = 0
            failed_jobs = 0
            job_types = set()
            
            job_patterns = COMPLETION_JOB_PATTERNS
            job_types_to_check = (
                expected_job_types if expected_job_types else set(job_patterns.keys())
            )
            
            # Check for each job type
            for job_type in job_types_to_check:
                patterns = job_patterns.get(job_type, [])
                job_found = False
                for pattern in patterns:
                    if "*" in pattern:
                        # Use glob pattern
                        matching_files = list(sample_dir.glob(pattern))
                        if matching_files:
                            job_found = True
                            break
                    else:
                        # Check for exact file
                        if (sample_dir / pattern).exists():
                            job_found = True
                            break
                
                if job_found:
                    total_jobs += 1
                    completed_jobs += 1  # Assume completed if files exist
                    job_types.add(job_type)

            if expected_job_types:
                total_jobs = len(expected_job_types)
            
            return {
                "total_jobs": total_jobs,
                "completed_jobs": completed_jobs,
                "failed_jobs": failed_jobs,
                "job_types": ",".join(sorted(job_types)) if job_types else ""
            }
            
        except Exception as e:
            logging.warning(f"Error calculating job counts from files for {sample_dir}: {e}")
            return {
                "total_jobs": 0,
                "completed_jobs": 0,
                "failed_jobs": 0,
                "job_types": ""
            }

    async def _load_samples_progressively(self, sample_dirs: List[Path], batch_size: int = 10) -> None:
        """Load samples in small batches to avoid blocking the UI"""
        try:
            total_dirs = len(sample_dirs)
            processed = 0
            
            # Show progress
            if hasattr(self, "samples_loading_container"):
                self.samples_loading_container.clear()
                with self.samples_loading_container:
                    ui.spinner(size="lg", color="primary")
                    ui.label("Loading samples...").classes("ml-2 text-lg")
                    progress_label = ui.label(f"Processing {processed}/{total_dirs} samples").classes("text-sm text-gray-500 mt-2")
            
            rows: List[Dict[str, Any]] = []
            
            for i in range(0, total_dirs, batch_size):
                batch = sample_dirs[i:i + batch_size]
                
                # Process batch asynchronously in background thread
                import asyncio
                batch_rows = await asyncio.to_thread(self._process_sample_batch, batch)
                rows.extend(batch_rows)
                
                processed += len(batch)
                
                # Update progress
                if hasattr(self, "samples_loading_container") and 'progress_label' in locals():
                    progress_label.set_text(f"Processing {processed}/{total_dirs} samples")
                
                # Small delay to keep UI responsive
                await asyncio.sleep(0.1)
            
            # Update table with all rows
            if rows:
                self._last_samples_rows = rows
                self._last_cache_time = time.time()
                
                if hasattr(self, "samples_table"):
                    try:
                        self.samples_table.rows = rows
                        self.samples_table.update()
                    except Exception:
                        pass
            
            # Hide loading container
            if hasattr(self, "samples_loading_container"):
                self.samples_loading_container.set_visibility(False)
                
        except Exception as e:
            logging.error(f"Error in progressive loading: {e}")
            if hasattr(self, "samples_loading_container"):
                self.samples_loading_container.set_visibility(False)

    def _process_sample_batch(self, sample_dirs: List[Path]) -> List[Dict[str, Any]]:
        """Process a batch of sample directories synchronously"""
        rows: List[Dict[str, Any]] = []
        expected_job_types = self._get_expected_completion_job_types()
        
        for sample_dir in sample_dirs:
            if not sample_dir.is_dir():
                continue
                
            master = sample_dir / "master.csv"
            if master.exists():
                sid = sample_dir.name
                
                # Determine last_seen from file mtime
                last_seen = master.stat().st_mtime
                
                # Try to load persisted overview and run info
                run_start = ""
                device = ""
                flowcell = ""
                ov_active = 0
                ov_pending = 0
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
                            saved_last = float(
                                first_row.get("samples_overview_last_seen", 0.0)
                            )
                            if saved_last:
                                last_seen = saved_last
                        except Exception:
                            pass
                            
                        # Try to get overview data from master.csv first
                        ov_active = int(
                            first_row.get("samples_overview_active_jobs", 0) or 0
                        )
                        ov_pending = int(
                            first_row.get("samples_overview_pending_jobs", 0) or 0
                        )
                        ov_total = int(
                            first_row.get("samples_overview_total_jobs", 0) or 0
                        )
                        ov_completed = int(
                            first_row.get("samples_overview_completed_jobs", 0) or 0
                        )
                        ov_failed = int(
                            first_row.get("samples_overview_failed_jobs", 0) or 0
                        )
                        ov_job_types = str(
                            first_row.get("samples_overview_job_types", "") or ""
                        )
                        
                        # If overview data is not available (0 values), calculate from actual analysis files
                        if ov_total == 0 and ov_completed == 0:
                            calculated_counts = self._calculate_job_counts_from_files(
                                sample_dir,
                                expected_job_types=expected_job_types or None,
                            )
                            ov_total = calculated_counts["total_jobs"]
                            ov_completed = calculated_counts["completed_jobs"]
                            ov_failed = calculated_counts["failed_jobs"]
                            ov_job_types = calculated_counts["job_types"]
                except Exception:
                    pass
                    
                origin_value = (
                    "Pre-existing"
                    if sid in self._preexisting_sample_ids and (time.time() - last_seen) >= self.completion_timeout_seconds
                    else "Live"
                )
                try:
                    # Only mark as Complete if timeout passed AND no active jobs
                    if origin_value == "Live" and (time.time() - last_seen) >= self.completion_timeout_seconds:
                        if ov_active == 0:
                            should_complete = True
                            if expected_job_types:
                                should_complete = self._expected_jobs_completed(
                                    sample_dir, expected_job_types
                                )
                            if should_complete:
                                origin_value = "Complete"
                        # If there are active jobs, keep as Live even if timeout passed
                except Exception:
                    pass
                    
                # File progress - use the same job counts that drive other columns
                # This ensures consistency with the Active, Total, Completed columns
                files_seen = ov_total  # Use total_jobs as files_seen
                files_processed = ov_completed  # Use completed_jobs as files_processed
                file_progress = files_processed / files_seen if files_seen > 0 else 0.0
                
                row_data = {
                    "sample_id": sid,
                    "origin": origin_value,
                    "run_start": self._format_timestamp_for_display(run_start),
                    "device": device,
                    "flowcell": flowcell,
                    "active_jobs": ov_active,
                    "pending_jobs": ov_pending,
                    "total_jobs": ov_total,
                    "completed_jobs": ov_completed,
                    "failed_jobs": ov_failed,
                    "job_types": ov_job_types,
                    "files_seen": files_seen,
                    "files_processed": files_processed,
                    "file_progress": file_progress,
                    "last_seen": time.strftime(
                        "%Y-%m-%d %H:%M:%S", time.localtime(last_seen)
                    ),
                    "_last_seen_raw": last_seen,
                }
                
                rows.append(row_data)
        
        return rows

    async def _scan_and_seed_samples_async(self, preexisting: bool = False) -> None:
        """Asynchronous version that runs file operations in background thread"""
        try:
            # Check if we have recent cached data
            cached_data = self._get_cached_samples()
            if cached_data and not preexisting:
                # Use cached data and update UI immediately
                if hasattr(self, "samples_table"):
                    try:
                        self.samples_table.rows = cached_data
                        self.samples_table.update()
                        # Hide loading indicator
                        if hasattr(self, "samples_loading_container"):
                            self.samples_loading_container.set_visibility(False)
                    except Exception:
                        pass
                return
            
            # Check if we need progressive loading for large directories
            base = Path(self.monitored_directory) if self.monitored_directory else None
            if base and base.exists():
                sample_dirs = [d for d in base.iterdir() if d.is_dir()]
                
                # Use progressive loading for directories with many samples
                if len(sample_dirs) > 20:
                    await self._load_samples_progressively(sample_dirs)
                    return
            
            # Show loading indicator
            if hasattr(self, "samples_loading_indicator"):
                self.samples_loading_indicator.set_visibility(True)
            
            # Run the synchronous file operations in a background thread using asyncio.to_thread
            # This prevents blocking the main thread and GUI updates
            import asyncio
            result = await asyncio.to_thread(self._scan_and_seed_samples, preexisting)
            
            # Update cache timestamp
            self._last_cache_time = time.time()
            
            # Hide loading indicator
            if hasattr(self, "samples_loading_indicator"):
                self.samples_loading_indicator.set_visibility(False)
            if hasattr(self, "samples_loading_container"):
                self.samples_loading_container.set_visibility(False)
            
        except Exception as e:
            logging.error(f"Error in async sample scanning: {e}")
            
            # Hide loading indicators on error
            if hasattr(self, "samples_loading_indicator"):
                self.samples_loading_indicator.set_visibility(False)
            if hasattr(self, "samples_loading_container"):
                self.samples_loading_container.set_visibility(False)

    async def _scan_for_new_samples_async(self) -> None:
        """Asynchronous version of new samples scanning"""
        try:
            # Always scan for new samples to ensure we get updates
            # The cache will be updated by _scan_for_new_samples if there are changes
            self._scan_for_new_samples()
        except Exception as e:
            logging.error(f"Error in async new samples scanning: {e}")

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
                        ov_pending = 0
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
                                saved_last = float(
                                    first_row.get("samples_overview_last_seen", 0.0)
                                )
                                if saved_last:
                                    last_seen = saved_last
                            except Exception:
                                pass
                            ov_active = int(
                                first_row.get("samples_overview_active_jobs", 0) or 0
                            )
                            ov_pending = int(
                                first_row.get("samples_overview_pending_jobs", 0) or 0
                            )
                            ov_total = int(
                                first_row.get("samples_overview_total_jobs", 0) or 0
                            )
                            ov_completed = int(
                                first_row.get("samples_overview_completed_jobs", 0) or 0
                            )
                            ov_failed = int(
                                first_row.get("samples_overview_failed_jobs", 0) or 0
                            )
                            ov_job_types = str(
                                first_row.get("samples_overview_job_types", "") or ""
                            )
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
                                self._format_timestamp_for_display(run_start)
                                if run_start
                                else ""
                            )
                            updated["run_start"] = (
                                formatted_run_start or existing_row.get("run_start", "")
                            )
                            updated["device"] = device or existing_row.get("device", "")
                            updated["flowcell"] = flowcell or existing_row.get(
                                "flowcell", ""
                            )
                            if ov_job_types:
                                updated["job_types"] = ov_job_types
                            updated["active_jobs"] = ov_active
                            updated["pending_jobs"] = ov_pending
                            updated["total_jobs"] = ov_total
                            updated["completed_jobs"] = ov_completed
                            updated["failed_jobs"] = ov_failed
                            updated["last_seen"] = time.strftime(
                                "%H:%M:%S", time.localtime(last_seen)
                            )
                            updated["_last_seen_raw"] = last_seen
                            # Determine origin based on inactivity threshold AND active jobs status
                            prev_origin = existing_row.get("origin")
                            try:
                                # Only mark as Complete if timeout passed AND no active jobs
                                if (time.time() - last_seen) >= self.completion_timeout_seconds:
                                    if ov_active == 0 and ov_pending == 0:
                                        updated["origin"] = "Complete"
                                        # Trigger finalization if transitioning from Live to Complete
                                        if prev_origin == "Live":
                                            self._trigger_target_bam_finalization(sid)
                                    else:
                                        # If there are active jobs, keep as Live even if timeout passed
                                        updated["origin"] = "Live"
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
                # Update cache timestamp when we update the samples data
                self._last_cache_time = time.time()
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
            if not hasattr(self, "workflow_runner") or self.workflow_runner is None:
                print(f"[GUI] No workflow runner available for {job_type} job")
                return False

            if sample_id is None:
                sample_id = Path(sample_dir).name

            # Check if it's a Simple workflow or Ray workflow
            if hasattr(self.workflow_runner, "submit_sample_job"):
                # Simple workflow
                return self.workflow_runner.submit_sample_job(
                    sample_dir, job_type, sample_id
                )
            elif hasattr(self.workflow_runner, "manager") and hasattr(
                self.workflow_runner.manager, "submit_sample_job"
            ):
                # Ray workflow - this is more complex, so we'll provide a fallback
                print(
                    f"[GUI] Ray workflow detected for {job_type} job - manual submission not yet supported"
                )
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
    center: str = None,
) -> GUILauncher:
    """Legacy launch function (kept for backward compatibility).

    Internally delegates to the refactored launcher in `robin.gui.app`.
    """
    from .gui.app import launch_gui as _launch  # type: ignore

    return _launch(
        host=host,
        port=port,
        show=show,
        reload=False,
        workflow_runner=workflow_runner,
        workflow_steps=workflow_steps,
        monitored_directory=monitored_directory,
        center=center,
    )


def get_gui_launcher() -> Optional[GUILauncher]:
    """Compatibility shim that delegates to robin.gui.app.get_gui_launcher."""
    try:
        from .gui.app import get_gui_launcher as _get  # type: ignore

        return _get()
    except Exception:
        return None


def send_gui_update(update_type: UpdateType, data: Dict[str, Any], priority: int = 0):
    """Compatibility shim that delegates to robin.gui.app.send_gui_update."""
    try:
        from .gui.app import send_gui_update as _send  # type: ignore

        _send(update_type, data, priority)
    except Exception:
        logging.info("[GUI] No GUI launcher available for update (update dropped)")


# Global reference retained for backward compatibility (now managed in gui.app)
_current_gui_launcher = None


if __name__ in {"__main__", "__mp_main__"}:
    """
    Command-line interface for direct GUI launching.

    Usage:
        # Launch GUI with default settings
        python gui_launcher.py

        # Launch GUI pointing to specific directory
        python gui_launcher.py test_out_priority_ray4

        # Launch GUI on specific port
        python gui_launcher.py test_out_priority_ray4 --port 8081

        # Launch GUI without opening browser
        python gui_launcher.py test_out_priority_ray4 --no-show

        # Launch GUI on different host
        python gui_launcher.py test_out_priority_ray4 --host 0.0.0.0

        # Launch GUI with auto-reload for development
        python gui_launcher.py test_out_priority_ray4 --reload
    """
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Direct GUI launcher for robin workflow monitoring",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Examples:
        python gui_launcher.py                                    # Launch with defaults
        python gui_launcher.py test_out_priority_ray4            # Monitor specific directory
        python gui_launcher.py test_out_priority_ray4 --port 8081 # Use different port
        python gui_launcher.py test_out_priority_ray4 --no-show   # Don't open browser
        python gui_launcher.py test_out_priority_ray4 --reload    # Enable auto-reload for development
                """,
    )

    parser.add_argument(
        "monitored_directory",
        nargs="?",
        default="",
        help="Directory containing sample output folders to monitor",
    )

    parser.add_argument(
        "--host",
        default="localhost",
        help="Host to bind the GUI to (default: localhost)",
    )

    parser.add_argument(
        "--port", type=int, default=8080, help="Port to run the GUI on (default: 8080)"
    )

    parser.add_argument(
        "--no-show", action="store_true", help="Don't automatically open the browser"
    )

    parser.add_argument(
        "--reload",
        action="store_true",
        help="Enable auto-reload for development (watches for file changes)",
    )

    args = parser.parse_args()

    # Validate monitored directory if provided
    if args.monitored_directory:
        monitored_path = Path(args.monitored_directory)
        if not monitored_path.exists():
            print(f"Error: Directory '{args.monitored_directory}' does not exist")
            sys.exit(1)
        if not monitored_path.is_dir():
            print(f"Error: '{args.monitored_directory}' is not a directory")
            sys.exit(1)
        print(f"Will monitor: {monitored_path.resolve()}")

    try:
        # Create and launch GUI directly to avoid relative import issues
        print("Creating GUI launcher...")

        # Create the GUI launcher directly
        launcher = GUILauncher(
            host=args.host,
            port=args.port,
            reload=args.reload,
        )

        # Set monitored directory if provided
        if args.monitored_directory:
            launcher.monitored_directory = str(Path(args.monitored_directory).resolve())

        print(f"Launching full GUI on http://{args.host}:{args.port}")
        print("Use Ctrl+C to stop the GUI")

        # Launch the GUI directly
        success = launcher.launch_gui(
            monitored_directory=args.monitored_directory,
            workflow_runner=None,
            workflow_steps=[],
        )

        if success:
            print("GUI launched successfully!")
            print("🌐 Open your browser to the URL above")
            print("Press Ctrl+C to stop the GUI")

            # Keep the main thread alive while GUI runs
            try:
                while launcher.is_gui_running():
                    time.sleep(1)
            except KeyboardInterrupt:
                print("\n🛑 Shutting down GUI...")
                launcher.stop_gui()
                print("GUI stopped")
        else:
            print("Failed to launch GUI")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\nGUI stopped by user")
    except Exception as e:
        print(f"GUI failed to start: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
