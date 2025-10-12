#!/usr/bin/env python3
"""
Master CSV Manager for robin

This module provides functionality to create and update master.csv files
for each sample, tracking comprehensive metadata across multiple BAM files.
"""

import os
import pandas as pd
from typing import Dict, Any
from dataclasses import dataclass
import logging
import time

# Import file locking (Unix/macOS)
try:
    import fcntl
    HAS_FCNTL = True
except ImportError:
    HAS_FCNTL = False


@dataclass
class MasterCSVManager:
    """Manages the creation and updating of master.csv files for each sample"""

    work_dir: str

    def _acquire_lock(self, file_handle, exclusive=True, timeout=10.0):
        """
        Acquire a file lock with timeout.
        
        Args:
            file_handle: Open file handle to lock
            exclusive: If True, acquire exclusive lock (for writing). If False, shared lock (for reading)
            timeout: Maximum time to wait for lock in seconds
        
        Returns:
            True if lock acquired, False otherwise
        """
        if not HAS_FCNTL:
            return True  # No locking available, proceed anyway
        
        lock_type = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
        start_time = time.time()
        
        while True:
            try:
                fcntl.flock(file_handle.fileno(), lock_type | fcntl.LOCK_NB)
                return True
            except (IOError, OSError):
                # Lock is held by another process
                if time.time() - start_time >= timeout:
                    return False
                time.sleep(0.01)  # Wait 10ms before retry
    
    def _release_lock(self, file_handle):
        """Release a file lock."""
        if not HAS_FCNTL:
            return
        
        try:
            fcntl.flock(file_handle.fileno(), fcntl.LOCK_UN)
        except (IOError, OSError):
            pass  # Lock already released or file closed

    def update_master_csv(
        self, sample_id: str, bam_data: Dict[str, Any], bam_info: Dict[str, Any]
    ) -> None:
        """Update or create master.csv file for a sample with new BAM data"""
        try:
            # Create sample directory if it doesn't exist
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)

            master_csv_path = os.path.join(sample_dir, "master.csv")

            # Load existing data or create new
            existing_data = self._load_existing_data(master_csv_path)

            # Update counters based on BAM state
            if bam_info.get("state") == "pass":
                existing_data = self._update_pass_counters(existing_data, bam_data)
            else:
                existing_data = self._update_fail_counters(existing_data, bam_data)

            # Update combined counters
            existing_data = self._update_combined_counters(existing_data, bam_data)

            # Update run information
            existing_data = self._update_run_info(existing_data, bam_info)

            # Update BAM tracking
            existing_data = self._update_bam_tracking(existing_data, bam_info)

            # Save to CSV
            self._save_to_csv(existing_data, master_csv_path)

            logger = logging.getLogger("robin.master_csv")
            logger.info(f"Updated master.csv for sample {sample_id}")

        except Exception as e:
            logger = logging.getLogger("robin.master_csv")
            logger.error(f"Error updating master.csv for {sample_id}: {e}")

    def _load_existing_data(self, csv_path: str) -> Dict[str, Any]:
        """
        Load existing data from CSV or create default structure.
        
        No read locking is used because:
        1. Writers use atomic write-and-rename, so readers never see partial files
        2. Worst case is reading slightly stale data (acceptable for monitoring)
        3. Keeps GUI responsive and non-blocking
        """
        if os.path.exists(csv_path):
            try:
                # Read without locking - atomic writes ensure we never see partial data
                df = pd.read_csv(csv_path)
                if not df.empty:
                    data = df.iloc[0].to_dict()
                    # Ensure string fields are properly converted to strings
                    # to prevent float objects from being passed to split() methods
                    string_fields = [
                        "devices", "basecall_models", "run_time", "flowcell_ids",
                        "run_info_run_time", "run_info_device", "run_info_model", 
                        "run_info_flow_cell", "samples_overview_job_types", "analysis_panel"
                    ]
                    for field in string_fields:
                        if field in data and data[field] is not None:
                            data[field] = str(data[field])
                    return data
            except Exception as e:
                print(f"Warning: Error reading existing CSV: {e}")

        # Return default structure
        return self._get_default_structure()
    
    def _get_default_structure(self) -> Dict[str, Any]:
        return {
            "counter_bam_passed": 0,
            "counter_bam_failed": 0,
            "counter_mapped_count": 0,
            "counter_pass_mapped_count": 0,
            "counter_fail_mapped_count": 0,
            "counter_unmapped_count": 0,
            "counter_pass_unmapped_count": 0,
            "counter_fail_unmapped_count": 0,
            "counter_pass_bases_count": 0,
            "counter_fail_bases_count": 0,
            "counter_bases_count": 0,
            "counter_mapped_reads_num": 0,
            "counter_unmapped_reads_num": 0,
            "counter_pass_mapped_reads_num": 0,
            "counter_fail_mapped_reads_num": 0,
            "counter_pass_unmapped_reads_num": 0,
            "counter_fail_unmapped_reads_num": 0,
            "counter_mapped_bases": 0,
            "counter_unmapped_bases": 0,
            "counter_pass_mapped_bases": 0,
            "counter_fail_mapped_bases": 0,
            "counter_pass_unmapped_bases": 0,
            "counter_fail_unmapped_bases": 0,
            "devices": "",
            "basecall_models": "",
            "run_time": "",
            "flowcell_ids": "",
            "run_info_run_time": "",
            "run_info_device": "",
            "run_info_model": "",
            "run_info_flow_cell": "",
            "bam_tracking_counter": 0,
            "bam_tracking_total_files": 0,
            # Analysis panel information
            "analysis_panel": "",
            # Samples overview (persisted GUI table aggregates)
            "samples_overview_active_jobs": 0,
            "samples_overview_total_jobs": 0,
            "samples_overview_completed_jobs": 0,
            "samples_overview_failed_jobs": 0,
            # Comma-separated unique job types seen for the sample
            "samples_overview_job_types": "",
            # Unix epoch seconds for last activity used by GUI
            "samples_overview_last_seen": 0.0,
        }

    def _update_pass_counters(
        self, data: Dict[str, Any], bam_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Update counters for pass BAM files"""
        data["counter_bam_passed"] += 1
        data["counter_pass_mapped_count"] += bam_data.get("mapped_reads", 0)
        data["counter_pass_unmapped_count"] += bam_data.get("unmapped_reads", 0)
        data["counter_pass_bases_count"] += bam_data.get("yield_tracking", 0)
        data["counter_pass_mapped_reads_num"] += bam_data.get("mapped_reads_num", 0)
        data["counter_pass_unmapped_reads_num"] += bam_data.get("unmapped_reads_num", 0)
        data["counter_pass_mapped_bases"] += bam_data.get("mapped_bases", 0)
        data["counter_pass_unmapped_bases"] += bam_data.get("unmapped_bases", 0)
        return data

    def _update_fail_counters(
        self, data: Dict[str, Any], bam_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Update counters for fail BAM files"""
        data["counter_bam_failed"] += 1
        data["counter_fail_mapped_count"] += bam_data.get("mapped_reads", 0)
        data["counter_fail_unmapped_count"] += bam_data.get("unmapped_reads", 0)
        data["counter_fail_bases_count"] += bam_data.get("yield_tracking", 0)
        data["counter_fail_mapped_reads_num"] += bam_data.get("mapped_reads_num", 0)
        data["counter_fail_unmapped_reads_num"] += bam_data.get("unmapped_reads_num", 0)
        data["counter_fail_mapped_bases"] += bam_data.get("mapped_bases", 0)
        data["counter_fail_unmapped_bases"] += bam_data.get("unmapped_bases", 0)
        return data

    def _update_combined_counters(
        self, data: Dict[str, Any], bam_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Update combined counters for all BAM files"""
        data["counter_mapped_count"] += bam_data.get("mapped_reads", 0)
        data["counter_unmapped_count"] += bam_data.get("unmapped_reads", 0)
        data["counter_bases_count"] += bam_data.get("yield_tracking", 0)
        data["counter_mapped_reads_num"] += bam_data.get("mapped_reads_num", 0)
        data["counter_unmapped_reads_num"] += bam_data.get("unmapped_reads_num", 0)
        data["counter_mapped_bases"] += bam_data.get("mapped_bases", 0)
        data["counter_unmapped_bases"] += bam_data.get("unmapped_bases", 0)
        return data

    def _update_run_info(
        self, data: Dict[str, Any], bam_info: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Update run information arrays"""
        # Update devices
        if bam_info.get("device_position"):
            device_str = str(bam_info["device_position"])
            devices_str = str(data.get("devices", ""))
            devices = devices_str.split(",") if devices_str else []
            if device_str not in devices:
                devices.append(device_str)
            data["devices"] = ",".join(devices)
            data["run_info_device"] = device_str

        # Update basecall models
        if bam_info.get("basecall_model"):
            model_str = str(bam_info["basecall_model"])
            models_str = str(data.get("basecall_models", ""))
            models = models_str.split(",") if models_str else []
            if model_str not in models:
                models.append(model_str)
            data["basecall_models"] = ",".join(models)
            data["run_info_model"] = model_str

        # Update flow cell IDs
        if bam_info.get("flow_cell_id"):
            flowcell_str = str(bam_info["flow_cell_id"])
            flowcells_str = str(data.get("flowcell_ids", ""))
            flowcells = flowcells_str.split(",") if flowcells_str else []
            if flowcell_str not in flowcells:
                flowcells.append(flowcell_str)
            data["flowcell_ids"] = ",".join(flowcells)
            data["run_info_flow_cell"] = flowcell_str

        # Update run time
        if bam_info.get("time_of_run"):
            # Convert to string if it's a float or other type
            time_str = str(bam_info["time_of_run"])
            run_times_str = str(data.get("run_time", ""))
            run_times = run_times_str.split(",") if run_times_str else []
            if time_str not in run_times:
                run_times.append(time_str)
            data["run_time"] = ",".join(run_times)
            data["run_info_run_time"] = time_str

        return data

    def _update_bam_tracking(
        self, data: Dict[str, Any], bam_info: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Update BAM tracking information"""
        data["bam_tracking_counter"] += 1
        data["bam_tracking_total_files"] += 1

        return data

    def update_analysis_panel(self, sample_id: str, panel: str) -> None:
        """
        Update the analysis panel information for a sample.
        
        Args:
            sample_id: The sample ID
            panel: The analysis panel used (rCNS2, AML, PanCan)
        """
        try:
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            master_csv_path = os.path.join(sample_dir, "master.csv")

            existing_data = self._load_existing_data(master_csv_path)
            
            # Update the analysis panel
            existing_data["analysis_panel"] = str(panel)
            
            # Save to CSV
            self._save_to_csv(existing_data, master_csv_path)
            
            logger = logging.getLogger("robin.master_csv")
            logger.info(f"Updated analysis panel to '{panel}' for sample {sample_id}")
            
        except Exception as e:
            logger = logging.getLogger("robin.master_csv")
            logger.error(f"Error updating analysis panel for {sample_id}: {e}")

    def _save_to_csv(self, data: Dict[str, Any], csv_path: str) -> None:
        """
        Save data to CSV file with exclusive locking to prevent concurrent writes.
        Uses atomic write-and-rename pattern to prevent partial reads.
        """
        import tempfile
        
        # Create dataframe
        df = pd.DataFrame([data])
        
        # Write to temporary file in the same directory (for atomic rename)
        csv_dir = os.path.dirname(csv_path)
        temp_fd, temp_path = tempfile.mkstemp(dir=csv_dir, prefix='.master_', suffix='.csv.tmp')
        
        try:
            # Write to temp file with exclusive lock
            with os.fdopen(temp_fd, 'w') as f:
                # Acquire exclusive lock
                if not self._acquire_lock(f, exclusive=True, timeout=10.0):
                    raise IOError(f"Could not acquire exclusive lock for writing {csv_path}")
                
                try:
                    # Write CSV content
                    df.to_csv(f, index=False)
                    # Ensure data is written to disk
                    f.flush()
                    os.fsync(f.fileno())
                finally:
                    self._release_lock(f)
            
            # Atomically rename temp file to target (this is atomic on Unix/macOS)
            os.replace(temp_path, csv_path)
            
        except Exception as e:
            # Clean up temp file if something went wrong
            try:
                if os.path.exists(temp_path):
                    os.unlink(temp_path)
            except Exception:
                pass
            raise e

    # ---------------------------------------------------------------------
    # Public helpers for GUI/workflow to persist overview stats
    # ---------------------------------------------------------------------
    def update_sample_overview(self, sample_id: str, overview: Dict[str, Any]) -> None:
        """
        Persist high-level per-sample overview stats used by the GUI table.

        Expected keys in ``overview`` (all optional, defaults applied when missing):
          - active_jobs: int
          - total_jobs: int
          - completed_jobs: int
          - failed_jobs: int
          - job_types: Iterable[str] | str (stored as comma-separated unique list)
          - last_seen: float | int (unix epoch seconds)
        """
        try:
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            master_csv_path = os.path.join(sample_dir, "master.csv")

            existing_data = self._load_existing_data(master_csv_path)

            # Update numeric counters
            existing_data["samples_overview_active_jobs"] = int(
                overview.get(
                    "active_jobs", existing_data.get("samples_overview_active_jobs", 0)
                )
            )
            existing_data["samples_overview_total_jobs"] = int(
                overview.get(
                    "total_jobs", existing_data.get("samples_overview_total_jobs", 0)
                )
            )
            existing_data["samples_overview_completed_jobs"] = int(
                overview.get(
                    "completed_jobs",
                    existing_data.get("samples_overview_completed_jobs", 0),
                )
            )
            existing_data["samples_overview_failed_jobs"] = int(
                overview.get(
                    "failed_jobs", existing_data.get("samples_overview_failed_jobs", 0)
                )
            )

            # Job types: normalize to comma-separated sorted unique values
            jt_value = overview.get("job_types")
            if isinstance(jt_value, str):
                new_types = {t.strip() for t in jt_value.split(",") if t.strip()}
            elif jt_value is None:
                new_types = set()
            else:
                try:
                    new_types = {str(t).strip() for t in jt_value if str(t).strip()}
                except Exception:
                    new_types = set()

            if existing_data.get("samples_overview_job_types"):
                old_types = {
                    t.strip()
                    for t in str(
                        existing_data.get("samples_overview_job_types", "")
                    ).split(",")
                    if t.strip()
                }
            else:
                old_types = set()
            all_types = sorted(old_types.union(new_types))
            existing_data["samples_overview_job_types"] = ",".join(all_types)

            # Last seen timestamp
            try:
                last_seen = float(overview.get("last_seen"))
            except Exception:
                last_seen = float(existing_data.get("samples_overview_last_seen", 0.0))
            existing_data["samples_overview_last_seen"] = last_seen

            # Save
            self._save_to_csv(existing_data, master_csv_path)
            logger = logging.getLogger("robin.master_csv")
            logger.debug(
                f"Updated samples overview in master.csv for sample {sample_id}"
            )
        except Exception as e:
            logger = logging.getLogger("robin.master_csv")
            logger.warning(f"Failed to update samples overview for {sample_id}: {e}")
