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


@dataclass
class MasterCSVManager:
    """Manages the creation and updating of master.csv files for each sample"""

    work_dir: str

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
        """Load existing data from CSV or create default structure"""
        if os.path.exists(csv_path):
            try:
                df = pd.read_csv(csv_path)
                if not df.empty:
                    return df.iloc[0].to_dict()
            except Exception as e:
                print(f"Warning: Error reading existing CSV: {e}")

        # Return default structure
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
            devices = data["devices"].split(",") if data["devices"] else []
            if device_str not in devices:
                devices.append(device_str)
            data["devices"] = ",".join(devices)
            data["run_info_device"] = device_str

        # Update basecall models
        if bam_info.get("basecall_model"):
            model_str = str(bam_info["basecall_model"])
            models = (
                data["basecall_models"].split(",") if data["basecall_models"] else []
            )
            if model_str not in models:
                models.append(model_str)
            data["basecall_models"] = ",".join(models)
            data["run_info_model"] = model_str

        # Update flow cell IDs
        if bam_info.get("flow_cell_id"):
            flowcell_str = str(bam_info["flow_cell_id"])
            flowcells = data["flowcell_ids"].split(",") if data["flowcell_ids"] else []
            if flowcell_str not in flowcells:
                flowcells.append(flowcell_str)
            data["flowcell_ids"] = ",".join(flowcells)
            data["run_info_flow_cell"] = flowcell_str

        # Update run time
        if bam_info.get("time_of_run"):
            # Convert to string if it's a float or other type
            time_str = str(bam_info["time_of_run"])
            run_times = data["run_time"].split(",") if data["run_time"] else []
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

    def _save_to_csv(self, data: Dict[str, Any], csv_path: str) -> None:
        """Save data to CSV file"""
        df = pd.DataFrame([data])
        df.to_csv(csv_path, index=False)

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
