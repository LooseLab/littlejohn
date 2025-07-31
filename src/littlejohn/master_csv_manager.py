#!/usr/bin/env python3
"""
Master CSV Manager for LittleJohn

This module provides functionality to create and update master.csv files
for each sample, tracking comprehensive metadata across multiple BAM files.
"""

import os
import pandas as pd
from typing import Dict, Any, List
from dataclasses import dataclass
import logging


@dataclass
class MasterCSVManager:
    """Manages the creation and updating of master.csv files for each sample"""
    
    work_dir: str
    
    def update_master_csv(self, sample_id: str, bam_data: Dict[str, Any], bam_info: Dict[str, Any]) -> None:
        """Update or create master.csv file for a sample with new BAM data"""
        try:
            # Create sample directory if it doesn't exist
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            
            master_csv_path = os.path.join(sample_dir, "master.csv")
            
            # Load existing data or create new
            existing_data = self._load_existing_data(master_csv_path)
            
            # Update counters based on BAM state
            if bam_info.get('state') == 'pass':
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
            
            logger = logging.getLogger("littlejohn.master_csv")
            logger.info(f"Updated master.csv for sample {sample_id}")
            
        except Exception as e:
            logger = logging.getLogger("littlejohn.master_csv")
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
            'counter_bam_passed': 0,
            'counter_bam_failed': 0,
            'counter_mapped_count': 0,
            'counter_pass_mapped_count': 0,
            'counter_fail_mapped_count': 0,
            'counter_unmapped_count': 0,
            'counter_pass_unmapped_count': 0,
            'counter_fail_unmapped_count': 0,
            'counter_pass_bases_count': 0,
            'counter_fail_bases_count': 0,
            'counter_bases_count': 0,
            'counter_mapped_reads_num': 0,
            'counter_unmapped_reads_num': 0,
            'counter_pass_mapped_reads_num': 0,
            'counter_fail_mapped_reads_num': 0,
            'counter_pass_unmapped_reads_num': 0,
            'counter_fail_unmapped_reads_num': 0,
            'counter_mapped_bases': 0,
            'counter_unmapped_bases': 0,
            'counter_pass_mapped_bases': 0,
            'counter_fail_mapped_bases': 0,
            'counter_pass_unmapped_bases': 0,
            'counter_fail_unmapped_bases': 0,
            'devices': '',
            'basecall_models': '',
            'run_time': '',
            'flowcell_ids': '',
            'run_info_run_time': '',
            'run_info_device': '',
            'run_info_model': '',
            'run_info_flow_cell': '',
            'bam_tracking_counter': 0,
            'bam_tracking_total_files': 0,
            'bam_files': ''
        }
    
    def _update_pass_counters(self, data: Dict[str, Any], bam_data: Dict[str, Any]) -> Dict[str, Any]:
        """Update counters for pass BAM files"""
        data['counter_bam_passed'] += 1
        data['counter_pass_mapped_count'] += bam_data.get('mapped_reads', 0)
        data['counter_pass_unmapped_count'] += bam_data.get('unmapped_reads', 0)
        data['counter_pass_bases_count'] += bam_data.get('yield_tracking', 0)
        data['counter_pass_mapped_reads_num'] += bam_data.get('mapped_reads_num', 0)
        data['counter_pass_unmapped_reads_num'] += bam_data.get('unmapped_reads_num', 0)
        data['counter_pass_mapped_bases'] += bam_data.get('mapped_bases', 0)
        data['counter_pass_unmapped_bases'] += bam_data.get('unmapped_bases', 0)
        return data
    
    def _update_fail_counters(self, data: Dict[str, Any], bam_data: Dict[str, Any]) -> Dict[str, Any]:
        """Update counters for fail BAM files"""
        data['counter_bam_failed'] += 1
        data['counter_fail_mapped_count'] += bam_data.get('mapped_reads', 0)
        data['counter_fail_unmapped_count'] += bam_data.get('unmapped_reads', 0)
        data['counter_fail_bases_count'] += bam_data.get('yield_tracking', 0)
        data['counter_fail_mapped_reads_num'] += bam_data.get('mapped_reads_num', 0)
        data['counter_fail_unmapped_reads_num'] += bam_data.get('unmapped_reads_num', 0)
        data['counter_fail_mapped_bases'] += bam_data.get('mapped_bases', 0)
        data['counter_fail_unmapped_bases'] += bam_data.get('unmapped_bases', 0)
        return data
    
    def _update_combined_counters(self, data: Dict[str, Any], bam_data: Dict[str, Any]) -> Dict[str, Any]:
        """Update combined counters for all BAM files"""
        data['counter_mapped_count'] += bam_data.get('mapped_reads', 0)
        data['counter_unmapped_count'] += bam_data.get('unmapped_reads', 0)
        data['counter_bases_count'] += bam_data.get('yield_tracking', 0)
        data['counter_mapped_reads_num'] += bam_data.get('mapped_reads_num', 0)
        data['counter_unmapped_reads_num'] += bam_data.get('unmapped_reads_num', 0)
        data['counter_mapped_bases'] += bam_data.get('mapped_bases', 0)
        data['counter_unmapped_bases'] += bam_data.get('unmapped_bases', 0)
        return data
    
    def _update_run_info(self, data: Dict[str, Any], bam_info: Dict[str, Any]) -> Dict[str, Any]:
        """Update run information arrays"""
        # Update devices
        if bam_info.get('device_position'):
            device_str = str(bam_info['device_position'])
            devices = data['devices'].split(',') if data['devices'] else []
            if device_str not in devices:
                devices.append(device_str)
            data['devices'] = ','.join(devices)
            data['run_info_device'] = device_str
        
        # Update basecall models
        if bam_info.get('basecall_model'):
            model_str = str(bam_info['basecall_model'])
            models = data['basecall_models'].split(',') if data['basecall_models'] else []
            if model_str not in models:
                models.append(model_str)
            data['basecall_models'] = ','.join(models)
            data['run_info_model'] = model_str
        
        # Update flow cell IDs
        if bam_info.get('flow_cell_id'):
            flowcell_str = str(bam_info['flow_cell_id'])
            flowcells = data['flowcell_ids'].split(',') if data['flowcell_ids'] else []
            if flowcell_str not in flowcells:
                flowcells.append(flowcell_str)
            data['flowcell_ids'] = ','.join(flowcells)
            data['run_info_flow_cell'] = flowcell_str
        
        # Update run time
        if bam_info.get('time_of_run'):
            # Convert to string if it's a float or other type
            time_str = str(bam_info['time_of_run'])
            run_times = data['run_time'].split(',') if data['run_time'] else []
            if time_str not in run_times:
                run_times.append(time_str)
            data['run_time'] = ','.join(run_times)
            data['run_info_run_time'] = time_str
        
        return data
    
    def _update_bam_tracking(self, data: Dict[str, Any], bam_info: Dict[str, Any]) -> Dict[str, Any]:
        """Update BAM tracking information"""
        data['bam_tracking_counter'] += 1
        data['bam_tracking_total_files'] += 1
        
        # Update bam_files list
        bam_files = data['bam_files'].split(',') if data['bam_files'] else []
        filename = os.path.basename(bam_info.get('file_path', ''))
        if filename and filename not in bam_files:
            bam_files.append(filename)
        data['bam_files'] = ','.join(bam_files)
        
        return data
    
    def _save_to_csv(self, data: Dict[str, Any], csv_path: str) -> None:
        """Save data to CSV file"""
        df = pd.DataFrame([data])
        df.to_csv(csv_path, index=False)