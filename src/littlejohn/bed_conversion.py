#!/usr/bin/env python3
"""
BAM to parquet file conversion module for LittleJohn.

This module provides automated conversion of BAM files to parquet format using
matkit from the robin package. The analysis integrates with LittleJohn's
preprocessing pipeline and provides comprehensive metadata output.

Features:
- Automated BAM processing with LittleJohn's preprocessing pipeline
- BAM to parquet conversion using matkit from robin package
- Parquet file management for incremental processing
- Integration with CPG master file for methylation analysis
- Comprehensive metadata extraction and logging
"""

import os
import time
import tempfile
import shutil
import itertools
import logging
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional
from pathlib import Path

# Import robin utilities
try:
    from robin.utils import run_matkit
    from robin import resources
except ImportError:
    run_matkit = None
    resources = None

# Import local utilities for merge_modkit_files
try:
    from littlejohn.temp_utilities import merge_modkit_files
except ImportError:
    merge_modkit_files = None

from littlejohn.logging_config import get_job_logger


@dataclass
class BedConversionMetadata:
    """Container for BAM to parquet conversion metadata and results"""
    sample_id: str
    bam_path: str
    analysis_timestamp: float
    parquet_path: Optional[str] = None
    processing_steps: List[str] = field(default_factory=list)
    error_message: Optional[str] = None
    results: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.results is None:
            self.results = {}


class BedConversionAnalysis:
    """BAM to parquet conversion analysis worker"""
    
    def __init__(self, work_dir=None, threads=4):
        self.work_dir = work_dir or os.getcwd()
        self.threads = threads
        
        # Find required files and paths
        self.cpgs_master_file = self._find_cpgs_master_file()
        
        # Initialize counter for incremental file naming
        self.file_counter = 1
        
        # Use logger for initialization
        logger = logging.getLogger("littlejohn.bed_conversion")
        logger.info(f"BAM to parquet conversion analysis initialized")
        logger.debug(f"Work directory: {self.work_dir}")
        logger.debug(f"CPGs master file: {self.cpgs_master_file}")
        logger.debug(f"Threads: {self.threads}")
        
    def _find_cpgs_master_file(self) -> str:
        """Find the CPGs master file from robin resources"""
        if resources is not None:
            try:
                cpgs_path = os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)),
                    "sturg_nanodx_cpgs_0125.bed.gz",
                )
                if os.path.exists(cpgs_path):
                    return cpgs_path
            except Exception:
                pass
        
        # Fallback paths
        possible_paths = [
            "sturg_nanodx_cpgs_0125.bed.gz",
            "data/sturg_nanodx_cpgs_0125.bed.gz",
            "/usr/local/share/sturg_nanodx_cpgs_0125.bed.gz"
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        
        # If not found, create a placeholder (this will cause an error later)
        logger = logging.getLogger("littlejohn.bed_conversion")
        logger.warning("CPGs master file not found, will use placeholder")
        return "sturg_nanodx_cpgs_0125.bed.gz"
    
    def _get_next_file_number(self) -> int:
        """Get the next file number for incremental naming"""
        current = self.file_counter
        self.file_counter += 1
        return current
    
    def process_bam_file(self, bam_path: str, metadata: Dict[str, Any]) -> BedConversionMetadata:
        """Process a single BAM file for parquet conversion"""
        logger = logging.getLogger("littlejohn.bed_conversion")
        
        logger.info(f"Processing BAM file: {bam_path} for parquet conversion")
        sample_id = metadata.get('sample_id', 'unknown')
        start_time = time.time()
        
        # Get next file number for incremental naming
        file_number = self._get_next_file_number()
        
        bed_result = BedConversionMetadata(
            sample_id=sample_id,
            bam_path=bam_path,
            analysis_timestamp=start_time
        )
        
        logger.info(f"Starting BAM to parquet conversion for sample: {sample_id}")
        logger.debug(f"BAM file: {os.path.basename(bam_path)}")
        logger.debug(f"Work directory: {self.work_dir}")
        
        try:
            # Step 1: Check if BAM file exists and is readable
            if not os.path.exists(bam_path):
                raise FileNotFoundError(f"BAM file not found: {bam_path}")
            
            bed_result.processing_steps.append("file_validation")
            
            # Step 2: Create sample-specific directory
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            logger.info(f"Created sample directory: {sample_dir}")
            
            bed_result.processing_steps.append("directory_created")
            
            # Step 3: Process BAM file with matkit
            logger.info(f"Converting BAM to parquet using matkit...")
            
            # Create output path specific to this sample
            parquet_path = os.path.join(
                self.work_dir,
                sample_id,
                f"{sample_id}.parquet",  # Use sampleID for the output filename
            )
            
            # Process the BAM file - this follows the exact same pattern as the working code
            data = self._process_bams([bam_path], sample_dir)
            
            # Update state with new data - this creates the parquet file using the same approach
            state = self._update_state(parquet_path, data, sample_id, file_number)
            
            bed_result.parquet_path = state
            bed_result.processing_steps.append("parquet_conversion_complete")
            
            logger.info(f"BAM to parquet conversion complete for {sample_id}")
            logger.info(f"   Parquet file: {parquet_path}")
            logger.info(f"   Processed files: {len(data)}")
            
            return bed_result
            
        except Exception as e:
            logger.error(f"Error in BAM to parquet conversion for {sample_id}: {e}")
            bed_result.error_message = str(e)
            bed_result.processing_steps.append('conversion_failed')
            return bed_result
    
    def _process_bams(self, bams: List[str], work_dir: str) -> List[str]:
        """Process BAM files using matkit and return list of processed file paths"""
        logger = logging.getLogger("littlejohn.bed_conversion")
        processed_files = []
        
        for bam in bams:
            logger.debug(f"Processing BAM file: {bam}")
            
            # Create temporary file for matkit output - exactly like the working code
            temp_file = tempfile.NamedTemporaryFile(
                mode='w',
                suffix='.txt',
                prefix='modkit_',
                dir=work_dir,
                delete=False
            )
            temp_file.close()
            
            try:
                # Run matkit on the BAM file - using the same function as working code
                if run_matkit is not None:
                    run_matkit(bam, temp_file.name)
                else:
                    # Fallback if robin is not available
                    with open(temp_file.name, 'w') as f:
                        f.write(f"# Dummy output for {bam}\n")
                
                processed_files.append(temp_file.name)
                logger.debug(f"Successfully processed: {temp_file.name}")
                
            except Exception as e:
                logger.error(f"Failed to process {bam}: {e}")
                # Clean up the temporary file
                try:
                    os.remove(temp_file.name)
                except:
                    pass
                raise
        
        return processed_files
    
    def _update_state(self, state: str, data: List[str], sample_id: str, file_number: int) -> str:
        """Update the state with new data from BAM processing - creates parquet file directly"""
        logger = logging.getLogger("littlejohn.bed_conversion")
        
        if merge_modkit_files is None:
            logger.error("littlejohn.temp_utilities.merge_modkit_files not available")
            raise ImportError("littlejohn.temp_utilities.merge_modkit_files is required for parquet creation")
        
        # Configuration for mnpflex - exactly like the working code
        mnpflex_config = {
            "mnpuser": None,
            "mnppass": None
        }
        
        logger.debug(f"Creating parquet file: {state}")
        logger.debug(f"Processing {len(data)} matkit files")
        
        # Use merge_modkit_files exactly like the working code
        # This function is designed to create binary parquet files, not text files
        merge_modkit_files(
            data, 
            state,  # Output parquet file path
            state,  # State parquet file path (same as output)
            self.cpgs_master_file, 
            sample_id, 
            self.work_dir, 
            mnpflex_config, 
            1  # Number of BAM files that contributed
        )
        
        # Verify the parquet file was created
        if not os.path.exists(state):
            raise FileNotFoundError(f"Parquet file was not created: {state}")
        
        parquet_size = os.path.getsize(state)
        logger.debug(f"Parquet file created successfully: {state} ({parquet_size} bytes)")
        
        # Clean up temporary files - exactly like the working code
        for file_path in data:
            try:
                os.remove(file_path)
                logger.debug(f"Deleted temporary file: {file_path}")
            except OSError as e:
                logger.warning(f"Failed to delete temporary file {file_path}: {e}")
        
        return state


def grouper(iterable, n):
    """Group items from iterable into chunks of size n"""
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk


def bed_conversion_handler(job, work_dir=None):
    """
    Handler function for BAM to parquet conversion jobs.
    This function processes BAM files for parquet conversion analysis.
    
    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to BAM file directory)
    """
    try:
        bam_path = job.context.filepath
        
        # Get job-specific logger
        logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
        logger.info(f"Starting BAM to parquet conversion for: {os.path.basename(bam_path)}")
        
        # Get metadata from preprocessing
        bam_metadata = job.context.metadata.get('bam_metadata', {})
        
        # Determine work directory
        if work_dir is None:
            # Default to BAM file directory
            work_dir = os.path.dirname(bam_path)
        else:
            # Use specified work directory, create if it doesn't exist
            os.makedirs(work_dir, exist_ok=True)
            logger.debug(f"Using specified work directory: {work_dir}")
        
        # Create BAM to parquet conversion analysis instance
        bed_analyzer = BedConversionAnalysis(work_dir=work_dir)
        
        # Process the BAM file
        bed_result = bed_analyzer.process_bam_file(bam_path, bam_metadata)
        
        # Store results in job context
        job.context.add_metadata('bed_conversion', bed_result.results)
        job.context.add_metadata('bed_processing_steps', bed_result.processing_steps)
        
        if bed_result.error_message:
            job.context.add_error('bed_conversion', bed_result.error_message)
            logger.error(f"BAM to parquet conversion failed: {bed_result.error_message}")
        else:
            job.context.add_result('bed_conversion', {
                'status': 'success',
                'sample_id': bed_result.sample_id,
                'analysis_time': bed_result.analysis_timestamp,
                'parquet_path': bed_result.parquet_path,
                'processing_steps': bed_result.processing_steps
            })
            logger.info(f"BAM to parquet conversion complete for {os.path.basename(bam_path)}")
            logger.info(f"Sample ID: {bed_result.sample_id}")
            logger.info(f"Parquet file: {bed_result.parquet_path}")
            logger.debug(f"Processing steps: {', '.join(bed_result.processing_steps)}")
            
    except Exception as e:
        job.context.add_error('bed_conversion', str(e))
        logger.error(f"Error in BAM to parquet conversion for {job.context.filepath}: {e}") 