#!/usr/bin/env python3
"""
Target Analysis Module for LittleJohn

This module provides automated target analysis following the exact architecture
specified in the Target Analysis Documentation. It integrates with LittleJohn's 
workflow system and processes files for target-specific analysis.

Features:
- Automated target analysis using various input file types
- Integration with LittleJohn's workflow system
- Sample-specific output directories
- Comprehensive metadata extraction and logging
- Error handling and result tracking
- State persistence and incremental processing

Classes
-------
TargetMetadata
    Container for target analysis metadata and results.

TargetAnalysis
    Main analysis class that processes files for target analysis.

Dependencies
-----------
- pandas: Data manipulation and analysis
- numpy: Numerical computations
- logging: Logging for debugging and monitoring
- typing: Type hints
- tempfile: Temporary file creation
- pathlib: File system paths
- os: Operating system interface
- time: Time utilities
- json: JSON serialization
- pickle: Python object serialization
- gc: Garbage collection
- pysam: BAM file processing
- asyncio: Asynchronous processing support
- subprocess: External command execution

Example Usage
-----------
.. code-block:: python

    from littlejohn.analysis.target_analysis import TargetAnalysis

    # Initialize analysis
    target_analysis = TargetAnalysis(
        work_dir="output/",
        config_path="target_config.json"
    )

    # Process files
    target_analysis.process_file("sample.bam")

Notes
-----
The module follows the LittleJohn framework patterns for:
- Integration with workflow system
- Worker process management
- State tracking and persistence
- Error handling and logging
- Output generation and file management

Authors
-------
Matt Loose
"""

import os
import tempfile
import logging
import time
import pickle
import gc
import json
import asyncio
import subprocess
import shutil
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass, field
from collections import Counter
from io import StringIO
import numpy as np
import pandas as pd
import pysam
from littlejohn.logging_config import get_job_logger


def json_serializable(obj):
    """Convert numpy types to JSON serializable types."""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, pd.Series):
        return obj.tolist()
    elif isinstance(obj, pd.DataFrame):
        return obj.to_dict('records')
    return obj

# Import robin resources for BED files
try:
    from robin import resources
except ImportError:
    resources = None


def run_bedtools(bamfile, bedfile, tempbamfile):
    """
    This function extracts the target sites from the bamfile.

    Parameters
    ----------
    bamfile : str
        Path to the input BAM file
    bedfile : str
        Path to the BED file defining regions
    tempbamfile : str
        Path where the output BAM file should be written
    """
    logger = logging.getLogger("littlejohn.target")
    
    try:
        # Use subprocess.run with shell=True for commands with redirection
        # Or open the output file and redirect stdout there
        with open(tempbamfile, "w") as outfile:
            result = subprocess.run(
                ["bedtools", "intersect", "-a", bamfile, "-b", bedfile],
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )

        if result.returncode != 0:
            logger.error(f"Error running bedtools: {result.stderr}")
            return

        pysam.index(tempbamfile)
    except Exception as e:
        logger.error(f"Error in run_bedtools: {e}")


def get_covdfs(bamfile):
    """
    Extract coverage information from a BAM file.

    This function calculates coverage statistics for both the entire genome
    and specific target regions defined in a BED file.

    Parameters
    ----------
    bamfile : str
        Path to the input BAM file.

    Returns
    -------
    tuple
        A tuple containing two pandas DataFrames:
        - Coverage statistics for the entire genome
        - Coverage statistics for target regions

    Notes
    -----
    Uses pysam for efficient BAM file processing.
    """
    logger = logging.getLogger("littlejohn.target")
    
    try:
        # Get genome-wide coverage using pysam.coverage
        coverage_output = pysam.coverage(f"{bamfile}")
        
        newcovdf = pd.read_csv(StringIO(coverage_output), sep="\t")
        
        logger.info(f"Raw pysam.coverage columns: {list(newcovdf.columns)}")
        logger.info(f"Sample raw coverage data: {newcovdf.head(2).to_dict('records')}")
        
        newcovdf.drop(
            columns=["coverage", "meanbaseq", "meanmapq"],
            inplace=True,
        )
        
        logger.info(f"After dropping columns: {list(newcovdf.columns)}")

        # Find unique_genes.bed file from robin resources
        unique_genes_bed = None
        if resources is not None:
            try:
                unique_genes_bed = os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)),
                    "unique_genes.bed",
                )
                if not os.path.exists(unique_genes_bed):
                    unique_genes_bed = None
            except Exception:
                pass
        
        # Fallback paths for unique_genes.bed
        if unique_genes_bed is None:
            possible_paths = [
                "unique_genes.bed",
                "data/unique_genes.bed",
                "/usr/local/share/unique_genes.bed"
            ]
            for path in possible_paths:
                if os.path.exists(path):
                    unique_genes_bed = path
                    break
        
        if unique_genes_bed is None:
            logger.warning("unique_genes.bed not found, skipping bedcov analysis")
            bedcovdf = pd.DataFrame(columns=["chrom", "startpos", "endpos", "name", "bases"])
        else:
            # Get target region coverage using pysam.bedcov
            bedcovdf = pd.read_csv(
                StringIO(
                    pysam.bedcov(
                        unique_genes_bed,
                        f"{bamfile}",
                    )
                ),
                names=["chrom", "startpos", "endpos", "name", "bases"],
                sep="\t",
            )
        
        logger.info(f"Successfully extracted coverage data from {bamfile}")
        logger.info(f"Genome coverage: {len(newcovdf)} regions")
        logger.info(f"Target coverage: {len(bedcovdf)} regions")
        
        return newcovdf, bedcovdf
        
    except Exception as e:
        logger.error(f"Error in get_covdfs: {str(e)}")
        return None, None


def run_bedmerge(newcovdf, cov_df_main, bedcovdf, bedcov_df_main):
    """
    Merge coverage dataframes for incremental processing.
    
    Parameters
    ----------
    newcovdf : pd.DataFrame
        New genome coverage data
    cov_df_main : pd.DataFrame
        Existing genome coverage data
    bedcovdf : pd.DataFrame
        New target coverage data
    bedcov_df_main : pd.DataFrame
        Existing target coverage data
        
    Returns
    -------
    tuple
        Merged coverage dataframes
    """
    logger = logging.getLogger("littlejohn.target")
    
    logger.info(f"Starting coverage merge:")
    logger.info(f"  New genome coverage: {len(newcovdf) if newcovdf is not None else 0} regions")
    logger.info(f"  Existing genome coverage: {len(cov_df_main) if cov_df_main is not None else 0} regions")
    logger.info(f"  New target coverage: {len(bedcovdf) if bedcovdf is not None else 0} targets")
    logger.info(f"  Existing target coverage: {len(bedcov_df_main) if bedcov_df_main is not None else 0} targets")
    
    # Merge genome coverage data
    if cov_df_main is not None and not cov_df_main.empty:
        merged_df = pd.merge(
            newcovdf,
            cov_df_main,
            on=["#rname", "startpos", "endpos"],
            suffixes=("_df1", "_df2"),
        )
        merged_df["numreads"] = merged_df["numreads_df1"] + merged_df["numreads_df2"]
        merged_df["covbases"] = merged_df["covbases_df1"] + merged_df["covbases_df2"]
        merged_df["meandepth"] = merged_df["meandepth_df1"] + merged_df["meandepth_df2"]

        merged_df.drop(
            columns=[
                "numreads_df1",
                "numreads_df2",
                "meandepth_df1",
                "meandepth_df2",
                "covbases_df1",
                "covbases_df2",
            ],
            inplace=True,
        )
        
        logger.info(f"Merged genome coverage: {len(cov_df_main)} + {len(newcovdf)} -> {len(merged_df)} regions")
        logger.info(f"Sample merged genome data: {merged_df.head(3).to_dict('records')}")
    else:
        merged_df = newcovdf
        logger.info(f"No existing genome coverage, using new data: {len(merged_df)} regions")
    
    # Merge target coverage data
    if bedcov_df_main is not None and not bedcov_df_main.empty:
        merged_bed_df = pd.merge(
            bedcovdf,
            bedcov_df_main,
            on=["chrom", "startpos", "endpos", "name"],
            suffixes=("_df1", "_df2"),
        )
        merged_bed_df["bases"] = merged_bed_df["bases_df1"] + merged_bed_df["bases_df2"]
        merged_bed_df.drop(columns=["bases_df1", "bases_df2"], inplace=True)
        
        logger.info(f"Merged target coverage: {len(bedcov_df_main)} + {len(bedcovdf)} -> {len(merged_bed_df)} targets")
        logger.info(f"Sample merged target data: {merged_bed_df.head(3).to_dict('records')}")
    else:
        merged_bed_df = bedcovdf
        logger.info(f"No existing target coverage, using new data: {len(merged_bed_df)} targets")
    
    return merged_df, merged_bed_df


@dataclass
class TargetMetadata:
    """Container for target analysis metadata and results"""
    sample_id: str
    file_path: str
    analysis_timestamp: float
    target_data_path: Optional[str] = None
    target_plot_path: Optional[str] = None
    analysis_results: Optional[Dict] = None
    processing_steps: List[str] = None
    error_message: Optional[str] = None
    coverage_data: Optional[Dict] = None
    target_bam_path: Optional[str] = None
    coverage_over_time: Optional[np.ndarray] = None
    
    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.analysis_results is None:
            self.analysis_results = {}
        if self.coverage_data is None:
            self.coverage_data = {}


class TargetAnalysis:
    """Target analysis worker"""
    
    def __init__(self, work_dir=None, config_path=None, threads=4):
        logger = logging.getLogger("littlejohn.target")
        
        self.work_dir = work_dir or os.getcwd()
        self.config_path = config_path
        self.threads = threads
        
        # Initialize counter for incremental file naming
        self.file_counter = 1
        
        # Load configuration if provided
        self.config = self._load_config()
        
        # Note: State is now persisted to files instead of in-memory dictionaries
        
        # Find BED file for target extraction
        self.bedfile = self._find_target_bed()
        
        # Configuration parameters
        self.callthreshold = self.config.get('call_threshold', 0.1)
        self.simtime = self.config.get('simtime', False)
        self.reference = self.config.get('reference', None)
        self.snp_calling = self.config.get('snp_calling', False)
        
        logger.info(f"Target Analysis initialized")
        
    def _find_target_bed(self) -> str:
        """Find the target BED file from robin resources"""
        if resources is not None:
            try:
                bed_path = os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)),
                    "rCNS2_panel_name_uniq.bed",
                )
                if os.path.exists(bed_path):
                    return bed_path
            except Exception:
                pass
        
        # Fallback paths
        possible_paths = [
            "rCNS2_panel_name_uniq.bed",
            "data/rCNS2_panel_name_uniq.bed",
            "/usr/local/share/rCNS2_panel_name_uniq.bed"
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        
        # If not found, create a placeholder (this will cause an error later)
        logger = logging.getLogger("littlejohn.target")
        logger.warning("Target BED file not found, will use placeholder")
        return "rCNS2_panel_name_uniq.bed"
    
    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from file or use defaults"""
        if self.config_path and os.path.exists(self.config_path):
            try:
                with open(self.config_path, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logging.warning(f"Error loading config from {self.config_path}: {e}")
        
        # Default configuration
        return {
            'analysis_type': 'standard',
            'output_format': 'json',
            'include_plots': True,
            'threshold': 0.5,
            'call_threshold': 0.1,
            'simtime': False,
            'snp_calling': False
        }
    
    def _get_next_file_number(self) -> int:
        """Get the next file number for incremental processing"""
        return self.file_counter
    
    def _check_and_create_folder(self, base_dir: str, sample_id: str) -> str:
        """Create sample-specific directory and return path"""
        sample_dir = os.path.join(base_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        return sample_dir
    
    def process_file(self, file_path: str, metadata: Dict[str, Any], timestamp: Optional[float] = None) -> TargetMetadata:
        """
        Process a file for target analysis.
        
        Args:
            file_path: Path to the input file
            metadata: File metadata from preprocessing
            timestamp: Optional timestamp for coverage tracking
            
        Returns:
            TargetMetadata object with analysis results
        """
        logger = logging.getLogger("littlejohn.target")
        
        logger.info(f"Processing file: {file_path} for target analysis")
        start_time = time.time()
        
        # Extract sample ID from metadata
        sample_id = metadata.get('sample_id', 'unknown')
        logger.debug(f"Extracted sample_id: {sample_id}")
        
        target_result = TargetMetadata(
            sample_id=sample_id,
            file_path=file_path,
            analysis_timestamp=start_time
        )
        logger.debug(f"Created TargetMetadata object: {target_result}")
        
        logger.info(f"Starting target analysis for sample: {sample_id}")
        logger.debug(f"File: {os.path.basename(file_path)}")
        logger.debug(f"Work directory: {self.work_dir}")
        
        try:
            # Step 1: Validate input file
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Input file not found: {file_path}")
            
            target_result.processing_steps.append("file_validation")
            
            # Step 2: Create sample-specific output directory
            sample_output_dir = self._check_and_create_folder(self.work_dir, sample_id)
            logger.info(f"Created sample output directory: {sample_output_dir}")
            
            target_result.processing_steps.append("directory_created")
            
            # Step 3: Load analysis counter from disk
            analysis_counter = self._load_analysis_counter(sample_id, self.work_dir, logger)
            logger.info(f"Loaded analysis counter for {sample_id}: {analysis_counter}")
            
            # Step 4: Load existing accumulated data from files
            existing_covdf = self._load_existing_coverage_data(sample_output_dir, "coverage_main.csv", logger)
            existing_bedcovdf = self._load_existing_coverage_data(sample_output_dir, "bed_coverage_main.csv", logger)
            existing_coverage_over_time = self._load_existing_coverage_over_time(sample_output_dir, logger)
            existing_target_bam_path = os.path.join(sample_output_dir, "target.bam")
            
            # Step 5: Create temporary directory for processing
            with tempfile.TemporaryDirectory(dir=sample_output_dir) as temp_dir:
                logger.info(f"Created temporary directory: {temp_dir}")
                
                # Step 6: Extract coverage data using get_covdfs
                logger.info(f"Extracting coverage data...")
                newcovdf, bedcovdf = get_covdfs(file_path)
                
                if newcovdf is None or bedcovdf is None:
                    raise RuntimeError("Failed to extract coverage data from BAM file")
                
                # Store coverage data in metadata
                target_result.coverage_data = {
                    'genome_coverage': newcovdf.to_dict('records') if not newcovdf.empty else [],
                    'target_coverage': bedcovdf.to_dict('records') if not bedcovdf.empty else [],
                    'genome_coverage_shape': newcovdf.shape,
                    'target_coverage_shape': bedcovdf.shape
                }
                
                target_result.processing_steps.append("coverage_extracted")
                logger.info(f"Coverage data extracted: genome={newcovdf.shape}, targets={bedcovdf.shape}")
                
                # Step 7: Extract target regions using bedtools
                logger.info(f"Extracting target regions using bedtools...")
                tempbamfile = tempfile.NamedTemporaryFile(
                    dir=sample_output_dir, suffix=".bam"
                )
                
                # Run bedtools intersection
                run_bedtools(file_path, self.bedfile, tempbamfile.name)
                
                # Check if target BAM has reads
                if pysam.AlignmentFile(tempbamfile.name, "rb").count(until_eof=True) > 0:
                    logger.info(f"Target regions found in {sample_id}")
                    
                    # Handle target BAM file accumulation
                    if os.path.exists(existing_target_bam_path):
                        # Merge with existing target BAM
                        tempbamholder = tempfile.NamedTemporaryFile(
                            dir=sample_output_dir, suffix=".bam"
                        )
                        pysam.cat(
                            "-o",
                            tempbamholder.name,
                            existing_target_bam_path,
                            tempbamfile.name,
                        )
                        shutil.copy2(tempbamholder.name, existing_target_bam_path)
                        target_result.target_bam_path = existing_target_bam_path
                        logger.info(f"Updated target BAM file: {existing_target_bam_path}")
                        
                        # Clean up temporary files
                        try:
                            os.remove(f"{tempbamholder.name}.bai")
                        except FileNotFoundError:
                            pass
                    else:
                        # Create new target BAM file
                        shutil.copy2(tempbamfile.name, existing_target_bam_path)
                        target_result.target_bam_path = existing_target_bam_path
                        logger.info(f"Created target BAM file: {existing_target_bam_path}")
                else:
                    logger.info(f"No target regions found in {sample_id}")
                
                # Clean up temporary BAM index
                try:
                    os.remove(f"{tempbamfile.name}.bai")
                except FileNotFoundError:
                    pass
                
                target_result.processing_steps.append("target_regions_extracted")
                
                # Step 8: Update coverage data
                if existing_covdf is not None:
                    updated_covdf, updated_bedcovdf = run_bedmerge(
                        newcovdf, existing_covdf, bedcovdf, existing_bedcovdf
                    )
                else:
                    updated_covdf, updated_bedcovdf = newcovdf, bedcovdf
                
                target_result.processing_steps.append("coverage_updated")
                
                # Step 9: Calculate coverage statistics
                bases = updated_covdf["covbases"].sum()
                genome = updated_covdf["endpos"].sum()
                coverage = bases / genome if genome > 0 else 0.0
                
                logger.info(f"Coverage calculated: {coverage:.4f} ({bases} bases / {genome} genome)")
                logger.info(f"Coverage data shape: {updated_covdf.shape}")
                logger.info(f"Coverage columns: {list(updated_covdf.columns)}")
                logger.info(f"Sample coverage data:")
                logger.info(f"  First few rows: {updated_covdf.head(3).to_dict('records')}")
                logger.info(f"  covbases sum: {bases}")
                logger.info(f"  endpos sum: {genome}")
                
                # Step 10: Track coverage over time
                if self.simtime and timestamp:
                    currenttime = timestamp * 1000
                else:
                    currenttime = time.time() * 1000
                
                if existing_coverage_over_time is not None:
                    updated_coverage_over_time = np.vstack(
                        [existing_coverage_over_time, [(currenttime, coverage)]]
                    )
                else:
                    updated_coverage_over_time = np.array([[currenttime, coverage]])
                
                target_result.coverage_over_time = updated_coverage_over_time
                target_result.processing_steps.append("coverage_tracked")
                
                # Step 11: Save updated coverage data to files
                np.save(
                    os.path.join(sample_output_dir, "coverage_time_chart.npy"),
                    updated_coverage_over_time,
                )
                
                updated_covdf.to_csv(
                    os.path.join(sample_output_dir, "coverage_main.csv"),
                    index=False,
                )
                
                # Calculate length and coverage for bed_coverage_main.csv (matching reference format)
                bed_coverage_main_df = updated_bedcovdf.copy()
                bed_coverage_main_df["length"] = (
                    bed_coverage_main_df["endpos"] - bed_coverage_main_df["startpos"] + 1
                )
                
                bed_coverage_main_df["coverage"] = (
                    bed_coverage_main_df["bases"] / bed_coverage_main_df["length"]
                )
                
                # Reorder columns to match reference format: chrom,startpos,endpos,name,length,coverage,bases
                bed_coverage_main_df = bed_coverage_main_df[["chrom", "startpos", "endpos", "name", "length", "coverage", "bases"]]
                
                bed_coverage_main_df.to_csv(
                    os.path.join(sample_output_dir, "bed_coverage_main.csv"),
                    index=False,
                )
                
                target_result.processing_steps.append("coverage_saved")
                
                # Step 12: Calculate target coverage
                target_coverage_df = updated_bedcovdf.copy()
                target_coverage_df["length"] = (
                    target_coverage_df["endpos"] - target_coverage_df["startpos"] + 1
                )
                
                target_coverage_df["coverage"] = (
                    target_coverage_df["bases"] / target_coverage_df["length"]
                )
                
                # Reorder columns to match reference format: chrom,startpos,endpos,name,length,coverage,bases
                target_coverage_df = target_coverage_df[["chrom", "startpos", "endpos", "name", "length", "coverage", "bases"]]
                
                target_coverage_df.to_csv(
                    os.path.join(sample_output_dir, "target_coverage.csv"),
                    index=False,
                )
                
                target_result.processing_steps.append("target_coverage_calculated")
                
                # Step 13: Identify targets exceeding threshold
                run_list = target_coverage_df[
                    target_coverage_df["coverage"].ge(self.callthreshold)
                ]
                
                # Load existing targets exceeding threshold count
                targets_exceeding_file = os.path.join(sample_output_dir, "targets_exceeding_threshold_count.txt")
                existing_targets_count = 0
                if os.path.exists(targets_exceeding_file):
                    try:
                        with open(targets_exceeding_file, 'r') as f:
                            existing_targets_count = int(f.read().strip())
                    except (ValueError, IOError):
                        existing_targets_count = 0
                
                if self.reference:
                    logger.info(f"Reference genome found: {self.reference}")
                    
                    if len(run_list) > 0:
                        logger.info(f"Found {len(run_list)} regions exceeding threshold")
                        
                        # Save updated count
                        with open(targets_exceeding_file, 'w') as f:
                            f.write(str(len(run_list)))
                        
                        run_list[["chrom", "startpos", "endpos"]].to_csv(
                            os.path.join(sample_output_dir, "targets_exceeding_threshold.bed"),
                            sep="\t",
                            header=None,
                            index=None,
                        )
                        
                        logger.info(f"Targets exceeding threshold saved: {len(run_list)} regions")
                    else:
                        logger.info("No regions exceed coverage threshold")
                else:
                    logger.info("No reference genome provided")
                
                target_result.processing_steps.append("threshold_analysis_complete")
                
                # Step 13: Perform target analysis
                logger.info(f"Performing target analysis...")
                analysis_results = self._perform_target_analysis(file_path, temp_dir, metadata, newcovdf, bedcovdf)
                target_result.analysis_results = analysis_results
                target_result.processing_steps.append("analysis_completed")
                
                # Step 14: Generate output files
                logger.info(f"Generating output files...")
                output_files = self._generate_output_files(
                    sample_output_dir, analysis_counter, analysis_results, newcovdf, bedcovdf, logger, metadata
                )
                # No JSON files created - data is accumulated in CSV files
                target_result.target_data_path = None
                target_result.target_plot_path = None
                target_result.processing_steps.append("output_generated")
                
                # Step 15: Update analysis counter
                analysis_counter += 1
                self._save_analysis_counter(sample_id, analysis_counter, self.work_dir, logger)
                target_result.processing_steps.append("counter_updated")
                
                # Step 16: Force garbage collection
                gc.collect()
                
                target_result.processing_steps.append("target_analysis_complete")
                logger.info(f"Target analysis complete for {sample_id}")
                logger.info(f"   Results: {len(analysis_results)} items")
                logger.info(f"   Output files: {len(output_files)} files")
                logger.info(f"   Targets exceeding threshold: {len(run_list) if 'run_list' in locals() else 0}")
                
                return target_result
                
        except Exception as e:
            error_details = f"Error in target analysis for {sample_id}: {str(e)}"
            logger.error(error_details)
            target_result.error_message = error_details
            target_result.processing_steps.append("analysis_failed")
            return target_result
    
    def _perform_target_analysis(self, file_path: str, temp_dir: str, metadata: Dict[str, Any], 
                                newcovdf: pd.DataFrame, bedcovdf: pd.DataFrame) -> Dict[str, Any]:
        """
        Perform the actual target analysis.
        
        Args:
            file_path: Path to input file
            temp_dir: Temporary directory for processing
            metadata: File metadata
            newcovdf: Genome coverage DataFrame
            bedcovdf: Target coverage DataFrame
            
        Returns:
            Dictionary containing analysis results
        """
        logger = logging.getLogger("littlejohn.target")
        
        # This is where the actual target analysis logic would go
        # For now, we'll create a placeholder implementation that uses the coverage data
        
        results = {
            'analysis_type': self.config.get('analysis_type', 'standard'),
            'input_file': os.path.basename(file_path),
            'sample_id': metadata.get('sample_id', 'unknown'),
            'timestamp': time.time(),
            'targets_found': 0,
            'targets_analyzed': [],
            'coverage_summary': {
                'genome_regions': len(newcovdf) if newcovdf is not None else 0,
                'target_regions': len(bedcovdf) if bedcovdf is not None else 0,
                'total_bases_covered': 0,
                'average_coverage': 0.0
            },
            'summary_stats': {
                'total_targets': 0,
                'significant_targets': 0,
                'average_score': 0.0
            }
        }
        
        # Calculate coverage statistics
        if newcovdf is not None and not newcovdf.empty:
            results['coverage_summary']['total_bases_covered'] = newcovdf.get('numreads', pd.Series([0])).sum()
            results['coverage_summary']['average_coverage'] = newcovdf.get('numreads', pd.Series([0])).mean()
        
        # Analyze target regions
        if bedcovdf is not None and not bedcovdf.empty:
            results['targets_found'] = len(bedcovdf)
            results['summary_stats']['total_targets'] = len(bedcovdf)
            
            # Calculate significance based on coverage
            significant_threshold = self.config.get('threshold', 0.5)
            significant_targets = bedcovdf[bedcovdf['bases'] > significant_threshold * bedcovdf['bases'].max()]
            results['summary_stats']['significant_targets'] = len(significant_targets)
            
            # Add target analysis results
            for idx, row in bedcovdf.iterrows():
                target_info = {
                    'target_id': row.get('name', f'target_{idx+1:03d}'),
                    'position': f"{row['chrom']}:{row['startpos']}-{row['endpos']}",
                    'coverage': row['bases'],
                    'significance': 'high' if row['bases'] > significant_threshold * bedcovdf['bases'].max() else 'medium'
                }
                results['targets_analyzed'].append(target_info)
            
            # Calculate average score
            if results['targets_analyzed']:
                scores = [t['coverage'] for t in results['targets_analyzed']]
                results['summary_stats']['average_score'] = sum(scores) / len(scores)
        
        logger.info(f"Analysis completed: {results['targets_found']} targets found")
        logger.info(f"Coverage summary: {results['coverage_summary']}")
        return results
    
    def _generate_output_files(self, sample_dir: str, analysis_counter: int, 
                              analysis_results: Dict[str, Any], newcovdf: pd.DataFrame, 
                              bedcovdf: pd.DataFrame, logger, metadata: Dict[str, Any]) -> Dict[str, str]:
        """
        Generate minimal output files for the analysis.
        
        Args:
            sample_dir: Sample output directory
            analysis_counter: Current analysis counter
            analysis_results: Analysis results
            newcovdf: Genome coverage DataFrame
            bedcovdf: Target coverage DataFrame
            logger: Logger instance
            metadata: File metadata
            
        Returns:
            Dictionary with paths to generated files
        """
        output_files = {}
        
        try:
            # No JSON files needed - accumulated data is already saved in CSV files
            # The coverage data is persisted in:
            # - coverage_main.csv (genome coverage)
            # - bed_coverage_main.csv (target coverage) 
            # - coverage_time_chart.npy (coverage over time)
            # - target_coverage.csv (calculated target coverage)
            
            logger.info(f"No additional output files needed - data accumulated in CSV files")
            return output_files
            
        except Exception as e:
            logger.error(f"Error in output file generation: {e}")
            return {}
    
    def _load_analysis_counter(self, sample_id: str, work_dir: str, logger) -> int:
        """Load the analysis counter for a sample from disk"""
        counter_file = os.path.join(work_dir, sample_id, "target_analysis_counter.txt")
        if os.path.exists(counter_file):
            try:
                with open(counter_file, 'r') as f:
                    return int(f.read().strip())
            except (ValueError, IOError) as e:
                logger.warning(f"Error loading counter for {sample_id}: {e}")
        return 0
    
    def _save_analysis_counter(self, sample_id: str, counter: int, work_dir: str, logger) -> None:
        """Save the analysis counter for a sample to disk"""
        counter_file = os.path.join(work_dir, sample_id, "target_analysis_counter.txt")
        try:
            # Create sample directory if it doesn't exist
            os.makedirs(os.path.dirname(counter_file), exist_ok=True)
            with open(counter_file, 'w') as f:
                f.write(str(counter))
        except IOError as e:
            logger.error(f"Error saving counter for {sample_id}: {e}")
    
    def _load_existing_coverage_data(self, sample_dir: str, filename: str, logger) -> Optional[pd.DataFrame]:
        """Load existing coverage data from CSV file if it exists"""
        file_path = os.path.join(sample_dir, filename)
        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path)
                logger.info(f"Loaded existing coverage data from {filename}: {df.shape}")
                
                # For bed_coverage_main.csv, extract only the raw data columns needed for merging
                if filename == "bed_coverage_main.csv" and "length" in df.columns and "coverage" in df.columns:
                    # Extract only the raw data columns: chrom,startpos,endpos,name,bases
                    df = df[["chrom", "startpos", "endpos", "name", "bases"]]
                    logger.info(f"Extracted raw data columns from {filename}: {df.shape}")
                
                return df
            except Exception as e:
                logger.warning(f"Error loading existing coverage data from {filename}: {e}")
                return None
        return None
    
    def _load_existing_coverage_over_time(self, sample_dir: str, logger) -> Optional[np.ndarray]:
        """Load existing coverage over time data from numpy file if it exists"""
        file_path = os.path.join(sample_dir, "coverage_time_chart.npy")
        if os.path.exists(file_path):
            try:
                data = np.load(file_path)
                logger.info(f"Loaded existing coverage over time data: {data.shape}")
                return data
            except Exception as e:
                logger.warning(f"Error loading existing coverage over time data: {e}")
                return None
        return None


def process_single_file(file_path: str, metadata: Dict[str, Any], work_dir: str, logger) -> Dict[str, Any]:
    """
    Process a single file for target analysis using the complete pipeline.
    
    Args:
        file_path: Path to input file
        metadata: File metadata
        work_dir: Working directory
        logger: Logger instance
            
    Returns:
        Dictionary with target analysis results
    """
    sample_id = metadata.get('sample_id', 'unknown')
    logger.info(f"🎯 Starting target analysis for sample: {sample_id}")
    
    # Create sample-specific output directory
    sample_output_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)
    logger.info(f"Created sample output directory: {sample_output_dir}")
    
    analysis_result = {
        'sample_id': sample_id,
        'file_path': file_path,
        'analysis_timestamp': time.time(),
        'target_data_path': None,
        'target_plot_path': None,
        'analysis_results': {},
        'processing_steps': [],
        'error_message': None,
        'coverage_data': {},
        'target_bam_path': None,
        'coverage_over_time': None
    }
    
    try:
        # Initialize target analysis
        target_analysis = TargetAnalysis(work_dir=work_dir)
        
        # Process the file
        target_metadata = target_analysis.process_file(file_path, metadata)
        
        # Update analysis result with target metadata
        analysis_result.update({
            'target_data_path': target_metadata.target_data_path,
            'target_plot_path': target_metadata.target_plot_path,
            'analysis_results': target_metadata.analysis_results,
            'processing_steps': target_metadata.processing_steps,
            'error_message': target_metadata.error_message,
            'coverage_data': target_metadata.coverage_data,
            'target_bam_path': target_metadata.target_bam_path,
            'coverage_over_time': target_metadata.coverage_over_time
        })
        
        if target_metadata.error_message:
            logger.error(f"Target analysis failed: {target_metadata.error_message}")
        else:
            logger.info(f"Target analysis complete for {sample_id}")
            logger.info(f"   Results: {len(target_metadata.analysis_results)} items")
            logger.info(f"   Processing steps: {', '.join(target_metadata.processing_steps)}")
            logger.info(f"   Coverage data: genome={target_metadata.coverage_data.get('genome_coverage_shape', (0,0))}, targets={target_metadata.coverage_data.get('target_coverage_shape', (0,0))}")
            if target_metadata.target_bam_path:
                logger.info(f"   Target BAM: {target_metadata.target_bam_path}")
        
        return analysis_result
        
    except Exception as e:
        error_details = f"Error in target analysis for {sample_id}: {str(e)}"
        logger.error(error_details)
        analysis_result['error_message'] = error_details
        analysis_result['processing_steps'].append('analysis_failed')
        return analysis_result


def target_handler(job, work_dir=None):
    """
    Handler function for target analysis jobs.
    This function processes files for target-specific analysis.
    
    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to file directory)
    """
    # Get job-specific logger
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    
    try:
        file_path = job.context.filepath
        
        logger.info(f"Starting target analysis for: {os.path.basename(file_path)}")
        
        # Get metadata from preprocessing
        file_metadata = job.context.metadata.get('bam_metadata', {})
        
        # Determine work directory
        if work_dir is None:
            # Default to file directory
            work_dir = os.path.dirname(file_path)
        else:
            # Use specified work directory, create if it doesn't exist
            os.makedirs(work_dir, exist_ok=True)
        
        # Process the file
        result = process_single_file(file_path, file_metadata, work_dir, logger)
        
        # Store results in job context
        job.context.add_metadata('target_analysis', result)
        
        if result.get('error_message'):
            job.context.add_error('target_analysis', result['error_message'])
            logger.error(f"Target analysis failed: {result['error_message']}")
        else:
            job.context.add_result('target_analysis', {
                'status': 'success',
                'sample_id': result.get('sample_id', 'unknown'),
                'analysis_time': result.get('analysis_timestamp', 0),
                'targets_found': result.get('analysis_results', {}).get('targets_found', 0),
                'processing_steps': result.get('processing_steps', []),
                'target_data_path': result.get('target_data_path', ''),
                'target_plot_path': result.get('target_plot_path', ''),
                'target_bam_path': result.get('target_bam_path', ''),
                'coverage_summary': result.get('analysis_results', {}).get('coverage_summary', {})
            })
            logger.info(f"Target analysis complete for {os.path.basename(file_path)}")
            logger.info(f"Sample ID: {result.get('sample_id', 'unknown')}")
            logger.info(f"Targets found: {result.get('analysis_results', {}).get('targets_found', 0)}")
            
    except Exception as e:
        error_details = f"Error in target analysis for {job.context.filepath}: {str(e)}"
        job.context.add_error('target_analysis', error_details)
        logger.error(error_details) 