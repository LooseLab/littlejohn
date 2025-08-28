#!/usr/bin/env python3
"""
Target Analysis Module for robin

This module provides automated target analysis following the exact architecture
specified in the Target Analysis Documentation. It integrates with robin's
workflow system and processes files for target-specific analysis.

Features:
- Automated target analysis using various input file types
- Integration with robin's workflow system
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

    from robin.analysis.target_analysis import TargetAnalysis

    # Initialize analysis
    target_analysis = TargetAnalysis(
        work_dir="output/",
        config_path="target_config.json"
    )

    # Process files
    target_analysis.process_file("sample.bam")

Notes
-----
The module follows the robin framework patterns for:
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
import gc
import json
import subprocess
import shutil
from typing import Dict, Any, Optional, List
from dataclasses import dataclass
from io import StringIO
import numpy as np
import pandas as pd
import pysam
from robin.logging_config import get_job_logger

# Optional import for Docker functionality
try:
    import docker
except ImportError:
    docker = None


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
        return obj.to_dict("records")
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
    logger = logging.getLogger("robin.target")

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
    logger = logging.getLogger("robin.target")

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
                "/usr/local/share/unique_genes.bed",
            ]
            for path in possible_paths:
                if os.path.exists(path):
                    unique_genes_bed = path
                    break

        if unique_genes_bed is None:
            logger.warning("unique_genes.bed not found, skipping bedcov analysis")
            bedcovdf = pd.DataFrame(
                columns=["chrom", "startpos", "endpos", "name", "bases"]
            )
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
    logger = logging.getLogger("robin.target")

    logger.info("Starting coverage merge:")
    logger.info(
        f"  New genome coverage: {len(newcovdf) if newcovdf is not None else 0} regions"
    )
    logger.info(
        f"  Existing genome coverage: {len(cov_df_main) if cov_df_main is not None else 0} regions"
    )
    logger.info(
        f"  New target coverage: {len(bedcovdf) if bedcovdf is not None else 0} targets"
    )
    logger.info(
        f"  Existing target coverage: {len(bedcov_df_main) if bedcov_df_main is not None else 0} targets"
    )

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

        logger.info(
            f"Merged genome coverage: {len(cov_df_main)} + {len(newcovdf)} -> {len(merged_df)} regions"
        )
        logger.info(
            f"Sample merged genome data: {merged_df.head(3).to_dict('records')}"
        )
    else:
        merged_df = newcovdf
        logger.info(
            f"No existing genome coverage, using new data: {len(merged_df)} regions"
        )

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

        logger.info(
            f"Merged target coverage: {len(bedcov_df_main)} + {len(bedcovdf)} -> {len(merged_bed_df)} targets"
        )
        logger.info(
            f"Sample merged target data: {merged_bed_df.head(3).to_dict('records')}"
        )
    else:
        merged_bed_df = bedcovdf
        logger.info(
            f"No existing target coverage, using new data: {len(merged_bed_df)} targets"
        )

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
        logger = logging.getLogger("robin.target")

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
        self.callthreshold = self.config.get("call_threshold", 0.1)
        self.simtime = self.config.get("simtime", False)
        self.reference = self.config.get("reference", None)
        self.snp_calling = self.config.get("snp_calling", False)
        
        # Check environment variable for reference genome
        if not self.reference:
            self.reference = os.environ.get("LITTLEJOHN_REFERENCE")
        
        # Log reference genome status
        if self.reference:
            logger.info(f"Reference genome configured: {self.reference}")
        else:
            logger.info("No reference genome configured - SNP calling will not be available")

        logger.info("Target Analysis initialized")

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
            "/usr/local/share/rCNS2_panel_name_uniq.bed",
        ]

        for path in possible_paths:
            if os.path.exists(path):
                return path

        # If not found, create a placeholder (this will cause an error later)
        logger = logging.getLogger("robin.target")
        logger.warning("Target BED file not found, will use placeholder")
        return "rCNS2_panel_name_uniq.bed"

    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from file or use defaults"""
        if self.config_path and os.path.exists(self.config_path):
            try:
                with open(self.config_path, "r") as f:
                    return json.load(f)
            except Exception as e:
                logging.warning(f"Error loading config from {self.config_path}: {e}")

        # Default configuration
        return {
            "analysis_type": "standard",
            "output_format": "json",
            "include_plots": True,
            "threshold": 0.5,
            "call_threshold": 0.1,
            "simtime": False,
            "snp_calling": False,
        }

    def _get_next_file_number(self) -> int:
        """Get the next file number for incremental processing"""
        return self.file_counter

    def _check_and_create_folder(self, base_dir: str, sample_id: str) -> str:
        """Create sample-specific directory and return path"""
        sample_dir = os.path.join(base_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        return sample_dir

    def process_file(
        self,
        file_path: str,
        metadata: Dict[str, Any],
        timestamp: Optional[float] = None,
    ) -> TargetMetadata:
        """
        Process a file for target analysis.

        Args:
            file_path: Path to the input file
            metadata: File metadata from preprocessing
            timestamp: Optional timestamp for coverage tracking

        Returns:
            TargetMetadata object with analysis results
        """
        logger = logging.getLogger("robin.target")

        logger.info(f"Processing file: {file_path} for target analysis")
        start_time = time.time()

        # Extract sample ID from metadata
        sample_id = metadata.get("sample_id", "unknown")
        logger.debug(f"Extracted sample_id: {sample_id}")

        target_result = TargetMetadata(
            sample_id=sample_id, file_path=file_path, analysis_timestamp=start_time
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
            analysis_counter = self._load_analysis_counter(
                sample_id, self.work_dir, logger
            )
            logger.info(f"Loaded analysis counter for {sample_id}: {analysis_counter}")

            # Step 4: Load existing accumulated data from files
            existing_covdf = self._load_existing_coverage_data(
                sample_output_dir, "coverage_main.csv", logger
            )
            existing_bedcovdf = self._load_existing_coverage_data(
                sample_output_dir, "bed_coverage_main.csv", logger
            )
            existing_coverage_over_time = self._load_existing_coverage_over_time(
                sample_output_dir, logger
            )
            existing_target_bam_path = os.path.join(sample_output_dir, "target.bam")

            # Step 5: Create temporary directory for processing
            with tempfile.TemporaryDirectory(dir=sample_output_dir) as temp_dir:
                logger.info(f"Created temporary directory: {temp_dir}")

                # Step 6: Extract coverage data using get_covdfs
                logger.info("Extracting coverage data...")
                newcovdf, bedcovdf = get_covdfs(file_path)

                if newcovdf is None or bedcovdf is None:
                    raise RuntimeError("Failed to extract coverage data from BAM file")

                # Store coverage data in metadata
                target_result.coverage_data = {
                    "genome_coverage": (
                        newcovdf.to_dict("records") if not newcovdf.empty else []
                    ),
                    "target_coverage": (
                        bedcovdf.to_dict("records") if not bedcovdf.empty else []
                    ),
                    "genome_coverage_shape": newcovdf.shape,
                    "target_coverage_shape": bedcovdf.shape,
                }

                target_result.processing_steps.append("coverage_extracted")
                logger.info(
                    f"Coverage data extracted: genome={newcovdf.shape}, targets={bedcovdf.shape}"
                )

                # Step 7: Extract target regions using bedtools
                logger.info("Extracting target regions using bedtools...")
                tempbamfile = tempfile.NamedTemporaryFile(
                    dir=sample_output_dir, suffix=".bam"
                )

                # Run bedtools intersection
                run_bedtools(file_path, self.bedfile, tempbamfile.name)

                # Check if target BAM has reads
                if (
                    pysam.AlignmentFile(tempbamfile.name, "rb").count(until_eof=True)
                    > 0
                ):
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
                        logger.info(
                            f"Updated target BAM file: {existing_target_bam_path}"
                        )

                        # Clean up temporary files
                        try:
                            os.remove(f"{tempbamholder.name}.bai")
                        except FileNotFoundError:
                            pass
                    else:
                        # Create new target BAM file
                        shutil.copy2(tempbamfile.name, existing_target_bam_path)
                        target_result.target_bam_path = existing_target_bam_path
                        logger.info(
                            f"Created target BAM file: {existing_target_bam_path}"
                        )
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

                logger.info(
                    f"Coverage calculated: {coverage:.4f} ({bases} bases / {genome} genome)"
                )
                logger.info(f"Coverage data shape: {updated_covdf.shape}")
                logger.info(f"Coverage columns: {list(updated_covdf.columns)}")
                logger.info("Sample coverage data:")
                logger.info(
                    f"  First few rows: {updated_covdf.head(3).to_dict('records')}"
                )
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
                    bed_coverage_main_df["endpos"]
                    - bed_coverage_main_df["startpos"]
                    + 1
                )

                bed_coverage_main_df["coverage"] = (
                    bed_coverage_main_df["bases"] / bed_coverage_main_df["length"]
                )

                # Reorder columns to match reference format: chrom,startpos,endpos,name,length,coverage,bases
                bed_coverage_main_df = bed_coverage_main_df[
                    [
                        "chrom",
                        "startpos",
                        "endpos",
                        "name",
                        "length",
                        "coverage",
                        "bases",
                    ]
                ]

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
                target_coverage_df = target_coverage_df[
                    [
                        "chrom",
                        "startpos",
                        "endpos",
                        "name",
                        "length",
                        "coverage",
                        "bases",
                    ]
                ]

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
                targets_exceeding_file = os.path.join(
                    sample_output_dir, "targets_exceeding_threshold_count.txt"
                )
                existing_targets_count = 0
                if os.path.exists(targets_exceeding_file):
                    try:
                        with open(targets_exceeding_file, "r") as f:
                            existing_targets_count = int(f.read().strip())
                    except (ValueError, IOError):
                        existing_targets_count = 0

                # Always generate targets_exceeding_threshold.bed for potential SNP analysis
                if len(run_list) > 0:
                    logger.info(
                        f"Found {len(run_list)} regions exceeding threshold"
                    )

                    # Save updated count
                    with open(targets_exceeding_file, "w") as f:
                        f.write(str(len(run_list)))

                    # Generate BED file for regions exceeding threshold
                    run_list[["chrom", "startpos", "endpos"]].to_csv(
                        os.path.join(
                            sample_output_dir, "targets_exceeding_threshold.bed"
                        ),
                        sep="\t",
                        header=None,
                        index=None,
                    )

                    logger.info(
                        f"Targets exceeding threshold saved: {len(run_list)} regions"
                    )
                    
                    # Log the results for debugging (no storage to avoid memory leaks)
                    logger.info(f"Generated BED file with {len(run_list)} regions exceeding threshold")
                else:
                    logger.info("No regions exceed coverage threshold")
                    
                    # Create empty BED file to avoid errors
                    empty_bed_file = os.path.join(
                        sample_output_dir, "targets_exceeding_threshold.bed"
                    )
                    with open(empty_bed_file, "w") as f:
                        # Write minimal BED header
                        pass
                    
                    logger.info("Created empty targets_exceeding_threshold.bed file")
                    
                    # Log the empty file creation (no storage to avoid memory leaks)
                    logger.info("Created empty BED file - no regions exceed threshold")
                
                # Log reference genome status for SNP calling
                if self.reference:
                    logger.info(f"Reference genome found: {self.reference} - SNP calling will be available")
                else:
                    logger.info("No reference genome provided - SNP calling will not be available")

                target_result.processing_steps.append("threshold_analysis_complete")

                # Step 13: Perform target analysis
                logger.info("Performing target analysis...")
                analysis_results = self._perform_target_analysis(
                    file_path, temp_dir, metadata, newcovdf, bedcovdf
                )
                target_result.analysis_results = analysis_results
                target_result.processing_steps.append("analysis_completed")

                # Step 14: Generate output files
                logger.info("Generating output files...")
                output_files = self._generate_output_files(
                    sample_output_dir,
                    analysis_counter,
                    analysis_results,
                    newcovdf,
                    bedcovdf,
                    logger,
                    metadata,
                )
                # No JSON files created - data is accumulated in CSV files
                target_result.target_data_path = None
                target_result.target_plot_path = None
                target_result.processing_steps.append("output_generated")

                # Step 15: Update analysis counter
                analysis_counter += 1
                self._save_analysis_counter(
                    sample_id, analysis_counter, self.work_dir, logger
                )
                target_result.processing_steps.append("counter_updated")

                # Step 16: Force garbage collection
                gc.collect()

                target_result.processing_steps.append("target_analysis_complete")
                logger.info(f"Target analysis complete for {sample_id}")
                logger.info(f"   Results: {len(analysis_results)} items")
                logger.info(f"   Output files: {len(output_files)} files")
                logger.info(
                    f"   Targets exceeding threshold: {len(run_list) if 'run_list' in locals() else 0}"
                )

                return target_result

        except Exception as e:
            error_details = f"Error in target analysis for {sample_id}: {str(e)}"
            logger.error(error_details)
            target_result.error_message = error_details
            target_result.processing_steps.append("analysis_failed")
            return target_result

    def _perform_target_analysis(
        self,
        file_path: str,
        temp_dir: str,
        metadata: Dict[str, Any],
        newcovdf: pd.DataFrame,
        bedcovdf: pd.DataFrame,
    ) -> Dict[str, Any]:
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
        logger = logging.getLogger("robin.target")

        # This is where the actual target analysis logic would go
        # For now, we'll create a placeholder implementation that uses the coverage data

        results = {
            "analysis_type": self.config.get("analysis_type", "standard"),
            "input_file": os.path.basename(file_path),
            "sample_id": metadata.get("sample_id", "unknown"),
            "timestamp": time.time(),
            "targets_found": 0,
            "targets_analyzed": [],
            "coverage_summary": {
                "genome_regions": len(newcovdf) if newcovdf is not None else 0,
                "target_regions": len(bedcovdf) if bedcovdf is not None else 0,
                "total_bases_covered": 0,
                "average_coverage": 0.0,
            },
            "summary_stats": {
                "total_targets": 0,
                "significant_targets": 0,
                "average_score": 0.0,
            },
        }

        # Calculate coverage statistics
        if newcovdf is not None and not newcovdf.empty:
            results["coverage_summary"]["total_bases_covered"] = newcovdf.get(
                "numreads", pd.Series([0])
            ).sum()
            results["coverage_summary"]["average_coverage"] = newcovdf.get(
                "numreads", pd.Series([0])
            ).mean()

        # Analyze target regions
        if bedcovdf is not None and not bedcovdf.empty:
            results["targets_found"] = len(bedcovdf)
            results["summary_stats"]["total_targets"] = len(bedcovdf)

            # Calculate significance based on coverage
            significant_threshold = self.config.get("threshold", 0.5)
            significant_targets = bedcovdf[
                bedcovdf["bases"] > significant_threshold * bedcovdf["bases"].max()
            ]
            results["summary_stats"]["significant_targets"] = len(significant_targets)

            # Add target analysis results
            for idx, row in bedcovdf.iterrows():
                target_info = {
                    "target_id": row.get("name", f"target_{idx+1:03d}"),
                    "position": f"{row['chrom']}:{row['startpos']}-{row['endpos']}",
                    "coverage": row["bases"],
                    "significance": (
                        "high"
                        if row["bases"]
                        > significant_threshold * bedcovdf["bases"].max()
                        else "medium"
                    ),
                }
                results["targets_analyzed"].append(target_info)

            # Calculate average score
            if results["targets_analyzed"]:
                scores = [t["coverage"] for t in results["targets_analyzed"]]
                results["summary_stats"]["average_score"] = sum(scores) / len(scores)

        logger.info(f"Analysis completed: {results['targets_found']} targets found")
        logger.info(f"Coverage summary: {results['coverage_summary']}")
        return results

    def _generate_output_files(
        self,
        sample_dir: str,
        analysis_counter: int,
        analysis_results: Dict[str, Any],
        newcovdf: pd.DataFrame,
        bedcovdf: pd.DataFrame,
        logger,
        metadata: Dict[str, Any],
    ) -> Dict[str, str]:
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

            logger.info(
                "No additional output files needed - data accumulated in CSV files"
            )
            return output_files

        except Exception as e:
            logger.error(f"Error in output file generation: {e}")
            return {}

    def _load_analysis_counter(self, sample_id: str, work_dir: str, logger) -> int:
        """Load the analysis counter for a sample from disk"""
        counter_file = os.path.join(work_dir, sample_id, "target_analysis_counter.txt")
        if os.path.exists(counter_file):
            try:
                with open(counter_file, "r") as f:
                    return int(f.read().strip())
            except (ValueError, IOError) as e:
                logger.warning(f"Error loading counter for {sample_id}: {e}")
        return 0

    def _save_analysis_counter(
        self, sample_id: str, counter: int, work_dir: str, logger
    ) -> None:
        """Save the analysis counter for a sample to disk"""
        counter_file = os.path.join(work_dir, sample_id, "target_analysis_counter.txt")
        try:
            # Create sample directory if it doesn't exist
            os.makedirs(os.path.dirname(counter_file), exist_ok=True)
            with open(counter_file, "w") as f:
                f.write(str(counter))
        except IOError as e:
            logger.error(f"Error saving counter for {sample_id}: {e}")

    def _load_existing_coverage_data(
        self, sample_dir: str, filename: str, logger
    ) -> Optional[pd.DataFrame]:
        """Load existing coverage data from CSV file if it exists"""
        file_path = os.path.join(sample_dir, filename)
        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path)
                logger.info(
                    f"Loaded existing coverage data from {filename}: {df.shape}"
                )

                # For bed_coverage_main.csv, extract only the raw data columns needed for merging
                if (
                    filename == "bed_coverage_main.csv"
                    and "length" in df.columns
                    and "coverage" in df.columns
                ):
                    # Extract only the raw data columns: chrom,startpos,endpos,name,bases
                    df = df[["chrom", "startpos", "endpos", "name", "bases"]]
                    logger.info(
                        f"Extracted raw data columns from {filename}: {df.shape}"
                    )

                return df
            except Exception as e:
                logger.warning(
                    f"Error loading existing coverage data from {filename}: {e}"
                )
                return None
        return None

    def _load_existing_coverage_over_time(
        self, sample_dir: str, logger
    ) -> Optional[np.ndarray]:
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


def process_single_file(
    file_path: str, metadata: Dict[str, Any], work_dir: str, logger, reference: Optional[str] = None
) -> Dict[str, Any]:
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
    sample_id = metadata.get("sample_id", "unknown")
    logger.info(f"🎯 Starting target analysis for sample: {sample_id}")

    # Create sample-specific output directory
    sample_output_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)
    logger.info(f"Created sample output directory: {sample_output_dir}")

    analysis_result = {
        "sample_id": sample_id,
        "file_path": file_path,
        "analysis_timestamp": time.time(),
        "target_data_path": None,
        "target_plot_path": None,
        "analysis_results": {},
        "processing_steps": [],
        "error_message": None,
        "coverage_data": {},
        "target_bam_path": None,
        "coverage_over_time": None,
    }

    try:
        # Initialize target analysis
        target_analysis = TargetAnalysis(work_dir=work_dir)
        
        # Set reference genome if provided
        if reference:
            target_analysis.reference = reference
            logger.info(f"Using reference genome: {reference}")

        # Process the file
        target_metadata = target_analysis.process_file(file_path, metadata)

        # Update analysis result with target metadata
        analysis_result.update(
            {
                "target_data_path": target_metadata.target_data_path,
                "target_plot_path": target_metadata.target_plot_path,
                "analysis_results": target_metadata.analysis_results,
                "processing_steps": target_metadata.processing_steps,
                "error_message": target_metadata.error_message,
                "coverage_data": target_metadata.coverage_data,
                "target_bam_path": target_metadata.target_bam_path,
                "coverage_over_time": target_metadata.coverage_over_time,
            }
        )

        if target_metadata.error_message:
            logger.error(f"Target analysis failed: {target_metadata.error_message}")
        else:
            logger.info(f"Target analysis complete for {sample_id}")
            logger.info(f"   Results: {len(target_metadata.analysis_results)} items")
            logger.info(
                f"   Processing steps: {', '.join(target_metadata.processing_steps)}"
            )
            logger.info(
                f"   Coverage data: genome={target_metadata.coverage_data.get('genome_coverage_shape', (0,0))}, targets={target_metadata.coverage_data.get('target_coverage_shape', (0,0))}"
            )
            if target_metadata.target_bam_path:
                logger.info(f"   Target BAM: {target_metadata.target_bam_path}")

        return analysis_result

    except Exception as e:
        error_details = f"Error in target analysis for {sample_id}: {str(e)}"
        logger.error(error_details)
        analysis_result["error_message"] = error_details
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def target_handler(job, work_dir=None, reference=None):
    """
    Handler function for target analysis jobs.
    This function processes files for target-specific analysis.

    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to file directory)
        reference: Optional path to reference genome for SNP calling
    """
    # Get job-specific logger
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)

    try:
        file_path = job.context.filepath

        logger.info(f"Starting target analysis for: {os.path.basename(file_path)}")
        
        # Debug: Log reference genome status
        if reference:
            logger.info(f"Reference genome provided to target_handler: {reference}")
        else:
            logger.info("No reference genome provided to target_handler")
            
        # Also check job metadata for reference
        job_reference = job.context.metadata.get("reference")
        if job_reference:
            logger.info(f"Reference genome found in job metadata: {job_reference}")
            # Use job metadata reference if not provided directly
            if not reference:
                reference = job_reference
                logger.info(f"Using reference from job metadata: {reference}")
        else:
            logger.info("No reference genome found in job metadata")

        # Get metadata from preprocessing
        file_metadata = job.context.metadata.get("bam_metadata", {})

        # Determine work directory
        if work_dir is None:
            # Default to file directory
            work_dir = os.path.dirname(file_path)
        else:
            # Use specified work directory, create if it doesn't exist
            os.makedirs(work_dir, exist_ok=True)

        # Process the file
        result = process_single_file(file_path, file_metadata, work_dir, logger, reference)

        # Store results in job context
        job.context.add_metadata("target_analysis", result)

        if result.get("error_message"):
            job.context.add_error("target_analysis", result["error_message"])
            logger.error(f"Target analysis failed: {result['error_message']}")
        else:
            job.context.add_result(
                "target_analysis",
                {
                    "status": "success",
                    "sample_id": result.get("sample_id", "unknown"),
                    "analysis_time": result.get("analysis_timestamp", 0),
                    "targets_found": result.get("analysis_results", {}).get(
                        "targets_found", 0
                    ),
                    "processing_steps": result.get("processing_steps", []),
                    "target_data_path": result.get("target_data_path", ""),
                    "target_plot_path": result.get("target_plot_path", ""),
                    "target_bam_path": result.get("target_bam_path", ""),
                    "coverage_summary": result.get("analysis_results", {}).get(
                        "coverage_summary", {}
                    ),
                },
            )
            logger.info(f"Target analysis complete for {os.path.basename(file_path)}")
            logger.info(f"Sample ID: {result.get('sample_id', 'unknown')}")
            logger.info(
                f"Targets found: {result.get('analysis_results', {}).get('targets_found', 0)}"
            )

    except Exception as e:
        error_details = f"Error in target analysis for {job.context.filepath}: {str(e)}"
        job.context.add_error("target_analysis", error_details)
        logger.error(error_details)


def ensure_sorted_igv_bam(sample_dir: str, threads: int = 4, force_regenerate: bool = False) -> str:
    """
    Ensure an IGV-ready coordinate-sorted and indexed BAM exists for a sample.

    Preferred output: ``sample_dir/igv/igv_ready.bam`` (+ .bai).
    Backward compat: if an existing Clair3-sorted BAM is present, reuse it.
    If neither exists, sort and index ``sample_dir/target.bam`` into
    ``sample_dir/igv/igv_ready.bam``.

    Args:
        sample_dir: Path to the sample output directory.
        threads: Number of threads to pass to pysam.sort.
        force_regenerate: If True, force recreation even if output exists.
        
    Returns:
        The path to the sorted BAM file (empty string on failure).
    """
    logger = logging.getLogger("robin.target")
    try:
        sample_dir = os.path.abspath(sample_dir)
        igv_dir = os.path.join(sample_dir, "igv")
        clair_dir = os.path.join(sample_dir, "clair3")
        os.makedirs(igv_dir, exist_ok=True)

        # Preferred IGV-ready path
        out_bam = os.path.join(igv_dir, "igv_ready.bam")
        out_bai = f"{out_bam}.bai"
        
        # Check if files exist and force_regenerate is not set
        if not force_regenerate and os.path.exists(out_bam) and os.path.exists(out_bai):
            logger.info(f"IGV BAM already present: {out_bam}")
            return out_bam

        # Backward compatibility: reuse Clair3-sorted BAM if present
        legacy_bams = [
            os.path.join(clair_dir, "sorted_targets_exceeding_rerun.bam"),
            os.path.join(clair_dir, "sorted_targets_exceeding.bam"),
        ]
        for legacy in legacy_bams:
            if os.path.exists(legacy) and os.path.exists(f"{legacy}.bai"):
                logger.info(f"Reusing legacy Clair3 BAM for IGV: {legacy}")
                # Optionally copy to IGV location to standardize
                try:
                    import shutil
                    shutil.copy2(legacy, out_bam)
                    shutil.copy2(f"{legacy}.bai", out_bai)
                    return out_bam
                except Exception:
                    # Fallback: just return legacy path
                    return legacy

        in_bam = os.path.join(sample_dir, "target.bam")
        if not os.path.exists(in_bam):
            logger.warning(f"target.bam not found in {sample_dir}; cannot build IGV BAM")
            return ""

        logger.info(f"Sorting BAM for IGV: {in_bam} -> {out_bam}")
        pysam.sort(f"-@{threads}", "-o", out_bam, in_bam)
        pysam.index(out_bam, out_bai)
        logger.info(f"IGV BAM ready: {out_bam}")
        return out_bam
    except Exception as e:
        logger.error(f"Failed to create IGV BAM in {sample_dir}: {e}")
        return ""


def igv_bam_handler(job, work_dir: Optional[str] = None) -> None:
    """
    Workflow handler to build an IGV-ready BAM for a sample.

    Expects a per-sample directory containing `target.bam`. Outputs to
    `sample_dir/igv/igv_ready.bam` (and .bai). On success, records the output
    path in the job context under results['igv_bam'].
    """
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    try:
        # Determine sample dir: prefer provided work_dir + sample_id
        sid = job.context.get_sample_id() if hasattr(job.context, "get_sample_id") else "unknown"
        base = work_dir or job.context.metadata.get("work_dir")
        sample_dir = None
        
        logger.info(f"IGV BAM job: sample_id={sid}, work_dir={work_dir}, base={base}")
        
        if base and sid and sid != "unknown":
            sample_dir = os.path.join(base, sid)
            logger.info(f"Using work_dir + sample_id path: {sample_dir}")
        else:
            # Fallback: directory of the current file
            sample_dir = os.path.dirname(job.context.filepath)
            logger.info(f"Using fallback path: {sample_dir}")
        
        # Additional fallback: if the filepath is already a sample directory, use it directly
        if os.path.isdir(job.context.filepath):
            sample_dir = job.context.filepath
            logger.info(f"Using filepath directly as sample directory: {sample_dir}")

        if not os.path.isdir(sample_dir):
            raise RuntimeError(f"Sample directory not found: {sample_dir}")
        
        logger.info(f"Final sample directory: {sample_dir}")
        logger.info(f"Looking for target.bam in: {os.path.join(sample_dir, 'target.bam')}")
        
        # Check if force_regenerate is requested
        force_regenerate = job.context.metadata.get("force_regenerate", False)
        if force_regenerate:
            logger.info("Force regenerate requested - will recreate IGV BAM even if it exists")
        
        out = ensure_sorted_igv_bam(sample_dir, force_regenerate=force_regenerate)
        if not out:
            raise RuntimeError("Failed to create IGV-ready BAM")

        # Record result
        try:
            job.context.add_result("igv_bam", {"bam": out, "bai": f"{out}.bai"})
        except Exception:
            pass
        logger.info(f"IGV-ready BAM created: {out}")
    except Exception as e:
        job.context.add_error("igv_bam", str(e))
        logger.error(f"igv_bam handler failed: {e}")


def run_snp_analysis(sample_dir: str, threads: int = 4, force_regenerate: bool = False, reference: Optional[str] = None) -> str:
    """
    Run SNP analysis using Clair3 for a sample.
    
    This function runs the complete SNP calling pipeline including:
    - BAM sorting and preparation
    - Clair3 variant calling
    - snpEff annotation
    - SnpSift annotation with ClinVar data
    
    Args:
        sample_dir: Path to the sample output directory
        threads: Number of threads to use for processing
        force_regenerate: If True, force recreation even if output exists
        
    Returns:
        The path to the output directory containing SNP results (empty string on failure)
    """
    logger = logging.getLogger("littlejohn.target")
    
    # CRITICAL: Immediate logging to see if we even get here
    logger.info("🚀 ENTERING run_snp_analysis function")
    logger.info(f"Parameters: sample_dir={sample_dir}, threads={threads}, force_regenerate={force_regenerate}, reference={reference}")
    
    # Also log to stderr to make sure we see it
    import sys
    
    try:
        logger.info("🔄 STEP 1: Setting up sample directory")
        sample_dir = os.path.abspath(sample_dir)
        logger.info(f"Absolute sample_dir: {sample_dir}")
        
        clair_dir = os.path.join(sample_dir, "clair3")
        logger.info(f"Clair3 output directory: {clair_dir}")
        
        os.makedirs(clair_dir, exist_ok=True)
        logger.info(f"✅ Clair3 directory created/verified: {clair_dir}")
        
        logger.info("🔄 STEP 2: Checking for existing SNP analysis")
        
        # Check if SNP analysis already exists and force_regenerate is not set
        snp_output_files = [
            os.path.join(clair_dir, "snpsift_output.vcf"),
            os.path.join(clair_dir, "snpsift_indel_output.vcf")
        ]
        
        if not force_regenerate and all(os.path.exists(f) for f in snp_output_files):
            logger.info(f"SNP analysis already present in {clair_dir}")
            return clair_dir
        
        logger.info("🔄 STEP 3: Checking for required input files")
        target_bam = os.path.join(sample_dir, "target.bam")
        targets_bed = os.path.join(sample_dir, "targets_exceeding_threshold.bed")
        
        logger.info(f"🔍 CHECKING INPUT FILES FOR CLAIR3")
        logger.info(f"  Sample directory: {sample_dir}")
        logger.info(f"  Target BAM path: {target_bam}")
        logger.info(f"  Targets BED path: {targets_bed}")
        logger.info(f"  Reference path: {reference}")
        
        
        if not os.path.exists(target_bam):
            logger.error(f"❌ target.bam not found: {target_bam}")
            return ""
            
        if not os.path.exists(targets_bed):
            logger.error(f"❌ targets_exceeding_threshold.bed not found: {targets_bed}")
            return ""
        
        
        # Check for reference genome
        if not reference:
            # Try to find reference from common locations
            possible_refs = [
                os.path.join(sample_dir, "reference.fasta"),
                os.path.join(sample_dir, "reference.fa"),
                os.path.join(os.path.dirname(sample_dir), "reference.fasta"),
                os.path.join(os.path.dirname(sample_dir), "reference.fa"),
                "/usr/local/share/reference.fasta",
                "/usr/local/share/reference.fa"
            ]
            
            for ref_path in possible_refs:
                if os.path.exists(ref_path):
                    reference = ref_path
                    break
        
        if not reference:
            logger.error("Reference genome not found. SNP calling requires a reference genome.")
            return ""
        
        logger.info(f"Starting SNP analysis for sample: {os.path.basename(sample_dir)}")
        logger.info(f"Reference: {reference}")
        logger.info(f"Target BAM: {target_bam}")
        logger.info(f"Targets BED: {targets_bed}")
        
        logger.info("🔄 STEP 4: Sorting target BAM for Clair3")
        sorted_bam = os.path.join(clair_dir, "sorted_targets_exceeding.bam")
        logger.info(f"Sorting BAM: {target_bam} -> {sorted_bam}")
        
        try:
            # Check if sorted BAM already exists and is valid
            if os.path.exists(sorted_bam) and os.path.exists(sorted_bam + ".bai"):
                logger.info(f"✅ Sorted BAM already exists: {sorted_bam}")
                # Verify the file is readable
                try:
                    test_bam = pysam.AlignmentFile(sorted_bam, "rb")
                    test_bam.close()
                    logger.info("✅ Sorted BAM file is valid and readable")
                except Exception as e:
                    logger.warning(f"⚠️ Existing sorted BAM appears corrupted, regenerating: {e}")
                    os.remove(sorted_bam)
                    if os.path.exists(sorted_bam + ".bai"):
                        os.remove(sorted_bam + ".bai")
                    raise Exception("Corrupted BAM file")
            else:
                logger.info("🔄 Creating new sorted BAM file...")
            
            # Sort the BAM file
            logger.info(f"🔄 Sorting BAM with {threads} threads...")
            pysam.sort(f"-@{threads}", "-o", sorted_bam, target_bam)
            
            # Verify the sorted file was created
            if not os.path.exists(sorted_bam):
                raise Exception(f"Sorted BAM file was not created: {sorted_bam}")
            
            # Create index
            logger.info("🔄 Creating BAM index...")
            pysam.index(sorted_bam)
            
            # Verify index was created
            if not os.path.exists(sorted_bam + ".bai"):
                raise Exception(f"BAM index was not created: {sorted_bam}.bai")
            
            # Final verification
            file_size = os.path.getsize(sorted_bam)
            logger.info(f"✅ BAM sorting completed successfully: {sorted_bam} ({file_size / (1024**2):.1f} MB)")
            
        except Exception as e:
            logger.error(f"❌ BAM sorting failed: {e}")
            # Clean up any partial files
            for file_path in [sorted_bam, sorted_bam + ".bai"]:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    logger.info(f"Cleaned up partial file: {file_path}")
            return ""
        
        logger.info("🔄 STEP 5: Running Clair3 pipeline")
        logger.info("Running Clair3 pipeline...")
        
        # Run Clair3 pipeline using Docker
        logger.info("🔍 CHECKING DOCKER AVAILABILITY")
        if docker is None:
            logger.error("❌ Docker module not available. Cannot run Clair3 pipeline.")
            return ""
        else:
            logger.info("✅ Docker module is available")
            
        try:
            client = docker.from_env()
            
            # Test Docker connectivity
            try:
                client.ping()
                logger.info("✅ Docker daemon is accessible")
            except Exception as docker_error:
                logger.error(f"❌ Docker daemon not accessible: {docker_error}")
                return ""
            
            # Check if Clair3 image exists
            try:
                client.images.get("hkubal/clairs-to:latest")
                logger.info("✅ Clair3 image found: hkubal/clairs-to:latest")
            except Exception as image_error:
                logger.warning(f"⚠️ Clair3 image not found, attempting to pull: hkubal/clairs-to:latest")
                try:
                    client.images.pull("hkubal/clairs-to:latest")
                    logger.info("✅ Clair3 image pulled successfully")
                except Exception as pull_error:
                    logger.error(f"❌ Failed to pull Clair3 image: {pull_error}")
                    return ""
            
            # Create container-specific paths using the original working pattern
            container_bamfile = f"/data/output/{os.path.basename(sorted_bam)}"  # BAM goes to output directory
            container_bedfile = f"/data/bed/{os.path.basename(targets_bed)}"   # BED goes to bed directory
            container_reference = f"/data/reference/{os.path.basename(reference)}"  # Reference goes to reference directory
            container_output = "/data/output"
                    
            # Clean up old debugging - no longer needed with the fixed approach
            
            # Use the original working pattern: bind each directory to its own container path
            volume_bindings = {
                os.path.abspath(os.path.dirname(sorted_bam)): {"bind": "/data/bam", "mode": "ro"},      # BAM directory (clair3)
                os.path.abspath(os.path.dirname(targets_bed)): {"bind": "/data/bed", "mode": "ro"},    # BED directory (Sample_103)
                os.path.abspath(os.path.dirname(reference)): {"bind": "/data/reference", "mode": "ro"}, # Reference directory
                os.path.abspath(clair_dir): {"bind": "/data/output", "mode": "rw"},                     # Output directory (clair3)
            }
            
            # Debug: Let's verify the BAM binding was added
            bam_dir = os.path.dirname(sorted_bam)
            
            logger.info("🔍 Verifying volume bindings...")
            for source_path, binding_info in volume_bindings.items():
                if not os.path.exists(source_path):
                    logger.error(f"❌ Volume source path does not exist: {source_path}")
                    return ""
                if not os.access(source_path, os.R_OK):
                    logger.error(f"❌ Volume source path not readable: {source_path}")
                    return ""
                logger.info(f"✅ Volume binding verified: {source_path} -> {binding_info['bind']}")
            
            # Verify the specific files exist in their directories
            logger.info("🔍 Verifying input files in volume directories...")
            sorted_bam_dir = os.path.dirname(sorted_bam)
            targets_bed_dir = os.path.dirname(targets_bed)
            reference_dir = os.path.dirname(reference)
            
            if not os.path.exists(sorted_bam):
                logger.error(f"❌ Sorted BAM file not found in directory: {sorted_bam}")
                return ""
            if not os.path.exists(targets_bed):
                logger.error(f"❌ Targets BED file not found in directory: {targets_bed}")
                return ""
            if not os.path.exists(reference):
                logger.error(f"❌ Reference file not found in directory: {reference}")
                return ""
            
            logger.info(f"✅ All input files verified in their directories")
            
            host_config = client.api.create_host_config(binds=volume_bindings)
            
            # Function to split BED file into manageable regions
            def split_bed_into_regions(bed_file, max_region_size=150000000):
                """Split BED file into efficient regions up to max_region_size bases"""
                regions = []
                try:
                    # Read all BED entries and sort them
                    bed_entries = []
                    with open(bed_file, 'r') as f:
                        for line_num, line in enumerate(f, 1):
                            line = line.strip()
                            if not line or line.startswith('#'):
                                continue
                            
                            parts = line.split('\t')
                            if len(parts) < 3:
                                logger.warning(f"Skipping invalid BED line {line_num}: {line}")
                                continue
                            
                            chrom = parts[0]
                            start = int(parts[1])
                            end = int(parts[2])
                            bed_entries.append((chrom, start, end))
                    
                    # Sort by chromosome and start position
                    bed_entries.sort(key=lambda x: (x[0], x[1]))
                    
                    # Group into efficient regions
                    current_chrom = None
                    current_start = None
                    current_end = None
                    
                    for chrom, start, end in bed_entries:
                        if chrom != current_chrom:
                            # New chromosome, start new region
                            if current_chrom is not None:
                                regions.append(f"{current_chrom}:{current_start+1}-{current_end}")
                            current_chrom = chrom
                            current_start = start
                            current_end = end
                        elif end - current_start <= max_region_size:
                            # Can extend current region
                            current_end = end
                        else:
                            # Current region would exceed max size, start new one
                            regions.append(f"{current_chrom}:{current_start+1}-{current_end}")
                            current_start = start
                            current_end = end
                    
                    # Add final region
                    if current_chrom is not None:
                        regions.append(f"{current_chrom}:{current_start+1}-{current_end}")
                    
                    logger.info(f"Combined BED entries into {len(regions)} efficient regions (max size: {max_region_size:,} bases)")
                    return regions
                except Exception as e:
                    logger.error(f"Error splitting BED file: {e}")
                    return []
            
            # Split BED file into regions
            logger.info("Splitting BED file into efficient regions...")
            regions = split_bed_into_regions(targets_bed, max_region_size=150000000)
            
            if not regions:
                logger.error("Failed to split BED file into regions")
                return ""
            
            logger.info(f"Processing {len(regions)} efficient regions: {regions[:5]}{'...' if len(regions) > 5 else ''}")
            
            # Log region sizes for efficiency analysis
            total_bases = 0
            for i, region in enumerate(regions):
                try:
                    chrom_part, pos_part = region.split(':')
                    start_end = pos_part.split('-')
                    start = int(start_end[0]) - 1  # Convert back to 0-based
                    end = int(start_end[1])
                    region_size = end - start
                    total_bases += region_size
                    logger.info(f"Region {i+1}: {region} (size: {region_size:,} bases)")
                except Exception as e:
                    logger.warning(f"Could not parse region size for {region}: {e}")
            
            logger.info(f"Total bases to process: {total_bases:,} across {len(regions)} regions")
            logger.info(f"Average region size: {total_bases // len(regions):,} bases")
            
            # Process each region separately
            all_snv_vcfs = []
            all_indel_vcfs = []
            
            for i, region in enumerate(regions):
                logger.info(f"Processing region {i+1}/{len(regions)}: {region}")
                print(f"Processing region {i+1}/{len(regions)}: {region} - PRINT STATEMENT")
                # Create region-specific output directory
                region_output = f"{container_output}/region_{i+1}"
                
                # Build Clair3 command for this region
                command = (
                    f"/opt/bin/run_clairs_to "
                    f"--tumor_bam_fn {container_bamfile} "
                    f"--ref_fn {container_reference} "
                    f"--threads 4 "
                    f"--remove_intermediate_dir "
                    f"--platform ont_r10_guppy_hac_5khz "
                    f"--output_dir {region_output} "
                    f"-b {container_bedfile} "
                    f"--chunk_size 5000000 "
                    f"--disable_intermediate_phasing "
                    f"--region {region}"
                )
            
            logger.info(f"Container input BAM: {container_bamfile}")
            logger.info(f"Container input BED: {container_bedfile}")
            logger.info(f"Container output dir: {container_output}")
            logger.info(f"Container reference: {container_reference}")
            logger.info(f"Volume bindings: {volume_bindings}")
            
            # Function to merge VCF files
            def merge_vcf_files(vcf_files, output_file, variant_type):
                """Merge multiple VCF files into a single output file"""
                try:
                    if not vcf_files:
                        logger.error(f"No {variant_type} VCF files to merge")
                        return False
                    
                    logger.info(f"Merging {len(vcf_files)} {variant_type} VCF files into {output_file}")
                    
                    # Read headers from first file
                    import gzip
                    headers = []
                    first_vcf = vcf_files[0]
                    
                    with gzip.open(first_vcf, 'rt') as f:
                        for line in f:
                            if line.startswith('#'):
                                headers.append(line)
                            else:
                                break
                    
                    # Write merged output
                    with gzip.open(output_file, 'wt') as out_f:
                        # Write headers
                        for header in headers:
                            out_f.write(header)
                        
                        # Write variants from each file
                        for vcf_file in vcf_files:
                            if os.path.exists(vcf_file):
                                with gzip.open(vcf_file, 'rt') as in_f:
                                    for line in in_f:
                                        if not line.startswith('#'):
                                            out_f.write(line)
                    
                    logger.info(f"Successfully merged {variant_type} VCF files into {output_file}")
                    return True
                    
                except Exception as e:
                    logger.error(f"Error merging {variant_type} VCF files: {e}")
                    return False
            
            # Log system resources before starting
            try:
                import psutil
                memory = psutil.virtual_memory()
                logger.info(f"System memory: {memory.available / (1024**3):.1f} GB available out of {memory.total / (1024**3):.1f} GB")
                logger.info(f"System CPU cores: {psutil.cpu_count()}")
            except ImportError:
                logger.info("psutil not available, skipping resource logging")
            
            
            
            
            # Process each region
            for i, region in enumerate(regions):
                print(f"Processing region {i+1}/{len(regions)}: {region} - PRINT STATEMENT")
                logger.info(f"Processing region {i+1}/{len(regions)}: {region}")
                
                # Create region-specific output directory
                region_output = f"{container_output}/region_{i+1}"
                
                # Build Clair3 command for this region
                command = (
                    f"/opt/bin/run_clairs_to "
                    f"--tumor_bam_fn {container_bamfile} "
                    f"--ref_fn {container_reference} "
                    f"--threads 4 "
                    f"--remove_intermediate_dir "
                    f"--platform ont_r10_guppy_hac_5khz "
                    f"--output_dir {region_output} "
                    f"-b {container_bedfile} "
                    f"--chunk_size 5000000 "
                    f"--disable_intermediate_phasing "
                    f"--region {region}"
                )
                
                logger.info(f"Running Clair3 command for region {i+1}: {command}")
                
                # Create and start container for this region
                logger.info(f"Creating Clair3 container for region {i+1}...")
                container = client.api.create_container(
                    image="hkubal/clairs-to:latest",
                    command=command,
                    volumes=list(volume_bindings.keys()),
                    host_config=host_config,
                )
                logger.info(f"Container created with ID: {container.get('Id')}")
                
                # Start container
                client.api.start(container=container.get("Id"))
                
                # Stream container logs
                for line in client.api.logs(container=container.get("Id"), stream=True, follow=True):
                    logger.info(f"Clair3 Region {i+1}: {line.decode().strip()}")
                
                # Wait for completion
                result = client.api.wait(container=container.get("Id"))
                if result["StatusCode"] != 0:
                    # Get detailed error logs before removing container
                    error_logs = client.api.logs(container=container.get("Id"), stderr=True, stdout=False)
                    error_details = error_logs.decode().strip()
                    if error_details:
                        logger.error(f"Clair3 region {i+1} error logs: {error_details}")
                    
                    # Also get stdout logs for context
                    stdout_logs = client.api.logs(container=container.get("Id"), stdout=True, stderr=False)
                    stdout_details = stdout_logs.decode().strip()
                    if stdout_details:
                        logger.info(f"Clair3 region {i+1} stdout logs: {stdout_details}")
                    
                    logger.warning(f"Clair3 region {i+1} failed with status code {result['StatusCode']}, skipping this region")
                    continue
                
                # Clean up container
                client.api.remove_container(container=container.get("Id"))
                
                # Check if output files were created for this region
                region_snv = f"{clair_dir}/region_{i+1}/snv.vcf.gz"
                region_indel = f"{clair_dir}/region_{i+1}/indel.vcf.gz"
                
                if os.path.exists(region_snv):
                    all_snv_vcfs.append(region_snv)
                    logger.info(f"Region {i+1} SNV output created: {region_snv}")
                else:
                    logger.warning(f"Region {i+1} SNV output not found: {region_snv}")
                
                if os.path.exists(region_indel):
                    all_indel_vcfs.append(region_indel)
                    logger.info(f"Region {i+1} INDEL output created: {region_indel}")
                else:
                    logger.warning(f"Region {i+1} INDEL output not found: {region_indel}")
            
            # Merge all region outputs
            logger.info(f"Merging {len(all_snv_vcfs)} SNV VCF files and {len(all_indel_vcfs)} INDEL VCF files...")
            
            if all_snv_vcfs:
                # Merge SNV VCFs
                merged_snv = f"{clair_dir}/merged_snv.vcf.gz"
                merge_vcf_files(all_snv_vcfs, merged_snv, "SNV")
                shutil.copy2(merged_snv, f"{clair_dir}/output_done.vcf.gz")
                logger.info(f"Merged SNV VCFs into: {clair_dir}/output_done.vcf.gz")
            else:
                logger.error("No SNV VCF files were successfully created")
                return ""
            
            if all_indel_vcfs:
                # Merge INDEL VCFs
                merged_indel = f"{clair_dir}/merged_indel.vcf.gz"
                merge_vcf_files(all_indel_vcfs, merged_indel, "INDEL")
                shutil.copy2(merged_indel, f"{clair_dir}/output_indel_done.vcf.gz")
                logger.info(f"Merged INDEL VCFs into: {clair_dir}/output_indel_done.vcf.gz")
            else:
                logger.warning("No INDEL VCF files were successfully created")
            
            logger.info("Clair3 pipeline completed successfully for all regions")
            
        except ImportError:
            logger.error("Docker not available. Cannot run Clair3 pipeline.")
            return ""
        except Exception as e:
            logger.error(f"Error running Clair3 pipeline: {e}")
            return ""
        
        logger.info("🔄 STEP 6: Running annotation pipeline")
        logger.info("Running annotation pipeline...")
        
        try:
            # SNP annotation with snpEff
            logger.info(f"Starting snpEff annotation for SNPs...")
            logger.info(f"Input file: {clair_dir}/output_done.vcf.gz")
            logger.info(f"Output file: {clair_dir}/snpeff_output.vcf")
            
            # Define output file paths
            snpeff_out = f"{clair_dir}/snpeff_output.vcf"
            snpsift_out = f"{clair_dir}/snpsift_output.vcf"
            snpeff_indel_out = f"{clair_dir}/snpeff_indel_output.vcf"
            snpsift_indel_out = f"{clair_dir}/snpsift_indel_output.vcf"
            
            # Check if input file exists and has content
            if os.path.exists(f"{clair_dir}/output_done.vcf.gz"):
                file_size = os.path.getsize(f"{clair_dir}/output_done.vcf.gz")
                logger.info(f"Input VCF file exists, size: {file_size} bytes")
                
                # Check first few lines to verify it's a valid VCF
                try:
                    import gzip
                    with gzip.open(f"{clair_dir}/output_done.vcf.gz", 'rt') as f:
                        first_lines = [next(f) for _ in range(5)]
                    logger.info(f"First 5 lines of input VCF:")
                    for i, line in enumerate(first_lines):
                        logger.info(f"  Line {i+1}: {line.strip()}")
                except Exception as e:
                    logger.warning(f"Could not read input VCF file: {e}")
            else:
                logger.error(f"Input VCF file does not exist: {clair_dir}/output_done.vcf.gz")
            
            snpeff_cmd = ["snpEff", "-q", "hg38", f"{clair_dir}/output_done.vcf.gz"]
            logger.info(f"Running snpEff command: {' '.join(snpeff_cmd)}")
            
            # Check if snpEff is available
            try:
                snpeff_version = subprocess.run(["snpEff", "-version"], capture_output=True, text=True)
                if snpeff_version.returncode == 0:
                    logger.info(f"snpEff version: {snpeff_version.stdout.strip()}")
                else:
                    logger.warning(f"Could not get snpEff version: {snpeff_version.stderr}")
            except FileNotFoundError:
                logger.error("snpEff command not found in PATH")
            
            with open(snpeff_out, "w") as fout:
                logger.info("Starting snpEff execution...")
                result = subprocess.run(snpeff_cmd, stdout=fout, stderr=subprocess.PIPE, text=True)
                logger.info(f"snpEff execution completed with return code: {result.returncode}")
            
            if result.returncode != 0:
                logger.warning(f"snpEff failed with return code: {result.returncode}")
                logger.warning(f"snpEff stderr output: {result.stderr}")
                logger.info("Using raw Clair3 output for SNP annotation")
                # Decompress the gzipped file when copying to .vcf extension
                import gzip
                with gzip.open(f"{clair_dir}/output_done.vcf.gz", 'rt') as f_in:
                    with open(snpeff_out, 'w') as f_out:
                        f_out.write(f_in.read())
                logger.info("Created uncompressed VCF from gzipped input")
            else:
                logger.info("snpEff completed successfully!")
                # Check the output file
                if os.path.exists(snpeff_out):
                    output_size = os.path.getsize(snpeff_out)
                    logger.info(f"snpEff output file created, size: {output_size} bytes")
                    
                    # Check first few lines of output to verify annotation
                    try:
                        with open(snpeff_out, 'r') as f:
                            first_lines = [next(f) for _ in range(5)]
                        logger.info(f"First 5 lines of snpEff output:")
                        for i, line in enumerate(first_lines):
                            logger.info(f"  Line {i+1}: {line.strip()}")
                    except Exception as e:
                        logger.warning(f"Could not read snpEff output: {e}")
                else:
                    logger.error("snpEff output file was not created!")
            
            # SnpSift annotation
            logger.info(f"Starting SnpSift annotation for SNPs...")
            logger.info(f"Input file: {snpeff_out}")
            logger.info(f"Output file: {snpsift_out}")
            
            # Check if snpEff output exists and has content
            if os.path.exists(snpeff_out):
                snpeff_size = os.path.getsize(snpeff_out)
                logger.info(f"snpEff output file exists, size: {snpeff_size} bytes")
            else:
                logger.error(f"snpEff output file does not exist: {snpeff_out}")
            
            try:
                from littlejohn import resources
                clinvar_path = os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)), "clinvar.vcf"
                )
                
                logger.info(f"Looking for ClinVar file at: {clinvar_path}")
                
                if os.path.exists(clinvar_path):
                    clinvar_size = os.path.getsize(clinvar_path)
                    logger.info(f"ClinVar file found, size: {clinvar_size} bytes")
                    
                    # Check if SnpSift is available
                    try:
                        snpsift_version = subprocess.run(["SnpSift", "-version"], capture_output=True, text=True)
                        if snpsift_version.returncode == 0:
                            logger.info(f"SnpSift version: {snpsift_version.stdout.strip()}")
                        else:
                            logger.warning(f"Could not get SnpSift version: {snpsift_version.stderr}")
                    except FileNotFoundError:
                        logger.error("SnpSift command not found in PATH")
                    
                    snpsift_cmd = ["SnpSift", "annotate", clinvar_path, snpeff_out]
                    logger.info(f"Running SnpSift command: {' '.join(snpsift_cmd)}")
                    
                    with open(snpsift_out, "w") as fout:
                        logger.info("Starting SnpSift execution...")
                        result = subprocess.run(snpsift_cmd, stdout=fout, stderr=subprocess.PIPE, text=True)
                        logger.info(f"SnpSift execution completed with return code: {result.returncode}")
                    
                    if result.returncode != 0:
                        logger.warning(f"SnpSift failed with return code: {result.returncode}")
                        logger.warning(f"SnpSift stderr output: {result.stderr}")
                        logger.info("Using snpEff output for final SNP annotation")
                        shutil.copy2(snpeff_out, snpsift_out)
                        logger.info("Copied snpEff output as final SNP annotation")
                    else:
                        logger.info("SnpSift completed successfully!")
                        # Check the output file
                        snpsift_out = f"{clair_dir}/snpsift_output.vcf"
                        if os.path.exists(snpsift_out):
                            output_size = os.path.getsize(snpsift_out)
                            logger.info(f"SnpSift output file created, size: {output_size} bytes")
                            
                            # Check first few lines of output to verify annotation
                            try:
                                with open(snpsift_out, 'r') as f:
                                    first_lines = [next(f) for _ in range(5)]
                                logger.info(f"First 5 lines of SnpSift output:")
                                for i, line in enumerate(first_lines):
                                    logger.info(f"  Line {i+1}: {line.strip()}")
                            except Exception as e:
                                logger.warning(f"Could not read SnpSift output: {e}")
                        else:
                            logger.error("SnpSift output file was not created!")
                else:
                    logger.warning(f"ClinVar file not found at {clinvar_path}, using snpEff output")
                    shutil.copy2(snpeff_out, snpsift_out)
                    logger.info("Copied snpEff output as final SNP annotation (no ClinVar)")
                    
            except Exception as e:
                logger.warning(f"SnpSift annotation failed with exception: {e}")
                logger.info("Using snpEff output for final SNP annotation")
                shutil.copy2(snpeff_out, snpsift_out)
                logger.info("Copied snpEff output as final SNP annotation due to exception")
            
            # INDEL annotation
            logger.info(f"Starting snpEff annotation for INDELs...")
            logger.info(f"Input file: {clair_dir}/output_indel_done.vcf.gz")
            logger.info(f"Output file: {snpeff_indel_out}")
            
            # Check if input file exists and has content
            if os.path.exists(f"{clair_dir}/output_indel_done.vcf.gz"):
                file_size = os.path.getsize(f"{clair_dir}/output_indel_done.vcf.gz")
                logger.info(f"INDEL input VCF file exists, size: {file_size} bytes")
                
                # Check first few lines to verify it's a valid VCF
                try:
                    import gzip
                    with gzip.open(f"{clair_dir}/output_indel_done.vcf.gz", 'rt') as f:
                        first_lines = [next(f) for _ in range(5)]
                    logger.info(f"First 5 lines of INDEL input VCF:")
                    for i, line in enumerate(first_lines):
                        logger.info(f"  Line {i+1}: {line.strip()}")
                except Exception as e:
                    logger.warning(f"Could not read INDEL input VCF file: {e}")
            else:
                logger.error(f"INDEL input VCF file does not exist: {clair_dir}/output_indel_done.vcf.gz")
            
            snpeff_indel_cmd = ["snpEff", "-q", "hg38", f"{clair_dir}/output_indel_done.vcf.gz"]
            logger.info(f"Running snpEff INDEL command: {' '.join(snpeff_indel_cmd)}")
            
            with open(snpeff_indel_out, "w") as fout:
                logger.info("Starting snpEff INDEL execution...")
                result = subprocess.run(snpeff_indel_cmd, stdout=fout, stderr=subprocess.PIPE, text=True)
                logger.info(f"snpEff INDEL execution completed with return code: {result.returncode}")
            
            if result.returncode != 0:
                logger.warning(f"snpEff (INDEL) failed with return code: {result.returncode}")
                logger.warning(f"snpEff (INDEL) stderr output: {result.stderr}")
                logger.info("Using raw Clair3 INDEL output for annotation")
                # Decompress the gzipped file when copying to .vcf extension
                import gzip
                with gzip.open(f"{clair_dir}/output_indel_done.vcf.gz", 'rt') as f_in:
                    with open(snpeff_indel_out, 'w') as f_out:
                        f_out.write(f_in.read())
                logger.info("Created uncompressed INDEL VCF from gzipped input")
            else:
                logger.info("snpEff INDEL completed successfully!")
                # Check the output file
                if os.path.exists(snpeff_indel_out):
                    output_size = os.path.getsize(snpeff_indel_out)
                    logger.info(f"snpEff INDEL output file created, size: {output_size} bytes")
                    
                    # Check first few lines of output to verify annotation
                    try:
                        with open(snpeff_indel_out, 'r') as f:
                            first_lines = [next(f) for _ in range(5)]
                        logger.info(f"First 5 lines of snpEff INDEL output:")
                        for i, line in enumerate(first_lines):
                            logger.info(f"  Line {i+1}: {line.strip()}")
                    except Exception as e:
                        logger.warning(f"Could not read snpEff INDEL output: {e}")
                else:
                    logger.error("snpEff INDEL output file was not created!")
            
            # SnpSift annotation for INDELs
            logger.info(f"Starting SnpSift annotation for INDELs...")
            logger.info(f"Input file: {snpeff_indel_out}")
            logger.info(f"Output file: {snpsift_indel_out}")
            
            # Check if snpEff INDEL output exists and has content
            if os.path.exists(snpeff_indel_out):
                snpeff_indel_size = os.path.getsize(snpeff_indel_out)
                logger.info(f"snpEff INDEL output file exists, size: {snpeff_indel_size} bytes")
            else:
                logger.error(f"snpEff INDEL output file does not exist: {snpeff_indel_out}")
            
            try:
                if os.path.exists(clinvar_path):
                    logger.info(f"Using ClinVar file for INDEL annotation: {clinvar_path}")
                    
                    snpsift_indel_cmd = ["SnpSift", "annotate", clinvar_path, snpeff_indel_out]
                    logger.info(f"Running SnpSift INDEL command: {' '.join(snpsift_indel_cmd)}")
                    
                    with open(snpsift_indel_out, "w") as fout:
                        logger.info("Starting SnpSift INDEL execution...")
                        result = subprocess.run(snpsift_indel_cmd, stdout=fout, stderr=subprocess.PIPE, text=True)
                        logger.info(f"SnpSift INDEL execution completed with return code: {result.returncode}")
                    
                    if result.returncode != 0:
                        logger.warning(f"SnpSift (INDEL) failed with return code: {result.returncode}")
                        logger.warning(f"SnpSift (INDEL) stderr output: {result.stderr}")
                        logger.info("Using snpEff INDEL output for final INDEL annotation")
                        shutil.copy2(snpeff_indel_out, snpsift_indel_out)
                        logger.info("Copied snpEff INDEL output as final INDEL annotation")
                    else:
                        logger.info("SnpSift INDEL completed successfully!")
                        # Check the output file
                        if os.path.exists(snpsift_indel_out):
                            output_size = os.path.getsize(snpsift_indel_out)
                            logger.info(f"SnpSift INDEL output file created, size: {output_size} bytes")
                            
                            # Check first few lines of output to verify annotation
                            try:
                                with open(snpsift_indel_out, 'r') as f:
                                    first_lines = [next(f) for _ in range(5)]
                                logger.info(f"First 5 lines of SnpSift INDEL output:")
                                for i, line in enumerate(first_lines):
                                    logger.info(f"  Line {i+1}: {line.strip()}")
                            except Exception as e:
                                logger.warning(f"Could not read SnpSift INDEL output: {e}")
                        else:
                            logger.error("SnpSift INDEL output file was not created!")
                else:
                    logger.warning(f"ClinVar file not found for INDELs, using snpEff INDEL output")
                    shutil.copy2(snpeff_indel_out, snpsift_indel_out)
                    logger.info("Copied snpEff INDEL output as final INDEL annotation (no ClinVar)")
                    
            except Exception as e:
                logger.warning(f"SnpSift INDEL annotation failed with exception: {e}")
                logger.info("Using snpEff INDEL output for final INDEL annotation")
                shutil.copy2(snpeff_indel_out, snpsift_indel_out)
                logger.info("Copied snpEff INDEL output as final INDEL annotation due to exception")
            
            logger.info("Annotation pipeline completed successfully")
            
        except Exception as e:
            logger.error(f"Error in annotation pipeline: {e}")
            return ""
        
        # Parse VCF files to CSV format
        logger.info("Parsing VCF files...")
        try:
            # Import parse_vcf function or implement basic parsing
            def basic_vcf_to_csv(vcf_file, csv_file):
                """Basic VCF to CSV conversion"""
                try:
                    with open(vcf_file, 'r') as f:
                        lines = [line.strip() for line in f if not line.startswith('#')]
                    
                    if not lines:
                        return
                    
                    # Parse header from first non-comment line
                    header = lines[0].split('\t')
                    
                    # Parse data
                    data = []
                    for line in lines[1:]:
                        fields = line.split('\t')
                        if len(fields) >= len(header):
                            data.append(fields[:len(header)])
                    
                    # Create DataFrame and save
                    df = pd.DataFrame(data, columns=header)
                    df.to_csv(csv_file, index=False)
                    logger.info(f"VCF parsed to CSV: {csv_file}")
                    
                except Exception as e:
                    logger.warning(f"Basic VCF parsing failed for {vcf_file}: {e}")
            
            # Parse both SNP and INDEL VCFs
            basic_vcf_to_csv(
                f"{clair_dir}/snpsift_output.vcf", 
                f"{clair_dir}/snpsift_output.vcf.csv"
            )
            basic_vcf_to_csv(
                f"{clair_dir}/snpsift_indel_output.vcf", 
                f"{clair_dir}/snpsift_indel_output.vcf.csv"
            )
            
            # Note: VCF files are now properly uncompressed when they have .vcf extensions
            # Compressed .vcf.gz files are only used as input to processing tools
            
        except Exception as e:
            logger.warning(f"VCF parsing failed: {e}")
        
        logger.info("🎉 FINAL SUCCESS: SNP analysis completed successfully!")
        logger.info(f"SNP analysis completed successfully for {os.path.basename(sample_dir)}")
        logger.info(f"Output directory: {clair_dir}")
        return clair_dir
        
    except Exception as e:
        print(f"❌ EXCEPTION in run_snp_analysis: {e} - PRINT STATEMENT")
        logger.error(f"Failed to run SNP analysis in {sample_dir}: {e}")
        import traceback
        print(f"❌ TRACEBACK: {traceback.format_exc()} - PRINT STATEMENT")
        return ""


def snp_analysis_handler(job, work_dir: Optional[str] = None) -> None:
    """
    Workflow handler to run SNP analysis using Clair3 for a sample.
    
    Expects a per-sample directory containing `target.bam` and `targets_exceeding_threshold.bed`.
    Outputs to `sample_dir/clair3/` with annotated VCF files and CSV summaries.
    On success, records the output directory in the job context under results['snp_analysis'].
    
    Required metadata:
    - reference: Path to reference genome (optional, will auto-detect)
    - threads: Number of threads to use (default: 4)
    - force_regenerate: Whether to force regeneration of existing results (default: False)
    """
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    
    try:
        # Determine sample dir: prefer provided work_dir + sample_id
        sid = job.context.get_sample_id() if hasattr(job.context, "get_sample_id") else "unknown"
        base = work_dir or job.context.metadata.get("work_dir")
        sample_dir = None
        
        logger.info(f"SNP analysis job: sample_id={sid}, work_dir={work_dir}, base={base}")
        
        if base and sid and sid != "unknown":
            sample_dir = os.path.join(base, sid)
            logger.info(f"Using work_dir + sample_id path: {sample_dir}")
        else:
            # Fallback: directory of the current file
            sample_dir = os.path.dirname(job.context.filepath)
            logger.info(f"Using fallback path: {sample_dir}")
        
        # Additional fallback: if the filepath is already a sample directory, use it directly
        if os.path.isdir(job.context.filepath):
            sample_dir = job.context.filepath
            logger.info(f"Using filepath directly as sample directory: {sample_dir}")

        if not os.path.isdir(sample_dir):
            raise RuntimeError(f"Sample directory not found: {sample_dir}")
        
        logger.info(f"Final sample directory: {sample_dir}")
        
        # Let's see what's actually in this directory
        try:
            if os.path.exists(sample_dir):
                dir_contents = os.listdir(sample_dir)
                logger.info(f"📁 Directory contents of {sample_dir}:")
                for item in dir_contents:
                    item_path = os.path.join(sample_dir, item)
                    if os.path.isfile(item_path):
                        size = os.path.getsize(item_path)
                        logger.info(f"  📄 {item} ({size / (1024**2):.1f} MB)")
                    else:
                        logger.info(f"  📁 {item}/")
            else:
                logger.error(f"❌ Sample directory does not exist: {sample_dir}")
        except Exception as e:
            logger.warning(f"⚠️ Could not list directory contents: {e}")
        
        # Check for required input files
        target_bam = os.path.join(sample_dir, "target.bam")
        targets_bed = os.path.join(sample_dir, "targets_exceeding_threshold.bed")
        
        if not os.path.exists(target_bam):
            raise RuntimeError(f"target.bam not found: {target_bam}")
        
        if not os.path.exists(targets_bed):
            raise RuntimeError(f"targets_exceeding_threshold.bed not found: {targets_bed}")
        
        # Get job parameters
        threads = job.context.metadata.get("threads", 4)
        force_regenerate = job.context.metadata.get("force_regenerate", False)
        reference = job.context.metadata.get("reference")
        
        # Debug logging for reference genome
        logger.info(f"=== SNP ANALYSIS HANDLER DEBUGGING ===")
        logger.info(f"Job metadata keys: {list(job.context.metadata.keys())}")
        logger.info(f"Job metadata: {job.context.metadata}")
        logger.info(f"Extracted reference: {reference}")
        logger.info(f"=== END SNP ANALYSIS HANDLER DEBUGGING ===")
        
        if force_regenerate:
            logger.info("Force regenerate requested - will recreate SNP analysis even if it exists")
        
        if reference:
            logger.info(f"✅ Using reference genome from job metadata: {reference}")
        else:
            logger.warning("❌ No reference genome in job metadata - will try to auto-detect")
        
        # Run SNP analysis
        logger.info("🚀 ABOUT TO CALL run_snp_analysis function")
        logger.info(f"Calling run_snp_analysis with: sample_dir={sample_dir}, threads={threads}, force_regenerate={force_regenerate}, reference={reference}")
        
        output_dir = run_snp_analysis(sample_dir, threads=threads, force_regenerate=force_regenerate, reference=reference)
        
        logger.info(f"run_snp_analysis returned: {output_dir}")
        
        if not output_dir:
            raise RuntimeError("Failed to complete SNP analysis")

        # Record result
        try:
            job.context.add_result("snp_analysis", {
                "output_dir": output_dir,
                "snps_vcf": os.path.join(output_dir, "snpsift_output.vcf"),
                "indels_vcf": os.path.join(output_dir, "snpsift_indel_output.vcf"),
                "snps_csv": os.path.join(output_dir, "snpsift_output.vcf.csv"),
                "indels_csv": os.path.join(output_dir, "snpsift_indel_output.vcf.csv"),
                "sample_id": sid,
                "analysis_timestamp": time.time()
            })
        except Exception:
            pass
        
        logger.info(f"SNP analysis completed successfully for {sid}")
        logger.info(f"Output directory: {output_dir}")
        
    except Exception as e:
        job.context.add_error("snp_analysis", str(e))
        logger.error(f"snp_analysis handler failed: {e}")
