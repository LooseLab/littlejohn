#!/usr/bin/env python3
"""
Copy Number Variation (CNV) Analysis Module for robin

This module provides comprehensive CNV analysis following the exact architecture
specified in the CNV BAM Processing Pipeline Documentation.

Key Features
-----------
- Two-pass CNV analysis (sample vs reference)
- Dynamic bin width calculation
- Breakpoint detection using Kernel Change Point Detection
- State persistence and incremental processing
- Comprehensive output file generation

Classes
-------
CNVMetadata
    Container for CNV analysis metadata and results.

CNVAnalysis
    Main analysis class that processes BAM files for CNV analysis.

Dependencies
-----------
- pandas: Data manipulation and analysis
- numpy: Numerical computations
- pysam: BAM file processing
- ruptures: Change point detection
- scipy: Scientific computing
- cnv_from_bam: CNV extraction from BAM files (required)
- logging: Logging for debugging and monitoring
- typing: Type hints
- tempfile: Temporary file creation
- pathlib: File system paths
- os: Operating system interface
- time: Time utilities
- json: JSON serialization
- pickle: Python object serialization
- gc: Garbage collection

Example Usage
-----------
.. code-block:: python

    from robin.analysis.cnv_analysis import CNVAnalysis

    # Initialize analysis
    cnv_analysis = CNVAnalysis(
        work_dir="output/",
        ref_cnv_dict_path="HG01280_control_new.pkl"
    )

    # Process BAM files
    cnv_analysis.process_bam("sample.bam")

Notes
-----
The module follows the robin framework patterns for:
- Integration with workflow system
- Worker process management
- State tracking and persistence
- Error handling and logging
- Output generation and file management

Requirements
-----------
- cnv_from_bam module must be available for CNV extraction
- Reference CNV data must be available in robin module resources

Authors
-------
Matt Loose
"""

import os
import logging
import time
import pickle
import gc
import subprocess
import sys
from typing import Dict, Optional, List, Tuple
from dataclasses import dataclass

import numpy as np
import pysam
from scipy.ndimage import uniform_filter1d

import ruptures as rpt
from robin.logging_config import get_job_logger
import robin.resources as resources

os.environ["CI"] = "1"


# Global configuration for performance tuning
CHUNK_SIZE = 1000  # For processing large arrays in chunks
# Removed MAX_WORKERS - no longer using parallel processing


def run_cnv_analysis_subprocess(
    bam_path,
    copy_numbers,
    ref_cnv_dict,
    temp_dir,
    logger,
    threads=1,
    mapq_filter=60,
    sample_id: str = None,
    update_cnv_dict_path: str = None,
):
    """
    Run CNV analysis using cnv_from_bam in a subprocess to avoid signal handling issues.

    Args:
        bam_path: Path to BAM file
        copy_numbers: Copy numbers dictionary
        ref_cnv_dict: Reference CNV dictionary
        temp_dir: Temporary directory for intermediate files
        logger: Logger instance
        threads: Number of threads to use
        mapq_filter: Mapping quality filter

    Returns:
        Dictionary with analysis results or None if failed
    """
    # try:
    # Ensure temp directory exists
    os.makedirs(temp_dir, exist_ok=True)

    
    # Determine reference CNV dict path: if provided a path, use it; otherwise serialize dict
    if isinstance(ref_cnv_dict, str) and os.path.exists(ref_cnv_dict):
        ref_cnv_dict_path = ref_cnv_dict
    else:
        ref_cnv_dict_path = os.path.join(temp_dir, "ref_cnv_dict.pkl")
        with open(ref_cnv_dict_path, "wb") as f:
            pickle.dump(ref_cnv_dict, f)

    # Get path to the subprocess script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    subprocess_script = os.path.join(script_dir, "cnv_subprocess.py")

    # Run the subprocess
    cmd = [
        sys.executable,  # Use the same Python interpreter
        subprocess_script,
        "--bam-path",
        bam_path,
        # Prefer multi-sample update_cnv_dict with sample id to avoid extra serialization
        "--update-cnv-dict-path",
        (
            update_cnv_dict_path
            if update_cnv_dict_path
            else os.path.join(temp_dir, "update_cnv_dict.pkl")
        ),
        "--sample-id",
        sample_id if sample_id else "unknown",
        "--ref-cnv-dict-path",
        ref_cnv_dict_path,
        "--output-dir",
        temp_dir,
        "--threads",
        str(threads),
        "--mapq-filter",
        str(mapq_filter),
    ]

    logger.debug(f"Running CNV analysis in subprocess: {' '.join(cmd)}")

    # Run subprocess with stdout/stderr redirected to files to avoid capture overhead
    stdout_log_path = os.path.join(temp_dir, "cnv_subprocess.stdout.log")
    stderr_log_path = os.path.join(temp_dir, "cnv_subprocess.stderr.log")
    with open(stdout_log_path, "w") as _out, open(stderr_log_path, "w") as _err:
        result = subprocess.run(
            cmd,
            stdout=_out,
            stderr=_err,
            timeout=3600,  # 1 hour timeout
        )

    if result.returncode != 0:
        logger.error(f"Subprocess failed with return code {result.returncode}")
        try:
            with open(stdout_log_path, "r") as f:
                stdout_content = f.read()
            if stdout_content.strip():
                logger.error(f"subprocess stdout:\n{stdout_content}")
        except Exception:
            pass
        try:
            with open(stderr_log_path, "r") as f:
                stderr_content = f.read()
            if stderr_content.strip():
                logger.error(f"subprocess stderr:\n{stderr_content}")
        except Exception:
            pass
        return None

    # Load results
    results_path = os.path.join(temp_dir, "cnv_analysis_results.pkl")
    if os.path.exists(results_path):
        with open(results_path, "rb") as f:
            return pickle.load(f)
    else:
        logger.error("Results file not found")
        return None

    # except subprocess.TimeoutExpired:
    #    logger.error("CNV analysis subprocess timed out")
    #    return None
    # except Exception as e:
    #    logger.error(f"Error running CNV analysis subprocess: {e}")
    #    return None


class Result:
    """
    A class to store CNV results.
    """

    def __init__(self, cnv_dict: dict) -> None:
        self.cnv = cnv_dict


@dataclass
class CNVMetadata:
    """Container for CNV analysis metadata and results"""

    sample_id: str
    bam_path: str
    analysis_timestamp: float
    cnv_data_path: Optional[str] = None
    cnv_plot_path: Optional[str] = None
    sex_estimate: Optional[str] = None
    breakpoints: Optional[List[Dict]] = None
    chromosome_stats: Optional[Dict] = None
    processing_steps: List[str] = None
    error_message: Optional[str] = None

    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []


def run_ruptures(
    r_cnv: list, penalty_value: int, bin_width: int
) -> List[Tuple[int, int]]:
    """
    Detect change points in a time series using the Kernel Change Point Detection algorithm.

    Args:
        r_cnv: A list containing the time series data.
        penalty_value (int): The penalty value for the change point detection algorithm.
        bin_width (int): The width of the bins used to scale the detected change points.

    Returns:
        List[Tuple[float, float]]: A list of tuples where each tuple represents a detected change point range
                                    as (start, end) with respect to the bin width.
    """
    # Convert to numpy array once for better performance
    r_cnv_array = np.array(r_cnv)

    # Initialize the Kernel Change Point Detection algorithm with the Radial Basis Function (RBF) kernel
    algo_c = rpt.KernelCPD(kernel="rbf").fit(r_cnv_array)

    # Predict the change points using the provided penalty value
    ruptures_result = algo_c.predict(pen=penalty_value)

    # Compute the ranges around each change point using vectorized operations
    return [
        (cp * bin_width - bin_width, cp * bin_width + bin_width)
        for cp in ruptures_result
    ]


class CNV_Difference:
    """
    A class to store CNV difference data.
    """

    def __init__(self, *args, **kwargs) -> None:
        self.cnv = {}


def moving_average(data: np.ndarray, n: int = 3) -> np.ndarray:
    """
    Calculate the moving average of a given data array using scipy's uniform_filter1d.

    Args:
        data (np.ndarray): The data array.
        n (int): The window size for the moving average.

    Returns:
        np.ndarray: The array of moving averages.
    """
    return uniform_filter1d(data, size=n, mode="nearest")


def pad_arrays(
    arr1: np.ndarray, arr2: np.ndarray, pad_value: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Pad two arrays to the same length with a specified value.

    Args:
        arr1 (np.ndarray): The first array.
        arr2 (np.ndarray): The second array.
        pad_value (int): The value to pad with.

    Returns:
        Tuple[np.ndarray, np.ndarray]: The padded arrays.
    """
    len1, len2 = len(arr1), len(arr2)
    max_len = max(len1, len2)

    # Use numpy's pad function more efficiently
    if len1 < max_len:
        arr1 = np.pad(
            arr1, (0, max_len - len1), mode="constant", constant_values=pad_value
        )
    if len2 < max_len:
        arr2 = np.pad(
            arr2, (0, max_len - len2), mode="constant", constant_values=pad_value
        )

    return arr1, arr2


def has_reads(bam_file: str) -> bool:
    """Quickly checks if a BAM file has any reads."""
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for _ in bam.fetch(until_eof=True):
                return True
        return False
    except Exception as e:
        # Replace print with logging
        logging.error(f"Error checking reads in {bam_file}: {e}")
        return False


def estimate_sex_from_cnv(cnv_data: Dict, logger) -> str:
    """
    Estimate genetic sex (XX or XY) based on CNV data.

    Args:
        cnv_data: Dictionary containing normalized CNV data
        logger: Logger instance

    Returns:
        String indicating estimated sex
    """
    try:
        if "chrX" not in cnv_data or "chrY" not in cnv_data:
            return "Unknown"

        # Calculate average CNV values for X and Y chromosomes more efficiently
        x_data = [i for i in cnv_data["chrX"] if i != 0]
        y_data = [i for i in cnv_data["chrY"] if i != 0]

        if not x_data or not y_data:
            return "Unknown"

        x_avg = np.mean(x_data)
        y_avg = np.mean(y_data)

        logger.debug(f"X chromosome average: {x_avg:.3f}")
        logger.debug(f"Y chromosome average: {y_avg:.3f}")

        # Sex estimation logic from original code
        if x_avg >= 0.1 and y_avg <= -0.1:
            return "Female"
        elif x_avg <= 0.1 and y_avg >= -0.2:
            return "Male"
        elif x_avg >= 0.1 and y_avg >= -0.1:
            return "Male (query X/Y copy number changes)"
        elif x_avg > 0.1 and y_avg > 0.1:
            return "Male (query X/Y copy number changes)"
        elif x_avg < 0.1 and y_avg < -0.2:
            return "Unknown (Query XY copy number changes)"
        else:
            return "Unknown"

    except Exception as e:
        logger.error(f"Error estimating sex: {e}")
        return "Unknown"


def detect_breakpoints_from_cnv(cnv_data: Dict, bin_width: int, logger) -> List[Dict]:
    """
    Detect breakpoints in CNV data using change point detection.

    Args:
        cnv_data: Dictionary containing CNV data
        bin_width: Width of bins
        logger: Logger instance

    Returns:
        List of detected breakpoints
    """
    breakpoints = []

    try:
        logger.debug("Detecting breakpoints in CNV data")

        # Simple sequential processing for all datasets
        for key, data in cnv_data.items():
            if key != "chrM" and len(data) > 3:
                try:
                    paired_changepoints = run_ruptures(
                        data, penalty_value=5, bin_width=bin_width
                    )

                    for start, end in paired_changepoints:
                        if start >= 0 and end > start:
                            breakpoints.append(
                                {
                                    "chromosome": key,
                                    "start": start,
                                    "end": end,
                                    "length": end - start,
                                }
                            )
                except Exception as e:
                    logger.debug(f"Error processing chromosome {key}: {e}")

        logger.debug(f"Detected {len(breakpoints)} breakpoints")
        return breakpoints

    except Exception as e:
        logger.error(f"Error in breakpoint detection: {e}")
        return []


def calculate_chromosome_stats_from_cnv(cnv_data: Dict, logger) -> Dict:
    """
    Calculate chromosome-wide statistics from normalized CNV data.

    Args:
        cnv_data: Dictionary containing normalized CNV data
        logger: Logger instance

    Returns:
        Dictionary of chromosome statistics
    """
    stats = {}

    try:
        logger.debug("Calculating chromosome statistics")

        for chrom, data in cnv_data.items():
            if len(data) == 0:
                continue

            # Calculate basic statistics more efficiently using numpy
            data_array = np.array(data)
            stats[chrom] = {
                "mean": float(np.mean(data_array)),
                "std": float(np.std(data_array)),
                "min": float(np.min(data_array)),
                "max": float(np.max(data_array)),
                "length": len(data_array),
            }

        logger.debug(f"Calculated statistics for {len(stats)} chromosomes")
        return stats

    except Exception as e:
        logger.error(f"Error calculating chromosome stats: {e}")
        return {}


def load_analysis_counter(sample_id: str, work_dir: str, logger) -> int:
    """Load the analysis counter for a sample from disk with caching"""
    counter_file = os.path.join(work_dir, sample_id, "cnv_analysis_counter.txt")
    if os.path.exists(counter_file):
        try:
            with open(counter_file, "r") as f:
                return int(f.read().strip())
        except (ValueError, IOError) as e:
            logger.warning(f"Error loading counter for {sample_id}: {e}")
    return 0


def save_analysis_counter(sample_id: str, counter: int, work_dir: str, logger) -> None:
    """Save the analysis counter for a sample to disk"""
    counter_file = os.path.join(work_dir, sample_id, "cnv_analysis_counter.txt")
    try:
        # Create sample directory if it doesn't exist
        os.makedirs(os.path.dirname(counter_file), exist_ok=True)
        with open(counter_file, "w") as f:
            f.write(str(counter))
    except IOError as e:
        logger.error(f"Error saving counter for {sample_id}: {e}")


def find_significant_regions(values, chrom):
    """
    Find regions with significant CNV changes.

    Args:
        values: CNV values for a chromosome
        chrom: Chromosome name

    Returns:
        List of significant regions
    """
    regions = []
    threshold = 0.5  # CNV change threshold

    # Use numpy for better performance
    values_array = np.array(values)
    significant_indices = np.where(np.abs(values_array) > threshold)[0]

    for i in significant_indices:
        regions.append(
            {
                "chromosome": chrom,
                "start": i * 1000000,  # Approximate position
                "end": (i + 1) * 1000000,
                "type": "gain" if values_array[i] > 0 else "loss",
            }
        )

    return regions


def generate_bed_files(bed_dir, analysis_counter, breakpoints, cnv_data, logger):
    """
    Generate BED files for CNV regions and breakpoints.

    Args:
        bed_dir: BED files directory
        analysis_counter: Current analysis counter
        breakpoints: Detected breakpoints
        cnv_data: CNV data
        logger: Logger instance
    """
    try:
        # Generate BED file for CNV regions
        bed_file = os.path.join(bed_dir, f"new_file_{analysis_counter:03d}.bed")
        with open(bed_file, "w") as f:
            for chrom, values in cnv_data.items():
                if len(values) > 0:
                    # Find regions with significant CNV changes
                    significant_regions = find_significant_regions(values, chrom)
                    for region in significant_regions:
                        f.write(
                            f"{chrom}\t{region['start']}\t{region['end']}\t{region['type']}\n"
                        )

        # Generate BED file for breakpoints
        if breakpoints:
            breakpoint_file = os.path.join(
                bed_dir, f"breakpoints_{analysis_counter:03d}.bed"
            )
            with open(breakpoint_file, "w") as f:
                for bp in breakpoints:
                    f.write(
                        f"{bp['chromosome']}\t{bp['start']}\t{bp['end']}\tbreakpoint\n"
                    )

        logger.debug(f"Generated BED files in {bed_dir}")

    except Exception as e:
        logger.error(f"Error generating BED files: {e}")


def save_cnv_files(
    sample_dir,
    analysis_counter,
    r_cnv,
    r2_cnv,
    result3_cnv,
    bin_width,
    variance,
    breakpoints,
    sex_estimate,
    update_cnv_dict,
    logger,
):
    """
    Save CNV data files in the specified format from the documentation.

    Args:
        sample_dir: Sample output directory
        analysis_counter: Current analysis counter
        r_cnv: Sample CNV data
        r2_cnv: Reference CNV data
        result3_cnv: Normalized CNV data
        bin_width: Bin width
        variance: Variance
        breakpoints: Detected breakpoints
        sex_estimate: Sex estimation result
        update_cnv_dict: Accumulated copy numbers dictionary for incremental processing
        logger: Logger instance
    """
    try:
        # Save CNV data files as specified in documentation
        np.save(os.path.join(sample_dir, "CNV.npy"), r_cnv)
        np.save(os.path.join(sample_dir, "CNV2.npy"), r2_cnv)
        np.save(os.path.join(sample_dir, "CNV3.npy"), result3_cnv)

        # Save processing metadata
        cnv_dict = {
            "bin_width": bin_width,
            "variance": variance,
            "analysis_counter": analysis_counter,
        }
        np.save(os.path.join(sample_dir, "CNV_dict.npy"), cnv_dict)

        # Save breakpoint data as structured array
        if breakpoints:
            breakpoint_array = np.array(
                [(bp["chromosome"], bp["start"], bp["end"]) for bp in breakpoints],
                dtype=[("name", "U10"), ("start", "i8"), ("end", "i8")],
            )
            np.save(os.path.join(sample_dir, "cnv_data_array.npy"), breakpoint_array)

        # Helper for atomic pickle writes to avoid readers seeing partial files
        def _atomic_pickle_dump(obj, target_path: str) -> None:
            tmp_path = target_path + ".tmp"
            with open(tmp_path, "wb") as tf:
                pickle.dump(obj, tf)
                tf.flush()
                os.fsync(tf.fileno())
            os.replace(tmp_path, target_path)

        # Save sex estimation (atomic)
        _atomic_pickle_dump(sex_estimate, os.path.join(sample_dir, "XYestimate.pkl"))

        # Save accumulated copy numbers for incremental processing (atomic)
        _atomic_pickle_dump(
            update_cnv_dict, os.path.join(sample_dir, "update_cnv_dict.pkl")
        )

        # Create state directory for CNV detector
        state_dir = os.path.join(sample_dir, "cnv_detector_state")
        os.makedirs(state_dir, exist_ok=True)

        # Save state metadata
        state_metadata = {
            "analysis_counter": analysis_counter,
            "timestamp": time.time(),
            "bin_width": bin_width,
            "variance": variance,
            "breakpoints_count": len(breakpoints),
            "sex_estimate": sex_estimate,
        }
        _atomic_pickle_dump(
            state_metadata, os.path.join(state_dir, "tracker_metadata.pkl")
        )

        # Create BED files directory
        bed_dir = os.path.join(sample_dir, "bed_files")
        os.makedirs(bed_dir, exist_ok=True)

        # Generate BED files for CNV regions and breakpoints
        generate_bed_files(bed_dir, analysis_counter, breakpoints, result3_cnv, logger)

        logger.debug(f"Saved all CNV files to {sample_dir}")

    except Exception as e:
        logger.error(f"Error saving CNV files: {e}")




def process_single_bam(bam_path, metadata, work_dir, logger):
    """
    Process a single BAM file for CNV analysis using the complete pipeline.

    Args:
        bam_path: Path to BAM file
        metadata: BamMetadata object
        work_dir: Working directory
        logger: Logger instance

    Returns:
        Dictionary with CNV analysis results
    """
    # print(f"Processing CNV for BAM file: {bam_path}")
    sample_id = metadata.extracted_data.get("sample_id", "unknown")
    logger.info(f"🧬 Starting CNV analysis for sample: {sample_id}")

    # Log essential metadata only
    logger.debug(f"BAM file: {metadata.file_path} ({metadata.file_size:,} bytes)")
    logger.debug(f"Sample ID: {sample_id}")

    # Create sample-specific output directory
    sample_output_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Load analysis counter from disk
    analysis_counter = load_analysis_counter(sample_id, work_dir, logger)
    logger.debug(f"Analysis counter: {analysis_counter}")

    analysis_result = {
        "sample_id": sample_id,
        "bam_path": bam_path,
        "analysis_timestamp": time.time(),
        "analysis_counter": analysis_counter,
        "cnv_data_path": None,
        "sex_estimate": None,
        "breakpoints": [],
        "chromosome_stats": {},
        "processing_steps": [],
        "error_message": None,
    }

    try:
        # Use the sample output directory directly for intermediate files
        temp_dir = sample_output_dir

        # Check if BAM has reads
        if not has_reads(bam_path):
            logger.warning(f"No reads found in BAM file for {sample_id}")
            analysis_result["error_message"] = "No reads found in BAM file"
            analysis_result["processing_steps"].append("no_reads_found")
            return analysis_result

        analysis_result["processing_steps"].append("reads_found")

        # Load accumulated copy numbers from previous iterations
        update_cnv_dict_path = os.path.join(sample_output_dir, "update_cnv_dict.pkl")
        if os.path.exists(update_cnv_dict_path):
            # Robust load with retries to avoid reading while a writer is replacing the file
            load_exc = None
            for attempt in range(5):
                try:
                    with open(update_cnv_dict_path, "rb") as f:
                        update_cnv_dict = pickle.load(f)
                    logger.debug(
                        "Loaded accumulated copy numbers from previous iterations"
                    )
                    if sample_id in update_cnv_dict:
                        copy_numbers = update_cnv_dict[sample_id]
                        logger.debug(f"Found existing data for {sample_id}")
                    else:
                        copy_numbers = {}
                        logger.debug(
                            f"No existing data for {sample_id}, starting fresh"
                        )
                    load_exc = None
                    break
                except Exception as e:
                    print("Cnv LOADING failed - trying again.")
                    load_exc = e
                    time.sleep(0.2 * (attempt + 1))
            if load_exc is not None:
                logger.warning(
                    f"Error loading accumulated copy numbers after retries: {load_exc}"
                )
                copy_numbers = {}
                update_cnv_dict = {}
        else:
            copy_numbers = {}
            update_cnv_dict = {}
            logger.debug("No previous copy numbers found, starting fresh")

        # Compute existing reference CNV dict path (avoid re-serializing; pass path)
        ref_cnv_path = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "HG01280_control_new.pkl",
        )

        # Process BAM file with cnv_from_bam using subprocess
        logger.debug("Processing BAM file with cnv_from_bam (subprocess)")
        try:
            # Run CNV analysis in subprocess
            subprocess_result = run_cnv_analysis_subprocess(
                bam_path,
                copy_numbers,
                ref_cnv_path,
                temp_dir,
                logger,
                threads=4,
                mapq_filter=60,
                sample_id=sample_id,
                update_cnv_dict_path=update_cnv_dict_path,
            )

            if subprocess_result is None or not subprocess_result.get("success", False):
                error_msg = (
                    subprocess_result.get("error", "Unknown error")
                    if subprocess_result
                    else "Subprocess failed"
                )
                raise RuntimeError(f"CNV analysis subprocess failed: {error_msg}")

            # Extract results from subprocess
            r_cnv = subprocess_result["r_cnv"]
            r_bin = subprocess_result["r_bin"]
            r_var = subprocess_result["r_var"]
            genome_length = subprocess_result["genome_length"]
            r2_cnv = subprocess_result["r2_cnv"]
            updated_copy_numbers = subprocess_result["updated_copy_numbers"]

            # Update the copy numbers dictionary
            if sample_id not in update_cnv_dict:
                update_cnv_dict[sample_id] = {}
            update_cnv_dict[sample_id] = updated_copy_numbers

            analysis_result["processing_steps"].append("cnv_data_extracted")
            logger.debug("CNV data extracted successfully")
            logger.debug(
                f"Bin width: {r_bin:,}, Variance: {r_var:.6f}, Genome length: {genome_length:,}"
            )

        except Exception as e:
            logger.error(f"Error extracting CNV data: {e}")
            analysis_result["error_message"] = f"CNV extraction failed: {str(e)}"
            analysis_result["processing_steps"].append("cnv_extraction_failed")
            return analysis_result

        # Calculate normalized CNV data (difference between sample and reference)
        logger.debug("Calculating normalized CNV data")
        result3_cnv = {}
        for key in r_cnv.keys():
            if key != "chrM" and key in r2_cnv:
                moving_avg_data1 = moving_average(r_cnv[key])
                moving_avg_data2 = moving_average(r2_cnv[key])
                moving_avg_data1, moving_avg_data2 = pad_arrays(
                    moving_avg_data1, moving_avg_data2
                )
                result3_cnv[key] = moving_avg_data1 - moving_avg_data2

        analysis_result["processing_steps"].append("normalized_cnv_calculated")

        # Estimate sex from CNV data
        sex_estimate = estimate_sex_from_cnv(result3_cnv, logger)
        analysis_result["sex_estimate"] = sex_estimate
        analysis_result["processing_steps"].append("sex_estimated")

        # Detect breakpoints
        breakpoints = detect_breakpoints_from_cnv(r_cnv, r_bin, logger)
        analysis_result["breakpoints"] = breakpoints
        analysis_result["processing_steps"].append("breakpoints_detected")

        # Calculate chromosome statistics
        chromosome_stats = calculate_chromosome_stats_from_cnv(result3_cnv, logger)
        analysis_result["chromosome_stats"] = chromosome_stats
        analysis_result["processing_steps"].append("chromosome_stats_calculated")

        # Save CNV data in the specified format
        analysis_counter += 1

        # Save CNV data files as specified in the documentation
        save_cnv_files(
            sample_output_dir,
            analysis_counter,
            r_cnv,
            r2_cnv,
            result3_cnv,
            r_bin,
            r_var,
            breakpoints,
            sex_estimate,
            update_cnv_dict,
            logger,
        )

        analysis_result["cnv_data_path"] = os.path.join(
            sample_output_dir, f"{analysis_counter}_cnv_data.json"
        )
        analysis_result["processing_steps"].append("cnv_data_saved")

        # Save updated analysis counter to disk
        save_analysis_counter(sample_id, analysis_counter, work_dir, logger)
        analysis_result["analysis_counter"] = analysis_counter
        logger.debug(f"Saved analysis counter: {analysis_counter}")

        # Force garbage collection
        gc.collect()

        analysis_result["processing_steps"].append("cnv_analysis_complete")

        # Replace print with logging
        logger.debug(f"Analysis result: {analysis_result}")
        logger.info(f"CNV analysis complete for {sample_id}")
        logger.info(f"Sex Estimate: {analysis_result['sex_estimate']}")
        logger.info(f"Breakpoints: {len(analysis_result['breakpoints'])} detected")

        return analysis_result

    except Exception as e:
        logger.error(f"Error in CNV analysis for {sample_id}: {e}")
        analysis_result["error_message"] = str(e)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def cnv_handler(job, work_dir=None, target_panel="rCNS2"):
    """
    Handler function for CNV analysis jobs.
    This function processes BAM files for copy number variation analysis.

    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to BAM file directory)
        target_panel: Target panel type for consistency with other analysis components
    """
    # Get job-specific logger
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)

    # try:
    bam_path = job.context.filepath

    logger.info(f"Starting CNV analysis for: {os.path.basename(bam_path)}")
    
    # Log and validate target panel
    job_panel = job.context.metadata.get("target_panel", target_panel)
    if job_panel != target_panel:
        logger.warning(f"Panel mismatch: job metadata has '{job_panel}' but handler received '{target_panel}'. Using '{job_panel}' from metadata.")
        target_panel = job_panel
    logger.info(f"Using target panel: {target_panel}")

    # Get metadata from preprocessing
    bam_metadata = job.context.metadata.get("bam_metadata", {})

    # Create BamMetadata object for compatibility
    from robin.analysis.bam_preprocessor import BamMetadata

    metadata = BamMetadata(
        file_path=bam_path,
        file_size=job.context.metadata.get("file_size", 0),
        creation_time=job.context.metadata.get("creation_time", time.time()),
        extracted_data=bam_metadata,
    )

    # Determine work directory
    if work_dir is None:
        # Default to BAM file directory
        work_dir = os.path.dirname(bam_path)
    else:
        # Use specified work directory, create if it doesn't exist
        os.makedirs(work_dir, exist_ok=True)
        logger.debug(f"Using specified work directory: {work_dir}")

    # Process the BAM file
    result = process_single_bam(bam_path, metadata, work_dir, logger)

    # Store results in job context
    job.context.add_metadata("cnv_analysis", result)

    if result.get("error_message"):
        job.context.add_error("cnv_analysis", result["error_message"])
        logger.error(f"CNV analysis failed: {result['error_message']}")
    else:
        job.context.add_result(
            "cnv_analysis",
            {
                "status": "success",
                "sample_id": result.get("sample_id", "unknown"),
                "analysis_time": result.get("analysis_timestamp", 0),
                "sex_estimate": result.get("sex_estimate", "Unknown"),
                "breakpoints_count": len(result.get("breakpoints", [])),
                "processing_steps": result.get("processing_steps", []),
                "cnv_data_path": result.get("cnv_data_path", ""),
                "analysis_counter": result.get("analysis_counter", 0),
            },
        )
        logger.info(f"CNV analysis complete for {os.path.basename(bam_path)}")
        logger.info(f"Sample ID: {result.get('sample_id', 'unknown')}")
        logger.info(f"Sex Estimate: {result.get('sex_estimate', 'Unknown')}")
        logger.info(f"Breakpoints: {len(result.get('breakpoints', []))}")
        logger.info(f"Analysis Counter: {result.get('analysis_counter', 0)}")
        logger.debug(
            f"Processing steps: {', '.join(result.get('processing_steps', []))}"
        )
        logger.debug(
            f"Output directory: {os.path.dirname(result.get('cnv_data_path', ''))}"
        )

    # except Exception as e:
    #    job.context.add_error('cnv_analysis', str(e))
    #    logger.error(f"Error in CNV analysis for {job.context.filepath}: {e}")
