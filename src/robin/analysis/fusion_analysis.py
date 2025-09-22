#!/usr/bin/env python3
"""
Simplified Fusion Analysis Module for robin

This module provides workflow integration and standalone analysis capabilities
for fusion detection. Core fusion detection logic is imported from fusion_work.py
to avoid code duplication.

Key functions:
- Workflow integration (fusion_handler)
- Standalone file processing (process_single_file)
- Metadata management (FusionMetadata)
"""

import os
import tempfile
import logging
import time
import json
from typing import Dict, Any, Optional, List, Tuple, Set
from dataclasses import dataclass
import numpy as np
import pandas as pd
import pysam
from robin.logging_config import get_job_logger

# Import core fusion detection logic from fusion_work.py
from robin.analysis.fusion_work import (
    process_bam_file,
    _generate_output_files,
    FusionMetadata,
    GeneRegion,
    _setup_file_paths,
    _load_bed_regions,
    _ensure_gene_regions_loaded,
    has_supplementary_alignments,
    find_reads_with_supplementary,
    _find_gene_intersections,
    _process_reads_for_fusions,
    _optimize_fusion_dataframe,
    _filter_fusion_candidates,
    process_bam_for_fusions_work,
    _merge_fusion_metadata_objects,
    _load_fusion_metadata,
    preprocess_fusion_data_standalone,
)

logger = logging.getLogger(__name__)


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
    else:
        return obj


def _merge_fusion_metadata(
    new_metadata: FusionMetadata, existing_metadata: Optional[Dict]
) -> FusionMetadata:
    """
    Merge new fusion metadata with existing metadata.

    Args:
        new_metadata: New fusion metadata object
        existing_metadata: Existing metadata dictionary

    Returns:
        Merged FusionMetadata object
    """
    if existing_metadata is None:
        return new_metadata

    # Merge processing steps
    existing_steps = existing_metadata.get("processing_steps", [])
    new_steps = new_metadata.processing_steps
    merged_steps = list(set(existing_steps + new_steps))  # Remove duplicates

    # Merge fusion data
    existing_fusion_data = existing_metadata.get("fusion_data", {})
    new_fusion_data = new_metadata.fusion_data or {}

    # Combine target candidates
    existing_target = existing_fusion_data.get("target_candidates", [])
    new_target = new_fusion_data.get("target_candidates", [])
    merged_target = existing_target + new_target

    # Combine genome-wide candidates
    existing_genome = existing_fusion_data.get("genome_wide_candidates", [])
    new_genome = new_fusion_data.get("genome_wide_candidates", [])
    merged_genome = existing_genome + new_genome

    # Update the new metadata object
    new_metadata.processing_steps = merged_steps
    new_metadata.fusion_data = {
        "target_candidates": merged_target,
        "genome_wide_candidates": merged_genome,
    }

    # Update analysis results counts
    if new_metadata.analysis_results:
        new_metadata.analysis_results["target_candidates_count"] = len(merged_target)
        new_metadata.analysis_results["genome_wide_candidates_count"] = len(
            merged_genome
        )

    logging.info(
        f"Merged fusion metadata: {len(existing_steps)} existing steps + {len(new_steps)} new steps"
    )

    return new_metadata


def _perform_fusion_analysis(
    file_path: str,
    temp_dir: str,
    metadata: Dict[str, Any],
    fusion_metadata: FusionMetadata,
    target_panel: str = "rCNS2",
    has_supplementary: bool = False,
    supplementary_read_ids: List[str] = [],
    work_dir: str = None,
) -> Dict[str, Any]:
    """
    Perform the actual fusion analysis using fusion_work.py.

    Args:
        file_path: Path to the file to process
        temp_dir: Temporary directory for intermediate files
        metadata: File metadata
        fusion_metadata: Fusion metadata object
        target_panel: Target panel to use
        has_supplementary: Whether file has supplementary reads
        supplementary_read_ids: List of supplementary read IDs
        work_dir: Working directory for saving metadata

    Returns:
        Dictionary with analysis results
    """
    fusion_metadata.processing_steps.append("analysis_started")
    results = process_bam_file(
        file_path,
        temp_dir,
        metadata,
        fusion_metadata,
        target_panel,
        has_supplementary,
        supplementary_read_ids,
        work_dir,
    )
    return results


def process_single_file(
    file_path: str,
    metadata: Dict[str, Any],
    work_dir: str,
    logger=None,
    target_panel: str = "rCNS2",
    has_supplementary: bool = False,
    supplementary_read_ids: List[str] = [],
) -> Dict[str, Any]:
    """
    Process a single file for fusion analysis (standalone function-based approach).

    Args:
        file_path: Path to the file to process
        metadata: Metadata about the file
        work_dir: Working directory for output
        logger: Logger instance
        target_panel: Target panel to use for analysis
        has_supplementary: Whether file has supplementary reads
        supplementary_read_ids: List of supplementary read IDs

    Returns:
        Dictionary with processing results
    """
    try:
        if logger:
            logger.info(f"Starting fusion analysis for {file_path}")
            logger.info(f"Using work directory: {work_dir}")
            logger.info(f"Target panel: {target_panel}")

        # Extract sample ID from metadata
        sample_id = metadata.get("sample_id", "unknown")

        # Create metadata object
        fusion_metadata = FusionMetadata(
            sample_id=sample_id,
            file_path=file_path,
            analysis_timestamp=time.time(),
            target_panel=target_panel,
        )

        fusion_metadata.processing_steps.append("started")

        # Check if file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        # Check if file is readable
        if not os.access(file_path, os.R_OK):
            raise PermissionError(f"File not readable: {file_path}")

        # Check if file is a BAM file
        if not file_path.lower().endswith(".bam"):
            if logger:
                logger.warning(f"File does not have .bam extension: {file_path}")

        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            fusion_metadata.processing_steps.append("temp_dir_created")

            # Perform fusion analysis using fusion_work.py
            analysis_results = _perform_fusion_analysis(
                file_path,
                temp_dir,
                metadata,
                fusion_metadata,
                target_panel,
                has_supplementary,
                supplementary_read_ids,
                work_dir,
            )

            # Generate output files using fusion_work.py
            output_files = _generate_output_files(
                sample_id, analysis_results, fusion_metadata, work_dir
            )

            # Update metadata with results (now includes merged data)
            fusion_metadata.analysis_results = analysis_results
            fusion_metadata.target_fusion_path = output_files.get("target_fusion")
            fusion_metadata.genome_wide_fusion_path = output_files.get(
                "genome_wide_fusion"
            )

            # Update analysis results with final merged counts
            if fusion_metadata.fusion_data:
                fusion_metadata.analysis_results["target_candidates_count"] = len(
                    fusion_metadata.fusion_data.get("target_candidates", [])
                )
                fusion_metadata.analysis_results["genome_wide_candidates_count"] = len(
                    fusion_metadata.fusion_data.get("genome_wide_candidates", [])
                )

            fusion_metadata.processing_steps.append("completed")

            if logger:
                logger.info(f"Fusion analysis completed for {file_path}")

            # Convert to dictionary
            result = {
                "success": fusion_metadata.error_message is None,
                "sample_id": fusion_metadata.sample_id,
                "file_path": fusion_metadata.file_path,
                "analysis_timestamp": fusion_metadata.analysis_timestamp,
                "target_fusion_path": fusion_metadata.target_fusion_path,
                "genome_wide_fusion_path": fusion_metadata.genome_wide_fusion_path,
                "processing_steps": fusion_metadata.processing_steps,
                "error_message": fusion_metadata.error_message,
                "analysis_results": fusion_metadata.analysis_results,
            }

            if fusion_metadata.error_message:
                if logger:
                    logger.error(
                        f"Fusion analysis failed: {fusion_metadata.error_message}"
                    )
                    logger.error(
                        f"Processing steps completed: {fusion_metadata.processing_steps}"
                    )
            else:
                if logger:
                    logger.info(
                        f"Fusion analysis completed successfully for {file_path}"
                    )
                    logger.info(f"Processing steps: {fusion_metadata.processing_steps}")

            return result

    except Exception as e:
        import traceback

        error_msg = f"Error in fusion analysis: {str(e)}"
        if logger:
            logger.error(error_msg)
            logger.error(f"Traceback: {traceback.format_exc()}")
        return {
            "success": False,
            "error_message": error_msg,
            "file_path": file_path,
        }


def fusion_handler(job, work_dir=None, target_panel="rCNS2"):
    """
    Handler function for fusion analysis jobs in the workflow system.

    Args:
        job: Job object from the workflow system
        work_dir: Working directory for output
        target_panel: Target gene panel for fusion analysis (rCNS2, AML, or PanCan)
    """
    try:
        # Get logger with proper parameters
        logger = get_job_logger(str(job.job_id), "fusion", job.context.filepath)

        # Extract file path and metadata from job
        file_path = job.context.filepath
        metadata = job.context.metadata.get("bam_metadata", {})

        logger.info(f"DEBUG: fusion_handler called with target_panel='{target_panel}'")
        logger.info(f"DEBUG: Job metadata: {job.context.metadata}")

        # Access supplementary read information
        has_supplementary = metadata.get("has_supplementary_reads", False)
        if has_supplementary:
            logger.info(f"Starting fusion analysis for {file_path}")
            logger.info(f"Metadata: {metadata}")

            # Prefer disk-based supplementary read IDs if available to avoid large in-memory lists
            supplementary_read_ids = metadata.get("supplementary_read_ids", [])
            supp_ids_path = metadata.get("supplementary_read_ids_path")
            if (
                (not supplementary_read_ids)
                and supp_ids_path
                and os.path.exists(supp_ids_path)
            ):
                try:
                    with open(supp_ids_path, "r") as f:
                        supplementary_read_ids = [
                            line.strip() for line in f if line.strip()
                        ]
                except Exception as e:
                    logger.warning(
                        f"Could not read supplementary_read_ids from {supp_ids_path}: {e}"
                    )

            # Set default work directory if not provided
            if work_dir is None:
                work_dir = "fusion_output"

            logger.info(f"Using work directory: {work_dir}")
            logger.info(f"Using target panel: {target_panel}")

            # Use the target_panel parameter passed to the handler
            # (overrides any value in metadata)

            # Process the file using fusion_work.py
            result = process_single_file(
                file_path,
                metadata,
                work_dir,
                logger,
                target_panel,
                has_supplementary=has_supplementary,
                supplementary_read_ids=supplementary_read_ids,
            )

            # Add result to job context
            job.context.add_result("fusion_analysis", result)

            if result["success"]:
                logger.info(f"Fusion analysis completed successfully for {file_path}")
                logger.info(f"Results: {result}")
                # Cleanup supplementary IDs file if present
                if supp_ids_path and os.path.exists(supp_ids_path):
                    try:
                        os.remove(supp_ids_path)
                    except Exception:
                        pass
            else:
                error_msg = result.get("error_message", "Unknown error")
                logger.error(f"Fusion analysis failed for {file_path}: {error_msg}")
                logger.error(f"Full result: {result}")
                job.context.add_error("fusion_analysis", error_msg)
        else:
            logger.info(f"No supplementary reads found for {file_path}")
            # This is not an error - just a normal skip condition
            job.context.add_result(
                "fusion_analysis",
                {
                    "success": True,
                    "skipped": True,
                    "reason": "No supplementary reads found",
                    "file_path": file_path,
                },
            )

    except Exception as e:
        import traceback

        error_msg = f"Error in fusion handler: {str(e)}"

        # Try to get a logger, but don't fail if we can't
        try:
            logger = get_job_logger(str(job.job_id), "fusion", job.context.filepath)
            logger.error(error_msg)
            logger.error(f"Traceback: {traceback.format_exc()}")
            logger.error(
                f"Job details - ID: {job.job_id}, Type: {job.job_type}, File: {job.context.filepath}"
            )
        except Exception as logger_error:
            # If we can't get a logger, print to stderr
            import sys

            print(f"FUSION HANDLER ERROR: {error_msg}", file=sys.stderr)
            print(f"LOGGER ERROR: {logger_error}", file=sys.stderr)
            print(f"TRACEBACK: {traceback.format_exc()}", file=sys.stderr)

        job.context.add_error("fusion_analysis", error_msg)