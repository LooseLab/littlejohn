#!/usr/bin/env python3
"""
Temporary utilities for LittleJohn BAM to parquet conversion.

This module provides utilities for merging modkit files and creating parquet files,
adapted from the marion package for use in LittleJohn.
"""

import warnings
from typing import List, Optional
import gc
import os
import logging
import subprocess
import time
import pickle
from datetime import datetime
from pathlib import Path
import pandas as pd
import polars as pl
from contextlib import contextmanager
import tempfile

# Suppress pkg_resources deprecation warnings from sorted_nearest
warnings.filterwarnings("ignore", message="pkg_resources is deprecated", category=UserWarning)

try:
    import pyranges as pr
    import pysam
    from alive_progress import alive_bar
    from tqdm import tqdm
    from robin.utilities.ReadBam import ReadBam
    from robin.utilities.mnp_flex import APIClient as MnpFlexClient
    from robin import resources
except ImportError as e:
    logging.warning(f"Some dependencies not available: {e}")

import json

# Simple cross-process file lock using POSIX flock when available (no-op on unsupported platforms)
try:
    import fcntl  # type: ignore
except Exception:
    fcntl = None  # type: ignore


@contextmanager
def _exclusive_file_lock(lock_path: str):
    """Provide an exclusive lock on a file path. Best-effort no-op when flock is unavailable."""
    lock_dir = os.path.dirname(lock_path) or "."
    os.makedirs(lock_dir, exist_ok=True)
    lock_file = open(lock_path, "a+")
    if fcntl is not None:
        try:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        except Exception:
            # Continue without locking if flock fails
            pass
    try:
        yield
    finally:
        if fcntl is not None:
            try:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
            except Exception:
                pass
        lock_file.close()


def merge_modkit_files(
    new_files: List[str],
    existing_file: str,
    output_file: str,
    filter_bed_file: str,
    sample_id: str,
    output_dir: str,
    mnpflex_config: Optional[dict],
    num_bam_files_seen: int,
) -> None:
    """
    Merge modkit files with optimized column set and improved caching.

    This function uses only essential columns to reduce memory usage and processing time:
    - chrom, chromStart: Required for all classifiers
    - mod_code: Required for Sturgeon filtering
    - strand: Required for proper aggregation
    - valid_cov: Required for coverage calculation
    - percent_modified: Primary methylation data (required)
    - n_mod, n_canonical: Required for modification counts

    Args:
        new_files (List[str]): List of new modkit files to merge
        existing_file (str): Path to existing parquet file
        output_file (str): Path to output parquet file
        filter_bed_file (str): Path to BED file for filtering
        sample_id (str): Sample ID for organizing output files
        output_dir (str): Base output directory
        mnpflex_config (Optional[dict]): Configuration for MNP-FLEX integration
        num_bam_files_seen (int): Number of BAM files being processed
    """
    # Create sample-specific output directory
    sample_output_dir = os.path.join(output_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Define optimized schema with only essential columns
    essential_cols = [
        "chrom",
        "chromStart",
        "mod_code",
        "strand",
        "valid_cov",
        "percent_modified",
        "n_mod",
        "n_canonical",
    ]

    categorical_cols = ["chrom", "mod_code", "strand"]
    unsigned_int_cols = ["chromStart", "valid_cov", "n_mod", "n_canonical"]
    float_cols = ["percent_modified"]

    try:
        # Track cumulative BAM file count
        cumulative_bam_file_count = 0

        # Check if existing file has metadata about BAM file count
        if os.path.exists(existing_file):
            try:
                # Read existing metadata
                metadata_file = existing_file.replace(".parquet", "_metadata.json")
                if os.path.exists(metadata_file):
                    with open(metadata_file, "r") as f:
                        metadata = json.load(f)
                        cumulative_bam_file_count = metadata.get("bam_file_count", 0)
                        logging.info(
                            f"Found existing metadata with {cumulative_bam_file_count} BAM files"
                        )
            except Exception as e:
                logging.warning(f"Could not read existing metadata: {str(e)}")
                cumulative_bam_file_count = 0

        # Add the number of new BAM files being processed
        cumulative_bam_file_count += num_bam_files_seen

        logging.info(
            f"Total cumulative BAM files contributing to parquet: {cumulative_bam_file_count} (added {num_bam_files_seen} new files)"
        )

        # Cache or build PyRanges filter with improved caching
        cache_path = os.path.join(
            sample_output_dir, f"{os.path.basename(filter_bed_file)}.pgr_cache"
        )
        if os.path.exists(cache_path):
            with open(cache_path, "rb") as f:
                filter_ranges = pickle.load(f)
        else:
            comp = "gzip" if filter_bed_file.endswith(".gz") else None
            bed_df = pd.read_csv(
                filter_bed_file,
                sep="\t",
                header=None,
                names=["Chromosome", "Start", "End", "cg_label"],
                compression=comp,
                dtype={"Chromosome": str},
            )
            filter_ranges = pr.PyRanges(bed_df[["Chromosome", "Start", "End"]])
            with open(cache_path, "wb") as f:
                pickle.dump(filter_ranges, f)

        # Process new files in chunks using Polars/pandas
        new_frames = []
        for bed in new_files:
            try:
                # Read with pandas first to handle regex separator
                # Read all 18 columns from the input file
                full_cols = [
                    "chrom",
                    "chromStart",
                    "chromEnd",
                    "mod_code",
                    "score_bed",
                    "strand",
                    "thickStart",
                    "thickEnd",
                    "color",
                    "valid_cov",
                    "percent_modified",
                    "n_mod",
                    "n_canonical",
                    "n_othermod",
                    "n_delete",
                    "n_fail",
                    "n_diff",
                    "n_nocall",
                ]

                df = pd.read_csv(
                    bed,
                    sep="\s+",
                    header=None,
                    names=full_cols,
                    dtype={c: str for c in ["chrom", "mod_code", "strand", "color"]},
                )

                # Validate required columns
                missing_cols = set(full_cols) - set(df.columns)
                if missing_cols:
                    raise ValueError(f"Missing required columns: {missing_cols}")

                # Extract only essential columns for processing
                df_essential = df[essential_cols].copy()

                # Convert to Polars for efficient processing
                pl_df = pl.from_pandas(df_essential)

                # Convert numeric columns with proper error handling
                for c in unsigned_int_cols:
                    pl_df = pl_df.with_columns(pl.col(c).cast(pl.UInt32, strict=False))
                for c in float_cols:
                    pl_df = pl_df.with_columns(pl.col(c).cast(pl.Float32, strict=False))

                # Filter using PyRanges
                # Rename columns using Polars syntax
                pr_df = pl_df.rename({"chrom": "Chromosome", "chromStart": "Start"})
                # Add chromEnd for filtering (calculate from chromStart)
                pr_df = pr_df.with_columns((pl.col("Start") + 1).alias("End"))

                # Convert to pandas for PyRanges
                pr_df_pandas = pr_df.to_pandas()
                gr = pr.PyRanges(pr_df_pandas[["Chromosome", "Start", "End"]])
                inter = gr.intersect(filter_ranges).df
                filt = pd.merge(
                    inter.rename(
                        columns={
                            "Chromosome": "chrom",
                            "Start": "chromStart",
                        }
                    ),
                    pl_df.to_pandas(),
                    on=["chrom", "chromStart"],
                    how="inner",
                )
                new_frames.append(filt)
            except Exception as e:
                logging.error(f"Error processing file {bed}: {str(e)}")
                continue

        if not new_frames:
            logging.warning("No valid data to merge after filtering")
            return

        # Combine new data
        if len(new_frames) > 1:
            new_df = pd.concat(new_frames, ignore_index=True)
        else:
            new_df = new_frames[0]
        # Critical section: read/merge/write parquet guarded by an exclusive lock
        lock_path = f"{output_file}.lock"
        with _exclusive_file_lock(lock_path):
            # Use a shared StringCache during the merge to avoid costly categorical re-encodings
            with pl.StringCache():
                # If no existing file, just save the new data
                if not os.path.exists(existing_file):
                    pl_df = pl.from_pandas(new_df)
                    # Atomic write: write to temp then replace
                    fd, tmp_path = tempfile.mkstemp(prefix="lj_parquet_", suffix=".parquet", dir=os.path.dirname(output_file) or ".")
                    os.close(fd)
                    try:
                        pl_df.write_parquet(tmp_path)
                        os.replace(tmp_path, output_file)
                    finally:
                        try:
                            if os.path.exists(tmp_path):
                                os.remove(tmp_path)
                        except Exception:
                            pass

                    # Save metadata with cumulative BAM file count (atomic)
                    metadata = {
                        "bam_file_count": cumulative_bam_file_count,
                        "last_updated": datetime.now().isoformat(),
                        "sample_id": sample_id,
                        "files_added_in_this_update": num_bam_files_seen,
                        "column_format": "optimized",  # Mark as optimized format
                    }
                    metadata_file = output_file.replace(".parquet", "_metadata.json")
                    fdm, tmp_meta = tempfile.mkstemp(prefix="lj_meta_", suffix=".json", dir=os.path.dirname(metadata_file) or ".")
                    os.close(fdm)
                    try:
                        with open(tmp_meta, "w") as f:
                            json.dump(metadata, f, indent=2)
                        os.replace(tmp_meta, metadata_file)
                    finally:
                        try:
                            if os.path.exists(tmp_meta):
                                os.remove(tmp_meta)
                        except Exception:
                            pass

                    logging.info(
                        f"Created new optimized parquet file with {cumulative_bam_file_count} cumulative BAM files"
                    )
                    return

                # Process existing data in chunks
                existing_df = pl.scan_parquet(existing_file)
                # Convert new data to Polars
                pl_new_df = pl.from_pandas(new_df)

                # Convert existing data to regular DataFrame for concatenation
                existing_df = existing_df.collect()

                # Check if existing file is in old format (18 columns) or new format (8 columns)
                is_old_format = (
                    len(existing_df.columns) > 10
                )  # More than 10 columns indicates old format

                if is_old_format:
                    # Convert old format to new format by selecting only essential columns
                    logging.info("Converting existing file from old format to optimized format")
                    existing_df = existing_df.select(essential_cols)
                
                # Ensure consistent data types and column order
                for c in categorical_cols:
                    if c in existing_df.columns and c in pl_new_df.columns:
                        existing_df = existing_df.with_columns(pl.col(c).cast(pl.Categorical))
                        pl_new_df = pl_new_df.with_columns(pl.col(c).cast(pl.Categorical))

                # Ensure columns are in the same order
                pl_new_df = pl_new_df.select(existing_df.columns)

                # Validate column names match
                if set(existing_df.columns) != set(pl_new_df.columns):
                    missing_cols = set(existing_df.columns) - set(pl_new_df.columns)
                    extra_cols = set(pl_new_df.columns) - set(existing_df.columns)
                    raise ValueError(
                        f"Column mismatch: missing {missing_cols}, extra {extra_cols}"
                    )

                # Combine existing and new data
                combined = pl.concat([existing_df, pl_new_df])

                # Define aggregation expressions for essential columns only
                exprs = [
                    pl.first("mod_code"),
                    pl.first("strand"),
                    pl.mean("percent_modified").alias("percent_modified"),
                    *[pl.sum(c).alias(c) for c in ["valid_cov", "n_mod", "n_canonical"]],
                ]

                # Perform groupby and aggregation
                grouped = combined.group_by(["chrom", "chromStart"]).agg(exprs)

                # Atomic write: write to temp then replace
                fd2, tmp_out = tempfile.mkstemp(prefix="lj_parquet_", suffix=".parquet", dir=os.path.dirname(output_file) or ".")
                os.close(fd2)
                try:
                    grouped.write_parquet(tmp_out)
                    os.replace(tmp_out, output_file)
                finally:
                    try:
                        if os.path.exists(tmp_out):
                            os.remove(tmp_out)
                    except Exception:
                        pass

                # Save metadata with updated cumulative BAM file count (atomic)
                metadata = {
                    "bam_file_count": cumulative_bam_file_count,
                    "last_updated": datetime.now().isoformat(),
                    "sample_id": sample_id,
                    "files_added_in_this_update": num_bam_files_seen,
                    "column_format": "optimized",  # Mark as optimized format
                }
                metadata_file = output_file.replace(".parquet", "_metadata.json")
                fd3, tmp_meta2 = tempfile.mkstemp(prefix="lj_meta_", suffix=".json", dir=os.path.dirname(metadata_file) or ".")
                os.close(fd3)
                try:
                    with open(tmp_meta2, "w") as f:
                        json.dump(metadata, f, indent=2)
                    os.replace(tmp_meta2, metadata_file)
                finally:
                    try:
                        if os.path.exists(tmp_meta2):
                            os.remove(tmp_meta2)
                    except Exception:
                        pass

                logging.info(
                    f"Updated optimized parquet file with {cumulative_bam_file_count} cumulative BAM files (added {num_bam_files_seen} in this update)"
                )

                logging.debug(
                    f"Merged with optimized Polars and cache saved to: {output_file}"
                )

    except Exception as e:
        logging.error(f"Error in merge_modkit_files: {str(e)}")
        raise
    finally:
        # Cleanup temporary files if needed
        if os.path.exists(cache_path) and not os.path.exists(filter_bed_file):
            try:
                os.remove(cache_path)
            except Exception as e:
                logging.error(f"Error removing cache file: {str(e)}")
        gc.collect() 