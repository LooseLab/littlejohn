#!/usr/bin/env python3
"""
Master BED File Generator for ROBIN

This module generates merged master BED files from CNV and fusion breakpoint BED files.
The master BED file combines:
- CNV regions (new_file_{counter}.bed) - split into start/end breakpoints, padded by +/- 10kb
- CNV breakpoints (breakpoints_{counter}.bed) - merged as-is
- Fusion breakpoints (fusion_breakpoints_{counter}.bed) - merged as-is

All regions are intersected with the target gene panel and merged to remove overlaps.
"""

import os
import sys
import logging
import glob
import fcntl
import time
import argparse
from typing import List, Tuple, Dict, Optional
from pathlib import Path

import pandas as pd
import numpy as np

logger = logging.getLogger("robin.analysis.master_bed")


class FileLock:
    """
    Simple file-based lock for coordinating access across processes/threads.
    Uses fcntl for POSIX systems.
    """
    
    def __init__(self, lock_file: str, timeout: float = 30.0):
        self.lock_file = lock_file
        self.timeout = timeout
        self.fd = None
    
    def __enter__(self):
        self.acquire()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.release()
    
    def acquire(self):
        """Acquire the lock with timeout"""
        os.makedirs(os.path.dirname(self.lock_file), exist_ok=True)
        self.fd = open(self.lock_file, 'w')
        
        start_time = time.time()
        while True:
            try:
                fcntl.flock(self.fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                return
            except IOError:
                if time.time() - start_time > self.timeout:
                    raise TimeoutError(f"Could not acquire lock {self.lock_file} within {self.timeout}s")
                time.sleep(0.1)
    
    def release(self):
        """Release the lock"""
        if self.fd:
            fcntl.flock(self.fd, fcntl.LOCK_UN)
            self.fd.close()
            self.fd = None


def _get_latest_bed_file(bed_dir: str, pattern: str) -> Optional[str]:
    """
    Find the latest BED file matching a pattern based on counter.
    
    Args:
        bed_dir: Directory containing BED files
        pattern: Glob pattern to match (e.g., "new_file_*.bed")
        
    Returns:
        Path to the latest BED file, or None if not found
    """
    bed_files = glob.glob(os.path.join(bed_dir, pattern))
    if not bed_files:
        return None
    
    # Extract counter from filename and find the latest
    latest_file = None
    latest_counter = -1
    
    for bed_file in bed_files:
        try:
            # Extract counter from filename (e.g., "new_file_001.bed" -> 1)
            basename = os.path.basename(bed_file)
            # Pattern: name_counter.bed
            parts = basename.rsplit("_", 1)
            if len(parts) == 2:
                counter_str = parts[1].replace(".bed", "")
                counter = int(counter_str)
                if counter > latest_counter:
                    latest_counter = counter
                    latest_file = bed_file
        except (ValueError, IndexError):
            # If we can't parse the counter, use modification time as fallback
            if latest_file is None or os.path.getmtime(bed_file) > os.path.getmtime(latest_file):
                latest_file = bed_file
    
    return latest_file


def _load_bed_file(bed_path: str, require_bed6: bool = False) -> pd.DataFrame:
    """
    Load a BED file into a DataFrame, preserving strand information.
    
    Args:
        bed_path: Path to BED file
        require_bed6: If True, expects BED6 format (chrom, start, end, name, score, strand)
        
    Returns:
        DataFrame with columns: chrom, start, end, name, score, strand
    """
    if not os.path.exists(bed_path):
        return pd.DataFrame()
    
    try:
        # Read BED file - handle variable number of columns (BED format can have 3-12 columns)
        # First, read without column names to detect actual number of columns
        df_temp = pd.read_csv(
            bed_path,
            sep="\t",
            header=None,
            nrows=1,  # Just read first row to check column count
        )
        num_cols = len(df_temp.columns)
        
        # Standard BED6 column names
        bed6_cols = ["chrom", "start", "end", "name", "score", "strand"]
        
        # Read full file with appropriate column names based on actual column count
        if num_cols >= 6:
            # File has at least 6 columns (BED6 format)
            col_names = bed6_cols + [f"col{i}" for i in range(6, num_cols)]
            df = pd.read_csv(
                bed_path,
                sep="\t",
                header=None,
                names=col_names[:num_cols],
                dtype={"chrom": str, "start": int, "end": int, "name": str, "score": str, "strand": str},
                na_values=["."],
            )
            # Keep BED6 columns
            df = df[bed6_cols].copy()
        elif num_cols >= 4:
            # File has at least 4 columns (chrom, start, end, name)
            col_names = ["chrom", "start", "end", "name"] + [f"col{i}" for i in range(4, num_cols)]
            df = pd.read_csv(
                bed_path,
                sep="\t",
                header=None,
                names=col_names[:num_cols],
                dtype={"chrom": str, "start": int, "end": int, "name": str},
                na_values=["."],
            )
            # Add missing BED6 columns
            df["score"] = "0"
            df["strand"] = "."
            df = df[bed6_cols].copy()
        elif num_cols == 3:
            # File has only 3 columns (chrom, start, end)
            df = pd.read_csv(
                bed_path,
                sep="\t",
                header=None,
                names=["chrom", "start", "end"],
                dtype={"chrom": str, "start": int, "end": int},
            )
            # Add missing BED6 columns
            df["name"] = "."
            df["score"] = "0"
            df["strand"] = "."
            df = df[bed6_cols].copy()
        else:
            logger.warning(f"BED file has too few columns ({num_cols}): {bed_path}")
            return pd.DataFrame()
        
        # Fill missing values
        df["name"] = df["name"].fillna(".")
        df["score"] = df["score"].fillna("0")
        df["strand"] = df["strand"].fillna(".")
        
        # Normalize strand values
        df["strand"] = df["strand"].replace({"": ".", "?": "."})
        
        return df[bed6_cols]
    except Exception as e:
        logger.warning(f"Error loading BED file {bed_path}: {e}")
        return pd.DataFrame()


def _expand_cnv_regions(df: pd.DataFrame, expand_size: int = 10000) -> pd.DataFrame:
    """
    Convert CNV regions into breakpoint regions at start and end positions.
    Each breakpoint is padded by +/- expand_size and duplicated for both strands.
    
    Processing:
    1. Split each CNV region into two breakpoints: one at start, one at end
    2. Pad each breakpoint by +/- expand_size (creates 2*expand_size wide regions)
    3. Duplicate each breakpoint for both + and - strands
    
    Args:
        df: DataFrame with chrom, start, end, name, score, strand columns
        expand_size: Size to pad around each breakpoint in base pairs (default 10kb)
        
    Returns:
        DataFrame with breakpoint regions, duplicated for both strands
    """
    if df.empty:
        return df
    
    breakpoint_regions = []
    
    for _, row in df.iterrows():
        chrom = row["chrom"]
        region_start = row["start"]
        region_end = row["end"]
        name = row.get("name", ".")
        score = row.get("score", "0")
        
        # Create start breakpoint: pad around region_start
        start_bp_start = max(0, region_start - expand_size)
        start_bp_end = region_start + expand_size
        
        # Create end breakpoint: pad around region_end
        end_bp_start = max(0, region_end - expand_size)
        end_bp_end = region_end + expand_size
        
        # Add start breakpoint for both strands
        breakpoint_regions.append({
            "chrom": chrom,
            "start": start_bp_start,
            "end": start_bp_end,
            "name": name,
            "score": score,
            "strand": "+"
        })
        breakpoint_regions.append({
            "chrom": chrom,
            "start": start_bp_start,
            "end": start_bp_end,
            "name": name,
            "score": score,
            "strand": "-"
        })
        
        # Add end breakpoint for both strands
        breakpoint_regions.append({
            "chrom": chrom,
            "start": end_bp_start,
            "end": end_bp_end,
            "name": name,
            "score": score,
            "strand": "+"
        })
        breakpoint_regions.append({
            "chrom": chrom,
            "start": end_bp_start,
            "end": end_bp_end,
            "name": name,
            "score": score,
            "strand": "-"
        })
    
    result = pd.DataFrame(breakpoint_regions)
    return result


def _merge_overlapping_regions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping regions in a BED DataFrame, preserving strand information.
    Regions are merged separately by chromosome and strand.
    
    Args:
        df: DataFrame with chrom, start, end, name, score, strand columns
        
    Returns:
        DataFrame with merged regions
    """
    if df.empty:
        return df
    
    merged_regions = []
    
    # Sort by chromosome, strand, then start position
    df_sorted = df.sort_values(by=["chrom", "strand", "start", "end"]).copy()
    
    # Group by chromosome and strand
    for (chrom, strand), group in df_sorted.groupby(["chrom", "strand"], observed=True):
        if group.empty:
            continue
        
        # Merge overlapping regions within chromosome and strand
        current_start = None
        current_end = None
        current_name = None
        
        for _, row in group.iterrows():
            start = row["start"]
            end = row["end"]
            name = row.get("name", "merged")
            
            if current_start is None:
                # First region
                current_start = start
                current_end = end
                current_name = name
            elif start <= current_end:
                # Overlaps with current region - extend it
                current_end = max(current_end, end)
                # Combine names if different
                if current_name != name and name != "merged":
                    if current_name == "merged":
                        current_name = name
                    elif name not in current_name:
                        current_name = f"{current_name},{name}"
            else:
                # No overlap - save current region and start new one
                merged_regions.append({
                    "chrom": chrom,
                    "start": current_start,
                    "end": current_end,
                    "name": current_name if current_name else "merged",
                    "score": row.get("score", "0"),
                    "strand": strand
                })
                current_start = start
                current_end = end
                current_name = name
        
        # Don't forget the last region
        if current_start is not None:
            merged_regions.append({
                "chrom": chrom,
                "start": current_start,
                "end": current_end,
                "name": current_name if current_name else "merged",
                "score": "0",
                "strand": strand
            })
    
    if not merged_regions:
        return pd.DataFrame()
    
    return pd.DataFrame(merged_regions)


def _duplicate_unstranded_regions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Duplicate regions with unstranded ('.') strand to both + and - strands.
    
    Args:
        df: DataFrame with chrom, start, end, name, score, strand columns
        
    Returns:
        DataFrame with unstranded regions duplicated for both strands
    """
    if df.empty:
        return df
    
    # Separate stranded and unstranded regions
    stranded = df[df["strand"].isin(["+", "-"])].copy()
    unstranded = df[df["strand"] == "."].copy()
    
    if unstranded.empty:
        return stranded
    
    # Duplicate unstranded regions for both strands
    plus_strand = unstranded.copy()
    plus_strand["strand"] = "+"
    
    minus_strand = unstranded.copy()
    minus_strand["strand"] = "-"
    
    # Combine stranded regions with duplicated unstranded regions
    result = pd.concat([stranded, plus_strand, minus_strand], ignore_index=True)
    
    return result


def _intersect_with_target_genes(
    regions_df: pd.DataFrame,
    target_bed_path: str
) -> pd.DataFrame:
    """
    Filter regions to only include those that overlap with target gene BED file.
    Keeps the original regions that overlap, preserving strand information.
    
    Args:
        regions_df: DataFrame with chrom, start, end, name, score, strand columns
        target_bed_path: Path to target gene BED file
        
    Returns:
        DataFrame with regions that overlap target genes
    """
    if regions_df.empty:
        return pd.DataFrame()
    
    if not os.path.exists(target_bed_path):
        logger.warning(f"Target BED file not found: {target_bed_path}, skipping intersection")
        return regions_df
    
    try:
        # Load target genes (BED6 format)
        target_df = _load_bed_file(target_bed_path, require_bed6=True)
        if target_df.empty:
            logger.warning(f"Target BED file is empty: {target_bed_path}")
            return pd.DataFrame()
        
        # Duplicate unstranded target regions for both strands
        target_df = _duplicate_unstranded_regions(target_df)
        
        # Filter regions to only those that overlap with target genes
        overlapping_regions = []
        
        for _, region_row in regions_df.iterrows():
            region_chrom = region_row["chrom"]
            region_start = region_row["start"]
            region_end = region_row["end"]
            region_strand = region_row.get("strand", ".")
            
            # Check if this region overlaps with any target gene (same chromosome)
            overlapping_targets = target_df[
                (target_df["chrom"] == region_chrom) &
                (target_df["start"] < region_end) &
                (target_df["end"] > region_start)
            ]
            
            if not overlapping_targets.empty:
                # Keep the original region (it overlaps with target genes)
                # Preserve the region's strand information
                overlapping_regions.append({
                    "chrom": region_chrom,
                    "start": region_start,
                    "end": region_end,
                    "name": region_row.get("name", "overlap"),
                    "score": region_row.get("score", "0"),
                    "strand": region_strand
                })
        
        if not overlapping_regions:
            return pd.DataFrame()
        
        return pd.DataFrame(overlapping_regions)
        
    except Exception as e:
        logger.error(f"Error intersecting with target genes: {e}")
        return regions_df


def _load_fai_file(fai_path: str) -> Dict[str, int]:
    """
    Load a FASTA index (.fai) file and return chromosome lengths.
    
    Args:
        fai_path: Path to .fai file
        
    Returns:
        Dictionary mapping chromosome names to lengths
    """
    chrom_lengths = {}
    
    if not os.path.exists(fai_path):
        logger.warning(f"FAI file not found: {fai_path}")
        return chrom_lengths
    
    try:
        with open(fai_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # FAI format: chrom_name    length    offset    linebases    linewidth
                parts = line.split('\t')
                if len(parts) >= 2:
                    chrom_name = parts[0]
                    chrom_length = int(parts[1])
                    chrom_lengths[chrom_name] = chrom_length
    except Exception as e:
        logger.error(f"Error loading FAI file {fai_path}: {e}")
    
    return chrom_lengths


def _calculate_genome_coverage(
    bed_df: pd.DataFrame,
    fai_path: Optional[str] = None
) -> Optional[Dict[str, float]]:
    """
    Calculate the proportion of the genome covered by BED regions, accounting for strand.
    
    Args:
        bed_df: DataFrame with chrom, start, end, strand columns
        fai_path: Optional path to .fai file for genome sizes
        
    Returns:
        Dictionary with coverage statistics, or None if calculation fails
    """
    if bed_df.empty:
        return None
    
    if not fai_path or not os.path.exists(fai_path):
        logger.warning("FAI file not provided or not found - cannot calculate genome coverage")
        return None
    
    # Load chromosome lengths
    chrom_lengths = _load_fai_file(fai_path)
    if not chrom_lengths:
        logger.warning("No chromosome lengths loaded from FAI file")
        return None
    
    # Calculate total genome size (sum of all chromosome lengths)
    total_genome_size = sum(chrom_lengths.values())
    if total_genome_size == 0:
        logger.warning("Total genome size is 0")
        return None
    
    # Check for chromosomes in BED file that aren't in FAI file
    bed_chroms = set(bed_df["chrom"].unique())
    fai_chroms = set(chrom_lengths.keys())
    missing_chroms = bed_chroms - fai_chroms
    if missing_chroms:
        logger.warning(f"Chromosomes in BED file not found in FAI file: {sorted(missing_chroms)}")
    
    # Calculate coverage for each strand separately
    strand_coverage = {"+": 0, "-": 0, ".": 0}
    strand_total_genome = {"+": total_genome_size, "-": total_genome_size, ".": total_genome_size}
    
    # Group by strand
    for strand in ["+", "-", "."]:
        strand_df = bed_df[bed_df["strand"] == strand].copy()
        if strand_df.empty:
            continue
        
        # Calculate total covered bases for this strand
        # Regions are already merged, so we can sum their lengths
        covered_bases = 0
        skipped_regions = 0
        for _, row in strand_df.iterrows():
            chrom = row["chrom"]
            start = row["start"]
            end = row["end"]
            
            # Only count if chromosome exists in FAI file
            if chrom in chrom_lengths:
                # Clip region to chromosome boundaries
                region_length = max(0, min(end, chrom_lengths[chrom]) - max(start, 0))
                covered_bases += region_length
            else:
                skipped_regions += 1
        
        if skipped_regions > 0:
            logger.debug(f"Skipped {skipped_regions} {strand} strand regions due to missing chromosomes in FAI file")
        
        strand_coverage[strand] = covered_bases
    
    # Calculate proportions
    results = {
        "total_genome_size": total_genome_size,
        "plus_strand_coverage": strand_coverage["+"],
        "minus_strand_coverage": strand_coverage["-"],
        "unstranded_coverage": strand_coverage["."],
        "plus_strand_proportion": strand_coverage["+"] / strand_total_genome["+"] if strand_total_genome["+"] > 0 else 0.0,
        "minus_strand_proportion": strand_coverage["-"] / strand_total_genome["-"] if strand_total_genome["-"] > 0 else 0.0,
        "unstranded_proportion": strand_coverage["."] / strand_total_genome["."] if strand_total_genome["."] > 0 else 0.0,
    }
    
    # Calculate combined coverage
    # + and - strands are treated separately (each can have different coverage)
    # Unstranded regions (if any) are counted once
    total_covered = strand_coverage["+"] + strand_coverage["-"] + strand_coverage["."]
    # Total possible coverage: 2 * genome_size (for + and - strands separately)
    # plus unstranded regions if they exist
    total_possible_coverage = 2 * total_genome_size
    results["total_covered_bases"] = total_covered
    results["total_proportion"] = total_covered / total_possible_coverage if total_possible_coverage > 0 else 0.0
    
    return results


def _sort_bed_regions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort BED regions by chromosome (numeric order), strand, then start, then end.
    
    Args:
        df: DataFrame with chrom, start, end, strand columns
        
    Returns:
        Sorted DataFrame
    """
    if df.empty:
        return df
    
    def chrom_sort_key(chrom):
        """Sort chromosomes: chr1-22, then chrX, chrY, chrM, then others"""
        if not chrom.startswith("chr"):
            return (300, chrom)  # Non-standard chromosomes at end
        
        chrom_suffix = chrom[3:]
        if chrom_suffix == "X":
            return (100, 0)
        elif chrom_suffix == "Y":
            return (101, 0)
        elif chrom_suffix == "M":
            return (102, 0)
        else:
            try:
                return (int(chrom_suffix), 0)
            except ValueError:
                return (300, chrom_suffix)
    
    def strand_sort_key(strand):
        """Sort strands: +, then -, then ."""
        if strand == "+":
            return 0
        elif strand == "-":
            return 1
        else:
            return 2
    
    df = df.copy()
    df["_chrom_sort_key"] = df["chrom"].apply(chrom_sort_key)
    df["_strand_sort_key"] = df.get("strand", ".").apply(strand_sort_key)
    df = df.sort_values(by=["_chrom_sort_key", "_strand_sort_key", "start", "end"])
    df = df.drop(columns=["_chrom_sort_key", "_strand_sort_key"])
    
    return df.reset_index(drop=True)


def _get_target_bed_path(target_panel: Optional[str] = None) -> Optional[str]:
    """
    Get the path to the target gene BED file.
    
    Args:
        target_panel: Target panel name (rCNS2, AML, PanCan, etc.)
        
    Returns:
        Path to target BED file, or None if not found
    """
    try:
        from robin import resources
        resources_dir = os.path.dirname(resources.__file__)
        
        if target_panel == "rCNS2":
            bed_path = os.path.join(resources_dir, "rCNS2_panel_name_uniq.bed")
        elif target_panel == "AML":
            bed_path = os.path.join(resources_dir, "AML_panel_name_uniq.bed")
        elif target_panel == "PanCan":
            bed_path = os.path.join(resources_dir, "2025-03_pan-cancer_merged_20kb.bed")
        else:
            # Try custom panel
            bed_path = os.path.join(resources_dir, f"{target_panel}_panel_name_uniq.bed")
        
        if os.path.exists(bed_path):
            return bed_path
        
        # Fallback paths
        if target_panel == "rCNS2":
            fallback = "rCNS2_panel_name_uniq.bed"
        elif target_panel == "AML":
            fallback = "AML_panel_name_uniq.bed"
        elif target_panel == "PanCan":
            fallback = "2025-03_pan-cancer_merged_20kb.bed"
        else:
            fallback = f"{target_panel}_panel_name_uniq.bed"
        
        if os.path.exists(fallback):
            return fallback
            
    except ImportError:
        pass
    
    return None


def _try_get_target_panel_from_fusion_metadata(
    sample_id: str,
    work_dir: str,
) -> Optional[str]:
    """
    Try to get target_panel from fusion metadata if available.
    
    Args:
        sample_id: Sample ID
        work_dir: Working directory
        
    Returns:
        Target panel name if found, None otherwise
    """
    try:
        from robin.analysis.fusion_work import _load_fusion_metadata
        
        fusion_metadata = _load_fusion_metadata(work_dir, sample_id)
        if fusion_metadata and hasattr(fusion_metadata, 'target_panel'):
            return fusion_metadata.target_panel
    except Exception:
        pass
    
    return None


def generate_master_bed(
    sample_id: str,
    work_dir: str,
    analysis_counter: int,
    target_panel: Optional[str] = None,
    logger_instance: Optional[logging.Logger] = None,
) -> Optional[str]:
    """
    Generate master BED file from all available BED file types.
    
    This function:
    1. Loads the latest versions of new_file, breakpoints, and fusion_breakpoints BED files
    2. Converts CNV regions (new_file) into breakpoints at start/end positions, padded by +/- 10kb
    3. Merges overlapping regions
    4. Intersects with target gene panel
    5. Writes master_{counter}.bed
    
    Args:
        sample_id: Sample ID
        work_dir: Working directory
        analysis_counter: Current analysis counter
        target_panel: Target panel name (for gene intersection)
        logger_instance: Optional logger instance
        
    Returns:
        Path to generated master BED file, or None if generation failed
    """
    if logger_instance:
        log = logger_instance
    else:
        log = logger
    
    sample_dir = os.path.join(work_dir, sample_id)
    bed_dir = os.path.join(sample_dir, "bed_files")
    
    if not os.path.exists(bed_dir):
        log.debug(f"BED directory does not exist: {bed_dir}")
        return None
    
    # Get lock file path
    lock_file = os.path.join(sample_dir, "_locks", "master_bed.lock")
    
    try:
        with FileLock(lock_file, timeout=60.0):
            log.debug(f"Generating master BED file for sample {sample_id} (counter: {analysis_counter})")
            
            all_regions = []
            
            # 1. Load and convert CNV regions (new_file_{counter}.bed) to breakpoints
            cnv_bed = _get_latest_bed_file(bed_dir, "new_file_*.bed")
            if cnv_bed:
                log.debug(f"Loading CNV regions from: {cnv_bed}")
                cnv_df = _load_bed_file(cnv_bed)
                if not cnv_df.empty:
                    # Convert each region to start/end breakpoints, padded by +/- 10kb
                    expanded_df = _expand_cnv_regions(cnv_df, expand_size=10000)
                    all_regions.append(expanded_df)
                    log.debug(f"Loaded {len(cnv_df)} CNV regions (converted to {len(expanded_df)} breakpoint regions)")
            
            # 2. Load CNV breakpoints (breakpoints_{counter}.bed)
            # CNV breakpoints should have both + and - strands
            breakpoints_bed = _get_latest_bed_file(bed_dir, "breakpoints_*.bed")
            if breakpoints_bed:
                log.debug(f"Loading CNV breakpoints from: {breakpoints_bed}")
                bp_df = _load_bed_file(breakpoints_bed)
                if not bp_df.empty:
                    # Ensure CNV breakpoints have both strands
                    # Duplicate unstranded breakpoints for both strands
                    bp_df = _duplicate_unstranded_regions(bp_df)
                    all_regions.append(bp_df)
                    log.debug(f"Loaded {len(bp_df)} CNV breakpoints (with both strands)")
            
            # 3. Load fusion breakpoints (fusion_breakpoints_{counter}.bed)
            # Fusion breakpoints preserve their strand information
            fusion_bed = _get_latest_bed_file(bed_dir, "fusion_breakpoints_*.bed")
            if fusion_bed:
                log.debug(f"Loading fusion breakpoints from: {fusion_bed}")
                fusion_df = _load_bed_file(fusion_bed)
                if not fusion_df.empty:
                    # Fusion breakpoints should preserve their strand
                    # If unstranded, duplicate for both strands
                    fusion_df = _duplicate_unstranded_regions(fusion_df)
                    all_regions.append(fusion_df)
                    log.debug(f"Loaded {len(fusion_df)} fusion breakpoints (preserving strand)")
            
            if not all_regions:
                log.debug("No BED regions found to merge")
                return None
            
            # Combine all regions
            combined_df = pd.concat(all_regions, ignore_index=True)
            log.debug(f"Combined {len(combined_df)} total regions")
            
            # Merge overlapping regions
            merged_df = _merge_overlapping_regions(combined_df)
            log.debug(f"Merged to {len(merged_df)} regions")
            
            # Add target gene panel regions if available
            if target_panel:
                target_bed_path = _get_target_bed_path(target_panel)
                if target_bed_path:
                    log.debug(f"Loading target panel regions: {target_panel}")
                    target_df = _load_bed_file(target_bed_path, require_bed6=True)
                    if not target_df.empty:
                        # Duplicate unstranded target regions for both strands
                        target_df = _duplicate_unstranded_regions(target_df)
                        
                        # Keep all breakpoint regions and add target gene regions
                        log.debug(f"Keeping all {len(merged_df)} breakpoint regions and adding {len(target_df)} target gene regions")
                        
                        # Combine all breakpoints with target gene regions
                        all_regions_with_target = pd.concat([merged_df, target_df], ignore_index=True)
                        
                        log.debug(f"Combined {len(merged_df)} breakpoints with {len(target_df)} target gene regions")
                        
                        # Merge overlapping regions (target genes may overlap with breakpoints)
                        # Merging preserves strand information
                        merged_df = _merge_overlapping_regions(all_regions_with_target)
                        log.debug(f"Final merged regions including target genes: {len(merged_df)}")
                    else:
                        log.warning(f"Target panel BED file is empty: {target_bed_path}")
                else:
                    log.warning(f"Target panel BED file not found for {target_panel}")
            else:
                log.debug("No target panel specified - master BED will include all breakpoint regions")
            
            if merged_df.empty:
                log.debug("No regions remaining after processing")
                return None
            
            # Sort regions (by chrom, strand, start, end)
            sorted_df = _sort_bed_regions(merged_df)
            
            # Ensure we have all BED6 columns
            bed6_cols = ["chrom", "start", "end", "name", "score", "strand"]
            for col in bed6_cols:
                if col not in sorted_df.columns:
                    if col == "score":
                        sorted_df[col] = "0"
                    elif col == "strand":
                        sorted_df[col] = "."
                    else:
                        sorted_df[col] = "."
            
            # Write master BED file in BED6 format
            master_bed_path = os.path.join(bed_dir, f"master_{analysis_counter:03d}.bed")
            sorted_df[bed6_cols].to_csv(
                master_bed_path,
                sep="\t",
                header=False,
                index=False,
            )
            
            log.info(f"Generated master BED file: {master_bed_path} with {len(sorted_df)} regions")
            return master_bed_path
            
    except TimeoutError as e:
        log.error(f"Could not acquire lock for master BED generation: {e}")
        return None
    except Exception as e:
        log.error(f"Error generating master BED file: {e}")
        import traceback
        log.debug(f"Traceback: {traceback.format_exc()}")
        return None


def generate_master_bed_from_files(
    cnv_regions_file: Optional[str] = None,
    cnv_breakpoints_file: Optional[str] = None,
    fusion_breakpoints_file: Optional[str] = None,
    target_bed_file: Optional[str] = None,
    output_file: str = "master.bed",
    expand_cnv_size: int = 10000,
    logger_instance: Optional[logging.Logger] = None,
) -> Optional[str]:
    """
    Generate master BED file from explicit file paths (for CLI/testing).
    
    Args:
        cnv_regions_file: Path to CNV regions BED file (will be expanded by +/- expand_cnv_size)
        cnv_breakpoints_file: Path to CNV breakpoints BED file
        fusion_breakpoints_file: Path to fusion breakpoints BED file
        target_bed_file: Path to target gene panel BED file (optional)
        output_file: Path to output master BED file
        expand_cnv_size: Size to expand CNV regions in base pairs (default 10kb)
        logger_instance: Optional logger instance
        
    Returns:
        Path to generated master BED file, or None if generation failed
    """
    if logger_instance:
        log = logger_instance
    else:
        log = logger
    
    try:
        all_regions = []
        
        # 1. Load and convert CNV regions to breakpoints
        if cnv_regions_file and os.path.exists(cnv_regions_file):
            log.debug(f"Loading CNV regions from: {cnv_regions_file}")
            cnv_df = _load_bed_file(cnv_regions_file)
            if not cnv_df.empty:
                expanded_df = _expand_cnv_regions(cnv_df, expand_size=expand_cnv_size)
                all_regions.append(expanded_df)
                log.debug(f"Loaded {len(cnv_df)} CNV regions (converted to {len(expanded_df)} breakpoint regions)")
        
        # 2. Load CNV breakpoints
        if cnv_breakpoints_file and os.path.exists(cnv_breakpoints_file):
            log.debug(f"Loading CNV breakpoints from: {cnv_breakpoints_file}")
            bp_df = _load_bed_file(cnv_breakpoints_file)
            if not bp_df.empty:
                # Ensure CNV breakpoints have both strands
                bp_df = _duplicate_unstranded_regions(bp_df)
                all_regions.append(bp_df)
                log.debug(f"Loaded {len(bp_df)} CNV breakpoints (with both strands)")
        
        # 3. Load fusion breakpoints
        if fusion_breakpoints_file and os.path.exists(fusion_breakpoints_file):
            log.debug(f"Loading fusion breakpoints from: {fusion_breakpoints_file}")
            fusion_df = _load_bed_file(fusion_breakpoints_file)
            if not fusion_df.empty:
                # Fusion breakpoints preserve their strand
                fusion_df = _duplicate_unstranded_regions(fusion_df)
                all_regions.append(fusion_df)
                log.debug(f"Loaded {len(fusion_df)} fusion breakpoints (preserving strand)")
        
        if not all_regions:
            log.error("No BED regions found to merge")
            return None
        
        # Combine all regions
        combined_df = pd.concat(all_regions, ignore_index=True)
        log.debug(f"Combined {len(combined_df)} total regions")
        
        # Merge overlapping regions
        merged_df = _merge_overlapping_regions(combined_df)
        log.debug(f"Merged to {len(merged_df)} regions")
        
        # Add target gene panel regions if available
        if target_bed_file and os.path.exists(target_bed_file):
            log.debug(f"Loading target panel regions from: {target_bed_file}")
            target_df = _load_bed_file(target_bed_file, require_bed6=True)
            if not target_df.empty:
                # Duplicate unstranded target regions for both strands
                target_df = _duplicate_unstranded_regions(target_df)
                
                # Keep all breakpoint regions and add target gene regions
                log.debug(f"Keeping all {len(merged_df)} breakpoint regions and adding {len(target_df)} target gene regions")
                
                # Combine all breakpoints with target gene regions
                all_regions_with_target = pd.concat([merged_df, target_df], ignore_index=True)
                
                log.debug(f"Combined {len(merged_df)} breakpoints with {len(target_df)} target gene regions")
                
                # Merge overlapping regions
                merged_df = _merge_overlapping_regions(all_regions_with_target)
                log.debug(f"Final merged regions including target genes: {len(merged_df)}")
        
        if merged_df.empty:
            log.error("No regions remaining after processing")
            return None
        
        # Sort regions
        sorted_df = _sort_bed_regions(merged_df)
        
        # Ensure we have all BED6 columns
        bed6_cols = ["chrom", "start", "end", "name", "score", "strand"]
        for col in bed6_cols:
            if col not in sorted_df.columns:
                if col == "score":
                    sorted_df[col] = "0"
                elif col == "strand":
                    sorted_df[col] = "."
                else:
                    sorted_df[col] = "."
        
        # Write master BED file in BED6 format
        sorted_df[bed6_cols].to_csv(
            output_file,
            sep="\t",
            header=False,
            index=False,
        )
        
        log.info(f"Generated master BED file: {output_file} with {len(sorted_df)} regions")
        return output_file
        
    except Exception as e:
        log.error(f"Error generating master BED file: {e}")
        import traceback
        log.debug(f"Traceback: {traceback.format_exc()}")
        return None


def main():
    """CLI entry point for standalone testing."""
    parser = argparse.ArgumentParser(
        description="Generate master BED file from CNV and fusion breakpoint BED files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Combine CNV regions and breakpoints:
  python master_bed_generator.py --cnv-regions new_file_001.bed \\
                                  --cnv-breakpoints breakpoints_001.bed \\
                                  --output master.bed

  # Include fusion breakpoints and target panel:
  python master_bed_generator.py --cnv-regions new_file_001.bed \\
                                  --cnv-breakpoints breakpoints_001.bed \\
                                  --fusion-breakpoints fusion_breakpoints_001.bed \\
                                  --target-panel rCNS2_panel_name_uniq.bed \\
                                  --output master.bed

  # Custom CNV expansion size:
  python master_bed_generator.py --cnv-regions new_file_001.bed \\
                                  --expand-size 5000 \\
                                  --output master.bed
        """
    )
    
    parser.add_argument(
        "--cnv-regions",
        type=str,
        help=(
            "Path to CNV regions BED file (new_file_*.bed). "
            "Processing: Each CNV region is split into two breakpoints (start and end positions). "
            "Each breakpoint is padded by +/- expand-size (default 10kb) around the breakpoint position. "
            "Each breakpoint region is duplicated for both + and - strands. "
            "Regions are then merged if they overlap (separately by strand)."
        )
    )
    parser.add_argument(
        "--cnv-breakpoints",
        type=str,
        help=(
            "Path to CNV breakpoints BED file (breakpoints_*.bed). "
            "Processing: Regions are loaded as-is. If regions have unstranded ('.') strand, "
            "they are duplicated for both + and - strands. Regions are merged if they overlap (separately by strand)."
        )
    )
    parser.add_argument(
        "--fusion-breakpoints",
        type=str,
        help=(
            "Path to fusion breakpoints BED file (fusion_breakpoints_*.bed). "
            "Processing: Regions preserve their original strand information. "
            "If regions have unstranded ('.') strand, they are duplicated for both + and - strands. "
            "Regions are merged if they overlap (separately by strand)."
        )
    )
    parser.add_argument(
        "--target-panel",
        type=str,
        help=(
            "Path to target gene panel BED file (optional, for intersection). "
            "Processing: Target regions are loaded and unstranded regions are duplicated for both + and - strands. "
            "All breakpoint regions are kept (not filtered). "
            "All target gene regions are then added to the final output. "
            "Final regions are merged if they overlap (separately by strand)."
        )
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to output master BED file"
    )
    parser.add_argument(
        "--expand-size",
        type=int,
        default=10000,
        help="Size to expand CNV regions in base pairs (default: 10000 = 10kb)"
    )
    parser.add_argument(
        "--fai-file",
        type=str,
        help=(
            "Path to FASTA index (.fai) file (optional). "
            "If provided, calculates the total proportion of the genome covered by the BED file, "
            "accounting for strand information (+ and - strands are calculated separately)."
        )
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    
    # Check that at least one input file is provided
    if not any([args.cnv_regions, args.cnv_breakpoints, args.fusion_breakpoints]):
        parser.error("At least one input BED file must be provided (--cnv-regions, --cnv-breakpoints, or --fusion-breakpoints)")
    
    # Generate master BED file
    result = generate_master_bed_from_files(
        cnv_regions_file=args.cnv_regions,
        cnv_breakpoints_file=args.cnv_breakpoints,
        fusion_breakpoints_file=args.fusion_breakpoints,
        target_bed_file=args.target_panel,
        output_file=args.output,
        expand_cnv_size=args.expand_size,
    )
    
    if result:
        print(f"Successfully generated master BED file: {result}")
        
        # Calculate genome coverage if FAI file is provided
        if args.fai_file:
            try:
                # Load the generated BED file
                bed_df = _load_bed_file(result)
                if not bed_df.empty:
                    coverage_stats = _calculate_genome_coverage(bed_df, args.fai_file)
                    if coverage_stats:
                        print("\n" + "="*60)
                        print("Genome Coverage Statistics")
                        print("="*60)
                        print(f"Total genome size: {coverage_stats['total_genome_size']:,} bp")
                        print(f"\nCoverage by strand:")
                        print(f"  + strand: {coverage_stats['plus_strand_coverage']:,} bp "
                              f"({coverage_stats['plus_strand_proportion']*100:.4f}%)")
                        print(f"  - strand: {coverage_stats['minus_strand_coverage']:,} bp "
                              f"({coverage_stats['minus_strand_proportion']*100:.4f}%)")
                        if coverage_stats['unstranded_coverage'] > 0:
                            print(f"  . strand: {coverage_stats['unstranded_coverage']:,} bp "
                                  f"({coverage_stats['unstranded_proportion']*100:.4f}%)")
                        print(f"\nTotal covered bases: {coverage_stats['total_covered_bases']:,} bp")
                        print(f"Total proportion covered: {coverage_stats['total_proportion']*100:.4f}%")
                        print("="*60)
                    else:
                        print("Warning: Could not calculate genome coverage statistics", file=sys.stderr)
                else:
                    print("Warning: Generated BED file is empty, cannot calculate coverage", file=sys.stderr)
            except Exception as e:
                logger.error(f"Error calculating genome coverage: {e}")
                print(f"Warning: Error calculating genome coverage: {e}", file=sys.stderr)
        
        sys.exit(0)
    else:
        print("Failed to generate master BED file", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
