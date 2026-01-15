"""
Fusion detection and processing module for BAM files.

This module provides functionality to detect gene fusions from BAM files by analyzing
supplementary alignments and identifying reads that map to multiple gene regions.

IMPORTANT: This module uses centralized configuration from classification_config.py:
1. Only processes reads with supplementary alignments (SA tag)
2. Filters out reads where the same genomic alignment is annotated with multiple overlapping genes
3. This prevents false positives from mapping artifacts where the same genomic region is annotated with multiple overlapping genes
4. True fusions require reads to map to multiple genomic locations (supplementary alignments)
5. Fusions must be supported by a configurable minimum number of reads (default: 3) to ensure reliability and reduce false positives

All thresholds and rules are centrally managed in classification_config.py for consistency across the application.
"""

# Standard library imports
import os
import random
import logging
import json
import glob
import fcntl
import time
import bisect
from datetime import datetime
from pathlib import Path
from collections import defaultdict
from typing import Dict, Any, Optional, List, Tuple, Set, Union
from dataclasses import dataclass, asdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# Third-party imports
import numpy as np
import pandas as pd
import pysam
import networkx as nx
from sklearn.cluster import DBSCAN
import pyarrow.parquet as pq
import shutil

# Local imports
from robin.classification_config import (
    get_fusion_threshold,
    get_fusion_rule,
    validate_fusion_candidate,
    are_coordinates_similar
)

# FusionMetadata is now defined in this file

# Module-level logger
logger = logging.getLogger("robin.analysis.fusion_work")


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


@dataclass
class FusionMetadata:
    """Container for fusion analysis metadata and results"""

    sample_id: str
    file_path: str
    analysis_timestamp: float
    target_fusion_path: Optional[str] = None
    genome_wide_fusion_path: Optional[str] = None
    analysis_results: Optional[Dict] = None
    processing_steps: List[str] = None
    error_message: Optional[str] = None
    fusion_data: Optional[Dict] = None
    target_panel: Optional[str] = None

    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.analysis_results is None:
            self.analysis_results = {}
        if self.fusion_data is None:
            self.fusion_data = {}


# Minimal helpers to avoid NameError if SV pre-processing is used without utilities
def safe_read_csv(
    csv_path: str, dtype: Optional[Dict[str, Any]] = None
) -> pd.DataFrame:
    try:
        return pd.read_csv(csv_path, dtype=dtype)
    except Exception:
        return pd.DataFrame()


def get_summary(df: pd.DataFrame, min_support: int = 2) -> pd.DataFrame:
    # Placeholder summary; return empty to gracefully skip downstream when unavailable
    try:
        return pd.DataFrame()
    except Exception:
        return pd.DataFrame()


# =============================================================================
# DATA STRUCTURES
# =============================================================================


@dataclass
class GeneRegion:
    """Represents a gene region with start, end, and name."""

    start: int
    end: int
    name: str

    def overlaps_with(
        self, other_start: int, other_end: int, min_overlap: int = None
    ) -> bool:
        """Check if this region overlaps with another region."""
        if min_overlap is None:
            min_overlap = get_fusion_threshold("gene_overlap")
        overlap_start = max(self.start, other_start)
        overlap_end = min(self.end, other_end)
        return (
            overlap_end > overlap_start and (overlap_end - overlap_start) > min_overlap
        )


@dataclass
class MasterBedRegion:
    """Represents a region from the master BED file (breakpoints, CNV regions, etc.).
    Simpler than GeneRegion - just tracks coordinates without gene names.
    """

    start: int
    end: int
    name: str  # Optional: may contain region identifier from BED file

    def overlaps_with(self, other_start: int, other_end: int) -> bool:
        """Check if this region overlaps with another region (any overlap counts)."""
        overlap_start = max(self.start, other_start)
        overlap_end = min(self.end, other_end)
        return overlap_end > overlap_start


# =============================================================================
# MODULE-LEVEL CACHES
# =============================================================================

# Module-level caches for gene regions (shared across all function calls)
_gene_regions_cache: Dict[str, Dict[str, List[GeneRegion]]] = {}
_all_gene_regions_cache: Dict[str, Dict[str, List[GeneRegion]]] = {}
_target_panel_cache: Dict[str, str] = {}

# Cache for start position lists (for binary search optimization)
_gene_region_starts_cache: Dict[str, Dict[str, List[int]]] = {}
_all_gene_region_starts_cache: Dict[str, Dict[str, List[int]]] = {}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================


def overlap_fraction(a_start, a_end, b_start, b_end):
    """Return fraction of overlap with respect to the smaller of the two spans."""
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    if overlap == 0:
        return 0
    len_a = a_end - a_start
    len_b = b_end - b_start
    return overlap / min(len_a, len_b)


def _generate_random_color() -> str:
    """Generates a random color for use in plotting."""
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))


def _optimize_dataframe_memory(df: pd.DataFrame) -> pd.DataFrame:
    """Optimize DataFrame memory usage through dtype optimization."""
    # Define optimal dtypes for common columns
    dtype_optimizations = {
        "read_id": "category",
        "col4": "category",  # Gene name
        "reference_id": "category",  # Chromosome
        "strand": "category",
        "mapping_quality": np.int8,
        "reference_start": np.int32,
        "reference_end": np.int32,
        "read_start": np.int32,
        "read_end": np.int32,
        "mapping_span": np.int32,
    }

    # Apply optimizations only for columns that exist
    existing_columns = {
        col: dtype for col, dtype in dtype_optimizations.items() if col in df.columns
    }

    return df.astype(existing_columns, errors="ignore")


# =============================================================================
# FILE PATH AND BED LOADING FUNCTIONS
# =============================================================================


def _setup_file_paths(target_panel: str) -> Tuple[str, str]:
    """Setup file paths for gene data."""
    # Try to import robin resources to get the correct path
    try:
        from robin import resources

        resources_dir = os.path.dirname(resources.__file__)

        # Look for BED files in robin resources first
        if target_panel == "rCNS2":
            rCNS2_path = os.path.join(resources_dir, "rCNS2_panel_name_uniq.bed")
            if os.path.exists(rCNS2_path):
                gene_bed = rCNS2_path
            else:
                gene_bed = "rCNS2_panel_name_uniq.bed"
        elif target_panel == "AML":
            a_ml_path = os.path.join(resources_dir, "AML_panel_name_uniq.bed")
            if os.path.exists(a_ml_path):
                gene_bed = a_ml_path
            else:
                gene_bed = "AML_panel_name_uniq.bed"
        elif target_panel == "PanCan":
            pancan_path = os.path.join(resources_dir, "2025-03_pan-cancer_merged_20kb.bed")
            if os.path.exists(pancan_path):
                gene_bed = pancan_path
            else:
                gene_bed = "2025-03_pan-cancer_merged_20kb.bed"
        else:
            # Check for custom panel
            custom_panel_path = os.path.join(resources_dir, f"{target_panel}_panel_name_uniq.bed")
            if os.path.exists(custom_panel_path):
                gene_bed = custom_panel_path
            else:
                # Custom panel not found - use placeholder but don't change target_panel
                gene_bed = f"{target_panel}_panel_name_uniq.bed"

        # Look for all genes BED file
        all_genes_path = os.path.join(resources_dir, "all_genes2.bed")
        if os.path.exists(all_genes_path):
            all_gene_bed = all_genes_path
        else:
            all_gene_bed = "all_genes2.bed"

    except ImportError:
        # Fallback to local files if robin is not available
        if target_panel == "rCNS2":
            gene_bed = "rCNS2_panel_name_uniq.bed"
        elif target_panel == "AML":
            gene_bed = "AML_panel_name_uniq.bed"
        elif target_panel == "PanCan":
            gene_bed = "2025-03_pan-cancer_merged_20kb.bed"
        else:
            # Custom panel - use as is, don't fallback to rCNS2
            gene_bed = f"{target_panel}_panel_name_uniq.bed"

        all_gene_bed = "all_genes2.bed"

    return gene_bed, all_gene_bed


def _load_bed_regions(bed_file: str) -> Dict[str, List[GeneRegion]]:
    """
    Load BED file regions into memory for efficient lookup.
    Regions are sorted by start position for fast binary search.
    """
    regions = defaultdict(list)

    if not os.path.exists(bed_file):
        logger.warning(f"BED file does not exist: {bed_file}")
        return dict(regions)
    
    logger.debug(f"Loading BED file: {bed_file}")

    try:
        with open(bed_file, "r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) < 4:
                    continue

                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    gene_name = parts[3]

                    regions[chrom].append(GeneRegion(start, end, gene_name))
                except ValueError:
                    continue

    except FileNotFoundError:
        # If file not found, return empty dict
        pass
    except Exception as e:
        # If any other error, return empty dict
        logger.error(f"Error loading BED file {bed_file}: {e}")
        pass

    # Sort regions by start position for efficient binary search
    result = {}
    for chrom, region_list in regions.items():
        # Sort by start position, then by end position for consistency
        sorted_regions = sorted(region_list, key=lambda r: (r.start, r.end))
        result[chrom] = sorted_regions
    
    total_regions = sum(len(regions) for regions in result.values())
    logger.debug(f"Loaded and sorted {total_regions} total regions from {bed_file}")
    return result


def _ensure_gene_regions_loaded(target_panel: str) -> None:
    """Ensure gene regions are loaded into cache for the given target panel."""
    logger.info(f"DEBUG: _ensure_gene_regions_loaded called with target_panel='{target_panel}'")
    logger.info(f"DEBUG: Current cache keys: {list(_gene_regions_cache.keys())}")
    
    if target_panel not in _gene_regions_cache:
        logger.info(f"DEBUG: Loading gene regions for target_panel='{target_panel}'")
        gene_bed, all_gene_bed = _setup_file_paths(target_panel)
        logger.info(f"DEBUG: Resolved gene_bed='{gene_bed}', all_gene_bed='{all_gene_bed}'")

        # Load target panel gene regions
        _gene_regions_cache[target_panel] = _load_bed_regions(gene_bed)
        total_regions = sum(len(regions) for regions in _gene_regions_cache[target_panel].values())
        logger.info(f"Loaded {len(_gene_regions_cache[target_panel])} chromosomes, {total_regions} total regions for target panel {target_panel}")
        
        # Pre-compute start position lists for binary search optimization
        _gene_region_starts_cache[target_panel] = {
            chrom: [region.start for region in regions]
            for chrom, regions in _gene_regions_cache[target_panel].items()
        }

        # Load genome-wide gene regions (shared across all panels)
        if not _all_gene_regions_cache:
            _all_gene_regions_cache["shared"] = _load_bed_regions(all_gene_bed)
            total_genome_regions = sum(len(regions) for regions in _all_gene_regions_cache["shared"].values())
            logger.info(f"Loaded {len(_all_gene_regions_cache['shared'])} chromosomes, {total_genome_regions} total regions for genome-wide genes")
            
            # Pre-compute start position lists for binary search optimization
            _all_gene_region_starts_cache["shared"] = {
                chrom: [region.start for region in regions]
                for chrom, regions in _all_gene_regions_cache["shared"].items()
            }


def _load_master_bed_regions(bed_file: str) -> Dict[str, List[MasterBedRegion]]:
    """
    Load master BED file regions into memory for efficient lookup.
    Master BED regions don't require gene overlap - breaks can occur anywhere.
    
    Args:
        bed_file: Path to master BED file
        
    Returns:
        Dictionary mapping chromosome names to lists of MasterBedRegion objects
    """
    regions = defaultdict(list)
    
    if not os.path.exists(bed_file):
        logger.debug(f"Master BED file does not exist: {bed_file}")
        return dict(regions)
    
    logger.info(f"Loading master BED file: {bed_file}")
    
    try:
        with open(bed_file, "r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                
                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3] if len(parts) > 3 else f"region_{line_num}"
                    
                    regions[chrom].append(MasterBedRegion(start, end, name))
                except (ValueError, IndexError) as e:
                    logger.debug(f"Skipping invalid line {line_num} in {bed_file}: {e}")
                    continue
    
    except Exception as e:
        logger.error(f"Error loading master BED file {bed_file}: {e}")
        return dict(regions)
    
    result = dict(regions)
    total_regions = sum(len(regions) for regions in result.values())
    logger.info(f"Loaded {total_regions} total regions from master BED file")
    return result


def _get_master_bed_path(work_dir: str, sample_id: str) -> Optional[str]:
    """
    Get the path to the master BED file for a sample.
    
    Args:
        work_dir: Working directory
        sample_id: Sample ID
        
    Returns:
        Path to master BED file, or None if not found
    """
    try:
        analysis_counter = _load_analysis_counter(sample_id, work_dir)
        sample_dir = os.path.join(work_dir, sample_id)
        bed_dir = os.path.join(sample_dir, "bed_files")
        master_bed_path = os.path.join(bed_dir, f"master_{analysis_counter:03d}.bed")
        
        if os.path.exists(master_bed_path):
            return master_bed_path
        
        # Try to find the latest master BED file if counter-based doesn't exist
        if os.path.exists(bed_dir):
            master_bed_files = glob.glob(os.path.join(bed_dir, "master_*.bed"))
            if master_bed_files:
                # Sort by modification time and return the latest
                latest = max(master_bed_files, key=os.path.getmtime)
                logger.debug(f"Using latest master BED file: {latest}")
                return latest
    except Exception as e:
        logger.debug(f"Error finding master BED file: {e}")
    
    return None


# =============================================================================
# BAM FILE ANALYSIS FUNCTIONS
# =============================================================================


def has_supplementary_alignments(bamfile: str) -> bool:
    """Quickly check if a BAM file has supplementary alignments."""
    try:
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            for read in bam:
                if read.has_tag("SA"):
                    return True
                if read.is_supplementary:
                    return True
        return False
    except Exception as e:
        logging.error(f"Error checking supplementary alignments in {bamfile}: {str(e)}")
        return False


def find_reads_mapping_to_master_bed(
    bamfile: str, master_bed_regions: Dict[str, List[MasterBedRegion]]
) -> Set[str]:
    """
    Find all read names that map to regions in the master BED file.
    
    Args:
        bamfile: Path to BAM file
        master_bed_regions: Dictionary of master BED regions by chromosome
        
    Returns:
        Set of read names that map to master BED regions
    """
    reads_mapping_to_master = set()
    
    if not master_bed_regions:
        return reads_mapping_to_master
    
    try:
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            for read in bam:
                if read.is_unmapped:
                    continue
                
                # Get reference information
                ref_name = (
                    bam.get_reference_name(read.reference_id)
                    if read.reference_id >= 0
                    else None
                )
                if not ref_name or ref_name == "chrM":
                    continue
                
                # Check if this read overlaps with any master BED region
                if ref_name in master_bed_regions:
                    ref_start = read.reference_start
                    ref_end = read.reference_end
                    
                    for region in master_bed_regions[ref_name]:
                        if region.overlaps_with(ref_start, ref_end):
                            reads_mapping_to_master.add(read.query_name)
                            break  # No need to check other regions for this read
    
    except Exception as e:
        logger.error(f"Error finding reads mapping to master BED in {bamfile}: {e}")
        raise
    
    logger.debug(f"Found {len(reads_mapping_to_master)} reads mapping to master BED regions")
    return reads_mapping_to_master


def find_reads_with_supplementary(bamfile: str) -> Set[str]:
    """
    Find all read names that have supplementary alignments.

    Args:
        bamfile: Path to BAM file

    Returns:
        Set of read names with supplementary alignments
    """
    reads_with_supp = set()

    try:
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            for read in bam:
                if read.is_unmapped:
                    continue
                # Check for supplementary alignments or SA tag
                if read.is_supplementary or read.has_tag("SA"):
                    reads_with_supp.add(read.query_name)

    except Exception as e:
        logger.error(
            f"Error finding reads with supplementary alignments in {bamfile}: {e}"
        )
        raise

    logger.debug(f"Found {len(reads_with_supp)} reads with supplementary alignments")
    return reads_with_supp


def _check_read_alignments_overlap(read_rows: List[Dict]) -> bool:
    """
    Check for false positives in read alignments:
    - Same genomic alignment annotated with multiple genes (overlapping gene regions)
    - Very similar (but not identical) alignments that are likely mapping artifacts
    
    Args:
        read_rows: List of alignment dictionaries for the same read
        
    Returns:
        True if any false positive pattern is detected, False otherwise
    """
    if len(read_rows) < 2:
        return False
    
    # Check: Same genomic alignment annotated with multiple genes
    # Group by genomic coordinates (chr:start-end)
    genomic_alignments = {}
    for row in read_rows:
        genomic_key = f"{row['reference_id']}:{row['reference_start']}-{row['reference_end']}"
        if genomic_key not in genomic_alignments:
            genomic_alignments[genomic_key] = []
        genomic_alignments[genomic_key].append(row)
    
    # Only filter out if we have the EXACT same genomic coordinates with different genes
    # This is more restrictive - only remove true mapping artifacts
    for genomic_key, alignments in genomic_alignments.items():
        if len(alignments) > 1:
            # Check if all alignments have the same gene annotation
            genes = [align['col4'] for align in alignments]
            if len(set(genes)) > 1:
                # Same genomic coordinates with different gene annotations = false positive
                return True
    
    # Check: Very similar alignments (coordinate similarity filter)
    if get_fusion_rule("coordinate_similarity_filter"):
        max_diff = get_fusion_threshold("coordinate_similarity")
        
        # Optimize: Group by chromosome first to reduce comparisons
        # Only compare alignments on the same chromosome
        alignments_by_chrom = defaultdict(list)
        for i, row in enumerate(read_rows):
            alignments_by_chrom[row['reference_id']].append((i, row))
        
        # Compare pairs within each chromosome (reduces comparisons significantly)
        for chrom, chrom_alignments in alignments_by_chrom.items():
            if len(chrom_alignments) < 2:
                continue
            
            # Only compare pairs on the same chromosome
            for i in range(len(chrom_alignments)):
                for j in range(i + 1, len(chrom_alignments)):
                    _, row1 = chrom_alignments[i]
                    _, row2 = chrom_alignments[j]
                    
                    if are_coordinates_similar(
                        row1['reference_start'], row1['reference_end'],
                        row2['reference_start'], row2['reference_end'],
                        max_diff
                    ):
                        # Very similar alignments detected - likely mapping artifact
                        logger.debug(f"Filtering similar alignments: {row1['reference_id']}:{row1['reference_start']}-{row1['reference_end']} vs {row2['reference_id']}:{row2['reference_start']}-{row2['reference_end']}")
                        return True
    
    return False


def _find_gene_intersections(
    read: pysam.AlignedSegment,
    ref_name: str,
    ref_start: int,
    ref_end: int,
    gene_regions: List[GeneRegion],
    region_starts: Optional[List[int]] = None,
) -> List[Dict]:
    """
    Find intersections between a read and gene regions using optimized binary search.
    Only processes reads that have supplementary alignments (true fusion candidates).
    
    Uses binary search on sorted gene regions for O(log n + k) complexity where:
    - n = number of gene regions
    - k = number of overlapping regions (typically small)

    Args:
        read: Pysam aligned segment
        ref_name: Reference chromosome name
        ref_start: Reference start position
        ref_end: Reference end position
        gene_regions: List of gene regions on this chromosome (must be sorted by start position)
        region_starts: Optional pre-computed list of start positions (for performance)

    Returns:
        List of intersection dictionaries with the exact column structure
    """
    read_rows = []

    # Only process reads that have supplementary alignments (SA tag)
    # This ensures we only look at reads that actually map to multiple locations
    if not read.has_tag("SA"):
        return read_rows

    # Early return if no gene regions
    if not gene_regions:
        return read_rows

    # Binary search optimization:
    # 1. Find the first region that could overlap (regions with start <= ref_end)
    # 2. Iterate backwards checking overlaps until we pass ref_start
    
    # Use cached start positions if provided, otherwise create on-the-fly
    if region_starts is None:
        region_starts = [region.start for region in gene_regions]
    
    # Find the rightmost region with start <= ref_end
    # This is the last region that could potentially overlap
    rightmost_idx = bisect.bisect_right(region_starts, ref_end)
    
    # Now iterate backwards from rightmost_idx to find all overlapping regions
    # We iterate backwards because regions are sorted by start, and we want to find
    # all regions that overlap with [ref_start, ref_end]
    for i in range(rightmost_idx - 1, -1, -1):
        gene_region = gene_regions[i]
        
        # If this region ends before ref_start, we've gone too far (no more overlaps)
        # Since regions are sorted by start, all previous regions will also end before ref_start
        if gene_region.end < ref_start:
            break
        
        # Check if this region overlaps with the read
        if gene_region.overlaps_with(ref_start, ref_end):
            # Create the EXACT column structure expected by the original code
            read_rows.append(
                {
                    "col1": ref_name,  # Chromosome
                    "col2": gene_region.start,  # Gene start
                    "col3": gene_region.end,  # Gene end
                    "col4": gene_region.name,  # Gene name
                    "reference_id": ref_name,  # Chromosome (duplicate for compatibility)
                    "reference_start": ref_start,  # Read start position
                    "reference_end": ref_end,  # Read end position
                    "read_id": read.query_name,  # Read ID
                    "mapping_quality": read.mapping_quality,  # Mapping quality
                    "strand": "-" if read.is_reverse else "+",  # Strand
                    "read_start": read.query_alignment_start,  # Read alignment start
                    "read_end": read.query_alignment_end,  # Read alignment end
                    "is_secondary": read.is_secondary,  # Is secondary alignment
                    "is_supplementary": read.is_supplementary,  # Is supplementary alignment
                    "mapping_span": ref_end - ref_start,  # Mapping span
                }
            )

    return read_rows


def _process_reads_for_fusions(
    bamfile: str,
    reads_with_supp: Set[str],
    gene_regions: Dict[str, List[GeneRegion]],
) -> Optional[pd.DataFrame]:
    """
    Process reads to find gene intersections and create fusion candidates.
    Now includes overlap filtering to remove false positives.

    Args:
        bamfile: Path to BAM file
        reads_with_supp: Set of read names with supplementary alignments
        gene_regions: Dictionary of gene regions by chromosome

    Returns:
        DataFrame with fusion candidates or None if no candidates found
    """
    # Collect all alignments per read first
    read_alignments = {}  # read_id -> list of alignment rows

    try:
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            for read in bam:
                # Skip if not in our target reads
                if read.query_name not in reads_with_supp:
                    continue

                # Skip secondary alignments and unmapped reads
                if read.is_secondary or read.is_unmapped:
                    continue

                # Get reference information
                ref_name = (
                    bam.get_reference_name(read.reference_id)
                    if read.reference_id >= 0
                    else None
                )
                if not ref_name or ref_name == "chrM":
                    continue

                ref_start = read.reference_start
                ref_end = read.reference_end

                # Check if this read intersects with any gene regions
                if ref_name in gene_regions:
                    read_rows = _find_gene_intersections(
                        read, ref_name, ref_start, ref_end, gene_regions[ref_name]
                    )
                    if read_rows:
                        read_id = read.query_name
                        if read_id not in read_alignments:
                            read_alignments[read_id] = []
                        read_alignments[read_id].extend(read_rows)

    except Exception as e:
        logger.error(f"Error processing reads for fusions: {str(e)}")
        raise

    # Now filter out reads with overlapping alignments
    filtered_rows = []
    for read_id, alignments in read_alignments.items():
        # Check for overlaps between alignments of the same read
        if not _check_read_alignments_overlap(alignments):
            filtered_rows.extend(alignments)

    if not filtered_rows:
        return None

    # Create DataFrame with the EXACT column structure expected by the original code
    df = pd.DataFrame(filtered_rows)

    # Apply memory optimizations
    df = _optimize_fusion_dataframe(df)

    # Apply filtering thresholds using centralized configuration
    min_mq = get_fusion_threshold("mapping_quality")
    min_span = get_fusion_threshold("mapping_span")
    df = df[(df["mapping_quality"] > min_mq) & (df["mapping_span"] > min_span)].reset_index(
        drop=True
    )

    return df if not df.empty else None


def _optimize_fusion_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Optimize fusion DataFrame memory usage.

    Args:
        df: Input DataFrame

    Returns:
        Memory-optimized DataFrame
    """
    # Use categorical dtypes for string columns
    string_columns = ["col1", "col4", "reference_id", "strand", "read_id"]
    for col in string_columns:
        if col in df.columns:
            df[col] = df[col].astype("category")

    # Use appropriate integer types
    int_columns = [
        "col2",
        "col3",
        "reference_start",
        "reference_end",
        "read_start",
        "read_end",
        "mapping_quality",
        "mapping_span",
    ]
    for col in int_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    return df


def _process_reads_for_master_bed_fusions(
    bamfile: str,
    reads_mapping_to_master: Set[str],
) -> Optional[pd.DataFrame]:
    """
    Process reads mapping to master BED regions and track ALL their supplementary alignments.
    Unlike regular fusion detection, we don't require gene overlap - breaks can occur anywhere.
    
    Args:
        bamfile: Path to BAM file
        reads_mapping_to_master: Set of read names that map to master BED regions
        
    Returns:
        DataFrame with fusion candidates or None if no candidates found
    """
    read_rows = []
    
    if not reads_mapping_to_master:
        return None
    
    try:
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            for read in bam:
                # Only process reads that map to master BED regions
                if read.query_name not in reads_mapping_to_master:
                    continue
                
                # Skip secondary alignments and unmapped reads
                if read.is_secondary or read.is_unmapped:
                    continue
                
                # Only process reads with supplementary alignments (SA tag)
                if not read.has_tag("SA"):
                    continue
                
                # Get reference information
                ref_name = (
                    bam.get_reference_name(read.reference_id)
                    if read.reference_id >= 0
                    else None
                )
                if not ref_name or ref_name == "chrM":
                    continue
                
                ref_start = read.reference_start
                ref_end = read.reference_end
                
                # For master BED fusions, we track ALL supplementary alignments
                # without requiring gene overlap. Create a row for each alignment.
                read_rows.append(
                    {
                        "col1": ref_name,  # Chromosome
                        "col2": ref_start,  # Region start (read start)
                        "col3": ref_end,  # Region end (read end)
                        "col4": "master_bed_region",  # Identifier
                        "reference_id": ref_name,  # Chromosome
                        "reference_start": ref_start,  # Read start position
                        "reference_end": ref_end,  # Read end position
                        "read_id": read.query_name,  # Read ID
                        "mapping_quality": read.mapping_quality,  # Mapping quality
                        "strand": "-" if read.is_reverse else "+",  # Strand
                        "read_start": read.query_alignment_start,  # Read alignment start
                        "read_end": read.query_alignment_end,  # Read alignment end
                        "is_secondary": read.is_secondary,  # Is secondary alignment
                        "is_supplementary": read.is_supplementary,  # Is supplementary alignment
                        "mapping_span": ref_end - ref_start,  # Mapping span
                    }
                )
                
                # Also parse SA tag to get all supplementary alignments
                try:
                    sa_tag = read.get_tag("SA")
                    # SA tag format: "chr,pos,strand,CIGAR,mapQ,NM;chr,pos,strand,CIGAR,mapQ,NM;..."
                    sa_entries = sa_tag.split(";")
                    for sa_entry in sa_entries:
                        if not sa_entry:
                            continue
                        sa_parts = sa_entry.split(",")
                        if len(sa_parts) >= 3:
                            sa_chrom = sa_parts[0]
                            sa_pos = int(sa_parts[1])
                            sa_strand = sa_parts[2]
                            sa_mapq = int(sa_parts[4]) if len(sa_parts) > 4 else 0
                            
                            # Estimate end position from CIGAR if available
                            if len(sa_parts) > 3:
                                cigar_str = sa_parts[3]
                                # Simple estimate: use read length as span (rough approximation)
                                sa_end = sa_pos + read.query_length
                            else:
                                sa_end = sa_pos + 100  # Default small span
                            
                            # Skip mitochondrial chromosomes
                            if sa_chrom == "chrM" or sa_chrom == "M":
                                continue
                            
                            # Add supplementary alignment as a row
                            read_rows.append(
                                {
                                    "col1": sa_chrom,  # Chromosome
                                    "col2": sa_pos,  # Region start
                                    "col3": sa_end,  # Region end
                                    "col4": "master_bed_supplementary",  # Identifier
                                    "reference_id": sa_chrom,  # Chromosome
                                    "reference_start": sa_pos,  # Alignment start position
                                    "reference_end": sa_end,  # Alignment end position
                                    "read_id": read.query_name,  # Read ID
                                    "mapping_quality": sa_mapq,  # Mapping quality from SA tag
                                    "strand": sa_strand,  # Strand
                                    "read_start": 0,  # Not available from SA tag
                                    "read_end": read.query_length,  # Not available from SA tag
                                    "is_secondary": False,  # SA tag entries are supplementary
                                    "is_supplementary": True,  # This is a supplementary alignment
                                    "mapping_span": sa_end - sa_pos,  # Mapping span
                                }
                            )
                except (ValueError, KeyError) as e:
                    logger.debug(f"Could not parse SA tag for read {read.query_name}: {e}")
                    continue
    
    except Exception as e:
        logger.error(f"Error processing reads for master BED fusions: {str(e)}")
        raise
    
    if not read_rows:
        return None
    
    # Create DataFrame
    df = pd.DataFrame(read_rows)
    
    # Apply memory optimizations
    df = _optimize_fusion_dataframe(df)
    
    # Apply basic filtering thresholds (mapping quality, mapping span)
    # but NOT gene overlap requirement
    min_mq = get_fusion_threshold("mapping_quality")
    min_span = get_fusion_threshold("mapping_span")
    df = df[(df["mapping_quality"] > min_mq) & (df["mapping_span"] > min_span)].reset_index(
        drop=True
    )
    
    if df.empty:
        return None
    
    return df


def _filter_fusion_candidates(df: pd.DataFrame) -> Optional[pd.DataFrame]:
    """
    Filter fusion candidates to only include reads with supplementary alignments.

    Args:
        df: DataFrame with fusion candidates

    Returns:
        Filtered DataFrame or None if no candidates
    """
    if df is None or df.empty:
        return None

    try:
        # Count unique genes per read_id to find fusion candidates
        gene_counts = df.groupby("read_id", observed=True)["col4"].nunique()
        
        # Filter for reads that map to more than 1 gene (fusion candidates)
        fusion_read_ids = gene_counts[gene_counts > 1].index
        
        if len(fusion_read_ids) == 0:
            return None

        # Return all rows for reads that map to multiple genes
        result = df[df["read_id"].isin(fusion_read_ids)]

        if result.empty:
            return None

        return result.reset_index(drop=True)

    except Exception as e:
        logger.error(f"Error filtering fusion candidates: {str(e)}")
        return None


# =============================================================================
# MAIN PROCESSING FUNCTIONS
# =============================================================================


def process_bam_for_master_bed_fusions(
    bamfile: str, work_dir: str, sample_id: str, supplementary_read_ids: Optional[List[str]] = None
) -> Optional[pd.DataFrame]:
    """
    Process BAM file to find fusion candidates by tracking ALL supplementary alignments.
    No longer requires primary alignment to overlap master BED regions.
    Tracks supplementary mappings from all reads with supplementary alignments,
    without requiring gene overlap (breaks can occur anywhere).
    
    Args:
        bamfile: Path to BAM file
        work_dir: Working directory (kept for API compatibility, not used for filtering)
        sample_id: Sample ID (kept for API compatibility, not used for filtering)
        supplementary_read_ids: Optional list of read IDs with supplementary alignments
                               (if available from metadata, avoids checking SA tag for all reads)
        
    Returns:
        DataFrame with master BED fusion candidates or None if no candidates found
    """
    try:
        # Convert supplementary_read_ids to set for fast lookup if provided
        supplementary_reads_set = set(supplementary_read_ids) if supplementary_read_ids else None
        
        # Process ALL reads with supplementary alignments (no master BED overlap requirement)
        read_rows = []
        reads_with_supplementary_count = 0
        
        try:
            with pysam.AlignmentFile(bamfile, "rb") as bam:
                for read in bam:
                    if read.is_unmapped:
                        continue
                    
                    # Skip secondary alignments (we only process primary alignments)
                    if read.is_secondary:
                        continue
                    
                    # Get reference information
                    ref_name = (
                        bam.get_reference_name(read.reference_id)
                        if read.reference_id >= 0
                        else None
                    )
                    if not ref_name or ref_name == "chrM":
                        continue
                    
                    # Check for supplementary alignments
                    # Always check SA tag directly to ensure we catch ALL reads with supplementary alignments
                    # The supplementary_read_ids list may be incomplete if preprocessing missed some reads
                    # or if reads were added after preprocessing
                    has_supplementary = read.has_tag("SA")
                    
                    # Optional optimization: if we have a supplementary_read_ids list and the read is not in it,
                    # we can skip the SA tag check (but this risks missing reads if the list is incomplete)
                    # For now, we always check SA tag to be safe
                    if not has_supplementary:
                        continue
                    
                    # Process this read: it has supplementary alignments
                    reads_with_supplementary_count += 1
                    
                    ref_start = read.reference_start
                    ref_end = read.reference_end
                    
                    # Process this read: add primary alignment
                    read_rows.append(
                        {
                            "col1": ref_name,  # Chromosome
                            "col2": ref_start,  # Region start (read start)
                            "col3": ref_end,  # Region end (read end)
                            "col4": "master_bed_region",  # Identifier
                            "reference_id": ref_name,  # Chromosome
                            "reference_start": ref_start,  # Read start position
                            "reference_end": ref_end,  # Read end position
                            "read_id": read.query_name,  # Read ID
                            "mapping_quality": read.mapping_quality,  # Mapping quality
                            "strand": "-" if read.is_reverse else "+",  # Strand
                            "read_start": read.query_alignment_start,  # Read alignment start
                            "read_end": read.query_alignment_end,  # Read alignment end
                            "is_secondary": read.is_secondary,  # Is secondary alignment
                            "is_supplementary": read.is_supplementary,  # Is supplementary alignment
                            "mapping_span": ref_end - ref_start,  # Mapping span
                        }
                    )
                    
                    # Parse SA tag to get all supplementary alignments
                    if read.has_tag("SA"):
                        try:
                            sa_tag = read.get_tag("SA")
                            # SA tag format: "chr,pos,strand,CIGAR,mapQ,NM;chr,pos,strand,CIGAR,mapQ,NM;..."
                            sa_entries = sa_tag.split(";")
                            for sa_entry in sa_entries:
                                if not sa_entry:
                                    continue
                                sa_parts = sa_entry.split(",")
                                if len(sa_parts) >= 3:
                                    sa_chrom = sa_parts[0]
                                    sa_pos = int(sa_parts[1])
                                    sa_strand = sa_parts[2]
                                    sa_mapq = int(sa_parts[4]) if len(sa_parts) > 4 else 0
                                    
                                    # Estimate end position from CIGAR if available
                                    if len(sa_parts) > 3:
                                        cigar_str = sa_parts[3]
                                        # Simple estimate: use read length as span (rough approximation)
                                        sa_end = sa_pos + read.query_length
                                    else:
                                        sa_end = sa_pos + 100  # Default small span
                                    
                                    # Skip mitochondrial chromosomes
                                    if sa_chrom == "chrM" or sa_chrom == "M":
                                        continue
                                    
                                    # Add supplementary alignment as a row
                                    read_rows.append(
                                        {
                                            "col1": sa_chrom,  # Chromosome
                                            "col2": sa_pos,  # Region start
                                            "col3": sa_end,  # Region end
                                            "col4": "master_bed_supplementary",  # Identifier
                                            "reference_id": sa_chrom,  # Chromosome
                                            "reference_start": sa_pos,  # Alignment start position
                                            "reference_end": sa_end,  # Alignment end position
                                            "read_id": read.query_name,  # Read ID
                                            "mapping_quality": sa_mapq,  # Mapping quality from SA tag
                                            "strand": sa_strand,  # Strand
                                            "read_start": 0,  # Not available from SA tag
                                            "read_end": read.query_length,  # Not available from SA tag
                                            "is_secondary": False,  # SA tag entries are supplementary
                                            "is_supplementary": True,  # This is a supplementary alignment
                                            "mapping_span": sa_end - sa_pos,  # Mapping span
                                        }
                                    )
                        except (ValueError, KeyError) as e:
                            logger.debug(f"Could not parse SA tag for read {read.query_name}: {e}")
                            continue
        
        except Exception as e:
            logger.error(f"Error processing BAM file for master BED fusions: {e}")
            raise
        
        if reads_with_supplementary_count == 0:
            logger.debug("No reads with supplementary alignments found")
            return None
        
        logger.info(f"Found {reads_with_supplementary_count} reads with supplementary alignments")
        
        if not read_rows:
            logger.debug("No master BED fusion candidates found (no supplementary alignments)")
            return None
        
        # Create DataFrame
        df = pd.DataFrame(read_rows)
        
        # Apply memory optimizations
        df = _optimize_fusion_dataframe(df)
        
        # Apply basic filtering thresholds (mapping quality, mapping span)
        # but NOT gene overlap requirement
        min_mq = get_fusion_threshold("mapping_quality")
        min_span = get_fusion_threshold("mapping_span")
        
        before_filter = len(df)
        df = df[(df["mapping_quality"] > min_mq) & (df["mapping_span"] > min_span)].reset_index(
            drop=True
        )
        
        if df.empty:
            return None
        
        logger.info(f"Found {len(df)} master BED fusion candidates")
        return df
        
    except Exception as e:
        logger.error(f"Error processing master BED fusions: {str(e)}")
        logger.error("Exception details:", exc_info=True)
        return None


def process_bam_single_pass(
    bamfile: str,
    target_panel: str,
    work_dir: Optional[str] = None,
    sample_id: Optional[str] = None,
    supplementary_read_ids: Optional[List[str]] = None,
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Process BAM file in a single pass to find all fusion candidates.
    
    This optimized function combines three separate BAM file passes into one:
    1. Finding reads with supplementary alignments
    2. Processing target panel and genome-wide fusions
    3. Processing master BED fusions
    
    This provides a 2-3x speedup by eliminating redundant BAM file I/O.

    Args:
        bamfile: Path to BAM file
        target_panel: Target panel to use (rCNS2, AML, or PanCan)
        work_dir: Working directory (optional, for master BED processing)
        sample_id: Sample ID (optional, for master BED processing)
        supplementary_read_ids: Optional list of read IDs with supplementary alignments
                               (if available, can skip SA tag check for some reads)

    Returns:
        Tuple of (target_panel_candidates, genome_wide_candidates, master_bed_candidates) DataFrames
    """
    try:
        logger.debug(f"Processing BAM file in single pass: {bamfile}")
        
        # Ensure gene regions are loaded (cached, so this is fast after first call)
        _ensure_gene_regions_loaded(target_panel)
        
        # Get gene region dictionaries
        target_regions = _gene_regions_cache.get(target_panel, {})
        if "shared" not in _all_gene_regions_cache:
            logger.error("'shared' key not found in _all_gene_regions_cache!")
            genome_regions = {}
        else:
            genome_regions = _all_gene_regions_cache["shared"]
        
        # Get cached start position lists for binary search optimization
        target_region_starts = _gene_region_starts_cache.get(target_panel, {})
        genome_region_starts = _all_gene_region_starts_cache.get("shared", {})
        
        # Convert supplementary_read_ids to set for fast lookup if provided
        supplementary_reads_set = set(supplementary_read_ids) if supplementary_read_ids else None
        
        # Collectors for all three fusion types
        target_read_alignments = {}  # read_id -> list of alignment rows
        genome_read_alignments = {}  # read_id -> list of alignment rows
        master_bed_rows = []  # List of master BED alignment rows
        
        reads_with_supplementary_count = 0
        
        # Single pass through BAM file
        try:
            with pysam.AlignmentFile(bamfile, "rb") as bam:
                for read in bam:
                    # Skip unmapped reads
                    if read.is_unmapped:
                        continue
                    
                    # Skip secondary alignments (we only process primary alignments)
                    if read.is_secondary:
                        continue
                    
                    # Get reference information
                    ref_name = (
                        bam.get_reference_name(read.reference_id)
                        if read.reference_id >= 0
                        else None
                    )
                    if not ref_name or ref_name == "chrM":
                        continue
                    
                    # Check for supplementary alignments
                    # Always check SA tag directly to ensure we catch ALL reads
                    has_supplementary = read.has_tag("SA")
                    
                    # Optional optimization: if we have a supplementary_read_ids list and the read is not in it,
                    # we can skip the SA tag check (but this risks missing reads if the list is incomplete)
                    if supplementary_reads_set is not None and read.query_name not in supplementary_reads_set:
                        # Skip if not in the list (but still check SA tag to be safe)
                        if not has_supplementary:
                            continue
                    
                    if not has_supplementary:
                        continue
                    
                    # This read has supplementary alignments - process it for all fusion types
                    reads_with_supplementary_count += 1
                    read_id = read.query_name
                    ref_start = read.reference_start
                    ref_end = read.reference_end
                    
                    # 1. Process for target panel fusions
                    if ref_name in target_regions:
                        target_starts = target_region_starts.get(ref_name)
                        target_rows = _find_gene_intersections(
                            read, ref_name, ref_start, ref_end, target_regions[ref_name], target_starts
                        )
                        if target_rows:
                            if read_id not in target_read_alignments:
                                target_read_alignments[read_id] = []
                            target_read_alignments[read_id].extend(target_rows)
                    
                    # 2. Process for genome-wide fusions
                    if ref_name in genome_regions:
                        genome_starts = genome_region_starts.get(ref_name)
                        genome_rows = _find_gene_intersections(
                            read, ref_name, ref_start, ref_end, genome_regions[ref_name], genome_starts
                        )
                        if genome_rows:
                            if read_id not in genome_read_alignments:
                                genome_read_alignments[read_id] = []
                            genome_read_alignments[read_id].extend(genome_rows)
                    
                    # 3. Process for master BED fusions
                    # Add primary alignment
                    master_bed_rows.append(
                        {
                            "col1": ref_name,
                            "col2": ref_start,
                            "col3": ref_end,
                            "col4": "master_bed_region",
                            "reference_id": ref_name,
                            "reference_start": ref_start,
                            "reference_end": ref_end,
                            "read_id": read_id,
                            "mapping_quality": read.mapping_quality,
                            "strand": "-" if read.is_reverse else "+",
                            "read_start": read.query_alignment_start,
                            "read_end": read.query_alignment_end,
                            "is_secondary": read.is_secondary,
                            "is_supplementary": read.is_supplementary,
                            "mapping_span": ref_end - ref_start,
                        }
                    )
                    
                    # Parse SA tag to get all supplementary alignments
                    if read.has_tag("SA"):
                        try:
                            sa_tag = read.get_tag("SA")
                            # SA tag format: "chr,pos,strand,CIGAR,mapQ,NM;chr,pos,strand,CIGAR,mapQ,NM;..."
                            sa_entries = sa_tag.split(";")
                            for sa_entry in sa_entries:
                                if not sa_entry:
                                    continue
                                sa_parts = sa_entry.split(",")
                                if len(sa_parts) >= 3:
                                    sa_chrom = sa_parts[0]
                                    sa_pos = int(sa_parts[1])
                                    sa_strand = sa_parts[2]
                                    sa_mapq = int(sa_parts[4]) if len(sa_parts) > 4 else 0
                                    
                                    # Estimate end position from CIGAR if available
                                    if len(sa_parts) > 3:
                                        # Simple estimate: use read length as span
                                        sa_end = sa_pos + read.query_length
                                    else:
                                        sa_end = sa_pos + 100  # Default small span
                                    
                                    # Skip mitochondrial chromosomes
                                    if sa_chrom == "chrM" or sa_chrom == "M":
                                        continue
                                    
                                    # Add supplementary alignment as a row
                                    master_bed_rows.append(
                                        {
                                            "col1": sa_chrom,
                                            "col2": sa_pos,
                                            "col3": sa_end,
                                            "col4": "master_bed_supplementary",
                                            "reference_id": sa_chrom,
                                            "reference_start": sa_pos,
                                            "reference_end": sa_end,
                                            "read_id": read_id,
                                            "mapping_quality": sa_mapq,
                                            "strand": sa_strand,
                                            "read_start": 0,
                                            "read_end": read.query_length,
                                            "is_secondary": False,
                                            "is_supplementary": True,
                                            "mapping_span": sa_end - sa_pos,
                                        }
                                    )
                        except (ValueError, KeyError) as e:
                            logger.debug(f"Could not parse SA tag for read {read.query_name}: {e}")
                            continue
        
        except Exception as e:
            logger.error(f"Error reading BAM file: {e}")
            raise
        
        logger.info(f"Found {reads_with_supplementary_count} reads with supplementary alignments")
        
        # Process target panel candidates
        target_candidates = None
        if target_read_alignments:
            # Filter out reads with overlapping alignments and apply quality thresholds early
            # This reduces memory usage by filtering before DataFrame creation
            min_mq = get_fusion_threshold("mapping_quality")
            min_span = get_fusion_threshold("mapping_span")
            filtered_target_rows = []
            for read_id, alignments in target_read_alignments.items():
                if not _check_read_alignments_overlap(alignments):
                    # Apply quality filters early (before DataFrame creation)
                    for align in alignments:
                        if (align.get("mapping_quality", 0) > min_mq and 
                            align.get("mapping_span", 0) > min_span):
                            filtered_target_rows.append(align)
            
            if filtered_target_rows:
                target_df = pd.DataFrame(filtered_target_rows)
                target_df = _optimize_fusion_dataframe(target_df)
                target_candidates = target_df if not target_df.empty else None
        
        # Process genome-wide candidates
        genome_wide_candidates = None
        if genome_read_alignments:
            # Filter out reads with overlapping alignments and apply quality thresholds early
            min_mq = get_fusion_threshold("mapping_quality")
            min_span = get_fusion_threshold("mapping_span")
            filtered_genome_rows = []
            for read_id, alignments in genome_read_alignments.items():
                if not _check_read_alignments_overlap(alignments):
                    # Apply quality filters early (before DataFrame creation)
                    for align in alignments:
                        if (align.get("mapping_quality", 0) > min_mq and 
                            align.get("mapping_span", 0) > min_span):
                            filtered_genome_rows.append(align)
            
            if filtered_genome_rows:
                genome_df = pd.DataFrame(filtered_genome_rows)
                genome_df = _optimize_fusion_dataframe(genome_df)
                genome_wide_candidates = genome_df if not genome_df.empty else None
        
        # Process master BED candidates
        master_bed_candidates = None
        if master_bed_rows:
            # Apply quality filters early (before DataFrame creation) to reduce memory
            min_mq = get_fusion_threshold("mapping_quality")
            min_span = get_fusion_threshold("mapping_span")
            filtered_master_bed_rows = [
                row for row in master_bed_rows
                if row.get("mapping_quality", 0) > min_mq and row.get("mapping_span", 0) > min_span
            ]
            
            if filtered_master_bed_rows:
                master_bed_df = pd.DataFrame(filtered_master_bed_rows)
                master_bed_df = _optimize_fusion_dataframe(master_bed_df)
                master_bed_candidates = master_bed_df if not master_bed_df.empty else None
        
        logger.info(
            f"Single-pass results: {len(target_candidates) if target_candidates is not None else 0} target, "
            f"{len(genome_wide_candidates) if genome_wide_candidates is not None else 0} genome-wide, "
            f"{len(master_bed_candidates) if master_bed_candidates is not None else 0} master BED candidates"
        )
        
        return target_candidates, genome_wide_candidates, master_bed_candidates

    except Exception as e:
        logger.error(f"Error processing BAM file in single pass: {str(e)}")
        logger.error("Exception details:", exc_info=True)
        raise


def process_bam_for_fusions_work(
    bamfile: str, target_panel: str
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Process BAM file to find fusion candidates with memory optimization.
    This method now uses the optimized single-pass approach.
    
    NOTE: This function is kept for backward compatibility but now uses
    the single-pass implementation internally.

    Args:
        bamfile: Path to BAM file
        target_panel: Target panel to use (rCNS2, AML, or PanCan)

    Returns:
        Tuple of (target_panel_candidates, genome_wide_candidates) DataFrames
    """
    try:
        logger.info(f"DEBUG: process_bam_for_fusions_work called with target_panel='{target_panel}'")
        
        # Use single-pass processing (master BED not needed here)
        target_candidates, genome_wide_candidates, _ = process_bam_single_pass(
            bamfile, target_panel
        )
        
        logger.info(f"Target panel candidates found: {len(target_candidates) if target_candidates is not None else 0}")
        logger.info(f"Genome-wide candidates found: {len(genome_wide_candidates) if genome_wide_candidates is not None else 0}")

        return target_candidates, genome_wide_candidates

    except Exception as e:
        logger.error(f"Error processing BAM file for fusions: {str(e)}")
        logger.error("Exception details:", exc_info=True)
        raise


# =============================================================================
# STAGING OPTIMIZATION FUNCTIONS
# =============================================================================


def _get_staging_dir(work_dir: str, sample_id: str) -> str:
    """Get staging directory for temporary per-file fusion results"""
    staging_dir = os.path.join(work_dir, sample_id, "_fusion_staging")
    os.makedirs(staging_dir, exist_ok=True)
    return staging_dir


def _get_lock_file(work_dir: str, sample_id: str, lock_type: str = "counter") -> str:
    """Get lock file path for coordinating concurrent access"""
    lock_dir = os.path.join(work_dir, sample_id, "_locks")
    os.makedirs(lock_dir, exist_ok=True)
    return os.path.join(lock_dir, f"fusion_{lock_type}.lock")


def _get_pending_count(work_dir: str, sample_id: str) -> int:
    """Get count of files pending accumulation (thread-safe)."""
    staging_dir = _get_staging_dir(work_dir, sample_id)
    count_file = os.path.join(staging_dir, "pending_count.txt")
    try:
        if os.path.exists(count_file):
            with open(count_file, "r") as f:
                return int(f.read().strip())
    except (ValueError, IOError) as e:
        logger.debug(f"Could not read pending_count for {sample_id}: {e}")
    # Fallback to filesystem scan if counter missing/corrupt
    staging_files = glob.glob(os.path.join(staging_dir, "target_*.parquet"))
    return len(staging_files)


def _set_pending_count(work_dir: str, sample_id: str, count: int) -> None:
    """Set the pending staging file count (thread-safe)."""
    lock_file = _get_lock_file(work_dir, sample_id, "pending")
    staging_dir = _get_staging_dir(work_dir, sample_id)
    count_file = os.path.join(staging_dir, "pending_count.txt")
    with FileLock(lock_file, timeout=30.0):
        with open(count_file, "w") as f:
            f.write(str(max(0, int(count))))


def _increment_pending_count(work_dir: str, sample_id: str, delta: int = 1) -> int:
    """Increment the pending staging count and return the new value."""
    lock_file = _get_lock_file(work_dir, sample_id, "pending")
    staging_dir = _get_staging_dir(work_dir, sample_id)
    count_file = os.path.join(staging_dir, "pending_count.txt")
    with FileLock(lock_file, timeout=30.0):
        current = _get_pending_count(work_dir, sample_id)
        new_count = max(0, current + int(delta))
        with open(count_file, "w") as f:
            f.write(str(new_count))
    return new_count


def _atomic_counter_increment(work_dir: str, sample_id: str) -> int:
    """
    Atomically increment and return the fusion file counter for a sample.
    Uses file locking to prevent race conditions.
    
    Returns:
        The counter value to use for this file
    """
    lock_file = _get_lock_file(work_dir, sample_id, "counter")
    counter_file = os.path.join(work_dir, sample_id, "fusion_analysis_counter.txt")
    
    with FileLock(lock_file, timeout=30.0):
        # Read current counter
        if os.path.exists(counter_file):
            try:
                with open(counter_file, "r") as f:
                    counter = int(f.read().strip())
            except (ValueError, IOError):
                counter = 0
        else:
            counter = 0
        
        # Write incremented counter
        os.makedirs(os.path.dirname(counter_file), exist_ok=True)
        with open(counter_file, "w") as f:
            f.write(str(counter + 1))
        
        return counter


def process_bam_with_staging(
    file_path: str,
    temp_dir: str,
    metadata: Dict[str, Any],
    fusion_metadata: FusionMetadata,
    target_panel: str,
    has_supplementary: bool = False,
    supplementary_read_ids: List[str] = [],
    work_dir: Optional[str] = None,
    batch_size: int = 10,
) -> Tuple[Dict[str, Any], bool]:
    """
    Fast per-file processing that saves results to staging area.
    Does NOT merge with accumulated data - much faster for large datasets.
    
    Args:
        file_path: Path to BAM file
        temp_dir: Temporary directory
        metadata: File metadata
        fusion_metadata: Current fusion metadata object
        target_panel: Target panel name
        has_supplementary: Whether file has supplementary reads
        supplementary_read_ids: List of supplementary read IDs
        work_dir: Working directory for staging
        batch_size: Number of files before accumulation triggers
    
    Returns:
        Tuple of (results_dict, should_accumulate)
    """
    sample_id = fusion_metadata.sample_id
    
    if not has_supplementary:
        return {
            "has_supplementary": False,
            "target_candidates_count": 0,
            "genome_wide_candidates_count": 0,
            "target_candidates": None,
            "genome_wide_candidates": None,
        }, False
    
    if not work_dir:
        # Fallback to non-staging mode
        logger.warning("No work_dir provided - falling back to non-staging mode")
        return process_bam_file(
            file_path, temp_dir, metadata, fusion_metadata,
            target_panel, has_supplementary, supplementary_read_ids, work_dir
        ), False
    
    try:
        # Get atomic counter (thread-safe)
        counter = _atomic_counter_increment(work_dir, sample_id)
        logger.debug(f"Assigned fusion file counter: {counter}")
        
        # Process BAM file in single pass (combines all three operations)
        target_candidates, genome_wide_candidates, master_bed_candidates = process_bam_single_pass(
            file_path, target_panel, work_dir, sample_id, supplementary_read_ids
        )
        
        # Apply fusion candidate filtering
        if target_candidates is not None and not target_candidates.empty:
            target_candidates = _filter_fusion_candidates(target_candidates)
            
        if genome_wide_candidates is not None and not genome_wide_candidates.empty:
            genome_wide_candidates = _filter_fusion_candidates(genome_wide_candidates)
        
        # Save to staging (Parquet is ~5-10x faster than JSON)
        # Note: staging_dir is pre-created at batch level, but this ensures it exists
        staging_dir = _get_staging_dir(work_dir, sample_id)
        
        target_staging = os.path.join(staging_dir, f"target_{counter:06d}.parquet")
        genome_staging = os.path.join(staging_dir, f"genome_{counter:06d}.parquet")
        master_bed_staging = os.path.join(staging_dir, f"master_bed_{counter:06d}.parquet")
        
        # Save candidates to staging
        if target_candidates is not None and not target_candidates.empty:
            target_candidates.to_parquet(target_staging, index=False)
            logger.info(f"Saved {len(target_candidates)} target candidates to staging")
        else:
            # Save empty marker file
            pd.DataFrame().to_parquet(target_staging, index=False)
        
        if genome_wide_candidates is not None and not genome_wide_candidates.empty:
            genome_wide_candidates.to_parquet(genome_staging, index=False)
            logger.debug(f"Saved {len(genome_wide_candidates)} genome-wide candidates to staging")
        else:
            # Save empty marker file
            pd.DataFrame().to_parquet(genome_staging, index=False)
        
        # Save master BED candidates to staging
        if master_bed_candidates is not None and not master_bed_candidates.empty:
            master_bed_candidates.to_parquet(master_bed_staging, index=False)
            logger.info(f"Saved {len(master_bed_candidates)} master BED candidates to staging")
        else:
            # Save empty marker file
            pd.DataFrame().to_parquet(master_bed_staging, index=False)
        
        # Check if accumulation should run
        # Note: This check is not atomic - multiple workers might see threshold reached
        # The actual accumulation function will re-check inside the lock to prevent duplicate work
        pending_count = _increment_pending_count(work_dir, sample_id, delta=1)
        should_accumulate = pending_count >= batch_size
        
        logger.info(
            f"Fusion staging complete. Pending files: {pending_count}/{batch_size}"
        )
        
        if should_accumulate:
            logger.info(
                f"Accumulation threshold reached ({pending_count} >= {batch_size}) - will attempt accumulation"
            )
        
        results = {
            "has_supplementary": True,
            "target_candidates_count": (
                len(target_candidates) if target_candidates is not None else 0
            ),
            "genome_wide_candidates_count": (
                len(genome_wide_candidates) if genome_wide_candidates is not None else 0
            ),
            "master_bed_candidates_count": (
                len(master_bed_candidates) if master_bed_candidates is not None else 0
            ),
            "target_candidates": target_candidates,
            "genome_wide_candidates": genome_wide_candidates,
            "master_bed_candidates": master_bed_candidates,
        }
        
        return results, should_accumulate
        
    except Exception as e:
        logger.error(f"Error in fusion staging for {sample_id}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        # Fall back to non-staging mode on error
        return process_bam_file(
            file_path, temp_dir, metadata, fusion_metadata,
            target_panel, has_supplementary, supplementary_read_ids, work_dir
        ), False


def accumulate_fusion_candidates(
    work_dir: str,
    sample_id: str,
    target_panel: str,
    force: bool = False,
    batch_size: int = 10,
    reference: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Batch accumulation of staged fusion candidates with proper locking.
    
    This method:
    1. Locks the accumulation process to prevent concurrent accumulations
    2. Loads all staged files
    3. Efficiently concatenates them (faster than iterative merge)
    4. Merges batch with existing accumulated data
    5. Saves updated accumulated data
    6. Cleans up staging files
    
    Args:
        work_dir: Working directory
        sample_id: Sample identifier
        target_panel: Target panel type
        force: If True, accumulate even if below threshold (for end-of-run)
        batch_size: Minimum number of files to accumulate
    
    Returns:
        Dictionary with accumulation results
    """
    lock_file = _get_lock_file(work_dir, sample_id, "accumulation")
    
    lock_acquire_start = time.time()
    try:
        # Use shorter timeout to avoid workers blocking each other
        # If lock can't be acquired quickly, another worker is likely already accumulating
        with FileLock(lock_file, timeout=5.0):
            logger.info(f"Starting fusion batch accumulation for {sample_id} (force={force})")
            start_time = time.time()
            
            staging_dir = _get_staging_dir(work_dir, sample_id)
            
            # Find all staging files (re-check inside lock to avoid race conditions)
            # Multiple workers might have triggered accumulation, but only one should proceed
            target_files = sorted(glob.glob(os.path.join(staging_dir, "target_*.parquet")))
            genome_files = sorted(glob.glob(os.path.join(staging_dir, "genome_*.parquet")))
            master_bed_files = sorted(glob.glob(os.path.join(staging_dir, "master_bed_*.parquet")))
            _set_pending_count(work_dir, sample_id, len(target_files))
            
            # Track which master_bed staging files are new (for incremental processing)
            # These are the files being accumulated in this batch - all of them are "new" for this accumulation
            new_master_bed_files = set(master_bed_files)
            
            if not target_files:
                logger.info(f"No staged fusion files to accumulate for {sample_id}")
                return {"status": "no_files", "files_processed": 0}
            
            # Re-check if we should accumulate based on count (inside lock to prevent race conditions)
            # This prevents multiple workers from all trying to accumulate when threshold is reached
            if not force and len(target_files) < batch_size:
                logger.info(
                    f"Skipping fusion accumulation - only {len(target_files)} files staged "
                    f"(threshold: {batch_size}, force={force}). Another worker may have already accumulated."
                )
                return {"status": "below_threshold", "files_pending": len(target_files)}
            
            logger.info(
                f"Accumulating {len(target_files)} staged fusion files for {sample_id}"
            )
            
            # Load all staged files efficiently
            # Use list comprehensions to load all files of each type
            # This is more efficient than sequential loading and doesn't require index matching
            
            def load_parquet_files(file_list, file_type="file", max_workers=4):
                """Load multiple parquet files in parallel, skipping empty or failed files."""
                if not file_list:
                    return []
                
                def load_file(file_path):
                    try:
                        df = pd.read_parquet(file_path)
                        if not df.empty:
                            return df
                        return None
                    except Exception as e:
                        logger.warning(f"Error loading {file_type} staging file {os.path.basename(file_path)}: {e}")
                        return None
                
                # Use parallel loading for better performance with many files
                dfs = []
                if len(file_list) > 1:
                    with ThreadPoolExecutor(max_workers=max_workers) as executor:
                        results = list(executor.map(load_file, file_list))
                        dfs = [df for df in results if df is not None]
                else:
                    # Single file - no need for threading overhead
                    df = load_file(file_list[0])
                    if df is not None:
                        dfs.append(df)
                
                return dfs
            
            # Load all files of each type in parallel (no need to match by index - just load all available)
            load_start = time.time()
            target_dfs = load_parquet_files(target_files, "target")
            genome_dfs = load_parquet_files(genome_files, "genome")
            master_bed_dfs = load_parquet_files(master_bed_files, "master_bed")
            
            # Track which master_bed files are new (for incremental processing)
            # These are the files being accumulated in this batch - all of them are "new" for this accumulation
            new_master_bed_files = set(master_bed_files)
            logger.info(
                f"Loaded {len(target_dfs)} target, {len(genome_dfs)} genome-wide, "
                f"and {len(master_bed_dfs)} master BED staging files"
            )
            
            # Efficient batch merge using concat
            # Optimize for common cases: single file (no concat needed) or empty list
            def concat_dataframes(df_list):
                """Efficiently concatenate a list of DataFrames, handling edge cases."""
                if not df_list:
                    return pd.DataFrame()
                elif len(df_list) == 1:
                    return df_list[0].copy()  # Return a copy to avoid modifying original
                else:
                    return pd.concat(df_list, ignore_index=True)
            
            batch_target = concat_dataframes(target_dfs)
            batch_genome = concat_dataframes(genome_dfs)
            batch_master_bed = concat_dataframes(master_bed_dfs)
            
            logger.info(
                f"Batch merged: {len(batch_target)} target, {len(batch_genome)} genome-wide, "
                f"{len(batch_master_bed)} master BED candidates"
            )
            
            # Debug logging for master BED candidates
            if len(master_bed_files) > 0:
                logger.debug(f"Found {len(master_bed_files)} master Bed staging files")
                logger.debug(f"Loaded {len(master_bed_dfs)} non-empty master BED DataFrames from staging")
            else:
                logger.debug("No master BED staging files found")
            
            counts = _load_fusion_counts(work_dir, sample_id)
            # Append batch to dataset (append-only, avoids re-reading full history)
            batch_id = _extract_staging_batch_id(target_files)
            _append_fusion_candidates_parquet(batch_target, "target_candidates", work_dir, sample_id, batch_id)
            _append_fusion_candidates_parquet(batch_genome, "genome_wide_candidates", work_dir, sample_id, batch_id)
            _append_fusion_candidates_parquet(batch_master_bed, "master_bed_candidates", work_dir, sample_id, batch_id)

            counts["target_candidates"] = counts.get("target_candidates", 0) + len(batch_target)
            counts["genome_wide_candidates"] = counts.get("genome_wide_candidates", 0) + len(batch_genome)
            counts["master_bed_candidates"] = counts.get("master_bed_candidates", 0) + len(batch_master_bed)
            _save_fusion_counts(work_dir, sample_id, counts)

            logger.info(
                f"Final accumulated: {counts.get('target_candidates', 0)} target, "
                f"{counts.get('genome_wide_candidates', 0)} genome-wide, "
                f"{counts.get('master_bed_candidates', 0)} master BED candidates"
            )
            
            # Create updated metadata (with empty lists - data is in Parquet files)
            # This keeps the metadata structure but avoids storing large lists in JSON
            fusion_metadata = FusionMetadata(
                sample_id=sample_id,
                file_path="accumulated",
                analysis_timestamp=time.time(),
                target_panel=target_panel,
                fusion_data={
                    "target_candidates": [],  # Data stored in Parquet, not JSON
                    "genome_wide_candidates": [],  # Data stored in Parquet, not JSON
                    "master_bed_candidates": [],  # Data stored in Parquet, not JSON
                },
                analysis_results={
                    "target_candidates_count": counts.get("target_candidates", 0),
                    "genome_wide_candidates_count": counts.get("genome_wide_candidates", 0),
                    "master_bed_candidates_count": counts.get("master_bed_candidates", 0),
                },
                processing_steps=["accumulated"],
            )
            
            # Save metadata (without large candidate lists - they're in Parquet)
            _save_fusion_metadata(fusion_metadata, work_dir, sample_id)
            
            # Generate output files only on final accumulation (force=True)
            # This avoids expensive groupby operations and CSV generation during intermediate accumulations
            # Output files are only needed at the end, not after every batch
            # However, master BED breakpoint extraction should still run incrementally to build up the BED file
            if force:
                logger.info("Final accumulation detected - generating output files (CSV, BED, etc.)")
                output_start = time.time()
                _generate_output_files(
                    sample_id, 
                    fusion_metadata.analysis_results, 
                    fusion_metadata, 
                    work_dir, 
                    reference=reference,
                    generate_master_bed=True,  # Generate master BED on final accumulation
                    new_master_bed_files=new_master_bed_files,  # Pass new files for incremental processing
                )
            else:
                logger.debug("Intermediate accumulation - skipping output file generation (will be generated on final accumulation)")
            
            # Clean up staging files
            logger.info("Cleaning up fusion staging files...")
            for f in target_files + genome_files + master_bed_files:
                try:
                    os.remove(f)
                except OSError as e:
                    logger.warning(f"Could not remove staging file {f}: {e}")
            _set_pending_count(work_dir, sample_id, 0)
            
            elapsed = time.time() - start_time
            logger.info(
                f"Fusion batch accumulation complete for {sample_id}: "
                f"{len(target_files)} files in {elapsed:.2f}s "
                f"({elapsed/len(target_files):.3f}s per file)"
            )
            
            return {
                "status": "success",
                "files_processed": len(target_files),
                "target_candidates": counts.get("target_candidates", 0),
                "genome_wide_candidates": counts.get("genome_wide_candidates", 0),
                "master_bed_candidates": counts.get("master_bed_candidates", 0),
                "elapsed_time": elapsed,
            }
    
    except TimeoutError as e:
        logger.error(f"Could not acquire fusion accumulation lock for {sample_id}: {e}")
        return {"status": "lock_timeout", "error": str(e)}
    except Exception as e:
        logger.error(f"Error during fusion batch accumulation for {sample_id}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {"status": "error", "error": str(e)}


def process_bam_file(
    file_path,
    temp_dir,
    metadata,
    fusion_metadata,
    target_panel,
    has_supplementary,
    supplementary_read_ids,
    work_dir=None,
):
    """
    Process a single BAM file for fusion detection with disk persistence.

    Args:
        file_path: Path to BAM file
        temp_dir: Temporary directory
        metadata: File metadata
        fusion_metadata: Current fusion metadata object
        target_panel: Target panel name
        has_supplementary: Whether file has supplementary reads
        supplementary_read_ids: List of supplementary read IDs
        work_dir: Working directory for saving metadata (optional)
    """
    has_sup = has_supplementary
    sample_id = fusion_metadata.sample_id

    # Load existing fusion metadata from disk if work_dir is provided
    if work_dir and sample_id:
        existing_metadata = _load_fusion_metadata(work_dir, sample_id)
        if existing_metadata:
            # Merge existing data with current metadata
            fusion_metadata = _merge_fusion_metadata_objects(
                existing_metadata, fusion_metadata
            )
            logger.info(
                f"Loaded and merged existing fusion metadata for sample {sample_id}"
            )

    if has_sup:
        # Process BAM file in single pass (combines all three operations)
        target_candidates, genome_wide_candidates, master_bed_candidates = process_bam_single_pass(
            file_path, target_panel, work_dir, sample_id
        )

        # Apply fusion candidate filtering to get the final results
        if target_candidates is not None and not target_candidates.empty:
            target_candidates = _filter_fusion_candidates(target_candidates)
            
        if genome_wide_candidates is not None and not genome_wide_candidates.empty:
            genome_wide_candidates = _filter_fusion_candidates(genome_wide_candidates)
            
        fusion_metadata.processing_steps.append("supplementary_found")
        results = {
            "has_supplementary": True,
            "target_candidates_count": (
                len(target_candidates) if target_candidates is not None else 0
            ),
            "genome_wide_candidates_count": (
                len(genome_wide_candidates) if genome_wide_candidates is not None else 0
            ),
            "master_bed_candidates_count": (
                len(master_bed_candidates) if master_bed_candidates is not None else 0
            ),
            "target_candidates": target_candidates,
            "genome_wide_candidates": genome_wide_candidates,
            "master_bed_candidates": master_bed_candidates,
        }

        # Persist candidates to append-only datasets (avoid JSON bloat)
        if work_dir and sample_id:
            counts = _load_fusion_counts(work_dir, sample_id)
            batch_id = _atomic_counter_increment(work_dir, sample_id)
            _append_fusion_candidates_parquet(
                target_candidates, "target_candidates", work_dir, sample_id, batch_id
            )
            _append_fusion_candidates_parquet(
                genome_wide_candidates, "genome_wide_candidates", work_dir, sample_id, batch_id
            )
            _append_fusion_candidates_parquet(
                master_bed_candidates, "master_bed_candidates", work_dir, sample_id, batch_id
            )
            counts["target_candidates"] = counts.get("target_candidates", 0) + results["target_candidates_count"]
            counts["genome_wide_candidates"] = counts.get("genome_wide_candidates", 0) + results["genome_wide_candidates_count"]
            counts["master_bed_candidates"] = counts.get("master_bed_candidates", 0) + results["master_bed_candidates_count"]
            _save_fusion_counts(work_dir, sample_id, counts)

        # Keep metadata structure without embedding large candidate lists
        fusion_metadata.fusion_data = {
            "target_candidates": [],
            "genome_wide_candidates": [],
            "master_bed_candidates": [],
        }

        # Save updated fusion metadata to disk
        if work_dir and sample_id:
            _save_fusion_metadata(fusion_metadata, work_dir, sample_id)
            logger.info(f"Saved updated fusion metadata for sample {sample_id}")

        return results
    else:
        # Even if no supplementary reads, save metadata to disk for consistency
        if work_dir and sample_id:
            _save_fusion_metadata(fusion_metadata, work_dir, sample_id)

        return {
            "has_supplementary": False,
            "target_candidates_count": 0,
            "genome_wide_candidates_count": 0,
            "target_candidates": None,
            "genome_wide_candidates": None,
        }


def _generate_output_files(
    sample_id: str,
    analysis_results: Dict[str, Any],
    fusion_metadata: FusionMetadata,
    work_dir: str,
    reference: Optional[str] = None,
    generate_master_bed: bool = False,
    new_master_bed_files: Optional[Set[str]] = None,
) -> Dict[str, str]:
    """
    Generate output files for the analysis, including visualization preprocessing.
    This function now works with accumulated data from fusion_metadata instead of
    just the current analysis_results.

    Args:
            sample_id: Sample ID
            analysis_results: Analysis results
            fusion_metadata: Fusion metadata object with accumulated data
            work_dir: Working directory

    Returns:
            Dictionary mapping file type to file path
    """
    output_start = time.time()
    # Use accumulated data from fusion_metadata instead of just current analysis_results
    # Load from Parquet files if available (fast path), otherwise use in-memory data
    target_candidates = None
    if (
        hasattr(fusion_metadata, "fusion_data")
        and fusion_metadata.fusion_data
    ):
        # Try loading from Parquet first (fast path)
        target_candidates = _load_fusion_candidates_parquet("target_candidates", work_dir, sample_id)
        
        # Fallback to in-memory data if Parquet not available
        if target_candidates is None or target_candidates.empty:
            target_candidates_data = fusion_metadata.fusion_data.get("target_candidates")
            if target_candidates_data:
                if isinstance(target_candidates_data, pd.DataFrame):
                    target_candidates = target_candidates_data
                elif isinstance(target_candidates_data, list):
                    target_candidates = pd.DataFrame(target_candidates_data)
    
    if target_candidates is not None and not target_candidates.empty:
            # Filter for fusion candidates using the corrected logic
            groupby_start = time.time()
            gene_counts = target_candidates.groupby("read_id", observed=True)["col4"].nunique()
            fusion_read_ids = gene_counts[gene_counts > 1].index
            
            if len(fusion_read_ids) > 0:
                result = target_candidates[target_candidates["read_id"].isin(fusion_read_ids)]
                
                # Apply minimum read support threshold (3 or more supporting reads per gene pair)
                if not result.empty:
                    tag_start = time.time()
                    # Create tag column by grouping genes per read_id (same logic as _annotate_results)
                    lookup = result.groupby("read_id", observed=True)["col4"].agg(
                        lambda x: ",".join(sorted(set(x)))
                    )
                    result["tag"] = result["read_id"].map(lookup)
                    
                    # Group by gene pair (tag) and count supporting reads
                    pair_count_start = time.time()
                    gene_pair_read_counts = result.groupby("tag", observed=True)["read_id"].nunique()
                    min_support = get_fusion_threshold("read_support")
                    valid_gene_pairs = gene_pair_read_counts[gene_pair_read_counts >= min_support].index
                    result = result[result["tag"].isin(valid_gene_pairs)]

                if not result.empty:
                    csv_start = time.time()
                    # Save the filtered fusion candidates to CSV
                    result.to_csv(
                        os.path.join(
                            work_dir, sample_id, "fusion_candidates_master.csv"
                        ),
                        index=False,
                    )

                    # Preprocess for visualization
                    preprocess_start = time.time()
                    preprocess_fusion_data_standalone(
                        result,
                        os.path.join(
                            work_dir,
                            sample_id,
                            "fusion_candidates_master_processed.pkl",
                        ),
                    )

    # Use accumulated data from fusion_metadata for genome-wide candidates
    # Load from Parquet files if available (fast path), otherwise use in-memory data
    genome_wide_candidates = None
    logger.info(f"Checking genome-wide candidates in fusion_metadata: has_fusion_data={hasattr(fusion_metadata, 'fusion_data')}, fusion_data_exists={fusion_metadata.fusion_data is not None if hasattr(fusion_metadata, 'fusion_data') else False}")
    if (
        hasattr(fusion_metadata, "fusion_data")
        and fusion_metadata.fusion_data
    ):
        # Try loading from Parquet first (fast path)
        genome_wide_candidates = _load_fusion_candidates_parquet("genome_wide_candidates", work_dir, sample_id)
        
        # Fallback to in-memory data if Parquet not available
        if genome_wide_candidates is None or genome_wide_candidates.empty:
            genome_wide_data = fusion_metadata.fusion_data.get("genome_wide_candidates")
            if genome_wide_data:
                if isinstance(genome_wide_data, pd.DataFrame):
                    genome_wide_candidates = genome_wide_data
                elif isinstance(genome_wide_data, list):
                    genome_wide_candidates = pd.DataFrame(genome_wide_data)
                    logger.info(f"Found {len(genome_wide_data)} genome-wide candidate records")
    
    if genome_wide_candidates is not None and not genome_wide_candidates.empty:
            # Filter for fusion candidates using the corrected logic
            groupby_start = time.time()
            gene_counts_all = genome_wide_candidates.groupby("read_id", observed=True)["col4"].nunique()
            fusion_read_ids_all = gene_counts_all[gene_counts_all > 1].index
            
            if len(fusion_read_ids_all) > 0:
                result_all = genome_wide_candidates[genome_wide_candidates["read_id"].isin(fusion_read_ids_all)]
                logger.info(f"Genome-wide fusion candidates after basic filtering: {len(result_all)} records")

                if not result_all.empty:
                    csv_start = time.time()
                    # Save the filtered genome-wide candidates to CSV
                    result_all.to_csv(
                        os.path.join(work_dir, sample_id, "fusion_candidates_all.csv"),
                        index=False,
                    )

                    # Use the same preprocessing pipeline as target candidates
                    # This ensures proper annotation with tags, colors, and gene group detection
                    preprocess_start = time.time()
                    preprocess_fusion_data_standalone(
                        result_all,
                        os.path.join(
                            work_dir, sample_id, "fusion_candidates_all_processed.pkl"
                        ),
                    )

    # Save master BED fusion candidates (no gene overlap filtering needed)
    # Load from Parquet files if available (fast path), otherwise use in-memory data
    master_bed_candidates = None
    if (
        hasattr(fusion_metadata, "fusion_data")
        and fusion_metadata.fusion_data
    ):
        # Try loading from Parquet first (fast path)
        master_bed_candidates = _load_fusion_candidates_parquet("master_bed_candidates", work_dir, sample_id)
        
        # Fallback to in-memory data if Parquet not available
        if master_bed_candidates is None or master_bed_candidates.empty:
            master_bed_data = fusion_metadata.fusion_data.get("master_bed_candidates")
            if master_bed_data:
                if isinstance(master_bed_data, pd.DataFrame):
                    master_bed_candidates = master_bed_data
                elif isinstance(master_bed_data, list):
                    master_bed_candidates = pd.DataFrame(master_bed_data)
                    logger.info(f"Found {len(master_bed_data)} master BED candidate records")
    
    if master_bed_candidates is not None and not master_bed_candidates.empty:
            # For master BED candidates, we track reads mapping to multiple locations
            # CSV generation for master BED candidates has been deprecated
            # The Parquet file (master_bed_candidates.parquet) is the source of truth
            # and is used for all BED generation. The CSV was only used by the GUI,
            # which has been updated to no longer display the master BED table.
            pass
    
    # Create sv_count.txt file with content "0" if it doesn't exist

    # Create sv_count.txt file with content "0" if it doesn't exist
    sv_count_file = os.path.join(work_dir, sample_id, "sv_count.txt")
    if not os.path.exists(sv_count_file):
        with open(sv_count_file, "w") as f:
            f.write("0")

    # Generate fusion breakpoint BED file
    _generate_fusion_breakpoint_bed(sample_id, fusion_metadata, work_dir)
    
    # Generate master BED breakpoint BED file (new target regions from supplementary alignments)
    # This is called incrementally as data accumulates. For large datasets, we use an incremental
    # approach: only process NEW staging files and merge with existing breakpoints.
    if generate_master_bed:
        master_bed_candidates = _load_fusion_candidates_parquet("master_bed_candidates", work_dir, sample_id)
        if master_bed_candidates is not None and not master_bed_candidates.empty:
            _generate_master_bed_breakpoint_bed(
                sample_id, 
                fusion_metadata, 
                work_dir,
                new_master_bed_files=new_master_bed_files,  # Pass new files for incremental processing
            )
    
    # Generate master BED file only if requested (should only be done once per batch at the end)
    # Use async (non-blocking) generation to avoid blocking the analysis pipeline
    if generate_master_bed:
        try:
            from robin.analysis.master_bed_generator import generate_master_bed_async
            
            # Get analysis counter
            analysis_counter = _load_analysis_counter(sample_id, work_dir)
            
            # Get target_panel from fusion_metadata
            target_panel = fusion_metadata.target_panel if hasattr(fusion_metadata, 'target_panel') else None
            
            # Generate asynchronously (non-blocking)
            generate_master_bed_async(
                sample_id=sample_id,
                work_dir=work_dir,
                analysis_counter=analysis_counter,
                target_panel=target_panel,
                logger_instance=logger,
                reference=reference,
            )
        except Exception as e:
            logger.warning(f"Could not start async master BED generation: {e}")

    output_paths = {
        "target_candidates_path": os.path.join(
            work_dir, sample_id, "fusion_candidates_master.csv"
        ),
        "genome_wide_candidates_path": os.path.join(
            work_dir, sample_id, "fusion_candidates_all.csv"
        ),
    }
    
    # Master BED CSV generation has been deprecated - Parquet file is used instead
    
    return output_paths


def _get_cnv_bin_width(work_dir: str, sample_id: str) -> int:
    """
    Try to load bin_width from CNV analysis results if available.
    
    Args:
        work_dir: Working directory
        sample_id: Sample ID
        
    Returns:
        bin_width in base pairs, or default 10000 (10kb) if not available
    """
    try:
        sample_dir = os.path.join(work_dir, sample_id)
        cnv_dict_path = os.path.join(sample_dir, "CNV_dict.npy")
        
        if os.path.exists(cnv_dict_path):
            cnv_dict = np.load(cnv_dict_path, allow_pickle=True).item()
            bin_width = cnv_dict.get("bin_width")
            if bin_width and bin_width > 0:
                logger.debug(f"Loaded bin_width={bin_width} from CNV analysis")
                return int(bin_width)
    except Exception as e:
        logger.debug(f"Could not load bin_width from CNV analysis: {e}")
    
    # Default to 10kb if CNV data not available
    default_bin_width = 10000
    logger.debug(f"Using default bin_width={default_bin_width}")
    return default_bin_width


def _extract_fusion_breakpoints(fusion_metadata: FusionMetadata) -> List[Dict[str, Any]]:
    """
    Extract fusion breakpoint coordinates from fusion metadata.
    Only includes fusions that meet the minimum read support threshold.
    
    Args:
        fusion_metadata: FusionMetadata object with fusion candidates
        
    Returns:
        List of dictionaries with fusion breakpoint information
    """
    breakpoints = []
    breakpoint_set = set()  # For deduplication
    
    fusion_data = fusion_metadata.fusion_data or {}
    min_support = get_fusion_threshold("read_support")
    
    # Process target candidates - load from in-memory data (Parquet loading handled by _load_fusion_metadata)
    target_candidates = fusion_data.get("target_candidates", [])
    target_df = None
    if target_candidates:
        if isinstance(target_candidates, list):
            target_df = pd.DataFrame(target_candidates)
        elif isinstance(target_candidates, pd.DataFrame):
            target_df = target_candidates
    
    if target_df is not None and not target_df.empty:
        
        if not target_df.empty and "reference_start" in target_df.columns:
            # Filter for fusion candidates (reads mapping to multiple genes)
            gene_counts = target_df.groupby("read_id", observed=True)["col4"].nunique()
            fusion_read_ids = gene_counts[gene_counts > 1].index
            
            if len(fusion_read_ids) > 0:
                fusion_df = target_df[target_df["read_id"].isin(fusion_read_ids)]
                
                # Apply minimum read support threshold
                if not fusion_df.empty:
                    # Create tag column by grouping genes per read_id
                    lookup = fusion_df.groupby("read_id", observed=True)["col4"].agg(
                        lambda x: ",".join(sorted(set(x)))
                    )
                    fusion_df["tag"] = fusion_df["read_id"].map(lookup)
                    
                    # Group by gene pair (tag) and count supporting reads
                    gene_pair_read_counts = fusion_df.groupby("tag", observed=True)["read_id"].nunique()
                    valid_gene_pairs = gene_pair_read_counts[gene_pair_read_counts >= min_support].index
                    filtered_df = fusion_df[fusion_df["tag"].isin(valid_gene_pairs)]
                    
                    # Extract breakpoint coordinates for valid fusions
                    for _, row in filtered_df.iterrows():
                        chrom = row.get("reference_id", "Unknown")
                        start = int(row.get("reference_start", 0))
                        end = int(row.get("reference_end", 0))
                        gene = row.get("col4", "Unknown")
                        
                        if start > 0 and end > start:
                            bp_key = (chrom, start, end)
                            if bp_key not in breakpoint_set:
                                breakpoint_set.add(bp_key)
                                breakpoints.append({
                                    "chromosome": chrom,
                                    "start": start,
                                    "end": end,
                                    "gene": gene,
                                    "read_id": row.get("read_id", "Unknown"),
                                    "source": "target"
                                })
    
    # Process genome-wide candidates - load from in-memory data (Parquet loading handled by _load_fusion_metadata)
    genome_wide_candidates = fusion_data.get("genome_wide_candidates", [])
    genome_df = None
    if genome_wide_candidates:
        if isinstance(genome_wide_candidates, list):
            genome_df = pd.DataFrame(genome_wide_candidates)
        elif isinstance(genome_wide_candidates, pd.DataFrame):
            genome_df = genome_wide_candidates
    
    if genome_df is not None and not genome_df.empty:
        
        if not genome_df.empty and "reference_start" in genome_df.columns:
            # Filter for fusion candidates (reads mapping to multiple genes)
            gene_counts_all = genome_df.groupby("read_id", observed=True)["col4"].nunique()
            fusion_read_ids_all = gene_counts_all[gene_counts_all > 1].index
            
            if len(fusion_read_ids_all) > 0:
                fusion_df_all = genome_df[genome_df["read_id"].isin(fusion_read_ids_all)]
                
                # Apply minimum read support threshold
                if not fusion_df_all.empty:
                    # Create tag column by grouping genes per read_id
                    lookup_all = fusion_df_all.groupby("read_id", observed=True)["col4"].agg(
                        lambda x: ",".join(sorted(set(x)))
                    )
                    fusion_df_all["tag"] = fusion_df_all["read_id"].map(lookup_all)
                    
                    # Group by gene pair (tag) and count supporting reads
                    gene_pair_read_counts_all = fusion_df_all.groupby("tag", observed=True)["read_id"].nunique()
                    valid_gene_pairs_all = gene_pair_read_counts_all[gene_pair_read_counts_all >= min_support].index
                    filtered_df_all = fusion_df_all[fusion_df_all["tag"].isin(valid_gene_pairs_all)]
                    
                    # Extract breakpoint coordinates for valid fusions
                    for _, row in filtered_df_all.iterrows():
                        chrom = row.get("reference_id", "Unknown")
                        start = int(row.get("reference_start", 0))
                        end = int(row.get("reference_end", 0))
                        gene = row.get("col4", "Unknown")
                        
                        if start > 0 and end > start:
                            bp_key = (chrom, start, end)
                            if bp_key not in breakpoint_set:
                                breakpoint_set.add(bp_key)
                                breakpoints.append({
                                    "chromosome": chrom,
                                    "start": start,
                                    "end": end,
                                    "gene": gene,
                                    "read_id": row.get("read_id", "Unknown"),
                                    "source": "genome_wide"
                                })
    
    return breakpoints


def _load_analysis_counter(sample_id: str, work_dir: str) -> int:
    """
    Load the analysis counter for a sample from disk.
    
    Args:
        sample_id: Sample ID
        work_dir: Working directory
        
    Returns:
        Analysis counter value, or 0 if not found
    """
    try:
        counter_file = os.path.join(work_dir, sample_id, "cnv_analysis_counter.txt")
        if os.path.exists(counter_file):
            with open(counter_file, "r") as f:
                return int(f.read().strip())
    except (ValueError, IOError) as e:
        logger.debug(f"Could not load analysis counter for {sample_id}: {e}")
    return 0


def _create_breakpoint_pairs_vectorized(
    primary_df: pd.DataFrame,
    supplementary_df: pd.DataFrame,
    read_id: str,
) -> List[Dict[str, Any]]:
    """
    Create breakpoint pairs using vectorized Pandas operations (much faster than nested loops).
    
    This function uses Pandas merge to create the cartesian product efficiently,
    which is significantly faster than nested Python loops.
    
    Args:
        primary_df: DataFrame with primary alignments for a single read
        supplementary_df: DataFrame with supplementary alignments for a single read
        read_id: Read ID
    
    Returns:
        List of breakpoint pair dictionaries
    """
    if primary_df.empty or supplementary_df.empty:
        return []
    
    # Prepare primary alignments
    primary_cols = ["reference_id", "reference_start", "reference_end"]
    if "mapping_quality" in primary_df.columns:
        primary_cols.append("mapping_quality")
    if "mapping_span" in primary_df.columns:
        primary_cols.append("mapping_span")
    
    # Prepare supplementary alignments
    supp_cols = ["reference_id", "reference_start", "reference_end"]
    if "mapping_quality" in supplementary_df.columns:
        supp_cols.append("mapping_quality")
    if "mapping_span" in supplementary_df.columns:
        supp_cols.append("mapping_span")
    
    # Select only the columns we need
    primaries = primary_df[primary_cols].copy()
    supplementaries = supplementary_df[supp_cols].copy()
    
    # Add a key column for cross join
    primaries["_key"] = 1
    supplementaries["_key"] = 1
    
    # Perform cross join using merge (much faster than nested loops)
    pairs_df = primaries.merge(
        supplementaries,
        on="_key",
        suffixes=("_primary", "_supplementary")
    ).drop("_key", axis=1)
    
    # Convert to list of dictionaries
    breakpoint_pairs = []
    for _, row in pairs_df.iterrows():
        pair = {
            "read_id": read_id,
            "primary_chrom": row["reference_id_primary"],
            "primary_start": int(row["reference_start_primary"]),
            "primary_end": int(row["reference_end_primary"]),
            "supp_chrom": row["reference_id_supplementary"],
            "supp_start": int(row["reference_start_supplementary"]),
            "supp_end": int(row["reference_end_supplementary"]),
        }
        
        # Add optional fields
        if "mapping_quality_primary" in row:
            pair["primary_mapq"] = row["mapping_quality_primary"] if pd.notna(row["mapping_quality_primary"]) else 0
        if "mapping_span_primary" in row:
            pair["primary_span"] = row["mapping_span_primary"] if pd.notna(row["mapping_span_primary"]) else 0
        if "mapping_quality_supplementary" in row:
            pair["supp_mapq"] = row["mapping_quality_supplementary"] if pd.notna(row["mapping_quality_supplementary"]) else 0
        if "mapping_span_supplementary" in row:
            pair["supp_span"] = row["mapping_span_supplementary"] if pd.notna(row["mapping_span_supplementary"]) else 0
        
        breakpoint_pairs.append(pair)
    
    return breakpoint_pairs


def _merge_nearby_events(
    events_df: pd.DataFrame,
    cluster_distance: int = 5000,
) -> pd.DataFrame:
    """
    Merge nearby events on the same chromosome to remove duplicates.
    
    Events are merged if they:
    - Are on the same chromosome
    - Have the same event_type
    - Overlap or are within cluster_distance of each other
    
    Args:
        events_df: DataFrame with events (chromosome, start, end, event_type, read_count, etc.)
        cluster_distance: Maximum distance for merging events (in base pairs)
    
    Returns:
        DataFrame with merged events
    """
    if events_df.empty:
        return events_df
    
    # Group by chromosome and event_type
    merged_events = []
    
    for (chrom, event_type), group in events_df.groupby(["chromosome", "event_type"], observed=True):
        if group.empty:
            continue
        
        # Sort by start position
        group_sorted = group.sort_values("start").copy()
        
        # Merge overlapping or nearby events
        current_events = []
        for _, row in group_sorted.iterrows():
            merged = False
            
            # Check if this event overlaps or is near any existing merged event
            for i, existing in enumerate(current_events):
                # Check if events overlap
                overlap_start = max(existing["start"], row["start"])
                overlap_end = min(existing["end"], row["end"])
                overlaps = overlap_end > overlap_start
                
                # Check if events are nearby (gap between them is <= cluster_distance)
                # Gap is the distance between the closest endpoints
                if row["start"] > existing["end"]:
                    gap = row["start"] - existing["end"]
                elif existing["start"] > row["end"]:
                    gap = existing["start"] - row["end"]
                else:
                    gap = 0  # Overlapping
                
                if overlaps or gap <= cluster_distance:
                    # Merge: take the union of coordinates and max read count
                    current_events[i] = {
                        "chromosome": chrom,
                        "start": min(existing["start"], row["start"]),
                        "end": max(existing["end"], row["end"]),
                        "event_type": event_type,
                        "read_count": max(existing["read_count"], row["read_count"]),
                        "avg_mapping_quality": max(
                            existing.get("avg_mapping_quality", 0),
                            row.get("avg_mapping_quality", 0)
                        ),
                        "avg_mapping_span": max(
                            existing.get("avg_mapping_span", 0),
                            row.get("avg_mapping_span", 0)
                        ),
                    }
                    merged = True
                    break
            
            if not merged:
                # Add as new event
                current_events.append({
                    "chromosome": chrom,
                    "start": row["start"],
                    "end": row["end"],
                    "event_type": event_type,
                    "read_count": row["read_count"],
                    "avg_mapping_quality": row.get("avg_mapping_quality", 0),
                    "avg_mapping_span": row.get("avg_mapping_span", 0),
                })
        
        merged_events.extend(current_events)
    
    if not merged_events:
        return pd.DataFrame()
    
    return pd.DataFrame(merged_events)


def _cluster_breakpoint_pairs_dbscan(
    breakpoint_pairs: List[Dict[str, Any]],
    cluster_distance: int = 5000,
    min_read_support: int = 3,
) -> List[Dict[str, Any]]:
    """
    Cluster breakpoint pairs using DBSCAN algorithm for efficient O(n log n) clustering.
    
    This replaces the O(n²) nested loop approach with DBSCAN, which is much faster
    for large datasets. Two breakpoint pairs are clustered if:
    - Primary locations are close (within cluster_distance)
    - Supplementary locations are close (within cluster_distance)
    
    Args:
        breakpoint_pairs: List of breakpoint pair dictionaries
        cluster_distance: Maximum distance for clustering (in base pairs)
        min_read_support: Minimum number of reads required per cluster
    
    Returns:
        List of clustered breakpoint pairs with aggregated information
    """
    if not breakpoint_pairs:
        return []
    
    # Group pairs by chromosome combination (required for clustering)
    pairs_by_chrom = {}
    for i, pair in enumerate(breakpoint_pairs):
        chrom_key = (pair["primary_chrom"], pair["supp_chrom"])
        if chrom_key not in pairs_by_chrom:
            pairs_by_chrom[chrom_key] = []
        pairs_by_chrom[chrom_key].append((i, pair))
    
    clustered_pairs = []
    
    # Cluster within each chromosome combination
    for chrom_key, chrom_pairs in pairs_by_chrom.items():
        if len(chrom_pairs) == 0:
            continue
        
        # Extract indices and pairs
        indices = [idx for idx, _ in chrom_pairs]
        pairs = [pair for _, pair in chrom_pairs]
        
        # Create feature matrix for DBSCAN: [primary_midpoint, supp_midpoint]
        # We use midpoints for clustering to reduce dimensionality
        features = np.array([
            [
                (pair["primary_start"] + pair["primary_end"]) / 2,  # Primary midpoint
                (pair["supp_start"] + pair["supp_end"]) / 2,  # Supplementary midpoint
            ]
            for pair in pairs
        ])
        
        # Use DBSCAN with Manhattan distance (L1 norm) for genomic coordinates
        # eps is the maximum distance between samples in the same cluster
        # We use cluster_distance * 1.5 to account for coordinate ranges
        eps = cluster_distance * 1.5
        min_samples = min_read_support  # Minimum samples in a cluster
        
        # Cluster using DBSCAN
        clustering = DBSCAN(eps=eps, min_samples=min_samples, metric='manhattan')
        labels = clustering.fit_predict(features)
        
        # Group pairs by cluster label
        clusters = {}
        for idx, label in enumerate(labels):
            if label == -1:  # Noise points (not in any cluster)
                continue
            
            if label not in clusters:
                clusters[label] = {
                    "indices": [],
                    "pairs": [],
                    "read_ids": set(),
                    "primary_starts": [],
                    "primary_ends": [],
                    "supp_starts": [],
                    "supp_ends": [],
                    "primary_mapqs": [],
                    "primary_spans": [],
                    "supp_mapqs": [],
                    "supp_spans": [],
                }
            
            clusters[label]["indices"].append(indices[idx])
            clusters[label]["pairs"].append(pairs[idx])
            clusters[label]["read_ids"].add(pairs[idx]["read_id"])
            clusters[label]["primary_starts"].append(pairs[idx]["primary_start"])
            clusters[label]["primary_ends"].append(pairs[idx]["primary_end"])
            clusters[label]["supp_starts"].append(pairs[idx].get("supp_start", 0))
            clusters[label]["supp_ends"].append(pairs[idx].get("supp_end", 0))
            
            # Optional fields
            if "primary_mapq" in pairs[idx]:
                clusters[label]["primary_mapqs"].append(pairs[idx]["primary_mapq"])
            if "primary_span" in pairs[idx]:
                clusters[label]["primary_spans"].append(pairs[idx]["primary_span"])
            if "supp_mapq" in pairs[idx]:
                clusters[label]["supp_mapqs"].append(pairs[idx]["supp_mapq"])
            if "supp_span" in pairs[idx]:
                clusters[label]["supp_spans"].append(pairs[idx]["supp_span"])
        
        # Create clustered breakpoint pairs
        for label, cluster_data in clusters.items():
            if len(cluster_data["read_ids"]) < min_read_support:
                continue
            
            clustered_pair = {
                "primary_chrom": pairs[0]["primary_chrom"],
                "primary_start": min(cluster_data["primary_starts"]),
                "primary_end": max(cluster_data["primary_ends"]),
                "supp_chrom": pairs[0]["supp_chrom"],
                "supp_start": min(cluster_data["supp_starts"]),
                "supp_end": max(cluster_data["supp_ends"]),
                "read_count": len(cluster_data["read_ids"]),
                "read_ids": cluster_data["read_ids"],
            }
            
            # Add optional fields if available
            if cluster_data["primary_mapqs"]:
                clustered_pair["primary_avg_mapq"] = sum(cluster_data["primary_mapqs"]) / len(cluster_data["primary_mapqs"])
            if cluster_data["primary_spans"]:
                clustered_pair["primary_avg_span"] = sum(cluster_data["primary_spans"]) / len(cluster_data["primary_spans"])
            if cluster_data["supp_mapqs"]:
                clustered_pair["supp_avg_mapq"] = sum(cluster_data["supp_mapqs"]) / len(cluster_data["supp_mapqs"])
            if cluster_data["supp_spans"]:
                clustered_pair["supp_avg_span"] = sum(cluster_data["supp_spans"]) / len(cluster_data["supp_spans"])
            
            clustered_pairs.append(clustered_pair)
    
    return clustered_pairs


def _extract_master_bed_breakpoints(fusion_metadata: FusionMetadata, work_dir: Optional[str] = None, min_read_support: int = 3) -> List[Dict[str, Any]]:
    """
    Extract breakpoint coordinates from master BED fusion candidates.
    Identifies genomic rearrangements (deletions, inversions, translocations) by finding
    breakpoint pairs (primary + supplementary alignments) supported by multiple reads.
    
    A rearrangement event is defined by a breakpoint pair:
    - Primary alignment location (where the read starts)
    - Supplementary alignment location (where the read continues)
    
    Only includes breakpoint pairs supported by at least min_read_support reads.
    The breakpoint pairs can be on different chromosomes (translocations).
    
    Args:
        fusion_metadata: FusionMetadata object with master BED candidates
        work_dir: Working directory (required to load from Parquet files)
        min_read_support: Minimum number of reads required to support a breakpoint pair (default: 3)
        
    Returns:
        List of dictionaries with master BED breakpoint information (both primary and supplementary regions)
    """
    breakpoints = []
    
    # Try loading from Parquet files first (fast path, new storage format)
    master_bed_df = None
    sample_id = fusion_metadata.sample_id
    
    if work_dir and sample_id:
        master_bed_df = _load_fusion_candidates_parquet("master_bed_candidates", work_dir, sample_id)
        if master_bed_df is not None and not master_bed_df.empty:
            logger.debug(f"Loaded {len(master_bed_df)} master BED candidates from Parquet for breakpoint extraction")
    
    # Fallback to in-memory data (for backward compatibility)
    if master_bed_df is None or master_bed_df.empty:
        fusion_data = fusion_metadata.fusion_data or {}
        master_bed_candidates = fusion_data.get("master_bed_candidates", [])
        if master_bed_candidates:
            if isinstance(master_bed_candidates, list):
                master_bed_df = pd.DataFrame(master_bed_candidates)
                logger.debug(f"Loaded {len(master_bed_candidates)} master BED candidates from in-memory data (fallback)")
            elif isinstance(master_bed_candidates, pd.DataFrame):
                master_bed_df = master_bed_candidates
                logger.debug(f"Loaded {len(master_bed_candidates)} master BED candidates from in-memory DataFrame (fallback)")
        else:
            logger.debug("No master BED candidates found in Parquet or in-memory data")
    
    # Process master BED candidates if we have data
    # TEMPORARY: Skip breakpoint extraction for performance testing
    # Set SKIP_MASTER_BED_BREAKPOINTS=True to disable this section
    SKIP_MASTER_BED_BREAKPOINTS = os.environ.get("SKIP_MASTER_BED_BREAKPOINTS", "False").lower() == "true"
    
    if SKIP_MASTER_BED_BREAKPOINTS:
        logger.info("SKIPPING master BED breakpoint extraction (SKIP_MASTER_BED_BREAKPOINTS=True)")
        return breakpoints
    
    if master_bed_df is not None and not master_bed_df.empty:
        extract_start = time.time()
        # Group reads by their breakpoint pairs (primary + supplementary alignments)
        # Each read with supplementary alignments represents a potential rearrangement
        if "read_id" not in master_bed_df.columns:
            logger.debug("Missing required column (read_id) in master BED candidates")
            return breakpoints
        
        # Separate primary and supplementary alignments using col4
        # col4 is more reliable: "master_bed_region" = primary alignment, "master_bed_supplementary" = supplementary
        # The is_supplementary flag can be inconsistent (some primary alignments may have is_supplementary=True
        # if they're part of a chimeric read set in the BAM file)
        filter_start = time.time()
        if "col4" in master_bed_df.columns:
            primary_df = master_bed_df[master_bed_df["col4"] == "master_bed_region"].copy()
            supplementary_df = master_bed_df[master_bed_df["col4"] == "master_bed_supplementary"].copy()
        else:
            # Fallback: use is_supplementary if col4 is not available
            logger.debug("col4 column not found, using is_supplementary as fallback")
            if "is_supplementary" not in master_bed_df.columns:
                logger.debug("Missing required columns (col4 or is_supplementary) in master BED candidates")
                return breakpoints
            primary_df = master_bed_df[master_bed_df["is_supplementary"] == False].copy()
            supplementary_df = master_bed_df[master_bed_df["is_supplementary"] == True].copy()
        
        # Debug: check for any misclassified alignments
        if "is_supplementary" in master_bed_df.columns and "col4" in master_bed_df.columns:
            misclassified = master_bed_df[
                (master_bed_df["col4"] == "master_bed_region") & (master_bed_df["is_supplementary"] == True)
            ]
            if not misclassified.empty:
                logger.debug(f"Found {len(misclassified)} primary alignments with is_supplementary=True (using col4 for classification)")
        
        if primary_df.empty or supplementary_df.empty:
            logger.debug("Need both primary and supplementary alignments to identify breakpoint pairs")
            return breakpoints
        
        # Find reads that have both primary and supplementary alignments
        read_id_start = time.time()
        primary_read_ids = set(primary_df["read_id"].unique())
        supplementary_read_ids = set(supplementary_df["read_id"].unique())
        reads_with_both = primary_read_ids & supplementary_read_ids
        
        logger.debug(
            f"Found {len(reads_with_both)} reads with both primary and supplementary alignments "
            f"(out of {len(primary_read_ids)} reads with primary, {len(supplementary_read_ids)} with supplementary)"
        )
        
        # Debug: log some example read IDs if available
        if reads_with_both:
            example_reads = list(reads_with_both)[:3]
            logger.debug(f"Example reads with both alignments: {example_reads}")
        elif primary_read_ids:
            example_primary = list(primary_read_ids)[:3]
            logger.debug(f"Example reads with only primary: {example_primary}")
        elif supplementary_read_ids:
            example_supp = list(supplementary_read_ids)[:3]
            logger.debug(f"Example reads with only supplementary: {example_supp}")
        
        if not reads_with_both:
            logger.debug("No reads have both primary and supplementary alignments")
            return breakpoints
        
        # Filter to only reads with both primary and supplementary alignments
        filter_both_start = time.time()
        primary_filtered = primary_df[primary_df["read_id"].isin(reads_with_both)].copy()
        supplementary_filtered = supplementary_df[supplementary_df["read_id"].isin(reads_with_both)].copy()
        
        # For each read, create breakpoint pairs (primary + supplementary)
        # A breakpoint pair represents a potential rearrangement event
        cluster_distance = 5000  # Maximum distance for clustering breakpoint pairs (5kb)
        
        # Build breakpoint pairs: for each read, pair its primary alignment with each supplementary alignment
        # Optimized: use vectorized operations instead of iterrows()
        breakpoint_pairs = []
        
        # Group by read_id for efficient lookup
        groupby_start = time.time()
        primary_by_read = primary_filtered.groupby("read_id", observed=True)
        supplementary_by_read = supplementary_filtered.groupby("read_id", observed=True)
        
        pair_creation_start = time.time()
        # Sort reads_with_both to ensure deterministic processing order
        # This prevents different results when DataFrames have different row orders
        for read_id in sorted(reads_with_both):
            # Get all primary and supplementary alignments for this read
            read_primaries = primary_by_read.get_group(read_id) if read_id in primary_by_read.groups else pd.DataFrame()
            read_supplementaries = supplementary_by_read.get_group(read_id) if read_id in supplementary_by_read.groups else pd.DataFrame()
            
            if read_primaries.empty or read_supplementaries.empty:
                continue
            
            # Create pairs using vectorized Pandas operations (much faster than nested loops)
            pairs = _create_breakpoint_pairs_vectorized(read_primaries, read_supplementaries, read_id)
            breakpoint_pairs.extend(pairs)
        
        if not breakpoint_pairs:
            logger.debug("No breakpoint pairs created")
            return breakpoints
        
        # Cluster similar breakpoint pairs using DBSCAN (much faster than O(n²) nested loops)
        clustering_start = time.time()
        clustered_pairs = _cluster_breakpoint_pairs_dbscan(
            breakpoint_pairs, cluster_distance=5000, min_read_support=min_read_support
        )
        
        # Filter for breakpoint pairs with sufficient read support (already filtered in DBSCAN, but keep for safety)
        supported_pairs = [
            p for p in clustered_pairs 
            if p["read_count"] >= min_read_support
        ]
        
        if len(supported_pairs) > 0:
            logger.info(
                f"Found {len(supported_pairs)} master BED breakpoint pairs "
                f"supported by >= {min_read_support} reads (out of {len(clustered_pairs)} clustered pairs, "
                f"{len(breakpoint_pairs)} original pairs)"
            )
            
            # Extract breakpoint coordinates for both primary and supplementary regions
            coord_extract_start = time.time()
            for pair in supported_pairs:
                # Add primary region breakpoint
                if pair["primary_start"] > 0 and pair["primary_end"] > pair["primary_start"]:
                    breakpoints.append({
                        "chromosome": pair["primary_chrom"],
                        "start": pair["primary_start"],
                        "end": pair["primary_end"],
                        "read_count": pair["read_count"],
                        "source": "master_bed"
                    })
                
                # Add supplementary region breakpoint (new target region)
                if pair["supp_start"] > 0 and pair["supp_end"] > pair["supp_start"]:
                    breakpoints.append({
                        "chromosome": pair["supp_chrom"],
                        "start": pair["supp_start"],
                        "end": pair["supp_end"],
                        "read_count": pair["read_count"],
                        "source": "master_bed"
                    })
        else:
            # Log details for debugging if no supported pairs found
            if clustered_pairs:
                max_read_count = max(p["read_count"] for p in clustered_pairs)
                logger.debug(
                    f"No master BED breakpoint pairs found with >= {min_read_support} read support. "
                    f"Found {len(clustered_pairs)} clustered pairs from {len(breakpoint_pairs)} original pairs. "
                    f"Maximum read count in any cluster: {max_read_count}"
                )
            else:
                logger.debug(
                    f"No clustered breakpoint pairs created from {len(breakpoint_pairs)} original pairs"
                )
    
    return breakpoints


def _extract_master_bed_breakpoints_incremental(
    fusion_metadata: FusionMetadata,
    work_dir: str,
    sample_id: str,
    new_master_bed_files: Set[str],
    min_read_support: int = 3,
) -> List[Dict[str, Any]]:
    """
    Incrementally extract master BED breakpoints by processing only NEW staging files
    and merging with existing breakpoints from the previous BED file.
    
    This avoids reprocessing the entire accumulated dataset (which can be 1M+ rows)
    by only processing the new data and merging results.
    
    Args:
        fusion_metadata: FusionMetadata object
        work_dir: Working directory
        sample_id: Sample ID
        new_master_bed_files: Set of new staging file paths to process
        min_read_support: Minimum read support threshold
        
    Returns:
        List of all breakpoint dictionaries (existing + new)
    """
    incremental_start = time.time()
    
    # Load existing breakpoints from the previous BED file (if it exists)
    existing_breakpoints = []
    analysis_counter = _load_analysis_counter(sample_id, work_dir)
    sample_dir = os.path.join(work_dir, sample_id)
    bed_dir = os.path.join(sample_dir, "bed_files")
    previous_bed_file = os.path.join(bed_dir, f"master_bed_breakpoints_{analysis_counter:03d}.bed")
    
    if os.path.exists(previous_bed_file):
        # Parse existing BED file to get breakpoint coordinates
        # BED format: chrom, start, end, name, score, strand
        try:
            with open(previous_bed_file, "r") as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split("\t")
                        if len(parts) >= 3:
                            # Extract coordinates (BED file has expanded regions, need to get midpoint)
                            chrom = parts[0]
                            bed_start = int(parts[1])
                            bed_end = int(parts[2])
                            # Calculate midpoint (original breakpoint location)
                            midpoint = (bed_start + bed_end) // 2
                            # Estimate original breakpoint (region is +/- bin_width around breakpoint)
                            bin_width = _get_cnv_bin_width(work_dir, sample_id)
                            original_start = max(0, midpoint - bin_width)
                            original_end = midpoint + bin_width
                            existing_breakpoints.append({
                                "chromosome": chrom,
                                "start": original_start,
                                "end": original_end,
                                "read_count": 1,  # Unknown from BED file, use minimum
                                "source": "master_bed"
                            })
        except Exception as e:
            logger.warning(f"Could not load existing breakpoints from {previous_bed_file}: {e}")
    
    # Process only NEW staging files
    new_data_start = time.time()
    new_master_bed_dfs = []
    for staging_file in new_master_bed_files:
        try:
            df = pd.read_parquet(staging_file)
            if not df.empty:
                new_master_bed_dfs.append(df)
        except Exception as e:
            logger.warning(f"Error loading staging file {os.path.basename(staging_file)}: {e}")
    
        if not new_master_bed_dfs:
            return existing_breakpoints
    
    # Combine new staging files
    new_master_bed_df = pd.concat(new_master_bed_dfs, ignore_index=True) if len(new_master_bed_dfs) > 1 else new_master_bed_dfs[0]
    
    # Create temporary FusionMetadata with only new data for extraction
    temp_metadata = FusionMetadata(
        sample_id=sample_id,
        file_path="incremental",
        analysis_timestamp=time.time(),
        target_panel=fusion_metadata.target_panel,
        fusion_data={"master_bed_candidates": new_master_bed_df},
    )
    
    # Extract breakpoints from new data only
    new_breakpoints = _extract_master_bed_breakpoints(temp_metadata, work_dir=work_dir, min_read_support=min_read_support)
    
    # Merge new breakpoints with existing ones
    # For now, simple merge (could be optimized to cluster nearby breakpoints)
    all_breakpoints = existing_breakpoints + new_breakpoints
    
    return all_breakpoints


def _generate_master_bed_breakpoint_bed(
    sample_id: str,
    fusion_metadata: FusionMetadata,
    work_dir: str,
    new_master_bed_files: Optional[Set[str]] = None,
) -> None:
    """
    Generate BED file for master BED breakpoints (new target regions from supplementary alignments).
    Creates regions with +/- 1 bin_width around breakpoints.
    Uses counter-based naming consistent with other BED files.
    
    Args:
        sample_id: Sample ID
        fusion_metadata: FusionMetadata object with master BED candidates
        work_dir: Working directory
    """
    try:
        # Extract master BED breakpoints
        # For large datasets, use incremental approach: only process new staging files
        # and merge with existing breakpoints from previous BED file
        master_bed_candidates = _load_fusion_candidates_parquet("master_bed_candidates", work_dir, sample_id)
        master_bed_size = len(master_bed_candidates) if master_bed_candidates is not None and not master_bed_candidates.empty else 0
        
        # For large datasets (>200k rows), use incremental extraction to avoid 2+ minute delays
        # Only process NEW staging files and merge with existing breakpoints
        if master_bed_size > 200000 and new_master_bed_files:
            master_bed_breakpoints = _extract_master_bed_breakpoints_incremental(
                fusion_metadata, work_dir, sample_id, new_master_bed_files
            )
        else:
            # For smaller datasets, process all data (faster for small datasets)
            master_bed_breakpoints = _extract_master_bed_breakpoints(fusion_metadata, work_dir=work_dir)
        
        if not master_bed_breakpoints:
            logger.debug("No master BED breakpoints found - skipping BED file generation")
            return
        
        # Get bin_width from CNV analysis if available, otherwise use default
        bin_width = _get_cnv_bin_width(work_dir, sample_id)
        
        # Load analysis counter for consistent naming
        analysis_counter = _load_analysis_counter(sample_id, work_dir)
        
        # Create bed_files directory if it doesn't exist
        sample_dir = os.path.join(work_dir, sample_id)
        bed_dir = os.path.join(sample_dir, "bed_files")
        os.makedirs(bed_dir, exist_ok=True)
        
        # Generate BED file for master BED breakpoints with counter-based naming
        master_bed_bp_file = os.path.join(bed_dir, f"master_bed_breakpoints_{analysis_counter:03d}.bed")
        
        # Sort breakpoints by chromosome, then start position, then end position
        def sort_key(bp):
            chrom = bp["chromosome"]
            start = bp["start"]
            end = bp["end"]
            
            # Extract chromosome number for numeric sorting
            chrom_num = None
            if chrom.startswith("chr"):
                chrom_suffix = chrom[3:]
                if chrom_suffix == "X":
                    chrom_num = 100
                elif chrom_suffix == "Y":
                    chrom_num = 101
                elif chrom_suffix == "M":
                    chrom_num = 102
                else:
                    try:
                        chrom_num = int(chrom_suffix)
                    except ValueError:
                        chrom_num = 200  # Put non-standard chromosomes at the end
            else:
                chrom_num = 200  # Put non-standard chromosomes at the end
            
            return (chrom_num, start, end)
        
        sorted_breakpoints = sorted(master_bed_breakpoints, key=sort_key)
        
        with open(master_bed_bp_file, "w") as f:
            for bp in sorted_breakpoints:
                chrom = bp["chromosome"]
                start = bp["start"]
                end = bp["end"]
                
                # Calculate midpoint for breakpoint
                midpoint = (start + end) // 2
                
                # Create region +/- 1 bin_width around breakpoint
                region_start = max(0, midpoint - bin_width)
                region_end = midpoint + bin_width
                
                # Write BED entry: chrom, start, end, name (master_bed-breakpoint)
                name = "master_bed-breakpoint"
                f.write(f"{chrom}\t{region_start}\t{region_end}\t{name}\t0\t.\n")
        
        # Log summary of read support
        if master_bed_breakpoints:
            read_counts = [bp.get("read_count", 0) for bp in master_bed_breakpoints]
            min_reads = min(read_counts) if read_counts else 0
            max_reads = max(read_counts) if read_counts else 0
            avg_reads = sum(read_counts) / len(read_counts) if read_counts else 0
            logger.info(
                f"Generated master BED breakpoint BED file: {master_bed_bp_file} with {len(master_bed_breakpoints)} breakpoints "
                f"(read support: min={min_reads}, max={max_reads}, avg={avg_reads:.1f})"
            )
        else:
            logger.info(f"Generated master BED breakpoint BED file: {master_bed_bp_file} with 0 breakpoints")
        
        # Generate master BED events summary CSV for GUI (pre-computed to avoid blocking UI)
        _generate_master_bed_events_summary(sample_id, fusion_metadata, work_dir)
        
    except Exception as e:
        logger.warning(f"Error generating master BED breakpoint BED file: {e}")
        import traceback
        logger.debug(f"Traceback: {traceback.format_exc()}")


def _generate_master_bed_events_summary(
    sample_id: str,
    fusion_metadata: FusionMetadata,
    work_dir: str,
    min_read_support: int = 3,
    min_mapq: int = 50,
    cluster_distance: int = 5000,
) -> None:
    """
    Generate master BED events summary CSV file for GUI display.
    
    This pre-computes the expensive breakpoint pair clustering so the GUI
    can just read the results instead of computing them on the UI thread.
    
    Uses the same breakpoint pair extraction logic as _extract_master_bed_breakpoints:
    - Identifies reads with both primary and supplementary alignments
    - Creates breakpoint pairs (primary + supplementary for each read)
    - Clusters similar breakpoint pairs (both primary and supplementary locations must be close)
    - Writes both primary and supplementary regions as separate events to CSV
    
    Args:
        sample_id: Sample ID
        fusion_metadata: FusionMetadata object with master BED candidates
        work_dir: Working directory
        min_read_support: Minimum number of reads required (default: 3)
        min_mapq: Minimum mapping quality for supplementary mappings (default: 50)
        cluster_distance: Maximum distance for clustering breakpoint pairs (default: 5000bp)
    """
    try:
        # Load master BED candidates from Parquet
        master_bed_df = _load_fusion_candidates_parquet("master_bed_candidates", work_dir, sample_id)
        
        if master_bed_df is None or master_bed_df.empty:
            logger.debug("No master BED candidates found - skipping events summary generation")
            # Write empty file so GUI knows there's no data
            summary_file = os.path.join(work_dir, sample_id, "master_bed_events_summary.csv")
            pd.DataFrame().to_csv(summary_file, index=False)
            return
        
        # Check required columns
        if "read_id" not in master_bed_df.columns:
            logger.debug("Missing required column (read_id) in master BED candidates")
            return
        
        # Filter to only high-quality supplementary mappings (MapQ >= min_mapq)
        if "mapping_quality" in master_bed_df.columns and "col4" in master_bed_df.columns:
            high_quality_mask = (
                (master_bed_df["col4"] == "master_bed_region") |
                (master_bed_df["mapping_quality"] >= min_mapq)
            )
            master_bed_df = master_bed_df[high_quality_mask].copy()
            
            if master_bed_df.empty:
                logger.debug(f"No high-quality mappings found (MapQ >= {min_mapq})")
                summary_file = os.path.join(work_dir, sample_id, "master_bed_events_summary.csv")
                pd.DataFrame().to_csv(summary_file, index=False)
                return
        
        # Separate primary and supplementary alignments using col4
        if "col4" in master_bed_df.columns:
            primary_df = master_bed_df[master_bed_df["col4"] == "master_bed_region"].copy()
            supplementary_df = master_bed_df[master_bed_df["col4"] == "master_bed_supplementary"].copy()
        else:
            if "is_supplementary" not in master_bed_df.columns:
                logger.debug("Missing required columns (col4 or is_supplementary)")
                return
            primary_df = master_bed_df[master_bed_df["is_supplementary"] == False].copy()
            supplementary_df = master_bed_df[master_bed_df["is_supplementary"] == True].copy()
        
        if primary_df.empty or supplementary_df.empty:
            logger.debug("Need both primary and supplementary alignments")
            summary_file = os.path.join(work_dir, sample_id, "master_bed_events_summary.csv")
            pd.DataFrame().to_csv(summary_file, index=False)
            return
        
        # Find reads with both primary and supplementary alignments
        primary_read_ids = set(primary_df["read_id"].unique())
        supplementary_read_ids = set(supplementary_df["read_id"].unique())
        reads_with_both = primary_read_ids & supplementary_read_ids
        
        if not reads_with_both:
            logger.debug("No reads have both primary and supplementary alignments")
            summary_file = os.path.join(work_dir, sample_id, "master_bed_events_summary.csv")
            pd.DataFrame().to_csv(summary_file, index=False)
            return
        
        # Filter to only reads with both
        primary_filtered = primary_df[primary_df["read_id"].isin(reads_with_both)].copy()
        supplementary_filtered = supplementary_df[supplementary_df["read_id"].isin(reads_with_both)].copy()
        
        # Build breakpoint pairs
        breakpoint_pairs = []
        primary_by_read = primary_filtered.groupby("read_id", observed=True)
        supplementary_by_read = supplementary_filtered.groupby("read_id", observed=True)
        
        for read_id in reads_with_both:
            read_primaries = primary_by_read.get_group(read_id) if read_id in primary_by_read.groups else pd.DataFrame()
            read_supplementaries = supplementary_by_read.get_group(read_id) if read_id in supplementary_by_read.groups else pd.DataFrame()
            
            if read_primaries.empty or read_supplementaries.empty:
                continue
            
            # Use vectorized breakpoint pair creation (much faster than nested loops)
            pairs = _create_breakpoint_pairs_vectorized(read_primaries, read_supplementaries, read_id)
            breakpoint_pairs.extend(pairs)
        
        if not breakpoint_pairs:
            logger.debug("No breakpoint pairs created")
            summary_file = os.path.join(work_dir, sample_id, "master_bed_events_summary.csv")
            pd.DataFrame().to_csv(summary_file, index=False)
            return
        
        # Cluster similar breakpoint pairs using DBSCAN (much faster than O(n²) nested loops)
        clustered_pairs = _cluster_breakpoint_pairs_dbscan(
            breakpoint_pairs, cluster_distance=cluster_distance, min_read_support=min_read_support
        )
        
        # Filter for breakpoint pairs with sufficient read support (already filtered in DBSCAN, but keep for safety)
        supported_pairs = [
            p for p in clustered_pairs 
            if p["read_count"] >= min_read_support
        ]
        
        if not supported_pairs:
            logger.debug(f"No breakpoint pairs found with >= {min_read_support} read support")
            summary_file = os.path.join(work_dir, sample_id, "master_bed_events_summary.csv")
            pd.DataFrame().to_csv(summary_file, index=False)
            return
        
        # Convert to DataFrame format - include both primary and supplementary regions
        events = []
        for pair in supported_pairs:
            if pair["primary_start"] > 0 and pair["primary_end"] > pair["primary_start"]:
                events.append({
                    "chromosome": pair["primary_chrom"],
                    "start": pair["primary_start"],
                    "end": pair["primary_end"],
                    "event_type": "breakpoint_pair_primary",
                    "read_count": pair["read_count"],
                    "avg_mapping_quality": round(pair["primary_avg_mapq"], 1),
                    "avg_mapping_span": round(pair["primary_avg_span"], 0),
                })
            
            if pair["supp_start"] > 0 and pair["supp_end"] > pair["supp_start"]:
                events.append({
                    "chromosome": pair["supp_chrom"],
                    "start": pair["supp_start"],
                    "end": pair["supp_end"],
                    "event_type": "breakpoint_pair_supplementary",
                    "read_count": pair["read_count"],
                    "avg_mapping_quality": round(pair["supp_avg_mapq"], 1),
                    "avg_mapping_span": round(pair["supp_avg_span"], 0),
                })
        
        if not events:
            summary_file = os.path.join(work_dir, sample_id, "master_bed_events_summary.csv")
            pd.DataFrame().to_csv(summary_file, index=False)
            return
        
        result_df = pd.DataFrame(events)
        
        # Merge nearby duplicate events on the same chromosome
        # This prevents showing multiple very similar events that are within cluster_distance
        result_df = _merge_nearby_events(result_df, cluster_distance=cluster_distance)
        
        result_df = result_df.sort_values(
            ["read_count", "chromosome", "start"],
            ascending=[False, True, True]
        )
        
        # Write summary CSV file
        summary_file = os.path.join(work_dir, sample_id, "master_bed_events_summary.csv")
        result_df.to_csv(summary_file, index=False)
        logger.info(f"Generated master BED events summary: {summary_file} with {len(result_df)} events")
        
    except Exception as e:
        logger.warning(f"Error generating master BED events summary: {e}")
        import traceback
        logger.debug(f"Traceback: {traceback.format_exc()}")


def _generate_fusion_breakpoint_bed(
    sample_id: str,
    fusion_metadata: FusionMetadata,
    work_dir: str,
) -> None:
    """
    Generate BED file for fusion breakpoints with +/- 1 bin_width regions.
    Only includes fusions that meet the minimum read support threshold.
    Uses counter-based naming consistent with other BED files.
    
    Args:
        sample_id: Sample ID
        fusion_metadata: FusionMetadata object with fusion candidates
        work_dir: Working directory
    """
    try:
        # Extract fusion breakpoints (already filtered by threshold)
        fusion_breakpoints = _extract_fusion_breakpoints(fusion_metadata)
        
        if not fusion_breakpoints:
            logger.debug("No fusion breakpoints found - skipping BED file generation")
            return
        
        # Get bin_width from CNV analysis if available, otherwise use default
        bin_width = _get_cnv_bin_width(work_dir, sample_id)
        
        # Load analysis counter for consistent naming
        analysis_counter = _load_analysis_counter(sample_id, work_dir)
        
        # Create bed_files directory if it doesn't exist
        sample_dir = os.path.join(work_dir, sample_id)
        bed_dir = os.path.join(sample_dir, "bed_files")
        os.makedirs(bed_dir, exist_ok=True)
        
        # Generate BED file for fusion breakpoints with counter-based naming
        fusion_bed_file = os.path.join(bed_dir, f"fusion_breakpoints_{analysis_counter:03d}.bed")
        
        # Sort breakpoints by chromosome, then start position, then end position
        # Handle chromosome sorting (chr1, chr2, ..., chr10, chr11, ..., chrX, chrY, chrM)
        def sort_key(bp):
            chrom = bp["chromosome"]
            start = bp["start"]
            end = bp["end"]
            
            # Extract chromosome number for numeric sorting
            chrom_num = None
            if chrom.startswith("chr"):
                chrom_suffix = chrom[3:]
                if chrom_suffix == "X":
                    chrom_num = 100
                elif chrom_suffix == "Y":
                    chrom_num = 101
                elif chrom_suffix == "M":
                    chrom_num = 102
                else:
                    try:
                        chrom_num = int(chrom_suffix)
                    except ValueError:
                        chrom_num = 200  # Put non-standard chromosomes at the end
            else:
                chrom_num = 200  # Put non-standard chromosomes at the end
            
            return (chrom_num, start, end)
        
        sorted_breakpoints = sorted(fusion_breakpoints, key=sort_key)
        
        with open(fusion_bed_file, "w") as f:
            for bp in sorted_breakpoints:
                chrom = bp["chromosome"]
                start = bp["start"]
                end = bp["end"]
                gene = bp.get("gene", "Unknown")
                source = bp.get("source", "unknown")
                
                # Calculate midpoint for breakpoint
                midpoint = (start + end) // 2
                
                # Create region +/- 1 bin_width around breakpoint
                region_start = max(0, midpoint - bin_width)
                region_end = midpoint + bin_width
                
                # Write BED entry: chrom, start, end, name (gene-source)
                name = f"{gene}-{source}"
                f.write(f"{chrom}\t{region_start}\t{region_end}\t{name}\n")
        
        logger.info(f"Generated fusion breakpoint BED file: {fusion_bed_file} with {len(fusion_breakpoints)} breakpoints")
        
    except Exception as e:
        logger.warning(f"Error generating fusion breakpoint BED file: {e}")
        import traceback
        logger.debug(f"Traceback: {traceback.format_exc()}")


def find_and_process_bam_files(root_dir):
    """Recursively walks through root_dir and processes each .bam file."""
    root_path = Path(root_dir)

    if not root_path.exists() or not root_path.is_dir():
        raise ValueError(f"{root_dir} is not a valid directory")

    for bam_file in root_path.rglob("*.bam"):
        process_bam_file(bam_file)


# =============================================================================
# FUSION METADATA PERSISTENCE FUNCTIONS
# =============================================================================


def _get_parquet_paths(work_dir: str, sample_id: str) -> Dict[str, str]:
    """
    Get file paths for Parquet storage of fusion candidates.
    
    Args:
        work_dir: Working directory
        sample_id: Sample identifier
        
    Returns:
        Dictionary mapping candidate type to file path
    """
    sample_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)
    
    return {
        "target_candidates": os.path.join(sample_dir, f"{sample_id}_target_candidates.parquet"),
        "genome_wide_candidates": os.path.join(sample_dir, f"{sample_id}_genome_wide_candidates.parquet"),
        "master_bed_candidates": os.path.join(sample_dir, f"{sample_id}_master_bed_candidates.parquet"),
    }


def _get_parquet_dataset_dir(work_dir: str, sample_id: str, candidate_type: str) -> str:
    """Get directory for append-only Parquet dataset storage."""
    sample_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)
    return os.path.join(sample_dir, f"{candidate_type}_dataset")


def _migrate_legacy_parquet_to_dataset(work_dir: str, sample_id: str) -> None:
    """Migrate legacy single-file Parquet into append-only dataset layout."""
    for candidate_type in ["target_candidates", "genome_wide_candidates", "master_bed_candidates"]:
        dataset_dir = _get_parquet_dataset_dir(work_dir, sample_id, candidate_type)
        dataset_files = glob.glob(os.path.join(dataset_dir, "*.parquet"))
        if dataset_files:
            continue

        legacy_path = _get_parquet_paths(work_dir, sample_id).get(candidate_type)
        if not legacy_path or not os.path.exists(legacy_path):
            continue

        os.makedirs(dataset_dir, exist_ok=True)
        part_path = os.path.join(dataset_dir, "part_legacy.parquet")
        if not os.path.exists(part_path):
            try:
                shutil.copy2(legacy_path, part_path)
                logger.info(
                    f"Migrated legacy {candidate_type} parquet to dataset: {part_path}"
                )
            except Exception as e:
                logger.warning(
                    f"Could not migrate legacy {candidate_type} parquet: {e}"
                )


def _extract_staging_batch_id(staging_files: List[str]) -> int:
    """Extract a reasonable batch id from staging filenames."""
    counters = []
    for path in staging_files:
        name = os.path.basename(path)
        try:
            counter_str = name.split("_", 1)[1].split(".", 1)[0]
            counters.append(int(counter_str))
        except (IndexError, ValueError):
            continue
    if counters:
        return max(counters)
    return int(time.time())


def _append_fusion_candidates_parquet(
    candidates_df: pd.DataFrame,
    candidate_type: str,
    work_dir: str,
    sample_id: str,
    batch_id: int,
) -> None:
    """Append candidates to an append-only Parquet dataset."""
    if candidates_df is None or candidates_df.empty:
        return
    dataset_dir = _get_parquet_dataset_dir(work_dir, sample_id, candidate_type)
    os.makedirs(dataset_dir, exist_ok=True)
    part_path = os.path.join(dataset_dir, f"part_{batch_id:06d}.parquet")
    if os.path.exists(part_path):
        part_path = os.path.join(dataset_dir, f"part_{batch_id:06d}_{int(time.time())}.parquet")
    candidates_df.to_parquet(part_path, index=False, engine="pyarrow")
    logger.debug(f"Appended {len(candidates_df)} {candidate_type} to {part_path}")


def _get_counts_path(work_dir: str, sample_id: str) -> str:
    sample_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)
    return os.path.join(sample_dir, "fusion_counts.json")


def _count_parquet_rows(paths: List[str]) -> int:
    total = 0
    for path in paths:
        try:
            total += pq.ParquetFile(path).metadata.num_rows
        except Exception as e:
            logger.debug(f"Could not count rows in {path}: {e}")
    return total


def _load_fusion_counts(work_dir: str, sample_id: str) -> Dict[str, int]:
    counts_path = _get_counts_path(work_dir, sample_id)
    if os.path.exists(counts_path):
        try:
            with open(counts_path, "r") as f:
                return json.load(f)
        except Exception as e:
            logger.debug(f"Could not read fusion_counts for {sample_id}: {e}")
    counts = {}
    for candidate_type in ["target_candidates", "genome_wide_candidates", "master_bed_candidates"]:
        dataset_dir = _get_parquet_dataset_dir(work_dir, sample_id, candidate_type)
        dataset_files = sorted(glob.glob(os.path.join(dataset_dir, "*.parquet")))
        if dataset_files:
            counts[candidate_type] = _count_parquet_rows(dataset_files)
            continue
        legacy_path = _get_parquet_paths(work_dir, sample_id).get(candidate_type)
        if legacy_path and os.path.exists(legacy_path):
            counts[candidate_type] = _count_parquet_rows([legacy_path])
        else:
            counts[candidate_type] = 0
    return counts


def _save_fusion_counts(work_dir: str, sample_id: str, counts: Dict[str, int]) -> None:
    counts_path = _get_counts_path(work_dir, sample_id)
    try:
        with open(counts_path, "w") as f:
            json.dump(counts, f, indent=2)
    except Exception as e:
        logger.debug(f"Could not write fusion_counts for {sample_id}: {e}")


def _save_fusion_candidates_parquet(
    candidates_df: pd.DataFrame,
    candidate_type: str,
    work_dir: str,
    sample_id: str,
) -> None:
    """
    Save fusion candidates DataFrame to Parquet file.
    
    Args:
        candidates_df: DataFrame to save
        candidate_type: Type of candidates ('target_candidates', 'genome_wide_candidates', 'master_bed_candidates')
        work_dir: Working directory
        sample_id: Sample identifier
    """
    if candidates_df is None or candidates_df.empty:
        return
    
    try:
        parquet_paths = _get_parquet_paths(work_dir, sample_id)
        parquet_path = parquet_paths.get(candidate_type)
        
        if parquet_path:
            candidates_df.to_parquet(parquet_path, index=False, engine='pyarrow')
            logger.debug(f"Saved {len(candidates_df)} {candidate_type} to {parquet_path}")
    except Exception as e:
        logger.error(f"Error saving {candidate_type} to Parquet: {e}")


def _load_fusion_candidates_parquet(
    candidate_type: str,
    work_dir: str,
    sample_id: str,
) -> Optional[pd.DataFrame]:
    """
    Load fusion candidates from Parquet file.
    
    Args:
        candidate_type: Type of candidates ('target_candidates', 'genome_wide_candidates', 'master_bed_candidates')
        work_dir: Working directory
        sample_id: Sample identifier
        
    Returns:
        DataFrame if file exists and is readable, None otherwise
    """
    try:
        _migrate_legacy_parquet_to_dataset(work_dir, sample_id)
        dataset_dir = _get_parquet_dataset_dir(work_dir, sample_id, candidate_type)
        dataset_files = glob.glob(os.path.join(dataset_dir, "*.parquet"))
        if dataset_files:
            df = pd.read_parquet(dataset_dir, engine="pyarrow")
            if df.empty:
                logger.debug(f"Parquet dataset {dataset_dir} exists but is empty")
                return None
            logger.debug(f"Loaded {len(df)} {candidate_type} from dataset {dataset_dir}")
            return df

        parquet_paths = _get_parquet_paths(work_dir, sample_id)
        parquet_path = parquet_paths.get(candidate_type)
        if parquet_path and os.path.exists(parquet_path):
            df = pd.read_parquet(parquet_path, engine="pyarrow")
            if df.empty:
                logger.debug(f"Parquet file {parquet_path} exists but is empty")
                return None  # Treat empty file same as missing file
            logger.debug(f"Loaded {len(df)} {candidate_type} from {parquet_path}")
            return df
    except Exception as e:
        logger.warning(f"Error loading {candidate_type} from Parquet: {e}")
    
    return None


def search_fusion_candidates_by_read_id(
    work_dir: str,
    sample_id: str,
    read_ids: Union[str, List[str]],
    candidate_types: Optional[List[str]] = None,
    report_totals: bool = True,
) -> Dict[str, pd.DataFrame]:
    """
    Search for one or more read IDs across all fusion candidate Parquet files.
    
    Args:
        work_dir: Working directory
        sample_id: Sample identifier
        read_ids: Single read ID (str) or list of read IDs to search for
        candidate_types: Optional list of candidate types to search. If None, searches all types.
                        Valid types: 'target_candidates', 'genome_wide_candidates', 'master_bed_candidates'
        report_totals: If True, print summary of total unique reads in each candidate type (default: True)
        
    Returns:
        Dictionary mapping candidate_type to DataFrame containing matching rows (empty DataFrame if no matches)
    """
    # Normalize input: convert single string to list
    if isinstance(read_ids, str):
        read_ids = [read_ids]
    elif not isinstance(read_ids, list):
        raise ValueError(f"read_ids must be a string or list of strings, got {type(read_ids)}")
    
    # Convert to set for efficient lookup
    read_ids_set = set(read_ids)
    
    if candidate_types is None:
        candidate_types = ["target_candidates", "genome_wide_candidates", "master_bed_candidates"]
    
    results = {}
    totals_summary = {}
    
    for candidate_type in candidate_types:
        try:
            # Load the Parquet file
            df = _load_fusion_candidates_parquet(candidate_type, work_dir, sample_id)
            
            if df is not None and not df.empty:
                if "read_id" not in df.columns:
                    logger.warning(f"read_id column not found in {candidate_type} DataFrame. Columns: {list(df.columns)}")
                    results[candidate_type] = pd.DataFrame()
                    totals_summary[candidate_type] = {"total_entries": len(df), "unique_reads": 0, "file_exists": True}
                    continue
                
                # Count unique reads
                unique_read_count = df["read_id"].nunique()
                total_entries = len(df)
                totals_summary[candidate_type] = {
                    "total_entries": total_entries,
                    "unique_reads": unique_read_count,
                    "file_exists": True
                }
                
                # Filter for any of the read IDs in the list
                matches = df[df["read_id"].isin(read_ids_set)]
                if not matches.empty:
                    results[candidate_type] = matches
                    found_reads = matches["read_id"].unique()
                    logger.info(
                        f"Found {len(matches)} entries for {len(found_reads)} read_id(s) in {candidate_type}: "
                        f"{', '.join(found_reads[:5])}{'...' if len(found_reads) > 5 else ''}"
                    )
                else:
                    # Check if any read IDs are similar (case/whitespace differences)
                    df_read_ids_set = set(df["read_id"].unique())
                    read_ids_lower = {rid.lower().strip() for rid in read_ids_set}
                    df_read_ids_lower = {rid.lower().strip() for rid in df_read_ids_set}
                    if read_ids_lower & df_read_ids_lower:
                        logger.warning(
                            f"Found similar read IDs in {candidate_type} (case/whitespace differences). "
                            f"Searching for: {read_ids_set}, found similar: {read_ids_lower & df_read_ids_lower}"
                        )
                    results[candidate_type] = pd.DataFrame()  # Empty DataFrame to indicate searched but no matches
            else:
                results[candidate_type] = pd.DataFrame()  # Empty DataFrame if file doesn't exist or is empty
                totals_summary[candidate_type] = {"total_entries": 0, "unique_reads": 0, "file_exists": False}
        except Exception as e:
            logger.warning(f"Error searching {candidate_type} for read_ids {read_ids}: {e}")
            import traceback
            logger.debug(f"Traceback: {traceback.format_exc()}")
            results[candidate_type] = pd.DataFrame()
            totals_summary[candidate_type] = {"total_entries": 0, "unique_reads": 0, "file_exists": False, "error": str(e)}
    
    # Log summary of totals if requested
    if report_totals:
        logger.info("Total unique reads in each candidate type:")
        for candidate_type, stats in totals_summary.items():
            if stats.get("file_exists", False):
                logger.info(
                    f"  {candidate_type}: {stats['unique_reads']:,} unique reads "
                    f"({stats['total_entries']:,} total entries)"
                )
            elif stats.get("error"):
                logger.warning(f"  {candidate_type}: Error - {stats['error']}")
            else:
                logger.debug(f"  {candidate_type}: File empty or doesn't exist")
    
    return results


def diagnose_master_bed_breakpoint_extraction(
    work_dir: str,
    sample_id: str,
    read_ids: Optional[List[str]] = None,
    min_read_support: int = 3,
    cluster_distance: int = 5000,
) -> None:
    """
    Diagnostic function to debug why master BED breakpoints aren't being detected.
    
    This function loads the master_bed_candidates data and traces through the entire
    breakpoint extraction logic to identify where the process fails.
    
    Args:
        work_dir: Working directory
        sample_id: Sample identifier
        read_ids: Optional list of specific read IDs to focus on (if None, analyzes all reads)
        min_read_support: Minimum read support threshold (default: 3)
        cluster_distance: Clustering distance threshold in bp (default: 5000)
    """
    logger.debug("\n" + "="*80)
    logger.debug(f"DIAGNOSTIC: Master BED Breakpoint Extraction")
    logger.debug(f"Sample: {sample_id}")
    logger.debug(f"Work dir: {work_dir}")
    logger.debug("="*80 + "\n")
    
    # Load master_bed_candidates
    logger.debug("[STEP 1] Loading master_bed_candidates from Parquet...")
    master_bed_df = _load_fusion_candidates_parquet("master_bed_candidates", work_dir, sample_id)
    
    if master_bed_df is None or master_bed_df.empty:
        logger.debug(f"  ERROR: No master_bed_candidates found or file is empty")
        logger.debug(f"  Parquet path: {_get_parquet_paths(work_dir, sample_id).get('master_bed_candidates', 'N/A')}")
        return
    
    logger.debug(f"  Loaded {len(master_bed_df)} total entries")
    logger.debug(f"  Columns: {list(master_bed_df.columns)}")
    
    # Filter to specific read IDs if provided
    if read_ids:
        read_ids_set = set(read_ids)
        master_bed_df = master_bed_df[master_bed_df["read_id"].isin(read_ids_set)].copy()
        logger.debug(f"\n[FILTER] Filtered to {len(master_bed_df)} entries for {len(read_ids)} read IDs")
        if master_bed_df.empty:
            logger.debug(f"  ERROR: None of the specified read IDs found in master_bed_candidates")
            logger.debug(f"  Looking for: {read_ids}")
            all_read_ids = _load_fusion_candidates_parquet("master_bed_candidates", work_dir, sample_id)
            if all_read_ids is not None and not all_read_ids.empty:
                unique_reads = all_read_ids["read_id"].unique()[:10]
                logger.debug(f"  Example read IDs in file: {list(unique_reads)}")
            return
    
    # Check required columns
    logger.debug("\n[STEP 2] Checking required columns...")
    required_cols = ["read_id", "reference_id", "reference_start", "reference_end"]
    missing_cols = [col for col in required_cols if col not in master_bed_df.columns]
    if missing_cols:
        logger.debug(f"  ERROR: Missing required columns: {missing_cols}")
        return
    logger.debug(f"  All required columns present")
    
    # Check classification columns
    logger.debug("\n[STEP 3] Checking alignment classification...")
    has_col4 = "col4" in master_bed_df.columns
    has_is_supp = "is_supplementary" in master_bed_df.columns
    
    if has_col4:
        logger.debug(f"  col4 column found")
        col4_values = master_bed_df["col4"].value_counts()
        logger.debug(f"  col4 value counts:\n{col4_values}")
    else:
        logger.debug(f"  WARNING: col4 column not found")
    
    if has_is_supp:
        logger.debug(f"  is_supplementary column found")
        is_supp_counts = master_bed_df["is_supplementary"].value_counts()
        logger.debug(f"  is_supplementary value counts:\n{is_supp_counts}")
    else:
        logger.debug(f"  WARNING: is_supplementary column not found")
    
    if not has_col4 and not has_is_supp:
        logger.debug(f"  ERROR: Neither col4 nor is_supplementary column found - cannot classify alignments")
        return
    
    # Separate primary and supplementary
    logger.debug("\n[STEP 4] Separating primary and supplementary alignments...")
    if has_col4:
        primary_df = master_bed_df[master_bed_df["col4"] == "master_bed_region"].copy()
        supplementary_df = master_bed_df[master_bed_df["col4"] == "master_bed_supplementary"].copy()
        logger.debug(f"  Using col4 for classification")
    else:
        primary_df = master_bed_df[master_bed_df["is_supplementary"] == False].copy()
        supplementary_df = master_bed_df[master_bed_df["is_supplementary"] == True].copy()
        logger.debug(f"  Using is_supplementary for classification (fallback)")
    
    logger.debug(f"  Primary alignments: {len(primary_df)}")
    logger.debug(f"  Supplementary alignments: {len(supplementary_df)}")
    
    if primary_df.empty:
        logger.debug(f"  ERROR: No primary alignments found")
        return
    
    if supplementary_df.empty:
        logger.debug(f"  ERROR: No supplementary alignments found")
        return
    
    # Check for reads with both
    logger.debug("\n[STEP 5] Finding reads with both primary and supplementary alignments...")
    primary_read_ids = set(primary_df["read_id"].unique())
    supplementary_read_ids = set(supplementary_df["read_id"].unique())
    reads_with_both = primary_read_ids & supplementary_read_ids
    
    logger.debug(f"  Reads with primary only: {len(primary_read_ids - supplementary_read_ids)}")
    logger.debug(f"  Reads with supplementary only: {len(supplementary_read_ids - primary_read_ids)}")
    logger.debug(f"  Reads with both: {len(reads_with_both)}")
    
    if not reads_with_both:
        logger.debug(f"  ERROR: No reads have both primary and supplementary alignments")
        if read_ids:
            logger.debug(f"\n  Checking specific read IDs:")
            for rid in read_ids:
                has_primary = rid in primary_read_ids
                has_supp = rid in supplementary_read_ids
                logger.debug(f"    {rid}: primary={has_primary}, supplementary={has_supp}")
        return
    
    # Show example reads
    if read_ids:
        logger.debug(f"\n  Checking specific read IDs:")
        for rid in read_ids:
            if rid in reads_with_both:
                primaries = primary_df[primary_df["read_id"] == rid]
                supplementaries = supplementary_df[supplementary_df["read_id"] == rid]
                logger.debug(f"    {rid}:")
                logger.debug(f"      Primary: {len(primaries)} alignment(s)")
                for _, row in primaries.iterrows():
                    logger.debug(f"        {row['reference_id']}:{row['reference_start']}-{row['reference_end']} (col4={row.get('col4', 'N/A')}, is_supp={row.get('is_supplementary', 'N/A')})")
                logger.debug(f"      Supplementary: {len(supplementaries)} alignment(s)")
                for _, row in supplementaries.iterrows():
                    logger.debug(f"        {row['reference_id']}:{row['reference_start']}-{row['reference_end']} (col4={row.get('col4', 'N/A')}, is_supp={row.get('is_supplementary', 'N/A')})")
    
    # Create breakpoint pairs
    logger.debug("\n[STEP 6] Creating breakpoint pairs...")
    primary_filtered = primary_df[primary_df["read_id"].isin(reads_with_both)].copy()
    supplementary_filtered = supplementary_df[supplementary_df["read_id"].isin(reads_with_both)].copy()
    
    breakpoint_pairs = []
    primary_by_read = primary_filtered.groupby("read_id", observed=True)
    supplementary_by_read = supplementary_filtered.groupby("read_id", observed=True)
    
    for read_id in reads_with_both:
        read_primaries = primary_by_read.get_group(read_id) if read_id in primary_by_read.groups else pd.DataFrame()
        read_supplementaries = supplementary_by_read.get_group(read_id) if read_id in supplementary_by_read.groups else pd.DataFrame()
        
        if read_primaries.empty or read_supplementaries.empty:
            continue
        
        # Create pairs using vectorized Pandas operations (much faster than nested loops)
        pairs = _create_breakpoint_pairs_vectorized(read_primaries, read_supplementaries, read_id)
        breakpoint_pairs.extend(pairs)
    
    logger.debug(f"  Created {len(breakpoint_pairs)} breakpoint pairs from {len(reads_with_both)} reads")
    
    if not breakpoint_pairs:
        logger.debug(f"  ERROR: No breakpoint pairs created")
        return
    
    # Show example pairs
    if read_ids:
        logger.debug(f"\n  Breakpoint pairs for specified read IDs:")
        for pair in breakpoint_pairs[:20]:  # Show first 20
            if pair["read_id"] in read_ids:
                logger.debug(f"    {pair['read_id']}: primary={pair['primary_chrom']}:{pair['primary_start']}-{pair['primary_end']}, "
                      f"supp={pair['supp_chrom']}:{pair['supp_start']}-{pair['supp_end']}")
    
    # Cluster pairs
    logger.debug(f"\n[STEP 7] Clustering breakpoint pairs (distance={cluster_distance}bp)...")
    
    # Group by chromosome combination
    pairs_by_chrom = {}
    for i, pair in enumerate(breakpoint_pairs):
        chrom_key = (pair["primary_chrom"], pair["supp_chrom"])
        if chrom_key not in pairs_by_chrom:
            pairs_by_chrom[chrom_key] = []
        pairs_by_chrom[chrom_key].append((i, pair))
    
    logger.debug(f"  Grouped into {len(pairs_by_chrom)} chromosome combinations")
    for chrom_key, pairs in pairs_by_chrom.items():
        logger.debug(f"    {chrom_key[0]} -> {chrom_key[1]}: {len(pairs)} pairs")
    
    # Cluster within each chromosome combination
    clustered_pairs = []
    used_indices = set()
    
    for chrom_key, chrom_pairs in pairs_by_chrom.items():
        if len(chrom_pairs) == 0:
            continue
        
        # Use DBSCAN clustering for this chromosome combination
        # Convert to list format expected by DBSCAN function
        chrom_pair_list = [pair for _, pair in chrom_pairs]
        chrom_clustered = _cluster_breakpoint_pairs_dbscan(
            chrom_pair_list, cluster_distance=cluster_distance, min_read_support=min_read_support
        )
        clustered_pairs.extend(chrom_clustered)
    
    logger.debug(f"  Clustered into {len(clustered_pairs)} clusters")
    
    # Show cluster details
    if clustered_pairs:
        logger.debug(f"\n  Cluster details:")
        sorted_clusters = sorted(clustered_pairs, key=lambda p: p["read_count"], reverse=True)
        for i, cluster in enumerate(sorted_clusters[:10]):  # Show top 10
            logger.debug(f"    Cluster {i+1}: {cluster['read_count']} reads")
            logger.debug(f"      Primary: {cluster['primary_chrom']}:{cluster['primary_start']}-{cluster['primary_end']}")
            logger.debug(f"      Supplementary: {cluster['supp_chrom']}:{cluster['supp_start']}-{cluster['supp_end']}")
            if read_ids:
                cluster_read_ids_list = list(cluster['read_ids'])
                matching_reads = [rid for rid in read_ids if rid in cluster_read_ids_list]
                if matching_reads:
                    logger.debug(f"      Matching specified reads: {matching_reads}")
    
    # Filter by min_read_support
    logger.debug(f"\n[STEP 8] Filtering by min_read_support (>= {min_read_support})...")
    supported_pairs = [
        p for p in clustered_pairs 
        if p["read_count"] >= min_read_support
    ]
    
    logger.debug(f"  Clusters with >= {min_read_support} reads: {len(supported_pairs)} (out of {len(clustered_pairs)} total)")
    
    if supported_pairs:
        logger.debug(f"\n  SUCCESS: Found {len(supported_pairs)} supported breakpoint pairs")
        for i, pair in enumerate(supported_pairs):
            logger.debug(f"    {i+1}. {pair['read_count']} reads: "
                  f"primary={pair['primary_chrom']}:{pair['primary_start']}-{pair['primary_end']}, "
                  f"supp={pair['supp_chrom']}:{pair['supp_start']}-{pair['supp_end']}")
    else:
        logger.debug(f"  ERROR: No clusters meet the min_read_support threshold of {min_read_support}")
        if clustered_pairs:
            max_read_count = max(p["read_count"] for p in clustered_pairs)
            logger.debug(f"  Maximum read count in any cluster: {max_read_count}")
            logger.debug(f"  Top clusters by read count:")
            sorted_clusters = sorted(clustered_pairs, key=lambda p: p["read_count"], reverse=True)
            for i, cluster in enumerate(sorted_clusters[:5]):
                logger.debug(f"    {i+1}. {cluster['read_count']} reads: "
                      f"primary={cluster['primary_chrom']}:{cluster['primary_start']}-{cluster['primary_end']}, "
                      f"supp={cluster['supp_chrom']}:{cluster['supp_start']}-{cluster['supp_end']}")
    
    logger.debug("\n" + "="*80)
    logger.debug("DIAGNOSTIC COMPLETE")
    logger.debug("="*80 + "\n")


def _migrate_json_to_parquet(work_dir: str, sample_id: str, fusion_data: Dict[str, Any]) -> None:
    """
    Migrate fusion candidates from JSON lists to Parquet files.
    This is called when loading old JSON-based metadata.
    
    Args:
        work_dir: Working directory
        sample_id: Sample identifier
        fusion_data: Dictionary containing lists of candidate dicts
    """
    for candidate_type in ["target_candidates", "genome_wide_candidates", "master_bed_candidates"]:
        candidates_list = fusion_data.get(candidate_type, [])
        if candidates_list and isinstance(candidates_list, list) and len(candidates_list) > 0:
            try:
                # Convert list of dicts to DataFrame
                df = pd.DataFrame(candidates_list)
                if not df.empty:
                    # Save to Parquet
                    _save_fusion_candidates_parquet(df, candidate_type, work_dir, sample_id)
                    logger.info(f"Migrated {len(candidates_list)} {candidate_type} from JSON to Parquet")
            except Exception as e:
                logger.warning(f"Error migrating {candidate_type} to Parquet: {e}")


def _save_fusion_metadata(
    fusion_metadata: FusionMetadata, work_dir: str, sample_id: str
) -> None:
    """
    Save fusion metadata to disk for persistence between BAM file processing.
    Now uses Parquet files for candidate data storage (much faster than JSON).

    Args:
        fusion_metadata: FusionMetadata object to save
        work_dir: Working directory
        sample_id: Sample ID for the metadata file name
    """
    try:
        # Ensure sample directory exists
        sample_dir = os.path.join(work_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)

        # Save fusion candidates to Parquet files (fast storage)
        if fusion_metadata.fusion_data:
            for candidate_type in ["target_candidates", "genome_wide_candidates", "master_bed_candidates"]:
                candidates = fusion_metadata.fusion_data.get(candidate_type)
                
                # Handle both DataFrame and list of dicts (for backward compatibility)
                if candidates is not None:
                    if isinstance(candidates, pd.DataFrame):
                        # Already a DataFrame - save directly
                        _save_fusion_candidates_parquet(candidates, candidate_type, work_dir, sample_id)
                    elif isinstance(candidates, list) and len(candidates) > 0:
                        # List of dicts - convert to DataFrame first
                        df = pd.DataFrame(candidates)
                        _save_fusion_candidates_parquet(df, candidate_type, work_dir, sample_id)
                    # Empty lists are skipped (no file created)

        # Create metadata file path (for non-candidate metadata only)
        metadata_path = os.path.join(sample_dir, f"{sample_id}_fusion_metadata.json")

        # Convert dataclass to dictionary for JSON serialization
        metadata_dict = asdict(fusion_metadata)
        
        # Remove large candidate lists from JSON (they're now in Parquet files)
        # Keep empty lists or None to indicate structure
        if "fusion_data" in metadata_dict and metadata_dict["fusion_data"]:
            # Replace large lists with empty lists or file indicators
            for key in ["target_candidates", "genome_wide_candidates", "master_bed_candidates"]:
                if key in metadata_dict["fusion_data"]:
                    candidates = metadata_dict["fusion_data"][key]
                    # Only store in JSON if it's a small list (for backward compatibility)
                    # Large lists are stored in Parquet
                    if isinstance(candidates, list) and len(candidates) > 1000:
                        # Store empty list in JSON, data is in Parquet
                        metadata_dict["fusion_data"][key] = []
                    elif isinstance(candidates, pd.DataFrame):
                        # DataFrame - don't store in JSON
                        metadata_dict["fusion_data"][key] = []

        # Save to JSON file (metadata only, candidates are in Parquet)
        with open(metadata_path, "w") as f:
            json.dump(metadata_dict, f, indent=2, default=str)  # default=str handles any non-serializable types

        logger.debug(f"Saved fusion metadata to {metadata_path}")

    except Exception as e:
        logger.error(f"Error saving fusion metadata: {str(e)}")
        # Don't raise - just log the error so processing can continue


def _migrate_old_fusion_metadata(metadata_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Migrate old format fusion metadata to new format for backward compatibility.
    
    This function handles data generated from the main branch and converts it
    to be compatible with the fusion_bug branch format.
    
    Args:
        metadata_dict: Dictionary loaded from JSON file (may be old format)
        
    Returns:
        Updated dictionary compatible with current FusionMetadata structure
    """
    migrated = metadata_dict.copy()
    
    # Ensure all required fields exist with defaults
    required_fields = {
        "sample_id": "unknown",
        "file_path": "",
        "analysis_timestamp": time.time(),
        "target_fusion_path": None,
        "genome_wide_fusion_path": None,
        "analysis_results": {},
        "processing_steps": [],
        "error_message": None,
        "fusion_data": {},
        "target_panel": None,
    }
    
    # Fill in missing fields
    for field, default_value in required_fields.items():
        if field not in migrated:
            migrated[field] = default_value
    
    # Ensure processing_steps is a list
    if not isinstance(migrated.get("processing_steps"), list):
        migrated["processing_steps"] = []
    
    # Ensure analysis_results is a dict
    if not isinstance(migrated.get("analysis_results"), dict):
        migrated["analysis_results"] = {}
    
    # Ensure fusion_data is a dict with required structure
    if not isinstance(migrated.get("fusion_data"), dict):
        migrated["fusion_data"] = {}
    
    fusion_data = migrated["fusion_data"]
    
    # Ensure required fusion_data keys exist
    if "target_candidates" not in fusion_data:
        fusion_data["target_candidates"] = []
    if "genome_wide_candidates" not in fusion_data:
        fusion_data["genome_wide_candidates"] = []
    
    # Migrate old CSV-based data to fusion_data structure if needed
    # Check if we have CSV files but no fusion_data entries
    # Try to determine sample directory from various paths
    sample_dir = None
    if metadata_dict.get("target_fusion_path"):
        sample_dir = os.path.dirname(metadata_dict.get("target_fusion_path"))
    elif metadata_dict.get("genome_wide_fusion_path"):
        sample_dir = os.path.dirname(metadata_dict.get("genome_wide_fusion_path"))
    elif metadata_dict.get("file_path"):
        # Try to infer from file_path
        file_path = metadata_dict.get("file_path")
        if file_path:
            # Assume work_dir structure: work_dir/sample_id/file.bam
            potential_dir = os.path.dirname(file_path)
            if os.path.basename(potential_dir):  # If there's a parent directory
                sample_dir = potential_dir
    
    if sample_dir and os.path.exists(sample_dir):
        target_csv = os.path.join(sample_dir, "target_fusion.csv")
        genome_csv = os.path.join(sample_dir, "genome_wide_fusion.csv")
        
        # If fusion_data is empty but CSV files exist, try to load from CSVs
        if (not fusion_data.get("target_candidates") and 
            not fusion_data.get("genome_wide_candidates") and
            (os.path.exists(target_csv) or os.path.exists(genome_csv))):
            
            logger.info(f"Migrating fusion data from CSV files for {migrated.get('sample_id', 'unknown')}")
            
            try:
                # Load target fusion CSV if it exists
                if os.path.exists(target_csv):
                    try:
                        target_df = pd.read_csv(target_csv)
                        if not target_df.empty:
                            fusion_data["target_candidates"] = target_df.to_dict("records")
                            logger.info(f"Migrated {len(fusion_data['target_candidates'])} target candidates from CSV")
                    except Exception as e:
                        logger.warning(f"Could not migrate target fusion CSV: {e}")
                
                # Load genome-wide fusion CSV if it exists
                if os.path.exists(genome_csv):
                    try:
                        genome_df = pd.read_csv(genome_csv)
                        if not genome_df.empty:
                            fusion_data["genome_wide_candidates"] = genome_df.to_dict("records")
                            logger.info(f"Migrated {len(fusion_data['genome_wide_candidates'])} genome-wide candidates from CSV")
                    except Exception as e:
                        logger.warning(f"Could not migrate genome-wide fusion CSV: {e}")
                
                # Update analysis results counts
                if fusion_data.get("target_candidates") or fusion_data.get("genome_wide_candidates"):
                    migrated["analysis_results"]["target_candidates_count"] = len(fusion_data.get("target_candidates", []))
                    migrated["analysis_results"]["genome_wide_candidates_count"] = len(fusion_data.get("genome_wide_candidates", []))
                    migrated["processing_steps"].append("migrated_from_csv")
                    
            except Exception as e:
                logger.warning(f"Error during CSV migration: {e}")
    
    return migrated


def _create_metadata_from_csv_files(work_dir: str, sample_id: str) -> Optional[FusionMetadata]:
    """
    Create FusionMetadata from CSV files when metadata JSON doesn't exist.
    This provides backward compatibility for data generated from main branch.
    
    Args:
        work_dir: Working directory
        sample_id: Sample ID
        
    Returns:
        FusionMetadata object if CSV files found, None otherwise
    """
    try:
        sample_dir = os.path.join(work_dir, sample_id)
        if not os.path.exists(sample_dir):
            return None
        
        target_csv = os.path.join(sample_dir, "target_fusion.csv")
        genome_csv = os.path.join(sample_dir, "genome_wide_fusion.csv")
        
        # Check if CSV files exist
        has_target = os.path.exists(target_csv)
        has_genome = os.path.exists(genome_csv)
        
        if not (has_target or has_genome):
            return None
        
        logger.info(f"Creating fusion metadata from CSV files for {sample_id}")
        
        # Create new metadata object
        fusion_metadata = FusionMetadata(
            sample_id=sample_id,
            file_path="",  # Unknown for old data
            analysis_timestamp=time.time(),
            target_fusion_path=target_csv if has_target else None,
            genome_wide_fusion_path=genome_csv if has_genome else None,
            target_panel=None,  # Unknown for old data
        )
        
        # Load CSV data into fusion_data
        fusion_data = {}
        
        if has_target:
            try:
                target_df = pd.read_csv(target_csv)
                if not target_df.empty:
                    fusion_data["target_candidates"] = target_df.to_dict("records")
                    logger.info(f"Loaded {len(fusion_data['target_candidates'])} target candidates from CSV")
            except Exception as e:
                logger.warning(f"Could not load target fusion CSV: {e}")
                fusion_data["target_candidates"] = []
        else:
            fusion_data["target_candidates"] = []
        
        if has_genome:
            try:
                genome_df = pd.read_csv(genome_csv)
                if not genome_df.empty:
                    fusion_data["genome_wide_candidates"] = genome_df.to_dict("records")
                    logger.info(f"Loaded {len(fusion_data['genome_wide_candidates'])} genome-wide candidates from CSV")
            except Exception as e:
                logger.warning(f"Could not load genome-wide fusion CSV: {e}")
                fusion_data["genome_wide_candidates"] = []
        else:
            fusion_data["genome_wide_candidates"] = []
        
        fusion_metadata.fusion_data = fusion_data
        fusion_metadata.analysis_results = {
            "target_candidates_count": len(fusion_data.get("target_candidates", [])),
            "genome_wide_candidates_count": len(fusion_data.get("genome_wide_candidates", [])),
        }
        fusion_metadata.processing_steps = ["created_from_csv_files"]
        
        # Save the created metadata for future use
        try:
            _save_fusion_metadata(fusion_metadata, work_dir, sample_id)
            logger.info(f"Saved created fusion metadata for {sample_id}")
        except Exception as e:
            logger.warning(f"Could not save created metadata: {e}")
        
        return fusion_metadata
        
    except Exception as e:
        logger.warning(f"Error creating metadata from CSV files: {e}")
        return None


def _load_fusion_metadata(work_dir: str, sample_id: str) -> Optional[FusionMetadata]:
    """
    Load existing fusion metadata from disk with backward compatibility.
    Now loads candidate data from Parquet files (much faster than JSON).

    Args:
        work_dir: Working directory
        sample_id: Sample ID for the metadata file name

    Returns:
        FusionMetadata object if found, None otherwise
    """
    try:
        # Create metadata file path
        sample_dir = os.path.join(work_dir, sample_id)
        metadata_path = os.path.join(sample_dir, f"{sample_id}_fusion_metadata.json")

        if not os.path.exists(metadata_path):
            # Try to create metadata from CSV files if they exist (backward compatibility)
            return _create_metadata_from_csv_files(work_dir, sample_id)

        # Load from JSON file
        with open(metadata_path, "r") as f:
            metadata_dict = json.load(f)

        # Migrate old format to new format if needed
        metadata_dict = _migrate_old_fusion_metadata(metadata_dict)

        # Check if we have Parquet files (new format) or JSON lists (old format)
        # Try loading from Parquet first (fast path)
        fusion_data = metadata_dict.get("fusion_data", {})
        parquet_loaded = False
        
        for candidate_type in ["target_candidates", "genome_wide_candidates", "master_bed_candidates"]:
            parquet_df = _load_fusion_candidates_parquet(candidate_type, work_dir, sample_id)
            if parquet_df is not None and not parquet_df.empty:
                # Parquet file exists - use it (convert to list for compatibility)
                fusion_data[candidate_type] = parquet_df.to_dict("records")
                parquet_loaded = True
            elif candidate_type in fusion_data:
                # Check if JSON has data (old format)
                candidates_list = fusion_data.get(candidate_type, [])
                if isinstance(candidates_list, list) and len(candidates_list) > 0:
                    # Migrate from JSON to Parquet
                    _migrate_json_to_parquet(work_dir, sample_id, fusion_data)
                    parquet_loaded = True
                    # Keep the list for now (will be replaced by Parquet on next save)
        
        # Update metadata_dict with loaded/migrated fusion_data
        metadata_dict["fusion_data"] = fusion_data

        # Convert dictionary back to FusionMetadata object
        try:
            fusion_metadata = FusionMetadata(**metadata_dict)
        except TypeError as e:
            # If there are extra fields in the old format, filter them out
            logger.warning(f"Extra fields in metadata, filtering: {e}")
            # Get only the fields that FusionMetadata expects
            from dataclasses import fields
            valid_fields = {f.name for f in fields(FusionMetadata)}
            filtered_dict = {k: v for k, v in metadata_dict.items() if k in valid_fields}
            fusion_metadata = FusionMetadata(**filtered_dict)

        # Ensure the loaded metadata has the correct structure
        if (
            not hasattr(fusion_metadata, "fusion_data")
            or fusion_metadata.fusion_data is None
        ):
            fusion_metadata.fusion_data = {}

        # Ensure required keys exist (as empty lists if not loaded)
        if "target_candidates" not in fusion_metadata.fusion_data:
            fusion_metadata.fusion_data["target_candidates"] = []
        if "genome_wide_candidates" not in fusion_metadata.fusion_data:
            fusion_metadata.fusion_data["genome_wide_candidates"] = []
        if "master_bed_candidates" not in fusion_metadata.fusion_data:
            fusion_metadata.fusion_data["master_bed_candidates"] = []

        # Save migrated metadata back to disk if migration occurred
        if "migrated_from_csv" in fusion_metadata.processing_steps or parquet_loaded:
            try:
                _save_fusion_metadata(fusion_metadata, work_dir, sample_id)
                logger.info(f"Saved migrated/updated fusion metadata for {sample_id}")
            except Exception as e:
                logger.warning(f"Could not save migrated metadata: {e}")

        logger.debug(f"Loaded fusion metadata from {metadata_path}")
        logger.debug(
            f"Existing data: {len(fusion_metadata.fusion_data.get('target_candidates', []))} target, {len(fusion_metadata.fusion_data.get('genome_wide_candidates', []))} genome-wide, {len(fusion_metadata.fusion_data.get('master_bed_candidates', []))} master BED"
        )

        return fusion_metadata

    except Exception as e:
        logger.warning(f"Error loading fusion metadata: {str(e)}")
        import traceback
        logger.debug(f"Traceback: {traceback.format_exc()}")
        return None


def _merge_fusion_metadata_objects(
    existing: FusionMetadata, new: FusionMetadata
) -> FusionMetadata:
    """
    Merge two FusionMetadata objects, combining their fusion data.

    Args:
        existing: Existing FusionMetadata from disk
        new: New FusionMetadata from current processing

    Returns:
        Merged FusionMetadata object
    """
    # Start with the new metadata object
    merged = new

    # Merge processing steps (remove duplicates)
    existing_steps = existing.processing_steps or []
    new_steps = new.processing_steps or []
    merged.processing_steps = list(set(existing_steps + new_steps))

    # Merge fusion data
    existing_fusion_data = existing.fusion_data or {}
    new_fusion_data = new.fusion_data or {}

    # Combine target candidates
    existing_target = existing_fusion_data.get("target_candidates", [])
    new_target = new_fusion_data.get("target_candidates", [])
    merged_target = existing_target + new_target

    # Combine genome-wide candidates
    existing_genome = existing_fusion_data.get("genome_wide_candidates", [])
    new_genome = new_fusion_data.get("genome_wide_candidates", [])
    merged_genome = existing_genome + new_genome

    # Update merged metadata - ensure complete structure
    merged.fusion_data = {
        "target_candidates": merged_target,
        "genome_wide_candidates": merged_genome,
    }

    # Ensure all required keys exist (defensive programming)
    if "target_candidates" not in merged.fusion_data:
        merged.fusion_data["target_candidates"] = []
    if "genome_wide_candidates" not in merged.fusion_data:
        merged.fusion_data["genome_wide_candidates"] = []

    # Update analysis results counts
    if merged.analysis_results:
        merged.analysis_results["target_candidates_count"] = len(merged_target)
        merged.analysis_results["genome_wide_candidates_count"] = len(merged_genome)

    logger.info(
        f"Merged fusion metadata: {len(existing_target)} + {len(new_target)} = {len(merged_target)} target candidates"
    )
    logger.info(
        f"Merged fusion metadata: {len(existing_genome)} + {len(new_genome)} = {len(merged_genome)} genome-wide candidates"
    )

    return merged


# =============================================================================
# FUSION DATA PROCESSING FUNCTIONS
# =============================================================================


def _annotate_results(result: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
    """Annotates the result DataFrame with tags and colors with memory optimization.
    
    Filters fusion candidates to require at least 3 supporting reads per gene pair
    to ensure reliable fusion detection and reduce false positives.
    """
    # Use inplace operations to avoid unnecessary copies
    result = result.copy()  # Only one copy needed

    # Group by read_id and aggregate col4 (Gene) values efficiently
    # Use sorted() to ensure deterministic tag generation regardless of input order
    lookup = result.groupby("read_id", observed=True)["col4"].agg(
        lambda x: ",".join(sorted(set(x)))
    )
    result["tag"] = result["read_id"].map(lookup)

    # Generate colors for each read_id group efficiently
    colors = result.groupby("read_id", observed=True)["col4"].apply(
        lambda x: _generate_random_color()
    )
    result["Color"] = result["read_id"].map(colors)

    # Clean string data efficiently
    result = result.apply(lambda x: x.strip() if isinstance(x, str) else x)

    # Find good pairs (gene pairs supported by minimum threshold for reliable fusion detection)
    min_support = get_fusion_threshold("read_support")
    goodpairs = result.groupby("tag", observed=True)["read_id"].transform("nunique") >= min_support

    return result, goodpairs


def get_gene_network(gene_pairs):
    """
    Build gene network from gene pairs with memory optimization.

    Args:
        gene_pairs: List of gene pairs, where each pair should be a tuple of 2 gene names

    Returns:
        List of connected components (gene groups)
    """
    if not gene_pairs:
        return []

    G = nx.Graph()

    for pair in gene_pairs:
        # Validate that pair is a tuple/list with exactly 2 elements
        if isinstance(pair, (tuple, list)) and len(pair) == 2:
            # Ensure both elements are strings and not empty
            if (
                pair[0]
                and pair[1]
                and isinstance(pair[0], str)
                and isinstance(pair[1], str)
            ):
                G.add_edge(pair[0], pair[1])

    if not G.nodes():
        return []

    connected_components = list(nx.connected_components(G))

    # Free the graph immediately
    del G

    return [list(component) for component in connected_components]


def collapse_ranges(df, max_distance):
    """Collapse ranges within a fixed distance with memory optimization."""
    collapsed = []
    current_range = None

    for _, row in df.iterrows():
        # Only unpack what we need
        start, end = row["start"], row["end"]

        if current_range is None:
            current_range = row
        else:
            if start <= current_range["end"] + max_distance:
                current_range["end"] = max(current_range["end"], end)
            else:
                collapsed.append(current_range)
                current_range = row

    if current_range is not None:
        collapsed.append(current_range)

    return pd.DataFrame(collapsed)


def _get_reads(reads: pd.DataFrame) -> pd.DataFrame:
    """Get reads for a specific gene with memory optimization."""
    df = reads.copy()  # Only one copy needed
    df.columns = [
        "chromosome",
        "start",
        "end",
        "gene",
        "chromosome2",
        "start2",
        "end2",
        "id",
        "quality",
        "strand",
        "read_start",
        "read_end",
        "secondary",
        "supplementary",
        "span",
        "tag",
        "color",
    ]

    try:
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
    except ValueError as e:
        logger.error(f"Error converting start/end to int: {str(e)}")
        problematic_start = df[df['start'].apply(lambda x: not str(x).isdigit())]['start'].tolist()
        problematic_end = df[df['end'].apply(lambda x: not str(x).isdigit())]['end'].tolist()
        logger.debug(f"Problematic values in start: {problematic_start}")
        logger.debug(f"Problematic values in end: {problematic_end}")
        raise

    # Sort the DataFrame by chromosome, start, and end positions
    df = df.sort_values(by=["chromosome", "start", "end"])

    df = df.drop_duplicates(subset=["start2", "end2", "id"])

    # Group by chromosome and collapse ranges within each group
    # Use a more explicit approach to avoid the FutureWarning
    result_list = []
    for chromosome, group in df.groupby("chromosome", observed=True):
        collapsed_group = collapse_ranges(group, 10000)
        if not collapsed_group.empty:
            collapsed_group["chromosome"] = chromosome
            result_list.append(collapsed_group)

    if result_list:
        result = pd.concat(result_list, ignore_index=True)
    else:
        result = pd.DataFrame()

    # Free the intermediate DataFrame
    del df

    return result


def _generate_fusion_summary_files(output_file: str, processed_data: dict) -> None:
    """Generate summary files for the summary component.
    
    Creates sv_count.txt and fusion_results.csv files that can be easily
    read by the summary component without needing to parse pickle files.
    """
    try:
        import os
        from pathlib import Path
        
        # Get the directory where the output file is located
        output_dir = Path(output_file).parent
        
        # Extract fusion counts
        candidate_count = processed_data.get("candidate_count", 0)
        
        # Determine if this is target panel or genome-wide based on filename
        is_target_panel = "master" in Path(output_file).name
        
        # Generate sv_count.txt file (simple count)
        sv_count_file = output_dir / "sv_count.txt"
        if is_target_panel:
            # For target panel, we still write to sv_count.txt but it represents target fusions
            with open(sv_count_file, "w") as f:
                f.write(str(candidate_count))
        else:
            # For genome-wide, write the genome-wide count
            with open(sv_count_file, "w") as f:
                f.write(str(candidate_count))
        
        # Generate fusion_results.csv file with detailed information
        fusion_results_file = output_dir / "fusion_results.csv"
        with open(fusion_results_file, "w", newline="") as f:
            import csv
            writer = csv.writer(f)
            
            # Write header
            if is_target_panel:
                writer.writerow(["target_fusions", "genome_fusions"])
                writer.writerow([candidate_count, 0])  # Target panel only has target fusions
            else:
                writer.writerow(["target_fusions", "genome_fusions"])
                writer.writerow([0, candidate_count])  # Genome-wide only has genome fusions
        
        # Generate a combined summary file that has both counts
        # This will be overwritten each time, with the final one containing the correct totals
        summary_file = output_dir / "fusion_summary.csv"
        
        # Try to read existing summary to get both counts
        target_count = 0
        genome_count = 0
        
        if summary_file.exists():
            try:
                with open(summary_file, "r", newline="") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        target_count = int(row.get("target_fusions", 0))
                        genome_count = int(row.get("genome_fusions", 0))
                        break
            except Exception:
                pass  # If we can't read it, we'll start fresh
        
        # Update the appropriate count - preserve existing counts from other processing
        if is_target_panel:
            target_count = candidate_count
        else:
            genome_count = candidate_count
        
        # Write the updated summary with both counts preserved
        with open(summary_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["target_fusions", "genome_fusions"])
            writer.writerow([target_count, genome_count])
        
        logger.info(f"[Fusion] Updated summary file: target={target_count}, genome={genome_count}")
        
    except Exception as e:
        logger.warning(f"Failed to generate fusion summary files: {e}")


def preprocess_fusion_data_standalone(
    fusion_data: pd.DataFrame, output_file: str
) -> None:
    """Standalone version of fusion data preprocessing for CPU-bound execution with memory optimization.
    
    Applies filtering to require at least 3 supporting reads per gene pair
    to ensure reliable fusion detection while maintaining quality control.
    """
    try:
        # Apply categorical data types for efficiency
        fusion_data = fusion_data.astype(
            {
                "read_id": "category",
                "col4": "category",  # Gene column
                "reference_id": "category",  # Chromosome column
                "strand": "category",
            }
        )

        # Annotate results (this is the heavy lifting)
        annotated_data, goodpairs = _annotate_results(fusion_data)

        # Free the original data immediately
        del fusion_data

        # Create processed data structure for FusionVis
        processed_data = {
            "annotated_data": annotated_data,
            "goodpairs": goodpairs,
            "gene_pairs": [],
            "gene_groups": [],
            "candidate_count": 0,
        }

        # Process gene pairs and groups if we have good pairs
        if not annotated_data.empty and goodpairs.any():
            gene_pairs = (
                annotated_data[goodpairs]
                .sort_values(by="reference_start")["tag"]
                .unique()
                .tolist()
            )
            
            # Filter out empty strings and create valid gene pairs
            valid_gene_pairs = []
            for pair in gene_pairs:
                if pair and isinstance(pair, str) and pair.strip():
                    # Split by comma and filter out empty strings
                    genes = [gene.strip() for gene in pair.split(",") if gene.strip()]
                    # Only create pairs if we have at least 2 genes
                    if len(genes) >= 2:
                        valid_gene_pairs.append(tuple(genes))

            gene_pairs = valid_gene_pairs

            gene_groups_test = get_gene_network(gene_pairs)
            gene_groups = []

            for gene_group in gene_groups_test:
                # Count unique reads for this gene group (simpler approach)
                gene_group_reads = annotated_data[goodpairs][
                    annotated_data[goodpairs]["col4"].isin(gene_group)
                ]
                unique_read_count = gene_group_reads["read_id"].nunique()
                min_support = get_fusion_threshold("read_support")
                if unique_read_count >= min_support:  # Require minimum supporting reads for reliable fusion detection
                    gene_groups.append(gene_group)

            processed_data.update(
                {
                    "gene_pairs": gene_pairs,
                    "gene_groups": gene_groups,
                    "candidate_count": len(gene_groups),
                }
            )

        # Save processed data as pickle for efficient loading
        # Use atomic writes to prevent truncation from read/write clashes
        import pickle
        
        # Write to temporary file first, then atomically rename to prevent truncation
        # This ensures readers never see a partially written file
        temp_file = output_file + ".tmp"
        try:
            with open(temp_file, "wb") as f:
                pickle.dump(processed_data, f, protocol=pickle.HIGHEST_PROTOCOL)
            # Atomic rename - this is the key to preventing truncation
            os.replace(temp_file, output_file)
        except Exception as e:
            # Clean up temp file on error
            if os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except:
                    pass
            raise
        
        # Generate summary files for the summary component
        _generate_fusion_summary_files(output_file, processed_data)
        
        # Free the processed data immediately
        del processed_data

    except Exception as e:
        logger.error(f"Error pre-processing fusion data: {str(e)}")
        import traceback
        logger.debug(f"Exception details: {traceback.format_exc()}")


# =============================================================================
# STRUCTURAL VARIANT PROCESSING FUNCTIONS
# =============================================================================


def preprocess_structural_variants_standalone(output_dir: str) -> None:
    """Standalone version of structural variant preprocessing for CPU-bound execution."""
    try:
        sv_links_file = os.path.join(output_dir, "structural_variant_links.csv")

        if os.path.exists(sv_links_file) and os.path.getsize(sv_links_file) > 0:
            # Read the structural variant links CSV using safe reader
            dtype_spec = {
                "QNAME": str,
                "RNAME.1": str,
                "RNAME.2": str,
                "coord_1": np.int64,
                "coord_2": np.int64,
                "genomic_gap": np.int64,
                "event": str,
            }
            sv_links_df = safe_read_csv(sv_links_file, dtype=dtype_spec)

            if not sv_links_df.empty:
                # Process the links data to create a summary for display
                sv_summary = get_summary(sv_links_df, min_support=2)

                if not sv_summary.empty:
                    # Convert the summary to the format expected by the UI
                    sv_df = pd.DataFrame(
                        {
                            "Event Type": sv_summary["predominant_event"],
                            "Primary Location": sv_summary.apply(
                                lambda row: f"{row['RNAME.1']}:{row['coord_1']:,}",
                                axis=1,
                            ),
                            "Partner Location": sv_summary.apply(
                                lambda row: f"{row['RNAME.2']}:{row['coord_2']:,}",
                                axis=1,
                            ),
                            "Size (bp)": sv_summary["median_genomic_gap"].apply(
                                lambda x: f"{x:,}" if x >= 0 else "N/A"
                            ),
                            "Strand": "Unknown",  # Not available in links data
                            "Full Location": sv_summary.apply(
                                lambda row: (
                                    f"{row['RNAME.1']}:{row['coord_1']:,}-{row['coord_2']:,}"
                                    if row["RNAME.1"] == row["RNAME.2"]
                                    else f"{row['RNAME.1']}:{row['coord_1']:,} ⟷ {row['RNAME.2']}:{row['coord_2']:,}"
                                ),
                                axis=1,
                            ),
                            "Support Count": sv_summary["support_count"],
                            "Supporting Reads": sv_summary["supporting_reads"].apply(
                                lambda x: ", ".join(x[:5])
                                + ("..." if len(x) > 5 else "")
                            ),
                        }
                    )

                    # Save processed structural variant data
                    sv_df.to_csv(
                        os.path.join(output_dir, "structural_variants_processed.csv"),
                        index=False,
                    )

                    # Save count for quick access
                    with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                        f.write(str(len(sv_df)))

                    logger.info(f"Pre-processed {len(sv_df)} structural variant events")
                else:
                    # Save empty count
                    with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                        f.write("0")
            else:
                # Save empty count
                with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                    f.write("0")

    except Exception as e:
        logger.error(f"Error pre-processing structural variants: {str(e)}")
        logger.error("Exception details:", exc_info=True)


def finalize_fusion_accumulation_for_sample(
    sample_id: str, work_dir: str, target_panel: str, reference: Optional[str] = None
) -> Dict[str, Any]:
    """
    Force final accumulation of any remaining staged fusion files for a sample.
    This will also generate the final master BED file once all files are processed.
    
    This should be called when:
    - All files for a sample have been processed
    - The workflow is completing
    - There are staged files that haven't been accumulated yet
    
    Args:
        sample_id: Sample identifier
        work_dir: Working directory containing sample data
        target_panel: Target panel type
        reference: Optional reference genome path
    
    Returns:
        Dictionary with accumulation results
    """
    try:
        logger.info(f"Finalizing fusion accumulation for sample {sample_id}")
        
        # Check if there are pending files
        pending_count = _get_pending_count(work_dir, sample_id)
        
        if pending_count == 0:
            logger.info(f"No pending fusion files for {sample_id} - checking if final master BED generation is needed")
            # Even if no pending files, we should generate master BED if it hasn't been generated yet
            # This handles the case where all files were accumulated but master BED wasn't generated
        else:
            logger.info(f"Found {pending_count} pending fusion files for {sample_id} - forcing accumulation")
        
        # Force accumulation of remaining files (batch_size=1 ensures accumulation runs)
        # force=True will trigger master BED generation in _generate_output_files
        result = accumulate_fusion_candidates(
            work_dir, sample_id, target_panel, force=True, batch_size=1, reference=reference
        )
        
        # If accumulation succeeded but master BED wasn't generated (e.g., no pending files),
        # generate it now as a final step
        if result.get("status") == "success" or pending_count == 0:
            try:
                from robin.analysis.master_bed_generator import generate_master_bed
                
                # Get analysis counter
                analysis_counter = _load_analysis_counter(sample_id, work_dir)
                
                # Generate master BED file (final generation after all files processed)
                logger.info(f"Generating final master BED file for sample {sample_id} (counter: {analysis_counter})")
                master_bed_path = generate_master_bed(
                    sample_id=sample_id,
                    work_dir=work_dir,
                    analysis_counter=analysis_counter,
                    target_panel=target_panel,
                    logger_instance=logger,
                    reference=reference,
                )
                if master_bed_path:
                    logger.info(f"Final master BED file generated: {master_bed_path}")
                    result["master_bed_generated"] = True
                    result["master_bed_path"] = master_bed_path
            except Exception as e:
                logger.warning(f"Could not generate final master BED file: {e}")
        
        logger.info(f"Final fusion accumulation complete for {sample_id}: {result}")
        return result
        
    except Exception as e:
        logger.error(f"Error during final fusion accumulation for {sample_id}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {"status": "error", "error": str(e), "sample_id": sample_id}
