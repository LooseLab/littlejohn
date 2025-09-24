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
from pathlib import Path
from collections import defaultdict
from typing import Dict, Any, Optional, List, Tuple, Set
from dataclasses import dataclass, asdict

# Third-party imports
import numpy as np
import pandas as pd
import pysam
import networkx as nx

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


# =============================================================================
# MODULE-LEVEL CACHES
# =============================================================================

# Module-level caches for gene regions (shared across all function calls)
_gene_regions_cache: Dict[str, Dict[str, List[GeneRegion]]] = {}
_all_gene_regions_cache: Dict[str, Dict[str, List[GeneRegion]]] = {}
_target_panel_cache: Dict[str, str] = {}

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


def _setup_file_paths(target_panel: str = "rCNS2") -> Tuple[str, str]:
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
                # Fallback to rCNS2 if custom panel not found
                target_panel = "rCNS2"
                rCNS2_path = os.path.join(resources_dir, "rCNS2_panel_name_uniq.bed")
                if os.path.exists(rCNS2_path):
                    gene_bed = rCNS2_path
                else:
                    gene_bed = "rCNS2_panel_name_uniq.bed"

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
            # Default to rCNS2
            target_panel = "rCNS2"
            gene_bed = "rCNS2_panel_name_uniq.bed"

        all_gene_bed = "all_genes2.bed"

    return gene_bed, all_gene_bed


def _load_bed_regions(bed_file: str) -> Dict[str, List[GeneRegion]]:
    """Load BED file regions into memory for efficient lookup."""
    regions = defaultdict(list)

    if not os.path.exists(bed_file):
        return dict(regions)

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
    except Exception:
        # If any other error, return empty dict
        pass

    return dict(regions)


def _ensure_gene_regions_loaded(target_panel: str = "rCNS2") -> None:
    """Ensure gene regions are loaded into cache for the given target panel."""
    logger.info(f"DEBUG: _ensure_gene_regions_loaded called with target_panel='{target_panel}'")
    logger.info(f"DEBUG: Current cache keys: {list(_gene_regions_cache.keys())}")
    
    if target_panel not in _gene_regions_cache:
        logger.info(f"DEBUG: Loading gene regions for target_panel='{target_panel}'")
        gene_bed, all_gene_bed = _setup_file_paths(target_panel)
        logger.info(f"DEBUG: Resolved gene_bed='{gene_bed}', all_gene_bed='{all_gene_bed}'")

        # Load target panel gene regions
        _gene_regions_cache[target_panel] = _load_bed_regions(gene_bed)
        logger.info(f"Loaded {len(_gene_regions_cache[target_panel])} chromosomes for target panel {target_panel} from {gene_bed}")
        
        # Debug: Log some details about loaded regions
        total_regions = sum(len(regions) for regions in _gene_regions_cache[target_panel].values())
        logger.info(f"DEBUG: Total gene regions loaded for {target_panel}: {total_regions}")

        # Load genome-wide gene regions (shared across all panels)
        if not _all_gene_regions_cache:
            _all_gene_regions_cache["shared"] = _load_bed_regions(all_gene_bed)
            logger.info(f"Loaded {len(_all_gene_regions_cache['shared'])} chromosomes for genome-wide genes from {all_gene_bed}")
    else:
        logger.info(f"DEBUG: Gene regions for target_panel='{target_panel}' already loaded in cache")


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
        
        # Compare all pairs of alignments for similarity
        for i in range(len(read_rows)):
            for j in range(i + 1, len(read_rows)):
                row1 = read_rows[i]
                row2 = read_rows[j]
                
                # Only check alignments on the same chromosome
                if row1['reference_id'] == row2['reference_id']:
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
) -> List[Dict]:
    """
    Find intersections between a read and gene regions.
    Only processes reads that have supplementary alignments (true fusion candidates).

    Args:
        read: Pysam aligned segment
        ref_name: Reference chromosome name
        ref_start: Reference start position
        ref_end: Reference end position
        gene_regions: List of gene regions on this chromosome

    Returns:
        List of intersection dictionaries with the exact column structure
    """
    read_rows = []

    # Only process reads that have supplementary alignments (SA tag)
    # This ensures we only look at reads that actually map to multiple locations
    if not read.has_tag("SA"):
        return read_rows

    for gene_region in gene_regions:
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


def process_bam_for_fusions_work(
    bamfile: str, target_panel: str = "rCNS2"
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Process BAM file to find fusion candidates with memory optimization.
    This method matches the original fusion code exactly.

    Args:
        bamfile: Path to BAM file
        target_panel: Target panel to use (rCNS2, AML, or PanCan)

    Returns:
        Tuple of (target_panel_candidates, genome_wide_candidates) DataFrames

    Memory optimizations:
    - Streaming processing
    - Efficient data structures
    - Immediate cleanup of intermediate data
    """
    try:
        logger.info(f"DEBUG: process_bam_for_fusions_work called with target_panel='{target_panel}'")
        
        # Ensure gene regions are loaded
        _ensure_gene_regions_loaded(target_panel)

        # Find reads with supplementary alignments
        reads_with_supp = find_reads_with_supplementary(bamfile)

        if not reads_with_supp:
            logger.info(f"DEBUG: No supplementary reads found in {bamfile}")
            return None, None

        logger.info(f"DEBUG: Found {len(reads_with_supp)} reads with supplementary alignments")

        # Process reads for target panel fusions
        logger.info(f"DEBUG: Processing target panel fusions using cache key '{target_panel}'")
        logger.info(f"DEBUG: Available cache keys: {list(_gene_regions_cache.keys())}")
        
        target_candidates = _process_reads_for_fusions(
            bamfile, reads_with_supp, _gene_regions_cache[target_panel]
        )
        logger.info(f"Target panel candidates found: {len(target_candidates) if target_candidates is not None else 0}")

        # Process reads for genome-wide fusions
        genome_wide_candidates = _process_reads_for_fusions(
            bamfile, reads_with_supp, _all_gene_regions_cache["shared"]
        )
        logger.info(f"Genome-wide candidates found: {len(genome_wide_candidates) if genome_wide_candidates is not None else 0}")

        # Clear the reads_with_supp set to free memory
        reads_with_supp.clear()
        del reads_with_supp

        return target_candidates, genome_wide_candidates

    except Exception as e:
        logger.error(f"Error processing BAM file for fusions: {str(e)}")
        logger.error("Exception details:", exc_info=True)
        raise


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
        target_candidates, genome_wide_candidates = process_bam_for_fusions_work(
            file_path, target_panel
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
            "target_candidates": target_candidates,
            "genome_wide_candidates": genome_wide_candidates,
        }

        # Store fusion data in metadata - APPEND to existing data, don't replace
        if (
            not hasattr(fusion_metadata, "fusion_data")
            or fusion_metadata.fusion_data is None
        ):
            fusion_metadata.fusion_data = {}

        # Ensure required keys exist
        if "target_candidates" not in fusion_metadata.fusion_data:
            fusion_metadata.fusion_data["target_candidates"] = []
        if "genome_wide_candidates" not in fusion_metadata.fusion_data:
            fusion_metadata.fusion_data["genome_wide_candidates"] = []

        # Append new target candidates to existing ones
        if target_candidates is not None:
            new_target_records = target_candidates.to_dict("records")
            fusion_metadata.fusion_data["target_candidates"].extend(new_target_records)

        # Append new genome-wide candidates to existing ones
        if genome_wide_candidates is not None:
            new_genome_records = genome_wide_candidates.to_dict("records")
            fusion_metadata.fusion_data["genome_wide_candidates"].extend(
                new_genome_records
            )

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
    # Use accumulated data from fusion_metadata instead of just current analysis_results
    if (
        hasattr(fusion_metadata, "fusion_data")
        and fusion_metadata.fusion_data
        and fusion_metadata.fusion_data.get("target_candidates")
    ):

        # Convert accumulated target candidates back to DataFrame
        target_candidates_data = fusion_metadata.fusion_data["target_candidates"]
        if target_candidates_data:
            target_candidates = pd.DataFrame(target_candidates_data)

            # Filter for fusion candidates using the corrected logic
            gene_counts = target_candidates.groupby("read_id", observed=True)["col4"].nunique()
            fusion_read_ids = gene_counts[gene_counts > 1].index
            
            if len(fusion_read_ids) > 0:
                result = target_candidates[target_candidates["read_id"].isin(fusion_read_ids)]
                
                # Apply minimum read support threshold (3 or more supporting reads per gene pair)
                if not result.empty:
                    # Create tag column by grouping genes per read_id (same logic as _annotate_results)
                    lookup = result.groupby("read_id", observed=True)["col4"].agg(
                        lambda x: ",".join(sorted(set(x)))
                    )
                    result["tag"] = result["read_id"].map(lookup)
                    
                    # Group by gene pair (tag) and count supporting reads
                    gene_pair_read_counts = result.groupby("tag", observed=True)["read_id"].nunique()
                    min_support = get_fusion_threshold("read_support")
                    valid_gene_pairs = gene_pair_read_counts[gene_pair_read_counts >= min_support].index
                    result = result[result["tag"].isin(valid_gene_pairs)]

                if not result.empty:
                    # Save the filtered fusion candidates to CSV
                    result.to_csv(
                        os.path.join(
                            work_dir, sample_id, "fusion_candidates_master.csv"
                        ),
                        index=False,
                    )

                    # Preprocess for visualization
                    preprocess_fusion_data_standalone(
                        result,
                        os.path.join(
                            work_dir,
                            sample_id,
                            "fusion_candidates_master_processed.csv",
                        ),
                    )

    # Use accumulated data from fusion_metadata for genome-wide candidates
    logger.info(f"Checking genome-wide candidates in fusion_metadata: has_fusion_data={hasattr(fusion_metadata, 'fusion_data')}, fusion_data_exists={fusion_metadata.fusion_data is not None if hasattr(fusion_metadata, 'fusion_data') else False}")
    if (
        hasattr(fusion_metadata, "fusion_data")
        and fusion_metadata.fusion_data
        and fusion_metadata.fusion_data.get("genome_wide_candidates")
    ):
        genome_wide_data = fusion_metadata.fusion_data["genome_wide_candidates"]
        logger.info(f"Found {len(genome_wide_data)} genome-wide candidate records")

        # Convert accumulated genome-wide candidates back to DataFrame
        if genome_wide_data:
            genome_wide_candidates = pd.DataFrame(genome_wide_data)

            # Filter for fusion candidates using the corrected logic
            gene_counts_all = genome_wide_candidates.groupby("read_id", observed=True)["col4"].nunique()
            fusion_read_ids_all = gene_counts_all[gene_counts_all > 1].index
            
            if len(fusion_read_ids_all) > 0:
                result_all = genome_wide_candidates[genome_wide_candidates["read_id"].isin(fusion_read_ids_all)]
                logger.info(f"Genome-wide fusion candidates after basic filtering: {len(result_all)} records")

                if not result_all.empty:
                    # Save the filtered genome-wide candidates to CSV
                    result_all.to_csv(
                        os.path.join(work_dir, sample_id, "fusion_candidates_all.csv"),
                        index=False,
                    )

                    # Use the same preprocessing pipeline as target candidates
                    # This ensures proper annotation with tags, colors, and gene group detection
                    preprocess_fusion_data_standalone(
                        result_all,
                        os.path.join(
                            work_dir, sample_id, "fusion_candidates_all_processed.csv"
                        ),
                    )

    # Create sv_count.txt file with content "0" if it doesn't exist
    sv_count_file = os.path.join(work_dir, sample_id, "sv_count.txt")
    if not os.path.exists(sv_count_file):
        with open(sv_count_file, "w") as f:
            f.write("0")

    return {
        "target_candidates_path": os.path.join(
            work_dir, sample_id, "fusion_candidates_master.csv"
        ),
        "genome_wide_candidates_path": os.path.join(
            work_dir, sample_id, "fusion_candidates_all.csv"
        ),
    }


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


def _save_fusion_metadata(
    fusion_metadata: FusionMetadata, work_dir: str, sample_id: str
) -> None:
    """
    Save fusion metadata to disk for persistence between BAM file processing.

    Args:
        fusion_metadata: FusionMetadata object to save
        work_dir: Working directory
        sample_id: Sample ID for the metadata file name
    """
    try:
        # Ensure sample directory exists
        sample_dir = os.path.join(work_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)

        # Create metadata file path
        metadata_path = os.path.join(sample_dir, f"{sample_id}_fusion_metadata.json")

        # Convert dataclass to dictionary for JSON serialization
        metadata_dict = asdict(fusion_metadata)

        # Save to JSON file
        with open(metadata_path, "w") as f:
            json.dump(metadata_dict, f, indent=2)

        logger.debug(f"Saved fusion metadata to {metadata_path}")

    except Exception as e:
        logger.error(f"Error saving fusion metadata: {str(e)}")
        # Don't raise - just log the error so processing can continue


def _load_fusion_metadata(work_dir: str, sample_id: str) -> Optional[FusionMetadata]:
    """
    Load existing fusion metadata from disk.

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
            return None

        # Load from JSON file
        with open(metadata_path, "r") as f:
            metadata_dict = json.load(f)

        # Convert dictionary back to FusionMetadata object
        fusion_metadata = FusionMetadata(**metadata_dict)

        # Ensure the loaded metadata has the correct structure
        if (
            not hasattr(fusion_metadata, "fusion_data")
            or fusion_metadata.fusion_data is None
        ):
            fusion_metadata.fusion_data = {}

        # Ensure required keys exist
        if "target_candidates" not in fusion_metadata.fusion_data:
            fusion_metadata.fusion_data["target_candidates"] = []
        if "genome_wide_candidates" not in fusion_metadata.fusion_data:
            fusion_metadata.fusion_data["genome_wide_candidates"] = []

        logger.debug(f"Loaded fusion metadata from {metadata_path}")
        logger.debug(
            f"Existing data: {len(fusion_metadata.fusion_data.get('target_candidates', []))} target, {len(fusion_metadata.fusion_data.get('genome_wide_candidates', []))} genome-wide"
        )

        return fusion_metadata

    except Exception as e:
        logger.warning(f"Error loading fusion metadata: {str(e)}")
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
    lookup = result.groupby("read_id", observed=True)["col4"].agg(
        lambda x: ",".join(set(x))
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
        print(f"Error converting start/end to int: {str(e)}")
        print(
            f"Problematic values in start: {df[df['start'].apply(lambda x: not str(x).isdigit())]['start'].tolist()}"
        )
        print(
            f"Problematic values in end: {df[df['end'].apply(lambda x: not str(x).isdigit())]['end'].tolist()}"
        )
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
        
        # Update the appropriate count
        if is_target_panel:
            target_count = candidate_count
        else:
            genome_count = candidate_count
        
        # Write the updated summary
        with open(summary_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["target_fusions", "genome_fusions"])
            writer.writerow([target_count, genome_count])
        
    except Exception as e:
        print(f"Warning: Failed to generate fusion summary files: {e}")


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

        # Save processed data as pickle for efficient loading (atomic write to prevent corruption)
        import pickle

        def _atomic_pickle_dump(obj, target_path: str) -> None:
            """Atomically write pickle file to prevent corruption during concurrent reads."""
            import tempfile
            import os
            tmp_path = target_path + ".tmp"
            try:
                with open(tmp_path, "wb") as tf:
                    pickle.dump(obj, tf)
                    tf.flush()
                    os.fsync(tf.fileno())  # Ensure data is written to disk
                os.replace(tmp_path, target_path)  # Atomic rename
            except Exception as e:
                # Clean up temp file on error
                try:
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)
                except:
                    pass
                raise e

        _atomic_pickle_dump(processed_data, output_file)
        
        # Generate summary files for the summary component
        _generate_fusion_summary_files(output_file, processed_data)
        
        # Free the processed data immediately
        del processed_data

    except Exception as e:
        print(f"Error pre-processing fusion data: {str(e)}")
        import traceback
        print(f"Exception details: {traceback.format_exc()}")


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
