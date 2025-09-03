#!/usr/bin/env python3
"""
Simplified Fusion Analysis Module for robin

This module provides a simplified, function-based approach to fusion analysis
that removes unnecessary class hierarchy while preserving important functionality.

Key improvements:
- Removes FusionAnalysis wrapper class
- Converts FusionProcessor to module-level functions with caching
- Maintains gene region caching for performance
- Simplifies the API while preserving all functionality
"""

import os
import tempfile
import logging
import time
import json
import random
from typing import Dict, Any, Optional, List, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict
import numpy as np
import pandas as pd
import pysam
import networkx as nx
from robin.logging_config import get_job_logger


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


@dataclass
class GeneRegion:
    """Represents a gene region with start, end, and name."""

    start: int
    end: int
    name: str

    def overlaps_with(
        self, other_start: int, other_end: int, min_overlap: int = 100
    ) -> bool:
        """Check if this region overlaps with another region."""
        overlap_start = max(self.start, other_start)
        overlap_end = min(self.end, other_end)
        return (
            overlap_end > overlap_start and (overlap_end - overlap_start) > min_overlap
        )


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


# Module-level caches for gene regions (shared across all function calls)
_gene_regions_cache: Dict[str, Dict[str, List[GeneRegion]]] = {}
_all_gene_regions_cache: Dict[str, Dict[str, List[GeneRegion]]] = {}
_target_panel_cache: Dict[str, str] = {}


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
            aml_path = os.path.join(resources_dir, "AML_panel_name_uniq.bed")
            if os.path.exists(aml_path):
                gene_bed = aml_path
            else:
                gene_bed = "AML_panel_name_uniq.bed"
        else:
            # Default to rCNS2
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
        else:
            # Default to rCNS2
            target_panel = "rCNS2"
            gene_bed = "rCNS2_panel_name_uniq.bed"

        all_gene_bed = "all_genes2.bed"

    # Log the file paths for debugging
    logging.info(f"Fusion analysis initialized with target_panel: {target_panel}")
    logging.info(f"Gene bed path: {gene_bed}")
    logging.info(f"All gene bed path: {all_gene_bed}")

    return gene_bed, all_gene_bed


def _load_bed_regions(bed_file: str) -> Dict[str, List[GeneRegion]]:
    """
    Load BED file regions into memory for efficient lookup.

    Args:
        bed_file: Path to BED file

    Returns:
        Dictionary mapping chromosome to list of GeneRegion objects
    """
    regions = defaultdict(list)

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
    if target_panel not in _gene_regions_cache:
        gene_bed, all_gene_bed = _setup_file_paths(target_panel)

        # Load target panel gene regions
        _gene_regions_cache[target_panel] = _load_bed_regions(gene_bed)

        # Load genome-wide gene regions (shared across all panels)
        if not _all_gene_regions_cache:
            _all_gene_regions_cache["shared"] = _load_bed_regions(all_gene_bed)

        logging.info(
            f"Gene regions loaded for {target_panel} - target cache: {len(_gene_regions_cache[target_panel])} regions"
        )
        logging.info(
            f"Gene regions loaded - genome-wide cache: {len(_all_gene_regions_cache['shared'])} regions"
        )


def has_supplementary_alignments(bamfile: str) -> bool:
    """
    Quickly check if a BAM file has supplementary alignments.

    Args:
        bamfile: Path to the BAM file

    Returns:
        True if supplementary alignment is found, else False
    """
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
    Find all reads that have supplementary alignments.

    Args:
        bamfile: Path to the BAM file

    Returns:
        Set of read names with supplementary alignments
    """
    reads_with_supp = set()

    try:
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            for read in bam:
                if read.is_unmapped:
                    continue
                # Check for supplementary alignments or SA tag (matching ROBIN implementation)
                if read.is_supplementary or read.has_tag("SA"):
                    reads_with_supp.add(read.query_name)
    except Exception:
        pass

    return reads_with_supp


def _find_gene_intersections(
    read: pysam.AlignedSegment,
    ref_name: str,
    ref_start: int,
    ref_end: int,
    gene_regions: List[GeneRegion],
) -> List[Dict]:
    """
    Find gene intersections for a read.

    Args:
        read: Pysam AlignedSegment object
        ref_name: Reference chromosome name
        ref_start: Reference start position
        ref_end: Reference end position
        gene_regions: List of gene regions for this chromosome

    Returns:
        List of dictionaries with fusion candidate data in the correct format
    """
    intersections = []

    for gene_region in gene_regions:
        if gene_region.overlaps_with(ref_start, ref_end, min_overlap=100):
            intersection = {
                # Gene information columns (col1-col4) - following reference format
                "col1": ref_name,  # Chromosome
                "col2": gene_region.start,  # Gene start position
                "col3": gene_region.end,  # Gene end position
                "col4": gene_region.name,  # Gene name
                # Reference alignment information
                "reference_id": ref_name,
                "reference_start": ref_start,
                "reference_end": ref_end,
                # Read information
                "read_id": read.query_name,
                "mapping_quality": read.mapping_quality,
                "strand": "-" if read.is_reverse else "+",
                "read_start": read.query_alignment_start,
                "read_end": read.query_alignment_end,
                # Alignment flags
                "is_secondary": read.is_secondary,
                "is_supplementary": read.is_supplementary,
                # Calculated values
                "mapping_span": ref_end - ref_start,
            }
            intersections.append(intersection)

    return intersections


def _process_reads_for_fusions(
    bamfile: str,
    reads_with_supp: Set[str],
    gene_regions: Dict[str, List[GeneRegion]],
) -> Optional[pd.DataFrame]:
    """
    Process reads to find gene intersections and create fusion candidates.

    Args:
        bamfile: Path to BAM file
        reads_with_supp: Set of read names with supplementary alignments
        gene_regions: Dictionary of gene regions by chromosome

    Returns:
        DataFrame with fusion candidates or None if no candidates found
    """
    rows = []

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
                    rows.extend(read_rows)

    except Exception as e:
        logging.error(f"Error processing reads for fusions: {str(e)}")
        return None

    if not rows:
        return None

    # Create DataFrame
    df = pd.DataFrame(rows)

    # Apply quality filters (matching ROBIN implementation)
    df = df[(df["mapping_quality"] > 50) & (df["mapping_span"] > 200)].reset_index(
        drop=True
    )

    if df.empty:
        return None

    # Apply fusion candidate filtering (matching ROBIN implementation)
    # Filter for duplicated reads (reads mapping to multiple locations)
    duplicated_mask = df["read_id"].duplicated(keep=False)
    if not duplicated_mask.any():
        return None

    doubles = df[duplicated_mask]

    # Count unique genes per read
    gene_counts = doubles.groupby("read_id")["col4"].transform("nunique")

    # Keep reads mapping to multiple genes
    result = doubles[gene_counts > 1]

    return result if not result.empty else None


def _optimize_fusion_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Optimize DataFrame memory usage.

    Args:
        df: Input DataFrame

    Returns:
        Optimized DataFrame
    """
    if df.empty:
        return df

    # Optimize dtypes
    for col in df.columns:
        if df[col].dtype == "object":
            if col in ["read_id", "gene_name", "chromosome", "strand"]:
                df[col] = df[col].astype("category")

    return df


def _annotate_results(result: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Annotates the result DataFrame with tags and colors with memory optimization.

    Memory optimizations:
    - Use inplace operations where possible
    - Combine groupby operations
    - Use categorical dtypes
    - Vectorize string operations
    - Minimize DataFrame copies
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

    # Find good pairs (reads that map to more than 2 genes)
    goodpairs = result.groupby("tag", observed=True)["read_id"].transform("nunique") > 2

    return result, goodpairs


def _generate_random_color() -> str:
    """
    Generates a random color for use in plotting.

    Returns:
        str: A random hex color code.
    """
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))


def get_gene_network(gene_pairs):
    """
    Build gene network from gene pairs with memory optimization.

    Args:
        gene_pairs: List of gene pairs

    Returns:
        List of connected components

    Memory optimizations:
    - Efficient graph construction
    - Minimal intermediate data storage
    """
    G = nx.Graph()
    for pair in gene_pairs:
        G.add_edge(pair[0], pair[1])
    connected_components = list(nx.connected_components(G))

    # Free the graph immediately
    del G

    return [list(component) for component in connected_components]


def collapse_ranges(df, max_distance):
    """
    Collapse ranges within a fixed distance with memory optimization.

    Args:
        df: DataFrame with ranges
        max_distance: Maximum distance for collapsing

    Returns:
        DataFrame with collapsed ranges

    Memory optimizations:
    - Efficient iteration
    - Minimal intermediate data storage
    """
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
    """
    Get reads for a specific gene with memory optimization.

    Args:
        reads (pd.DataFrame): DataFrame with reads

    Returns:
        pd.DataFrame: DataFrame with reads

    Memory optimizations:
    - Efficient column operations
    - Minimal data copying
    - Immediate cleanup of intermediate data
    """
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

    # Add logging for start and end columns before conversion
    logging.debug(
        f"Converting start column to int. First few values: {df['start'].head()}"
    )
    logging.debug(f"Converting end column to int. First few values: {df['end'].head()}")

    try:
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
    except ValueError as e:
        logging.error(f"Error converting start/end to int: {str(e)}")
        logging.error(
            f"Problematic values in start: {df[df['start'].apply(lambda x: not str(x).isdigit())]['start'].tolist()}"
        )
        logging.error(
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


def preprocess_fusion_data_standalone(
    fusion_data: pd.DataFrame, output_file: str
) -> None:
    """
    Standalone version of fusion data preprocessing for CPU-bound execution with memory optimization.

    Args:
        fusion_data: Raw fusion candidate data
        output_file: Path to save processed data

    Memory optimizations:
    - Efficient data type optimization
    - Minimal intermediate data storage
    - Immediate cleanup of large DataFrames
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
            gene_pairs = [tuple(pair.split(",")) for pair in gene_pairs]
            gene_groups_test = get_gene_network(gene_pairs)
            gene_groups = []

            for gene_group in gene_groups_test:
                reads = _get_reads(
                    annotated_data[goodpairs][
                        annotated_data[goodpairs]["col4"].isin(gene_group)
                    ]
                )
                if len(reads) > 1:
                    gene_groups.append(gene_group)

            processed_data.update(
                {
                    "gene_pairs": gene_pairs,
                    "gene_groups": gene_groups,
                    "candidate_count": len(gene_groups),
                }
            )

        # Save processed data as pickle for efficient loading
        import pickle

        with open(output_file, "wb") as f:
            pickle.dump(processed_data, f)

        logging.info(f"Pre-processed fusion data saved to {output_file}")

        # Free the processed data immediately
        del processed_data

    except Exception as e:
        logging.error(f"Error pre-processing fusion data: {str(e)}")
        logging.error("Exception details:", exc_info=True)
    finally:
        # Force garbage collection
        import gc

        gc.collect()


def process_bam_for_fusions(
    bamfile: str, target_panel: str = "rCNS2"
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Process BAM file to find fusion candidates.

    Args:
        bamfile: Path to BAM file
        target_panel: Target panel to use for analysis

    Returns:
        Tuple of (target_panel_candidates, genome_wide_candidates) DataFrames
    """
    # Ensure gene regions are loaded
    _ensure_gene_regions_loaded(target_panel)

    try:
        # Find reads with supplementary alignments
        reads_with_supp = find_reads_with_supplementary(bamfile)

        if not reads_with_supp:
            return None, None

        # Process reads for target panel fusions
        target_candidates = _process_reads_for_fusions(
            bamfile, reads_with_supp, _gene_regions_cache[target_panel]
        )

        # Process reads for genome-wide fusions
        genome_wide_candidates = _process_reads_for_fusions(
            bamfile, reads_with_supp, _all_gene_regions_cache["shared"]
        )

        return target_candidates, genome_wide_candidates

    except Exception as e:
        import traceback

        logging.error(f"Error processing BAM file for fusions: {str(e)}")
        logging.error(f"Traceback: {traceback.format_exc()}")
        return None, None


def _check_and_create_folder(base_dir: str, sample_id: str) -> str:
    """Check and create folder for sample output."""
    sample_dir = os.path.join(base_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)
    return sample_dir


def _generate_output_files(
    sample_id: str,
    analysis_results: Dict[str, Any],
    fusion_metadata: FusionMetadata,
    work_dir: str,
) -> Dict[str, str]:
    """Generate output files for fusion analysis."""
    output_files = {}

    # Process target panel candidates
    if analysis_results["target_candidates"] is not None:
        target_candidates = analysis_results["target_candidates"]

        if not target_candidates.empty:
            # Filter for fusion candidates using the original logic
            uniques = target_candidates["read_id"].duplicated(keep=False)
            if uniques.any():
                doubles = target_candidates[uniques]
                counts = doubles.groupby("read_id", observed=True)["col4"].transform(
                    "nunique"
                )
                result = doubles[counts > 1]

                if not result.empty:
                    # Save raw fusion candidates
                    target_output = os.path.join(
                        work_dir, sample_id, "fusion_candidates_master.csv"
                    )
                    result.to_csv(target_output, index=False)
                    output_files["fusion_candidates_master"] = target_output

                    # Generate processed data for visualization
                    processed_output = os.path.join(
                        work_dir, sample_id, "fusion_candidates_master_processed.csv"
                    )
                    preprocess_fusion_data_standalone(result, processed_output)
                    output_files["fusion_candidates_master_processed"] = (
                        processed_output
                    )

    # Process genome-wide candidates
    if analysis_results["genome_wide_candidates"] is not None:
        genome_wide_candidates = analysis_results["genome_wide_candidates"]

        if not genome_wide_candidates.empty:
            # Filter for fusion candidates using the original logic
            uniques_all = genome_wide_candidates["read_id"].duplicated(keep=False)
            if uniques_all.any():
                doubles_all = genome_wide_candidates[uniques_all]
                counts_all = doubles_all.groupby("read_id", observed=True)[
                    "col4"
                ].transform("nunique")
                result_all = doubles_all[counts_all > 1]

                if not result_all.empty:
                    # Save raw fusion candidates
                    all_output = os.path.join(
                        work_dir, sample_id, "fusion_candidates_all.csv"
                    )
                    result_all.to_csv(all_output, index=False)
                    output_files["fusion_candidates_all"] = all_output

                    # Generate processed data for visualization
                    processed_all_output = os.path.join(
                        work_dir, sample_id, "fusion_candidates_all_processed.csv"
                    )
                    preprocess_fusion_data_standalone(result_all, processed_all_output)
                    output_files["fusion_candidates_all_processed"] = (
                        processed_all_output
                    )

    return output_files


def _load_existing_fusion_results(
    sample_dir: str, sample_id: str
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], Optional[Dict]]:
    """
    Load existing fusion results if they exist.

    Args:
        sample_dir: Sample-specific directory
        sample_id: Sample ID

    Returns:
        Tuple of (existing_target_candidates, existing_genome_wide_candidates, existing_metadata)
    """
    existing_target = None
    existing_genome_wide = None
    existing_metadata = None

    # Check for existing target fusion candidates
    target_path = os.path.join(sample_dir, "fusion_candidates_master.csv")
    if os.path.exists(target_path):
        try:
            existing_target = pd.read_csv(target_path)
            logging.info(
                f"Loaded existing target fusion candidates: {len(existing_target)} entries"
            )
        except Exception as e:
            logging.warning(f"Error loading existing target fusion candidates: {e}")

    # Check for existing genome-wide fusion candidates
    genome_path = os.path.join(sample_dir, "fusion_candidates_all.csv")
    if os.path.exists(genome_path):
        try:
            existing_genome_wide = pd.read_csv(genome_path)
            logging.info(
                f"Loaded existing genome-wide fusion candidates: {len(existing_genome_wide)} entries"
            )
        except Exception as e:
            logging.warning(
                f"Error loading existing genome-wide fusion candidates: {e}"
            )

    # Check for existing metadata
    metadata_path = os.path.join(sample_dir, f"{sample_id}_fusion_metadata.json")
    if os.path.exists(metadata_path):
        try:
            with open(metadata_path, "r") as f:
                existing_metadata = json.load(f)
            logging.info("Loaded existing fusion metadata")
        except Exception as e:
            logging.warning(f"Error loading existing fusion metadata: {e}")

    return existing_target, existing_genome_wide, existing_metadata


def _merge_fusion_candidates(
    new_candidates: Optional[pd.DataFrame], existing_candidates: Optional[pd.DataFrame]
) -> Optional[pd.DataFrame]:
    """
    Merge new fusion candidates with existing ones.

    Args:
        new_candidates: New fusion candidates DataFrame
        existing_candidates: Existing fusion candidates DataFrame

    Returns:
        Merged DataFrame or None if no candidates
    """
    if new_candidates is None and existing_candidates is None:
        return None
    elif new_candidates is None:
        return existing_candidates
    elif existing_candidates is None:
        return new_candidates

    # Combine the DataFrames
    combined = pd.concat([existing_candidates, new_candidates], ignore_index=True)

    # Remove duplicates based on key columns (read_id, col1, reference_id, reference_start, reference_end)
    # Keep the first occurrence (existing data takes precedence)
    key_columns = [
        "read_id",
        "col1",
        "reference_id",
        "reference_start",
        "reference_end",
    ]
    if all(col in combined.columns for col in key_columns):
        combined = combined.drop_duplicates(subset=key_columns, keep="first")
        logging.info(
            f"Merged fusion candidates: {len(existing_candidates)} existing + {len(new_candidates)} new -> {len(combined)} total (after deduplication)"
        )
    else:
        logging.warning(
            f"Could not deduplicate fusion candidates - missing key columns: {key_columns}"
        )
        logging.info(
            f"Merged fusion candidates: {len(existing_candidates)} existing + {len(new_candidates)} new -> {len(combined)} total"
        )

    return combined


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
    Perform the actual fusion analysis.

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
    from robin.analysis.fusion_work import process_bam_file

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
    Process a single file for fusion analysis (simplified function-based approach).

    Args:
        file_path: Path to the file to process
        metadata: Metadata about the file
        work_dir: Working directory for output
        logger: Logger instance
        target_panel: Target panel to use for analysis

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

        # Try to open BAM file to check if it's valid
        # try:
        #    with pysam.AlignmentFile(file_path, "rb") as bam:
        #        # Just check if we can read the header
        #        header = bam.header
        #        if logger:
        #            logger.info(f"BAM file header read successfully, {len(header.get('SQ', []))} sequences")
        # except Exception as bam_error:
        #    raise ValueError(f"Invalid or corrupted BAM file {file_path}: {str(bam_error)}")

        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            fusion_metadata.processing_steps.append("temp_dir_created")

            # Perform fusion analysis
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

            # Generate output files
            from robin.analysis.fusion_work import _generate_output_files

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


def fusion_handler(job, work_dir=None):
    """
    Handler function for fusion analysis jobs in the workflow system.

    Args:
        job: Job object from the workflow system
        work_dir: Working directory for output
    """
    try:
        # Get logger with proper parameters
        logger = get_job_logger(str(job.job_id), "fusion", job.context.filepath)

        # Extract file path and metadata from job
        file_path = job.context.filepath
        metadata = job.context.metadata.get("bam_metadata", {})

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

            # Get target panel from metadata or use default
            target_panel = metadata.get("target_panel", "rCNS2")

            # Process the file
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
