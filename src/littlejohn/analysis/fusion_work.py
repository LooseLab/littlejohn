"""
Fusion detection and processing module for BAM files.

This module provides functionality to detect gene fusions from BAM files by analyzing
supplementary alignments and identifying reads that map to multiple gene regions.
"""

# Standard library imports
import os
import random
import logging
from pathlib import Path
from collections import defaultdict
from typing import Dict, Any, Optional, List, Tuple, Set
from dataclasses import dataclass

# Third-party imports
import numpy as np
import pandas as pd
import pysam
import networkx as nx

from littlejohn.analysis.fusion_analysis import FusionMetadata

# Module-level logger
logger = logging.getLogger("littlejohn.analysis.fusion_work")


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
        self, other_start: int, other_end: int, min_overlap: int = 100
    ) -> bool:
        """Check if this region overlaps with another region."""
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
        from littlejohn import resources

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

    return gene_bed, all_gene_bed


def _load_bed_regions(bed_file: str) -> Dict[str, List[GeneRegion]]:
    """Load BED file regions into memory for efficient lookup."""
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
    """Find all reads that have supplementary alignments."""
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
    """Find gene intersections for a read."""
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
    gene_regions_all: Dict[str, List[GeneRegion]],
) -> Optional[pd.DataFrame]:
    """Process reads to find gene intersections and create fusion candidates."""
    rows = []
    rows_all = []

    try:
        # Open the BAM file for reading in binary mode
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            # Iterate through the BAM file but only process reads we need
            # This is still more efficient than the original approach because we can
            # break early once we've found all the reads we need
            for read in bam:
                # Only process reads that are in our target set
                if (
                    read.reference_name == "chrM"
                    or read.is_secondary
                    or read.is_unmapped
                    or read.query_name not in reads_with_supp
                ):  # or read.mapping_quality < 50:
                    continue

                if read.reference_name in gene_regions:
                    # Find intersections between this read and gene regions on this chromosome
                    # This identifies which genes this read overlaps with
                    read_rows = _find_gene_intersections(
                        read,
                        read.reference_name,
                        read.reference_start,
                        read.reference_end,
                        gene_regions[read.reference_name],
                    )
                    # Add any found gene intersections to our results list
                    rows.extend(read_rows)

                if read.reference_name in gene_regions_all:
                    read_rows_all = _find_gene_intersections(
                        read,
                        read.reference_name,
                        read.reference_start,
                        read.reference_end,
                        gene_regions_all[read.reference_name],
                    )
                    rows_all.extend(read_rows_all)

    except Exception as e:
        # If any error occurs during BAM file processing, print the error and return None
        # This prevents the entire analysis from crashing due to file format issues
        print(
            f"\n\n\\nfrom django.utils.translation import ungettextError processing reads for fusions: {str(e)}\n\n\n"
        )
        return None, None

    if rows:
        df = pd.DataFrame(rows)
        df = df[(df["mapping_quality"] > 50) & (df["mapping_span"] > 150)].reset_index(
            drop=True
        )
    else:
        df = None

    if rows_all:
        df_all = pd.DataFrame(rows_all)
        df_all = df_all[
            (df_all["mapping_quality"] > 50) & (df_all["mapping_span"] > 150)
        ].reset_index(drop=True)
    else:
        df_all = None

    return df, df_all


def _filter_fusion_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """Filter fusion candidates based on mapping quality and span."""
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


# =============================================================================
# MAIN PROCESSING FUNCTIONS
# =============================================================================


def process_bam_for_fusions(
    bamfile: str, target_panel: str = "rCNS2", supplementary_read_ids: List[str] = []
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """Process BAM file to find fusion candidates."""
    # Ensure gene regions are loaded
    _ensure_gene_regions_loaded(target_panel)

    try:
        # Find reads with supplementary alignments
        # reads_with_supp = find_reads_with_supplementary(bamfile)
        reads_with_supp = supplementary_read_ids
        if not reads_with_supp:
            return None, None

        df, df_all = _process_reads_for_fusions(
            bamfile,
            reads_with_supp,
            _gene_regions_cache[target_panel],
            _all_gene_regions_cache["shared"],
        )

        # Process reads for target panel fusions
        # target_candidates = _process_reads_for_fusions(
        #    bamfile, reads_with_supp, _gene_regions_cache[target_panel]
        # )
        if df is not None:
            target_candidates = _filter_fusion_candidates(df)
        else:
            target_candidates = None
        if df_all is not None:
            genome_wide_candidates = _filter_fusion_candidates(df_all)
        else:
            genome_wide_candidates = None

        # Process reads for genome-wide fusions
        # genome_wide_candidates = _process_reads_for_fusions(
        #    bamfile, reads_with_supp, _all_gene_regions_cache['shared']
        # )

        return target_candidates, genome_wide_candidates

    except Exception as e:

        print(
            f"\n\n\\nfrom django.utils.translation import ungettextError processing BAM file for fusions: {str(e)}\n\n\n"
        )
        return None, None


def process_bam_file(
    file_path,
    temp_dir,
    metadata,
    fusion_metadata,
    target_panel,
    has_supplementary,
    supplementary_read_ids,
):
    """Process a single BAM file for fusion detection."""
    has_sup = has_supplementary

    if has_sup:
        target_candidates, genome_wide_candidates = process_bam_for_fusions(
            file_path, target_panel, supplementary_read_ids
        )
        if target_candidates is not None and not target_candidates.empty:
            candidates = _optimize_dataframe_memory(target_candidates)

        if genome_wide_candidates is not None and not genome_wide_candidates.empty:
            all_candidates = _optimize_dataframe_memory(genome_wide_candidates)

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

        # Store fusion data in metadata
        fusion_metadata.fusion_data = {
            "target_candidates": (
                target_candidates.to_dict("records")
                if target_candidates is not None
                else []
            ),
            "genome_wide_candidates": (
                genome_wide_candidates.to_dict("records")
                if genome_wide_candidates is not None
                else []
            ),
        }
        return results
    else:
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

    Args:
            sample_id: Sample ID
            analysis_results: Analysis results
            fusion_metadata: Fusion metadata object
            work_dir: Working directory

    Returns:
            Dictionary mapping file type to file path
    """
    if analysis_results["target_candidates"] is not None:
        target_candidates = analysis_results["target_candidates"]
        feather_file = os.path.join(work_dir, sample_id, "target_candidates.feather")
        if os.path.exists(feather_file):
            old_candidates = pd.read_feather(feather_file)
            target_candidates = pd.concat(
                [old_candidates, target_candidates], ignore_index=True
            )
            target_candidates.to_feather(feather_file)
        target_candidates.to_csv(
            os.path.join(work_dir, sample_id, "fusion_candidates_master.csv")
        )
        target_candidates.to_feather(
            os.path.join(work_dir, sample_id, "target_candidates.feather")
        )
        preprocess_fusion_data_standalone(
            target_candidates,
            os.path.join(work_dir, sample_id, "fusion_candidates_master_processed.csv"),
        )

    if analysis_results["genome_wide_candidates"] is not None:
        genome_wide_candidates = analysis_results["genome_wide_candidates"]
        feather_file = os.path.join(
            work_dir, sample_id, "all_target_candidates.feather"
        )
        if os.path.exists(feather_file):
            old_candidates = pd.read_feather(feather_file)
            genome_wide_candidates = pd.concat(
                [old_candidates, genome_wide_candidates], ignore_index=True
            )
            genome_wide_candidates.to_feather(feather_file)
        genome_wide_candidates.to_csv(
            os.path.join(work_dir, sample_id, "fusion_candidates_all.csv")
        )
        genome_wide_candidates.to_feather(
            os.path.join(work_dir, sample_id, "all_target_candidates.feather")
        )
        preprocess_fusion_data_standalone(
            genome_wide_candidates,
            os.path.join(work_dir, sample_id, "fusion_candidates_all_processed.csv"),
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
# FUSION DATA PROCESSING FUNCTIONS
# =============================================================================


def _annotate_results(result: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
    """Annotates the result DataFrame with tags and colors with memory optimization."""
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


def get_gene_network(gene_pairs):
    """This is very slow, but it's the only way to get the gene network."""
    """Build gene network from gene pairs with memory optimization."""
    G = nx.Graph()
    for pair in gene_pairs:
        G.add_edge(pair[0], pair[1])
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


def preprocess_fusion_data_standalone(
    fusion_data: pd.DataFrame, output_file: str
) -> None:
    """Standalone version of fusion data preprocessing for CPU-bound execution with memory optimization."""
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
        # Free the processed data immediately
        del processed_data

    except Exception as e:
        print(f"Error pre-processing fusion data: {str(e)}")
        print("Exception details:", exc_info=True)


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
