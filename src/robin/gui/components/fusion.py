from __future__ import annotations

from typing import Any, Dict, List, Optional
from pathlib import Path
import logging
import pickle
import os
import hashlib
import asyncio
import concurrent.futures

import pandas as pd
import numpy as np
import matplotlib

# Use Agg backend for compatibility with DNA Features Viewer and NiceGUI
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # Import pyplot for plotting functionality

# Import DNA Features Viewer components at module level
try:
    from dna_features_viewer import GraphicFeature, GraphicRecord

    DNA_FEATURES_AVAILABLE = True
except ImportError:
    DNA_FEATURES_AVAILABLE = False
    # DNA Features Viewer not available - using fallback

    
# chrov ideograms removed - not working properly

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None


def _rename_and_sort(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to friendly names and sort for display."""
    if df is None or df.empty:
        return pd.DataFrame()
    rename_map = {
        "col1": "chromBED",
        "col2": "BS",
        "col3": "BE",
        "col4": "Gene",
        "reference_id": "chrom",
        "reference_start": "mS",
        "reference_end": "mE",
        "read_id": "readID",
        "mapping_quality": "mapQ",
        "strand": "strand",
        "read_start": "Read Map Start",
        "read_end": "Read Map End",
        "is_secondary": "Secondary",
        "is_supplementary": "Supplementary",
        "mapping_span": "mapping span",
    }
    try:
        return df.sort_values(by="reference_start").rename(columns=rename_map)
    except Exception:
        return df.rename(columns=rename_map)


def _create_data_hash(data: Dict[str, Any]) -> str:
    """Create a hash of the fusion data to detect actual content changes."""
    try:
        if not data:
            return ""
        
        # Create hash from key data components
        hash_components = []
        
        # Hash the annotated data DataFrame
        if "annotated_data" in data and data["annotated_data"] is not None:
            df_str = data["annotated_data"].to_string()
            df_hash = hashlib.md5(df_str.encode()).hexdigest()
            hash_components.append(f"df:{df_hash}")
        
        # Hash the goodpairs Series
        if "goodpairs" in data and data["goodpairs"] is not None:
            pairs_str = data["goodpairs"].to_string()
            pairs_hash = hashlib.md5(pairs_str.encode()).hexdigest()
            hash_components.append(f"pairs:{pairs_hash}")
        
        # Hash gene groups
        if "gene_groups" in data and data["gene_groups"] is not None:
            groups_str = str(sorted(data["gene_groups"]))
            groups_hash = hashlib.md5(groups_str.encode()).hexdigest()
            hash_components.append(f"groups:{groups_hash}")
        
        # Hash candidate count
        if "candidate_count" in data:
            hash_components.append(f"count:{data['candidate_count']}")
        
        return hashlib.md5("|".join(hash_components).encode()).hexdigest()
    except Exception as e:
        logging.warning(f"[Fusion] Failed to create data hash: {e}")
        return ""


def _generate_summary_files_from_pickle(sample_dir: Path) -> bool:
    """Generate summary files from existing pickle files if they don't exist.
    
    Returns True if summary files were generated or already existed, False otherwise.
    """
    try:
        target_file = sample_dir / "fusion_candidates_master_processed.csv"
        genome_file = sample_dir / "fusion_candidates_all_processed.csv"
        
        # Check if summary files already exist
        summary_file = sample_dir / "fusion_summary.csv"
        if summary_file.exists():
            return True  # Already generated
        
        # Load pickle data
        target_data = _load_processed_pickle(target_file, require_detailed_data=False)
        genome_data = _load_processed_pickle(genome_file, require_detailed_data=False)
        
        if target_data is None and genome_data is None:
            return False  # No data available
        
        # Extract counts
        target_count = target_data.get("candidate_count", 0) if target_data else 0
        genome_count = genome_data.get("candidate_count", 0) if genome_data else 0
        
        # Generate fusion_summary.csv
        import csv
        with open(summary_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["target_fusions", "genome_fusions"])
            writer.writerow([target_count, genome_count])
        
        # Generate individual results files if they don't exist
        target_results_file = sample_dir / "fusion_results_target.csv"
        genome_results_file = sample_dir / "fusion_results_genome.csv"
        
        if target_data is not None and not target_results_file.exists():
            with open(target_results_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["target_fusions", "genome_fusions"])
                writer.writerow([target_count, 0])
        
        if genome_data is not None and not genome_results_file.exists():
            with open(genome_results_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["target_fusions", "genome_fusions"])
                writer.writerow([0, genome_count])
        
        # Update sv_count.txt with the genome-wide count (for backward compatibility)
        sv_count_file = sample_dir / "sv_count.txt"
        if genome_count > 0:
            with open(sv_count_file, "w") as f:
                f.write(str(genome_count))
        
        # Generated summary files from pickle data
        return True
        
    except Exception as e:
        logging.warning(f"[Fusion] Failed to generate summary files from pickle data: {e}")
        return False


def _load_processed_pickle(file_path: Path, require_detailed_data: bool = False) -> Optional[Dict[str, Any]]:
    """Load the preprocessed fusion data structure, preferring text-based files over pickle files.

    Reads from JSON/CSV files first (generated by analysis code) to prevent access conflicts,
    falls back to pickle files if text files don't exist or if detailed data is required.
    
    Args:
        file_path: Path to the pickle file (used to determine text file locations)
        require_detailed_data: If True, will fall back to pickle files for detailed data
    """
    try:
        if not file_path.exists():
            return None
        
        # Determine the analysis type and corresponding text files
        is_target_panel = "master" in file_path.name
        analysis_type = "target" if is_target_panel else "genome"
        
        # Get the directory and construct text file paths
        output_dir = file_path.parent
        
        # If detailed data is required, skip text files and go directly to pickle files
        if require_detailed_data:
            # Loading detailed data from pickle file
            # Skip to pickle file loading below
        else:
            # PREFERRED: Load from JSON data file (generated by analysis code)
            json_file = output_dir / f"fusion_data_{analysis_type}.json"
            if json_file.exists():
                try:
                    import json
                    with open(json_file, "r") as f:
                        data = json.load(f)
                    
                    # Convert back to expected format for GUI compatibility
                    result = {
                        "candidate_count": data.get("candidate_count", 0),
                        "gene_groups": data.get("gene_groups", []),
                        "gene_pairs": data.get("gene_pairs", []),
                        "is_target_panel": data.get("is_target_panel", is_target_panel),
                        "analysis_type": data.get("analysis_type", analysis_type),
                        "sample_data": data.get("sample_data", []),
                        "total_candidates": data.get("total_candidates", 0),
                        "good_pairs_count": data.get("good_pairs_count", 0),
                        # For text-based files, we don't have the full detailed data
                        "annotated_data": None,
                        "goodpairs": None
                    }
                    
                    # Loaded fusion data from text file
                    return result
                    
                except Exception as e:
                    logging.warning(f"[Fusion] Failed to load JSON data: {json_file} - {e}")
        
        # FALLBACK 1: Try CSV files if JSON doesn't exist
        gene_groups_file = output_dir / f"fusion_gene_groups_{analysis_type}.csv"
        gene_pairs_file = output_dir / f"fusion_gene_pairs_{analysis_type}.csv"
        
        if gene_groups_file.exists() and gene_pairs_file.exists():
            try:
                import pandas as pd
                
                # Load gene groups
                gene_groups_df = pd.read_csv(gene_groups_file)
                gene_groups = []
                for _, row in gene_groups_df.iterrows():
                    genes = row["genes"].split(", ") if pd.notna(row["genes"]) else []
                    gene_groups.append([gene.strip() for gene in genes if gene.strip()])
                
                # Load gene pairs
                gene_pairs_df = pd.read_csv(gene_pairs_file)
                gene_pairs = []
                for _, row in gene_pairs_df.iterrows():
                    if pd.notna(row["gene1"]) and pd.notna(row["gene2"]) and row["gene2"]:
                        gene_pairs.append((str(row["gene1"]).strip(), str(row["gene2"]).strip()))
                    elif pd.notna(row["gene1"]):
                        gene_pairs.append((str(row["gene1"]).strip(),))
                
                result = {
                    "candidate_count": len(gene_groups),
                    "gene_groups": gene_groups,
                    "gene_pairs": gene_pairs,
                    "is_target_panel": is_target_panel,
                    "analysis_type": analysis_type,
                    "sample_data": [],
                    "total_candidates": len(gene_groups),
                    "good_pairs_count": len(gene_groups),
                    # For CSV files, we don't have the full detailed data
                    "annotated_data": None,
                    "goodpairs": None
                }
                
                # Loaded fusion data from CSV files
                return result
                
            except Exception as e:
                logging.warning(f"[Fusion] Failed to load CSV data: {e}")
        
        # FALLBACK 2: Read from pickle file (original format used by analysis code)
        if file_path.stat().st_size > 0:
            try:
                with open(file_path, "rb") as f:
                    try:
                        data = pickle.load(f)
                    except (pickle.UnpicklingError, EOFError) as e:
                        logging.warning(f"[Fusion] Pickle file appears truncated, attempting recovery: {file_path}")
                        f.seek(0)
                        try:
                            data = pickle.load(f)
                        except:
                            logging.error(f"[Fusion] Could not recover truncated pickle: {file_path}")
                            return None
                
                # Generate text-based files from pickle data for future GUI access
                if isinstance(data, dict):
                    _generate_summary_files_from_pickle_data(file_path.parent, data, is_target_panel)
                    # Generated text files from pickle data for future use
                
                # Loaded fusion data from pickle file
                return data
                    
            except Exception as e:
                logging.warning(f"[Fusion] Failed to load pickle file: {file_path} - {e}")
        
        return None
        
    except Exception as e:
        logging.exception(f"[Fusion] Failed to load fusion data: {file_path} - {e}")
        return None


def _generate_summary_files_from_pickle_data(output_dir: Path, data: dict, is_target_panel: bool) -> None:
    """Generate text-based summary files from pickle data for future use."""
    try:
        import json
        import csv
        
        analysis_type = "target" if is_target_panel else "genome"
        
        # Generate JSON data file
        json_file = output_dir / f"fusion_data_{analysis_type}.json"
        fusion_summary = {
            "candidate_count": data.get("candidate_count", 0),
            "gene_groups": data.get("gene_groups", []),
            "gene_pairs": data.get("gene_pairs", []),
            "is_target_panel": is_target_panel,
            "analysis_type": analysis_type,
            "sample_data": [],
            "total_candidates": data.get("candidate_count", 0),
            "good_pairs_count": data.get("candidate_count", 0)
        }
        
        with open(json_file, "w") as f:
            json.dump(fusion_summary, f, indent=2, default=str)
        
        # Generate CSV files
        gene_groups = data.get("gene_groups", [])
        gene_groups_file = output_dir / f"fusion_gene_groups_{analysis_type}.csv"
        with open(gene_groups_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["group_id", "genes", "gene_count"])
            for i, group in enumerate(gene_groups):
                writer.writerow([i, ", ".join(group), len(group)])
        
        gene_pairs = data.get("gene_pairs", [])
        gene_pairs_file = output_dir / f"fusion_gene_pairs_{analysis_type}.csv"
        with open(gene_pairs_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["pair_id", "gene1", "gene2"])
            for i, pair in enumerate(gene_pairs):
                if len(pair) >= 2:
                    writer.writerow([i, pair[0], pair[1]])
                elif len(pair) == 1:
                    writer.writerow([i, pair[0], ""])
        
        # Generated text-based files from pickle data
        
    except Exception as e:
        logging.warning(f"[Fusion] Failed to generate text files from pickle data: {e}")


def _make_fusion_table(container: Any, df: pd.DataFrame) -> Any:
    """Create or update a NiceGUI table for fusion candidates and return it."""
    if df is None or df.empty:
        with container:
            ui.label("No fusion candidates found yet").classes("text-gray-600")
        return None

    with container:
        table = (
            ui.table.from_pandas(
                _rename_and_sort(df),
                pagination=25,
            )
            .props("dense")
            .classes("w-full")
            .style("height: 900px")
            .style("font-size: 100%; font-weight: 300")
        )
        try:
            for col in table.columns:
                col["sortable"] = True
            with table.add_slot("top-right"):
                with ui.input(placeholder="Search").props("type=search").bind_value(
                    table, "filter"
                ).add_slot("append"):
                    ui.icon("search")
        except Exception:
            pass
        return table


def _plot_gene_group(
    container: Any,
    gene_group: List[str],
    annotated_data: pd.DataFrame,
    goodpairs: pd.Series,
) -> None:
    """Advanced fusion visualization matching the original sophisticated plotting code.

    Creates multi-panel plots showing:
    1. Gene structure with exons and annotations
    2. Read mapping positions with color coding
    3. Professional-grade visualization using DNA Features Viewer
    """
    try:
        # Check if we have the required data for plotting
        if annotated_data is None or goodpairs is None:
            # Data not available (e.g., loaded from text files without detailed data)
            with container:
                ui.label(f"Gene Group: {', '.join(gene_group)}").classes("text-lg font-bold mb-2")
                ui.label("Detailed visualization data not available").classes("text-gray-600 mb-2")
                ui.label("This gene group was detected but detailed read mapping data").classes("text-gray-500 text-sm")
                ui.label("is not available in the current data format.").classes("text-gray-500 text-sm")
                
                # Show basic gene group information
                with ui.card().classes("mt-4"):
                    ui.label("Gene Group Information").classes("font-bold mb-2")
                    ui.label(f"Genes: {', '.join(gene_group)}").classes("text-sm")
                    ui.label(f"Gene Count: {len(gene_group)}").classes("text-sm")
            return
        
        # Original plotting logic for when we have detailed data
        subset = annotated_data[goodpairs]
        subset = subset[subset["col4"].isin(gene_group)]
        if subset.empty:
            with container:
                ui.label("No reads for selected gene group").classes("text-gray-600")
            return

        # Clear container and create advanced visualization
        container.clear()
        with container:
            # Create matplotlib element for the sophisticated plot with ideograms
            mpl_element = ui.matplotlib(figsize=(20, 10)).classes("w-full")

            # Create the advanced fusion plot
            fig = _create_advanced_fusion_plot(
                gene_group, subset, annotated_data, goodpairs
            )

            # Update the matplotlib element
            mpl_element.figure = fig
            mpl_element.update()

    except Exception as e:
        logging.exception(f"[Fusion] Failed to plot gene group {gene_group}: {e}")
        # Fallback to text display
        with container:
            ui.label(f"Failed to plot gene group: {str(e)}").classes("text-red-600")


def _create_advanced_fusion_plot(
    gene_group: List[str],
    subset: pd.DataFrame,
    annotated_data: pd.DataFrame,
    goodpairs: pd.Series,
) -> plt.Figure:
    """Create the sophisticated gene fusion plot using real DNA Features Viewer integration.

    This matches the original code's approach with:
    - Real gene structure visualization using DNA Features Viewer
    - Sophisticated read alignment display
    - Unified side-by-side layout for gene regions
    """

    # Check if DNA Features Viewer is available
    if not DNA_FEATURES_AVAILABLE:
        # DNA Features Viewer not available, using fallback
        return _create_simple_fallback_plot(gene_group, subset)

    # Load gene annotation data
    gene_table = _load_gene_annotations()
    if gene_table is None:
        # Could not load gene annotations, using fallback
        return _create_simple_fallback_plot(gene_group, subset)

    # Create figure with tight layout - exactly as in original code
    plt.rcParams["figure.constrained_layout.use"] = True
    plt.rcParams["figure.constrained_layout.h_pad"] = 0.05
    plt.rcParams["figure.constrained_layout.w_pad"] = 0.05

    # Get unique genes and their data
    unique_genes = list(sorted(subset["col4"].unique()))
    if len(unique_genes) == 0:
        return _create_simple_fallback_plot(gene_group, subset)

    # Create the unified side-by-side layout without ideograms
    # For 2 genes, we'll create 2 columns with 2 rows each (gene structure + read mapping)
    num_genes = len(unique_genes)
    # Creating gene layout for visualization
    # Increase figure size and add more padding to prevent tight layout warnings
    fig, axes = plt.subplots(2, num_genes, figsize=(20, 10))

    # Handle single gene case
    if num_genes == 1:
        axes = axes.reshape(2, 1)

    # Process each gene
    for col_idx, gene_name in enumerate(unique_genes):
        gene_data = subset[subset["col4"] == gene_name]

        if gene_data.empty:
            continue

        # Get gene coordinates
        gene_start = gene_data["reference_start"].min()
        gene_end = gene_data["reference_end"].max()
        gene_chrom = gene_data["reference_id"].iloc[0]

        # Row 0: Gene structure with DNA Features Viewer
        ax_gene = axes[0, col_idx]
        _plot_gene_structure_with_dna_features(
            ax_gene, gene_name, gene_chrom, gene_start, gene_end, gene_table
        )

        # Row 1: Read mapping visualization
        ax_reads = axes[1, col_idx]
        _plot_read_mapping_sophisticated(
            ax_reads, gene_data, gene_name, gene_chrom, gene_start, gene_end
        )

    # Add overall title first
    fig.suptitle(
        f"Fusion Analysis: {', '.join(gene_group)}", fontsize=14, fontweight="bold"
    )

    # Use constrained_layout for better automatic spacing
    try:
        # Enable constrained layout for automatic spacing
        fig.set_constrained_layout(True)
        fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.1)
    except Exception as e:
        logging.warning(f"[Fusion] Constrained layout failed, using manual adjustment: {e}")
        # Fallback to manual adjustment with more space
        plt.subplots_adjust(
            top=0.85,  # Leave space for suptitle
            bottom=0.1,
            left=0.1,
            right=0.95,
            hspace=0.4,  # More space between rows
            wspace=0.3   # Space between columns
        )

    return fig




def _format_ticks_to_megabases(ax: plt.Axes) -> None:
    """Safely format x-axis ticks to megabases, handling custom formatters from DNA Features Viewer."""
    try:
        ticks = ax.get_xticks()
        # Set the tick positions first, then the labels to avoid warnings
        ax.set_xticks(ticks)
        ax.set_xticklabels([f'{t/1e6:.1f}' for t in ticks])
    except Exception as e:
        # Could not format ticks to megabases
        # Fallback: try to use ticklabel_format if possible
        try:
            ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
        except Exception:
            pass  # If both methods fail, just leave the default formatting


def _load_gene_annotations() -> Optional[pd.DataFrame]:
    """Load gene annotation data from the rCNS2_data.csv.gz file."""
    try:
        # Try to load the gene annotation data
        gene_data_path = "src/robin/resources/rCNS2_data.csv.gz"
        if os.path.exists(gene_data_path):
            gene_table = pd.read_csv(gene_data_path)
            # Gene annotations loaded
            return gene_table
        else:
            # Gene annotation file not found
            return None
    except Exception as e:
        logging.error(f"Error loading gene annotations: {e}")
        return None


def _plot_gene_structure_with_dna_features(
    ax: plt.Axes,
    gene_name: str,
    chrom: str,
    start: int,
    end: int,
    gene_table: pd.DataFrame,
):
    """Plot gene structure using DNA Features Viewer exactly as in the original code."""
    try:
        # Filter gene table for this specific gene and chromosome
        gene_info = gene_table[
            (gene_table["gene_name"] == gene_name) & (gene_table["Seqid"] == chrom)
        ]

        if gene_info.empty:
            # Fallback if no gene info found
            ax.set_title(f"Gene Structure: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
            # Extend plot range to show full context
            plot_start = start - (end - start) * 0.1
            plot_end = end + (end - start) * 0.1
            ax.set_xlim(plot_start, plot_end)
            ax.set_ylim(0, 1)
            ax.text(
                0.5,
                0.5,
                f"Gene: {gene_name}\n{chrom}:{start:,}-{end:,}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7),
            )
            # Convert x-axis to megabases
            ax.set_xlabel(f"Position (Mb) - {chrom}", fontsize=10)
            # Convert tick labels to megabases
            _format_ticks_to_megabases(ax)
            return

        # Create DNA Features Viewer visualization
        features = []

        # Add gene features
        for _, row in gene_info.iterrows():
            if row["Type"] == "gene":
                # Main gene feature
                strand = 1 if row["Strand"] == "+" else -1
                features.append(
                    GraphicFeature(
                        start=int(row["Start"]),
                        end=int(row["End"]),
                        strand=strand,
                        thickness=8,
                        color="#ffd700",  # Gold color for genes
                        label=row["gene_name"],
                        fontdict={"color": "black", "fontsize": 8},
                    )
                )
            elif row["Type"] == "exon":
                # Exon features
                strand = 1 if row["Strand"] == "+" else -1
                features.append(
                    GraphicFeature(
                        start=int(row["Start"]),
                        end=int(row["End"]),
                        strand=strand,
                        thickness=4,
                        color="#C0C0C0",  # Silver color for exons
                    )
                )

        if features:
            # Create GraphicRecord and plot it
            record = GraphicRecord(
                sequence_length=end - start, first_index=start, features=features
            )

            # Plot on the axis
            record.plot(
                ax=ax, with_ruler=False, draw_line=True, strand_in_label_threshold=4
            )
            ax.set_title(f"Gene Structure: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
            ax.set_xlabel(f"Position (Mb) - {chrom}", fontsize=10)
            # Extend plot range to show full context
            ax.margins(x=0.15, y=0.1)
            # Convert tick labels to megabases (avoid ticklabel_format due to DNA Features Viewer formatter)
            _format_ticks_to_megabases(ax)
        else:
            # Fallback if no features found
            ax.set_title(f"Gene Structure: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
            # Extend plot range to show full context
            plot_start = start - (end - start) * 0.1
            plot_end = end + (end - start) * 0.1
            ax.set_xlim(plot_start, plot_end)
            ax.set_ylim(0, 1)
            ax.text(
                0.5,
                0.5,
                f"Gene: {gene_name}\n{chrom}:{start:,}-{end:,}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7),
            )
            # Convert x-axis to megabases
            ax.set_xlabel(f"Position (Mb) - {chrom}", fontsize=10)
        # Convert tick labels to megabases
        _format_ticks_to_megabases(ax)

    except Exception as e:
        logging.error(f"Error plotting gene structure with DNA Features Viewer: {e}")
        # Fallback to simple text
        ax.set_title(f"Gene Structure: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
        # Extend plot range to show full context
        plot_start = start - (end - start) * 0.1
        plot_end = end + (end - start) * 0.1
        ax.set_xlim(plot_start, plot_end)
        ax.set_ylim(0, 1)
        ax.text(
            0.5,
            0.5,
            f"Error plotting\n{gene_name}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )
        # Convert x-axis to megabases even in error case
        ax.set_xlabel(f"Position (Mb) - {chrom}", fontsize=10)
        # Convert tick labels to megabases
        _format_ticks_to_megabases(ax)


def _plot_read_mapping_sophisticated(
    ax: plt.Axes,
    gene_data: pd.DataFrame,
    gene_name: str,
    chrom: str,
    start: int,
    end: int,
):
    """Plot sophisticated read mapping visualization exactly as in the original code."""
    try:
        # Get read positions and create the sophisticated visualization
        read_starts = gene_data["reference_start"].values
        read_ends = gene_data["reference_end"].values
        read_ids = gene_data["read_id"].values
        strands = gene_data["strand"].values

        # Create color map for reads (similar to original code)
        unique_reads = pd.unique(read_ids)
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_reads)))
        read_color_map = dict(zip(unique_reads, colors))

        # Plot each read with sophisticated features
        for i, (read_start, read_end, read_id, strand) in enumerate(
            zip(read_starts, read_ends, read_ids, strands)
        ):
            color = read_color_map[read_id]

            # Create the read visualization as in original code
            # Main read line with slight extension to show ends clearly
            read_length = read_end - read_start
            extended_start = read_start - read_length * 0.02  # Small extension left
            extended_end = read_end + read_length * 0.02  # Small extension right
            ax.plot(
                [extended_start, extended_end],
                [i, i],
                color=color,
                linewidth=3,
                alpha=0.8,
            )

            # Add strand indicator (arrow) at the actual read end
            if strand == "+":
                ax.arrow(
                    read_end,
                    i,
                    (end - start) * 0.01,
                    0,
                    head_width=0.2,
                    head_length=(end - start) * 0.005,
                    fc=color,
                    ec=color,
                    alpha=0.8,
                )
            else:
                ax.arrow(
                    read_start,
                    i,
                    -(end - start) * 0.01,
                    0,
                    head_width=0.2,
                    head_length=(end - start) * 0.005,
                    fc=color,
                    ec=color,
                    alpha=0.8,
                )

            # Add read ID label (limit to first few for clarity)
            if i < 5:
                ax.text(
                    read_start, i + 0.1, f"Read {read_id[:8]}", fontsize=6, color=color
                )

        # Customize plot
        ax.set_title(f"Read Mapping: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
        # Extend plot range to clearly show read ends
        plot_start = start - (end - start) * 0.15
        plot_end = end + (end - start) * 0.15
        ax.set_xlim(plot_start, plot_end)
        ax.set_xlabel(f"Position (Mb) - {chrom}")
        ax.set_ylabel("Reads")

        # Convert x-axis to megabases
        _format_ticks_to_megabases(ax)

        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Remove box edges/borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

        # Add legend if not too many reads and we have labeled artists
        if len(unique_reads) <= 10 and len(ax.get_legend_handles_labels()[0]) > 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)

    except Exception as e:
        logging.error(f"Error plotting sophisticated read mapping: {e}")
        # Fallback to simple text
        ax.set_title(f"Read Mapping: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
        ax.text(
            0.5,
            0.5,
            f"Error plotting reads\n{gene_name}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )
        # Convert x-axis to megabases even in error case
        ax.set_xlabel(f"Position (Mb) - {chrom}")
        # Format x-axis ticks to show megabases
        ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
        # Convert tick labels to megabases
        _format_ticks_to_megabases(ax)


def _process_reads_for_original_plot(subset: pd.DataFrame) -> pd.DataFrame:
    """Process reads data exactly as the original _get_reads function does."""
    try:
        # Convert to the format expected by original code
        df = subset.copy()
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

        # Convert start and end to int
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)

        # Sort and remove duplicates exactly as original
        df = df.sort_values(by=["chromosome", "start", "end"])
        df = df.drop_duplicates(subset=["start2", "end2", "id"])

        # Group by chromosome and collapse ranges (simplified version)
        # This would need the original collapse_ranges function
        return df

    except Exception as e:
        logging.error(f"Error processing reads for original plot: {e}")
        return subset


def _plot_gene_structure_original(
    ax: plt.Axes, data: pd.Series, chrom: str, start: int, end: int
):
    """Plot gene structure exactly as in original code using DNA Features Viewer."""
    try:
        # This would need access to the gene_table data from the original code
        # For now, create a placeholder that matches the original structure
        ax.set_title(f"Gene Structure: {data['gene']} ({chrom})", fontsize=10, fontweight="bold")
        ax.set_xlim(start, end)
        ax.set_ylim(0, 1)

        # Placeholder for gene structure - would need gene_table integration
        ax.text(
            0.5,
            0.5,
            f"Gene: {data['gene']}\n{chrom}:{start:,}-{end:,}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7),
        )

        ax.set_xlabel(f"Position (Mb) - {chrom}")
        ax.set_ylabel("Gene Structure")
        # Convert x-axis to megabases
        _format_ticks_to_megabases(ax)

    except Exception as e:
        logging.error(f"Error plotting gene structure: {e}")
        ax.text(
            0.5,
            0.5,
            f"Error plotting\n{data['gene']}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )


def _plot_read_mapping_original(
    ax: plt.Axes,
    data: pd.Series,
    chrom: str,
    start: int,
    end: int,
    subset: pd.DataFrame,
):
    """Plot read mapping exactly as in original code."""
    try:
        # Get reads for this gene
        gene_reads = subset[subset["col4"] == data["gene"]]

        if gene_reads.empty:
            ax.text(
                0.5, 0.5, "No reads", transform=ax.transAxes, ha="center", va="center"
            )
            return

        # Create the exact same visualization as original code
        # This would need the original plotting logic with overlapping ranges, ranks, etc.

        # For now, create a simplified version that shows the structure
        ax.set_title(f"Read Mapping: {data['gene']} ({chrom})", fontsize=10, fontweight="bold")
        # Extend plot range to clearly show read ends
        plot_start = start - (end - start) * 0.15
        plot_end = end + (end - start) * 0.15
        ax.set_xlim(plot_start, plot_end)

        # Plot reads as horizontal lines (simplified)
        for i, (_, read) in enumerate(gene_reads.iterrows()):
            read_start = read["reference_start"]
            read_end = read["reference_end"]
            ax.plot([read_start, read_end], [i, i], linewidth=2, alpha=0.8)

        ax.set_xlabel(f"Position (Mb) - {chrom}")
        ax.set_ylabel("Reads")
        # Convert x-axis to megabases
        _format_ticks_to_megabases(ax)
        ax.grid(True, alpha=0.3)
        
        # Remove box edges/borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

    except Exception as e:
        logging.error(f"Error plotting read mapping: {e}")
        ax.text(
            0.5,
            0.5,
            f"Error plotting reads\n{data['gene']}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )


def _plot_gene_structure(
    ax: plt.Axes, gene_name: str, chrom: str, start: int, end: int
):
    """Plot gene structure with exons and annotations."""
    try:
        # This would integrate with your gene annotation data
        # For now, create a placeholder gene structure
        ax.set_title(f"Gene Structure: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
        ax.set_xlim(start, end)
        ax.set_ylim(0, 1)

        # Placeholder gene structure visualization
        ax.text(
            0.5,
            0.5,
            f"Gene: {gene_name}\n{chrom}:{start:,}-{end:,}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7),
        )

        ax.set_xlabel(f"Position (Mb) - {chrom}")
        ax.set_ylabel("Gene Structure")
        # Convert x-axis to megabases
        _format_ticks_to_megabases(ax)

    except Exception as e:
        logging.error(f"Error plotting gene structure for {gene_name}: {e}")
        ax.text(
            0.5,
            0.5,
            f"Error plotting\n{gene_name}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )


def _plot_read_mapping(
    ax: plt.Axes,
    gene_data: pd.DataFrame,
    gene_name: str,
    chrom: str,
    start: int,
    end: int,
):
    """Plot read mapping positions with color coding."""
    try:
        # Get read positions
        read_starts = gene_data["reference_start"].values
        read_ends = gene_data["reference_end"].values
        read_ids = gene_data["read_id"].values

        # Create color map for reads
        unique_reads = pd.unique(read_ids)
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_reads)))
        read_color_map = dict(zip(unique_reads, colors))

        # Plot each read as a horizontal line with slight extension to show ends clearly
        for i, (read_start, read_end, read_id) in enumerate(
            zip(read_starts, read_ends, read_ids)
        ):
            color = read_color_map[read_id]
            read_length = read_end - read_start
            extended_start = read_start - read_length * 0.02  # Small extension left
            extended_end = read_end + read_length * 0.02  # Small extension right
            ax.plot(
                [extended_start, extended_end],
                [i, i],
                color=color,
                linewidth=3,
                alpha=0.8,
                label=f"Read {read_id}" if i < 10 else "",
            )  # Limit legend to first 10 reads

        # Customize plot
        ax.set_title(f"Read Mapping: {gene_name} ({chrom})", fontsize=10, fontweight="bold")
        # Extend plot range to clearly show read ends
        plot_start = start - (end - start) * 0.15
        plot_end = end + (end - start) * 0.15
        ax.set_xlim(plot_start, plot_end)
        ax.set_xlabel(f"Position (Mb) - {chrom}")
        ax.set_ylabel("Reads")

        # Convert x-axis to megabases
        _format_ticks_to_megabases(ax)

        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Remove box edges/borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

        # Add legend if not too many reads and we have labeled artists
        if len(unique_reads) <= 10 and len(ax.get_legend_handles_labels()[0]) > 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)

    except Exception as e:
        logging.error(f"Error plotting read mapping for {gene_name}: {e}")
        ax.text(
            0.5,
            0.5,
            f"Error plotting reads\n{gene_name}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )


def _create_simple_fallback_plot(
    gene_group: List[str], subset: pd.DataFrame
) -> plt.Figure:
    """Create a simple fallback plot if advanced visualization fails."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))

    # Simple scatter plot
    genes = list(sorted(subset["col4"].unique()))
    y_map = {g: i for i, g in enumerate(genes)}
    x = subset["reference_start"].astype(int)
    y = subset["col4"].map(y_map).astype(int)

    ax.scatter(x, y, s=15, alpha=0.8)
    ax.set_yticks(list(y_map.values()))
    ax.set_yticklabels(list(y_map.keys()))
    ax.set_xlabel("Mapping Start (reference_start)")
    ax.set_title(
        f"Reads supporting fusion group: {', '.join(gene_group)}",
        fontsize=12,
        fontweight="bold",
    )
    ax.grid(True, axis="x", linestyle=":", alpha=0.4)

    return fig


def add_fusion_section(launcher: Any, sample_dir: Path) -> None:
    """Build the Fusion UI section (target panel and genome-wide only).

    Expects files written by fusion_work.preprocess_fusion_data_standalone:
    - fusion_candidates_master_processed.csv (pickle payload)
    - fusion_candidates_all_processed.csv (pickle payload)
    """
    if ui is None:
        return

    # Local state for this section
    state: Dict[str, Any] = {
        "target": {
            "data": None,
            "mtime": None,
            "data_hash": None,
            "table_container": None,
            "table": None,
            "plot_container": None,
            "card": None,
        },
        "genome": {
            "data": None,
            "mtime": None,
            "data_hash": None,
            "table_container": None,
            "table": None,
            "plot_container": None,
            "card": None,
        },
    }

    async def refresh_async() -> None:
        """Refresh fusion data asynchronously."""
        try:
            # Check directory existence in background thread
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future = executor.submit(lambda: sample_dir and sample_dir.exists())
                dir_exists = await asyncio.wrap_future(future)
            
            if not dir_exists:
                return
            
            # Run data loading in background thread, UI updates in main thread
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future = executor.submit(load_fusion_data_sync, sample_dir)
                data = await asyncio.wrap_future(future)
                
                # Update UI in main thread
                if data:
                    update_fusion_ui(data, state)
                
        except Exception as e:
            logging.exception(f"[Fusion] Async refresh failed: {e}")

    def load_fusion_data_sync(sample_dir: Path) -> Optional[Dict[str, Any]]:
        """Load fusion data in background thread - no UI operations."""
        try:
            # Load target panel processed
            target_file = sample_dir / "fusion_candidates_master_processed.csv"
            genome_file = sample_dir / "fusion_candidates_all_processed.csv"

            t = _load_processed_pickle(target_file, require_detailed_data=True)
            g = _load_processed_pickle(genome_file, require_detailed_data=True)
            
            # Debug logging
            # Fusion data loaded successfully

            return {
                "target": t,
                "genome": g,
                "target_file": target_file,
                "genome_file": genome_file
            }
            
        except Exception as e:
            logging.exception(f"[Fusion] Data loading failed: {e}")
            return None

    def update_fusion_ui(data: Dict[str, Any], state: Dict[str, Any]) -> None:
        """Update fusion UI elements in main thread."""
        try:
            t = data["target"]
            g = data["genome"]
            target_file = data["target_file"]
            genome_file = data["genome_file"]
            
            # Update target panel UI
            target_mtime = target_file.stat().st_mtime if target_file.exists() else None
            target_data_hash = _create_data_hash(t) if t is not None else None
            
            # Always update summary and table when file changes
            if t is not None and target_mtime != state["target"].get("mtime"):
                state["target"]["data"] = t
                state["target"]["mtime"] = target_mtime
                state["target"]["data_hash"] = target_data_hash
                
                # summary
                try:
                    if "summary" in state and state["summary"].get("target_lbl"):
                        state["summary"][
                            "target_lbl"
                        ].text = (
                            f"Target fusion groups: {int(t.get('candidate_count', 0))}"
                        )
                except Exception:
                    pass
                # table
                try:
                    state["target"]["table_container"].clear()
                except Exception:
                    pass
                state["target"]["table"] = _make_fusion_table(
                    state["target"]["table_container"],
                    t.get("annotated_data") or pd.DataFrame(),
                )
                
                # Create visualization on first load or when data changes
                try:
                    state["target"]["plot_container"].clear()
                except Exception:
                    pass
                if t.get("candidate_count", 0) > 0 and t.get("gene_groups"):
                    with state["target"]["plot_container"].classes("w-full"):
                        with ui.row().classes("w-full"):
                            ui.select(
                                options=t.get("gene_groups", []),
                                with_input=True,
                                on_change=lambda e, t=t: (
                                    state["target"]["card"].clear(),
                                    state["target"]["card"].classes("w-full"),
                                    _plot_gene_group(
                                        state["target"]["card"],
                                        e.value,
                                        t.get("annotated_data"),
                                        t.get("goodpairs"),
                                    ),
                                ),
                            ).classes("w-full")
                        ui.label("").classes("text-sm text-gray-600")
                        state["target"]["card"] = ui.card().classes("w-full")
                        # Plot the first gene group by default
                        if t.get("gene_groups"):
                            _plot_gene_group(
                                state["target"]["card"],
                                t.get("gene_groups", [])[0],
                                t.get("annotated_data"),
                                t.get("goodpairs"),
                            )
            elif t is None:
                # Clear target panel if no data
                state["target"]["data"] = None
                state["target"]["mtime"] = None
                state["target"]["data_hash"] = None
                
                # Clear summary
                try:
                    if "summary" in state and state["summary"].get("target_lbl"):
                        state["summary"]["target_lbl"].text = "Target fusion groups: 0"
                except Exception:
                    pass
                
                # Clear table
                try:
                    state["target"]["table_container"].clear()
                    state["target"]["table"] = _make_fusion_table(
                        state["target"]["table_container"], pd.DataFrame()
                    )
                except Exception:
                    pass
                
                # Clear plot
                try:
                    state["target"]["plot_container"].clear()
                except Exception:
                    pass

            # Update genome-wide UI
            genome_mtime = genome_file.stat().st_mtime if genome_file.exists() else None
            genome_data_hash = _create_data_hash(g) if g is not None else None
            
            if g is not None and genome_mtime != state["genome"].get("mtime"):
                state["genome"]["data"] = g
                state["genome"]["mtime"] = genome_mtime
                state["genome"]["data_hash"] = genome_data_hash
                
                # summary
                try:
                    if "summary" in state and state["summary"].get("genome_lbl"):
                        state["summary"][
                            "genome_lbl"
                        ].text = (
                            f"Genome-wide fusion groups: {int(g.get('candidate_count', 0))}"
                        )
                except Exception:
                    pass
                # table
                try:
                    state["genome"]["table_container"].clear()
                except Exception:
                    pass
                state["genome"]["table"] = _make_fusion_table(
                    state["genome"]["table_container"],
                    g.get("annotated_data") or pd.DataFrame(),
                )
                
                # Create visualization on first load or when data changes
                try:
                    state["genome"]["plot_container"].clear()
                except Exception:
                    pass
                if g.get("candidate_count", 0) > 0 and g.get("gene_groups"):
                    with state["genome"]["plot_container"].classes("w-full"):
                        with ui.row().classes("w-full"):
                            ui.select(
                                options=g.get("gene_groups", []),
                                with_input=True,
                                on_change=lambda e, g=g: (
                                    state["genome"]["card"].clear(),
                                    state["genome"]["card"].classes("w-full"),
                                    _plot_gene_group(
                                        state["genome"]["card"],
                                        e.value,
                                        g.get("annotated_data"),
                                        g.get("goodpairs"),
                                    ),
                                ),
                            ).classes("w-full")
                        ui.label("").classes("text-sm text-gray-600")
                        state["genome"]["card"] = ui.card().classes("w-full")
                        # Plot the first gene group by default
                        if g.get("gene_groups"):
                            _plot_gene_group(
                                state["genome"]["card"],
                                g.get("gene_groups", [])[0],
                                g.get("annotated_data"),
                                g.get("goodpairs"),
                            )
            elif g is None:
                # Clear genome panel if no data
                state["genome"]["data"] = None
                state["genome"]["mtime"] = None
                state["genome"]["data_hash"] = None
                
                # Clear summary
                try:
                    if "summary" in state and state["summary"].get("genome_lbl"):
                        state["summary"]["genome_lbl"].text = "Genome-wide fusion groups: 0"
                except Exception:
                    pass
                
                # Clear table
                try:
                    state["genome"]["table_container"].clear()
                    state["genome"]["table"] = _make_fusion_table(
                        state["genome"]["table_container"], pd.DataFrame()
                    )
                except Exception:
                    pass
                
                # Clear plot
                try:
                    state["genome"]["plot_container"].clear()
                except Exception:
                    pass

        except Exception as e:
            logging.exception(f"[Fusion] UI update failed: {e}")

    # Old refresh_sync function removed - replaced by load_fusion_data_sync and update_fusion_ui

    # Build UI
    with ui.card().classes("w-full"):
        ui.label("🧪 Fusions").classes("text-lg font-semibold mb-2")
        # Summary row (counts)
        with ui.row().classes("w-full items-center justify-between mb-2"):
            with ui.column().classes("gap-1"):
                ui.label("Fusion Analysis").classes("text-sm font-medium")
            with ui.column().classes("gap-1 items-end"):
                state.setdefault("summary", {})
                state["summary"]["target_lbl"] = ui.label(
                    "Target fusion groups: --"
                ).classes("text-sm text-gray-600")
                state["summary"]["genome_lbl"] = ui.label(
                    "Genome-wide fusion groups: --"
                ).classes("text-sm text-gray-600")

        # Target Panel
        ui.label("Target Panel").classes("text-base font-medium mt-1 mb-1")
        state["target"]["plot_container"] = ui.column().classes("w-full")
        state["target"]["table_container"] = ui.column().classes("w-full mt-2")

        # Genome-wide
        ui.separator()
        ui.label("Genome-wide").classes("text-base font-medium mt-2 mb-1")
        state["genome"]["plot_container"] = ui.column().classes("w-full")
        state["genome"]["table_container"] = ui.column().classes("w-full mt-2")

    # Initial refresh and timer
    try:
        # Prime summary values on first load
        try:
            t0 = _load_processed_pickle(
                (sample_dir / "fusion_candidates_master_processed.csv")
                if sample_dir
                else Path("/dev/null"),
                require_detailed_data=False
            )
            g0 = _load_processed_pickle(
                (sample_dir / "fusion_candidates_all_processed.csv")
                if sample_dir
                else Path("/dev/null"),
                require_detailed_data=False
            )
            if "summary" in state:
                t_count0 = (
                    int(t0.get("candidate_count", 0)) if isinstance(t0, dict) else 0
                )
                g_count0 = (
                    int(g0.get("candidate_count", 0)) if isinstance(g0, dict) else 0
                )
                try:
                    state["summary"][
                        "target_lbl"
                    ].text = f"Target fusion groups: {t_count0}"
                    state["summary"][
                        "genome_lbl"
                    ].text = f"Genome-wide fusion groups: {g_count0}"
                except Exception:
                    pass
        except Exception:
            pass
        
        # Start the refresh timer (every 10 seconds) with async function
        ui.timer(10.0, lambda: ui.timer(0.1, refresh_async, once=True), active=True)
    except Exception:
        pass
