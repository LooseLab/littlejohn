from __future__ import annotations

from typing import Any, Dict, List, Optional
from pathlib import Path
import logging
import pickle
import os
import hashlib

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
    logging.warning("DNA Features Viewer not available - using fallback visualization")

    
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


def _generate_summary_files_from_pickle(sample_dir: Path, force_regenerate: bool = False) -> bool:
    """Generate summary files from existing pickle files if they don't exist.
    
    This provides backward compatibility for existing analyses that were run
    before the summary file generation was implemented.
    
    Returns True if summary files were generated, False otherwise.
    """
    try:
        import csv
        from pathlib import Path
        
        # Check if summary files already exist
        summary_file = sample_dir / "fusion_summary.csv"
        if summary_file.exists() and not force_regenerate:
            return True  # Already have summary files
        elif summary_file.exists() and force_regenerate:
            pass  # Will overwrite existing file
        
        # Load data from pickle files
        target_file = sample_dir / "fusion_candidates_master_processed.csv"
        genome_file = sample_dir / "fusion_candidates_all_processed.csv"
        
        logging.info(f"[Fusion] Checking for target file: {target_file}")
        logging.info(f"[Fusion] Target file exists: {target_file.exists()}")
        logging.info(f"[Fusion] Checking for genome-wide file: {genome_file}")
        logging.info(f"[Fusion] Genome-wide file exists: {genome_file.exists()}")
        
        # Debug: List all fusion-related files in the directory
        fusion_files = list(sample_dir.glob("*fusion*"))
        logging.info(f"[Fusion] All fusion files in directory: {[f.name for f in fusion_files]}")
        
        # Debug: Check file sizes
        if target_file.exists():
            logging.info(f"[Fusion] Target file size: {target_file.stat().st_size} bytes")
        if genome_file.exists():
            logging.info(f"[Fusion] Genome-wide file size: {genome_file.stat().st_size} bytes")
        
        target_data = _load_processed_pickle(target_file)
        genome_data = _load_processed_pickle(genome_file)
        
        if target_data is not None:
            logging.info(f"[Fusion] Target data loaded: candidate_count={target_data.get('candidate_count', 0)}")
        else:
            logging.warning(f"[Fusion] Failed to load target data from: {target_file}")
            
        if genome_data is not None:
            logging.info(f"[Fusion] Genome-wide data loaded: candidate_count={genome_data.get('candidate_count', 0)}")
        else:
            logging.warning(f"[Fusion] Failed to load genome-wide data from: {genome_file}")
        
        # Extract counts
        target_count = 0
        genome_count = 0
        
        if target_data is not None and isinstance(target_data, dict):
            target_count = target_data.get("candidate_count", 0)
        
        if genome_data is not None and isinstance(genome_data, dict):
            # Use gene_pairs count if candidate_count is 0 (like reporting code does)
            genome_count = genome_data.get("candidate_count", 0)
            if genome_count == 0 and genome_data.get("gene_pairs"):
                genome_count = len(genome_data.get("gene_pairs", []))
                logging.info(f"[Fusion] Summary: Using gene_pairs count for genome-wide: {genome_count}")
        
        # Generate fusion_summary.csv
        with open(summary_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["target_fusions", "genome_fusions"])
            writer.writerow([target_count, genome_count])
        
        logging.info(f"[Fusion] Generated summary file from pickle: target={target_count}, genome={genome_count}")
        
        # Generate fusion_results.csv (use the one with more data)
        results_file = sample_dir / "fusion_results.csv"
        if target_count > 0:
            # Use target panel data
            with open(results_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["target_fusions", "genome_fusions"])
                writer.writerow([target_count, 0])
        elif genome_count > 0:
            # Use genome-wide data
            with open(results_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["target_fusions", "genome_fusions"])
                writer.writerow([0, genome_count])
        
        # Generate sv_count.txt (use genome count for backward compatibility)
        sv_count_file = sample_dir / "sv_count.txt"
        with open(sv_count_file, "w") as f:
            f.write(str(genome_count))
        
        logging.info(f"[Fusion] Generated summary files from pickle data - target: {target_count}, genome: {genome_count}")
        return True
        
    except Exception as e:
        logging.warning(f"[Fusion] Failed to generate summary files from pickle: {e}")
        return False


def _load_processed_pickle(file_path: Path) -> Optional[Dict[str, Any]]:
    """Load the preprocessed fusion data structure written by preprocess_fusion_data_standalone.

    Note: Although the extension is .csv, the file is a pickle per current pipeline.
    """
    try:
        logging.info(f"[Fusion] Attempting to load pickle from: {file_path}")
        if not file_path.exists():
            logging.warning(f"[Fusion] File does not exist: {file_path}")
            return None
        if file_path.stat().st_size == 0:
            logging.warning(f"[Fusion] File is empty: {file_path}")
            return None
        
        # Try to load the pickle with better error handling
        with open(file_path, "rb") as f:
            try:
                data = pickle.load(f)
            except (pickle.UnpicklingError, EOFError) as e:
                # If pickle is truncated, try to load what we can
                logging.warning(f"[Fusion] Pickle file appears truncated, attempting recovery: {file_path}")
                f.seek(0)
                try:
                    # Try loading with protocol 0 which is more forgiving
                    data = pickle.load(f)
                except:
                    # If all else fails, return None and let the system regenerate
                    logging.error(f"[Fusion] Could not recover truncated pickle: {file_path}")
                    return None
        
        # Expected keys: annotated_data (DataFrame), goodpairs (Series), gene_groups (list), candidate_count (int)
        if isinstance(data, dict):
            logging.info(f"[Fusion] Loaded pickle data with keys: {list(data.keys())}")
            logging.info(f"[Fusion] Raw candidate_count: {data.get('candidate_count', 'not found')}")
            logging.info(f"[Fusion] Raw gene_groups count: {len(data.get('gene_groups', []))}")
            
            # Apply the same filtering logic as the reporting code
            annotated_data = data.get("annotated_data", pd.DataFrame())
            goodpairs = data.get("goodpairs", pd.Series())
            
            logging.info(f"[Fusion] Raw annotated_data shape: {annotated_data.shape}")
            logging.info(f"[Fusion] Raw goodpairs length: {len(goodpairs)}")
            
            if not annotated_data.empty and not goodpairs.empty:
                # Only keep the good pairs (same as reporting code does)
                data["annotated_data"] = annotated_data[goodpairs]
                logging.info(f"[Fusion] Filtered data: {len(annotated_data)} -> {len(data['annotated_data'])} good pairs")
            else:
                logging.info(f"[Fusion] No filtering applied - annotated_data empty: {annotated_data.empty}, goodpairs empty: {goodpairs.empty}")
            
            logging.info(f"[Fusion] Final candidate_count: {data.get('candidate_count', 'not found')}")
            logging.info(f"[Fusion] Final gene_groups count: {len(data.get('gene_groups', []))}")
            
            return data
        else:
            logging.warning(f"[Fusion] Loaded data is not a dict, type: {type(data)}")
            return None
    except Exception as e:
        logging.exception(
            f"[Fusion] Failed to load processed fusion pickle: {file_path} - {e}"
        )
        return None


def _make_fusion_table(container: Any, df: pd.DataFrame) -> Any:
    """Create or update a NiceGUI table for fusion candidates and return it."""
    if df is None or df.empty:
        with container:
            ui.label("No fusion candidates found yet").classes("text-gray-600")
        return None

    with container:
        # Use styled_table for consistent styling
        from robin.gui.theme import styled_table
        
        # Convert DataFrame to rows format for styled_table
        processed_df = _rename_and_sort(df)
        
        # Create columns definition from DataFrame
        columns = []
        for col in processed_df.columns:
            columns.append({
                "name": col,
                "label": col.replace("_", " ").title(),
                "field": col,
                "sortable": True
            })
        
        # Create rows from DataFrame
        rows = processed_df.to_dict('records')
        
        # Create styled table
        table_container, table = styled_table(
            columns=columns,
            rows=rows,
            pagination=25,
            class_size="table-xs"
        )
        
        # Add search functionality
        try:
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
        # For genome-wide fusions, use raw annotated_data if goodpairs filters out everything
        if goodpairs.sum() == 0:
            # No good pairs - use raw data (like reporting code does)
            # Strip whitespace from gene names to handle leading/trailing spaces
            subset = annotated_data[annotated_data["col4"].str.strip().isin(gene_group)]
            logging.info(f"[Fusion] Using raw data for genome-wide gene group {gene_group}: {len(subset)} rows")
        else:
            # Target panel fusions - use filtered data
            # Fix pandas reindexing warning by ensuring indices match
            try:
                # Align indices to avoid reindexing warning
                aligned_goodpairs = goodpairs.reindex(annotated_data.index, fill_value=False)
                subset = annotated_data[aligned_goodpairs]
                
                # Strip whitespace from gene names to handle leading/trailing spaces
                subset = subset[subset["col4"].str.strip().isin(gene_group)]
                
                logging.info(f"[Fusion] Using filtered data for target gene group {gene_group}: {len(subset)} rows")
            except Exception as e:
                # Fallback to raw data if indexing fails
                logging.warning(f"[Fusion] Indexing failed, using raw data: {e}")
                # Strip whitespace from gene names to handle leading/trailing spaces
                subset = annotated_data[annotated_data["col4"].str.strip().isin(gene_group)]
                logging.info(f"[Fusion] Using raw data (fallback) for gene group {gene_group}: {len(subset)} rows")
        
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
        logging.warning("DNA Features Viewer not available, using fallback plot")
        return _create_simple_fallback_plot(gene_group, subset)

    # Load gene annotation data
    gene_table = _load_gene_annotations()
    if gene_table is None:
        logging.warning("Could not load gene annotations, using fallback plot")
        return _create_simple_fallback_plot(gene_group, subset)

    # Create figure with tight layout - disable constrained layout to avoid conflicts
    plt.rcParams["figure.constrained_layout.use"] = False

    # Get unique genes and their data
    unique_genes = list(sorted(subset["col4"].unique()))
    if len(unique_genes) == 0:
        logging.warning(f"[Fusion] No genes found for group {gene_group}, using fallback plot")
        return _create_simple_fallback_plot(gene_group, subset)
    
    # Ensure we have valid data to prevent zero-size axes
    if subset.empty:
        logging.warning(f"[Fusion] Empty subset for group {gene_group}, using fallback plot")
        return _create_simple_fallback_plot(gene_group, subset)

    # Create the unified side-by-side layout without ideograms
    # For 2 genes, we'll create 2 columns with 2 rows each (gene structure + read mapping)
    num_genes = len(unique_genes)
    logging.info(f"[Fusion] Creating {num_genes} gene layout with 2 rows (gene structure + read mapping)")
    # Increase figure size and add more padding to prevent tight layout warnings
    # Ensure minimum figure size to prevent zero-size axes warnings
    fig, axes = plt.subplots(2, num_genes, figsize=(max(20, 5 * num_genes), 10))

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

    # Use manual layout adjustment instead of tight_layout to avoid warnings
    # tight_layout often fails when there are many subplots or complex decorations
    plt.subplots_adjust(
        top=0.90,  # Leave space for suptitle
        bottom=0.1,
        left=0.08,
        right=0.95,
        hspace=0.4,  # More space between rows
        wspace=0.3,  # Space between columns
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
        logging.warning(f"Could not format ticks to megabases: {e}")
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
            logging.info(f"Loaded gene annotations: {len(gene_table)} entries")
            return gene_table
        else:
            logging.warning(f"Gene annotation file not found: {gene_data_path}")
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
    logging.info(f"[Fusion] add_fusion_section() called with sample_dir: {sample_dir}")
    
    if ui is None:
        logging.warning("[Fusion] ui is None, returning early")
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
            "status_container": None,
        },
    }

    def refresh_fusion() -> None:
        """Refresh fusion data."""
        try:
            logging.info(f"[Fusion] refresh_fusion() called for sample_dir: {sample_dir}")
            
            # Simple directory check
            if not sample_dir or not sample_dir.exists():
                logging.warning(f"[Fusion] Sample directory not found: {sample_dir}")
                return
            
            logging.info(f"[Fusion] Loading fusion data from: {sample_dir}")
            
            # Load fusion data directly (already in background)
            fusion_data = _load_fusion_data(sample_dir)
            
            logging.info(f"[Fusion] Loaded fusion data: {fusion_data}")
            
            # Update UI directly
            _update_fusion_ui(fusion_data, state)
                
        except Exception as e:
            logging.exception(f"[Fusion] Refresh failed: {e}")

    def _load_fusion_data(sample_dir: Path) -> Dict[str, Any]:
        """Load fusion data from files."""
        try:
            logging.info(f"[Fusion] _load_fusion_data() called with sample_dir: {sample_dir}")
            
            # Load target panel processed
            target_file = sample_dir / "fusion_candidates_master_processed.csv"
            genome_file = sample_dir / "fusion_candidates_all_processed.csv"
            
            logging.info(f"[Fusion] Target file: {target_file}")
            logging.info(f"[Fusion] Genome file: {genome_file}")

            t = _load_processed_pickle(target_file)
            g = _load_processed_pickle(genome_file)
            
            # Debug logging
            logging.info(f"[Fusion] Loaded target data: {t is not None}")
            logging.info(f"[Fusion] Loaded genome-wide data: {g is not None}")
            if g is not None:
                logging.info(f"[Fusion] Genome-wide candidate count: {g.get('candidate_count', 0)}")
                logging.info(f"[Fusion] Genome-wide gene groups: {len(g.get('gene_groups', []))}")
                logging.info(f"[Fusion] Genome-wide gene groups content: {g.get('gene_groups', [])}")
                logging.info(f"[Fusion] Genome-wide annotated_data shape: {g.get('annotated_data', pd.DataFrame()).shape}")
                logging.info(f"[Fusion] Genome-wide goodpairs shape: {g.get('goodpairs', pd.Series()).shape}")
                logging.info(f"[Fusion] Genome-wide gene_pairs count: {len(g.get('gene_pairs', []))}")
                logging.info(f"[Fusion] Genome-wide gene_pairs: {g.get('gene_pairs', [])[:10]}...")  # Show first 10
                
                # Use the same logic as reporting code - count gene_pairs instead of relying on candidate_count
                if g.get('gene_pairs') and len(g.get('gene_pairs', [])) > 0:
                    # Override candidate_count with the actual number of gene pairs (like reporting code does)
                    g['candidate_count'] = len(g.get('gene_pairs', []))
                    logging.info(f"[Fusion] Override genome-wide candidate_count to: {g['candidate_count']}")
                    
                    # Generate gene_groups from gene_pairs if missing (like reporting code does)
                    if not g.get('gene_groups') or len(g.get('gene_groups', [])) == 0:
                        # Convert gene_pairs to gene_groups format
                        gene_groups = []
                        for gene_pair in g.get('gene_pairs', []):
                            if isinstance(gene_pair, (tuple, list)) and len(gene_pair) >= 2:
                                gene_groups.append(list(gene_pair))
                        g['gene_groups'] = gene_groups
                        logging.info(f"[Fusion] Generated {len(gene_groups)} gene_groups from gene_pairs")
            else:
                logging.info(f"[Fusion] Genome-wide file exists: {genome_file.exists()}")
                if genome_file.exists():
                    logging.info(f"[Fusion] Genome-wide file size: {genome_file.stat().st_size} bytes")

            # Get file modification times
            target_mtime = target_file.stat().st_mtime if target_file.exists() else None
            genome_mtime = genome_file.stat().st_mtime if genome_file.exists() else None
            
            # Create data hashes
            target_data_hash = _create_data_hash(t) if t is not None else None
            genome_data_hash = _create_data_hash(g) if g is not None else None
            
            return {
                "target": {
                    "data": t,
                    "mtime": target_mtime,
                    "data_hash": target_data_hash
                },
                "genome": {
                    "data": g,
                    "mtime": genome_mtime,
                    "data_hash": genome_data_hash
                }
            }
        except Exception as e:
            logging.exception(f"[Fusion] Failed to load fusion data: {e}")
            return {"target": {"data": None, "mtime": None, "data_hash": None}, 
                   "genome": {"data": None, "mtime": None, "data_hash": None}}

    def _update_fusion_ui(fusion_data: Dict[str, Any], state: Dict[str, Any]) -> None:
        """Update fusion UI elements - runs on main UI thread."""
        try:
            target_data = fusion_data.get("target", {})
            genome_data = fusion_data.get("genome", {})
            
            t = target_data.get("data")
            g = genome_data.get("data")
            
            # Update target panel UI
            target_mtime = target_data.get("mtime")
            target_data_hash = target_data.get("data_hash")
            
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
                    t.get("annotated_data", pd.DataFrame()),
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
                                with_input=False,
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
                            ).classes("w-40")
                        with ui.row().classes("w-full"):
                            state["target"]["card"] = ui.card()
                            with state["target"]["card"]:
                                ui.label("Select gene pair to see results.").classes(
                                    "drop-shadow font-bold"
                                )
            
            # Only update visualization when data content actually changes (for background refreshes)
            elif t is not None and target_data_hash != state["target"].get("data_hash"):
                state["target"]["data_hash"] = target_data_hash
                # plot area
                try:
                    state["target"]["plot_container"].clear()
                except Exception:
                    pass
                if t.get("candidate_count", 0) > 0 and t.get("gene_groups"):
                    with state["target"]["plot_container"].classes("w-full"):
                        with ui.row().classes("w-full"):
                            ui.select(
                                options=t.get("gene_groups", []),
                                with_input=False,
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
                            ).classes("w-40")
                        with ui.row().classes("w-full"):
                            state["target"]["card"] = ui.card()
                            with state["target"]["card"]:
                                ui.label("Select gene pair to see results.").classes(
                                    "drop-shadow font-bold"
                                )

            # Update genome-wide UI
            genome_mtime = genome_data.get("mtime")
            genome_data_hash = genome_data.get("data_hash")
            logging.info(f"[Fusion] Genome-wide update check: g={g is not None}, mtime_changed={genome_mtime != state['genome'].get('mtime')}")
            if g is not None:
                logging.info(f"[Fusion] Genome-wide UI update: candidate_count={g.get('candidate_count', 0)}, gene_groups_count={len(g.get('gene_groups', []))}")
                logging.info(f"[Fusion] Genome-wide UI update: has_gene_groups={bool(g.get('gene_groups'))}")
                logging.info(f"[Fusion] Genome-wide UI update: will_show_dropdown={g.get('candidate_count', 0) > 0 and bool(g.get('gene_groups'))}")
            
            # Always update summary and table when file changes
            if g is not None and genome_mtime != state["genome"].get("mtime"):
                state["genome"]["data"] = g
                state["genome"]["mtime"] = genome_mtime
                state["genome"]["data_hash"] = genome_data_hash
                
                # summary
                try:
                    if "summary" in state and state["summary"].get("genome_lbl"):
                        state["summary"][
                            "genome_lbl"
                        ].text = f"Genome-wide fusion groups: {int(g.get('candidate_count', 0))}"
                except Exception:
                    pass
                # table
                try:
                    state["genome"]["table_container"].clear()
                except Exception:
                    pass
                state["genome"]["table"] = _make_fusion_table(
                    state["genome"]["table_container"],
                    g.get("annotated_data", pd.DataFrame()),
                )
                
                # Create visualization on first load or when data changes
                try:
                    state["genome"]["plot_container"].clear()
                    state["genome"]["status_container"].clear()
                except Exception:
                    pass
                logging.info(f"[Fusion] Genome-wide plotting check: candidate_count={g.get('candidate_count', 0)}, gene_groups={len(g.get('gene_groups', []))}")
                # More permissive condition - show dropdown if we have any gene groups, even if candidate_count is 0
                if g.get("gene_groups") and len(g.get("gene_groups", [])) > 0:
                    with state["genome"]["plot_container"].classes("w-full"):
                        with ui.row().classes("w-full"):
                            ui.select(
                                options=g.get("gene_groups", []),
                                with_input=False,
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
                            ).classes("w-40")
                        with ui.row().classes("w-full"):
                            state["genome"]["card"] = ui.card()
                            with state["genome"]["card"]:
                                ui.label("Select gene pair to see results.").classes(
                                    "drop-shadow font-bold"
                                )
                else:
                    # Show status message when genome-wide data is not available
                    with state["genome"]["status_container"].classes("w-full"):
                        ui.label("Genome-wide fusion analysis not available").classes("text-gray-600 text-sm")
                        ui.label("(Requires supplementary reads in BAM file)").classes("text-gray-500 text-xs")
            
            # Only update visualization when data content actually changes (for background refreshes)
            elif g is not None and genome_data_hash != state["genome"].get("data_hash"):
                state["genome"]["data_hash"] = genome_data_hash
                # plot area
                try:
                    state["genome"]["plot_container"].clear()
                    state["genome"]["status_container"].clear()
                except Exception:
                    pass
                logging.info(f"[Fusion] Genome-wide plotting check: candidate_count={g.get('candidate_count', 0)}, gene_groups={len(g.get('gene_groups', []))}")
                # More permissive condition - show dropdown if we have any gene groups, even if candidate_count is 0
                if g.get("gene_groups") and len(g.get("gene_groups", [])) > 0:
                    with state["genome"]["plot_container"].classes("w-full"):
                        with ui.row().classes("w-full"):
                            ui.select(
                                options=g.get("gene_groups", []),
                                with_input=False,
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
                            ).classes("w-40")
                        with ui.row().classes("w-full"):
                            state["genome"]["card"] = ui.card()
                            with state["genome"]["card"]:
                                ui.label("Select gene pair to see results.").classes(
                                    "drop-shadow font-bold"
                                )
                else:
                    # Show status message when genome-wide data is not available
                    with state["genome"]["status_container"].classes("w-full"):
                        ui.label("Genome-wide fusion analysis not available").classes("text-gray-600 text-sm")
                        ui.label("(Requires supplementary reads in BAM file)").classes("text-gray-500 text-xs")
        except Exception as e:
            logging.exception(f"[Fusion] Refresh failed: {e}")

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
        state["genome"]["status_container"] = ui.column().classes("w-full mt-2")

    # Initial refresh and timer
    try:
        # Prime summary values on first load
        try:
            t0 = _load_processed_pickle(
                (sample_dir / "fusion_candidates_master_processed.csv")
                if sample_dir
                else Path("/dev/null")
            )
            g0 = _load_processed_pickle(
                (sample_dir / "fusion_candidates_all_processed.csv")
                if sample_dir
                else Path("/dev/null")
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
        
        # Start the refresh timer (every 30 seconds)
        logging.info("[Fusion] Setting up refresh timer with immediate=True")
        ui.timer(30.0, refresh_fusion, active=True, immediate=True)
        logging.info("[Fusion] Timer set up successfully")
    except Exception as e:
        logging.exception(f"[Fusion] Exception in timer setup: {e}")
