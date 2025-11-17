from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List
import logging
import time

import pandas as pd


try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

from robin.gui.theme import styled_table


def add_mgmt_section(launcher: Any, sample_dir: Path) -> None:
    """Build the MGMT section and attach refresh timer. Uses launcher._mgmt_state."""
    with ui.card().classes("w-full"):
        with ui.card().classes(
            "w-full bg-gradient-to-r from-blue-50 to-indigo-50 mb-2 p-2"
        ):
            with ui.row().classes("gap-6 items-center"):
                mgmt_status = ui.label("MGMT status: Unknown").classes(
                    "text-sm font-semibold text-blue-800"
                )
                mgmt_avg = ui.label("Average: --%").classes("text-sm text-gray-700")
                mgmt_pred = ui.label("Prediction: --%").classes("text-sm text-gray-700")
        ui.label("MGMT methylation plot").classes("text-sm text-gray-700")
        # mgmt_img = ui.image('')
        # Matplotlib element for displaying the methylation plot
        # Using ui.matplotlib() as it works with pre-existing figure objects
        # Width is 2x relative to height (24:12 = 2:1 aspect ratio)
        mgmt_mpl = ui.matplotlib(figsize=(24, 12)).classes("w-full")
        ui.separator()
        ui.label("MGMT results (latest)").classes("text-sm text-gray-700 mt-2")
        _, mgmt_results_table = styled_table(
            columns=[
                {"name": "average", "label": "Average %", "field": "average"},
                {"name": "pred", "label": "Prediction %", "field": "pred"},
                {"name": "status", "label": "Status", "field": "status"},
            ],
            rows=[],
            pagination=0,
            class_size="table-xs",
        )
        ui.separator()
        ui.label("MGMT CpG Site Methylation Data").classes("text-sm text-gray-700 mt-2")
        _, mgmt_sites_table = styled_table(
            columns=[
                {"name": "site", "label": "Site", "field": "site"},
                {"name": "chr", "label": "Chr", "field": "chr"},
                {"name": "pos", "label": "CpG Position", "field": "pos"},
                {"name": "cov_fwd", "label": "Forward Cov", "field": "cov_fwd"},
                {"name": "cov_rev", "label": "Reverse Cov", "field": "cov_rev"},
                {"name": "cov_total", "label": "Total Cov", "field": "cov_total"},
                {"name": "meth", "label": "% Methylation", "field": "meth"},
                {"name": "meth_fwd", "label": "Forward Methylated", "field": "meth_fwd"},
                {"name": "meth_rev", "label": "Reverse Methylated", "field": "meth_rev"},
                {"name": "notes", "label": "Notes", "field": "notes"},
            ],
            rows=[],
            pagination=0,
            class_size="table-xs",
        )

    def _extract_mgmt_specific_sites(bed_path: Path) -> List[Dict[str, Any]]:
        try:
            import pandas as _pd

            df = _pd.read_csv(bed_path, sep="\t", header=None)
            
            # Check if column 10 contains space-separated values (old format)
            # Even if file has 18 columns, column 10 might still be space-separated
            has_space_separated_col10 = False
            if df.shape[1] > 10 and len(df) > 0:
                # Check if column 10 (index 9) contains space-separated values
                sample_val = str(df.iloc[0, 9])
                has_space_separated_col10 = ' ' in sample_val or '\t' in sample_val
            
            # Check if this is the new bedmethyl format (separate columns) or old format (space-separated column 10)
            if df.shape[1] >= 12 and not has_space_separated_col10:
                # New bedmethyl format with separate columns
                cols = [
                    "Chromosome",
                    "Start",
                    "End",
                    "Modified_Base_Code",
                    "Score",
                    "Strand",
                    "Start2",
                    "End2",
                    "RGB",
                    "Nvalid_cov",      # Column 10: Valid coverage (absolute count)
                    "Fraction_Modified",  # Column 11: Nmod / Nvalid_cov (fraction 0-1, not percentage)
                    "Nmod",            # Column 12: Absolute count of modified reads
                ]
                # Read at least the first 12 columns
                num_cols_to_read = min(len(cols), df.shape[1])
                df = df.iloc[:, :num_cols_to_read]
                df.columns = cols[:num_cols_to_read]
                
                # Convert to proper types
                df["Nvalid_cov"] = df["Nvalid_cov"].astype(float)
                df["Fraction_Modified"] = df["Fraction_Modified"].astype(float)
                if "Nmod" in df.columns:
                    df["Nmod"] = df["Nmod"].astype(float)
                else:
                    # If Nmod column doesn't exist, calculate it from fraction * coverage
                    df["Nmod"] = df["Nvalid_cov"] * df["Fraction_Modified"]
                
                # Ensure Start is integer for proper comparison
                df["Start"] = df["Start"].astype(int)
                
                # For backward compatibility, keep Coverage and Modified_Fraction columns
                df["Coverage"] = df["Nvalid_cov"]
                df["Modified_Fraction"] = df["Fraction_Modified"] * 100.0
            elif df.shape[1] >= 10:
                # Old format with space-separated Coverage_Info in column 10
                cols = [
                    "Chromosome",
                    "Start",
                    "End",
                    "Name",
                    "Score",
                    "Strand",
                    "Start2",
                    "End2",
                    "RGB",
                    "Coverage_Info",
                ]
                # Only use first 10 columns for old format, even if file has more columns
                df = df.iloc[:, : len(cols)]
                df.columns = cols
                
                # Parse Coverage_Info (space-separated: coverage fraction/percentage)
                # Format: "coverage fraction" or "coverage percentage"
                cov_split = df["Coverage_Info"].astype(str).str.split()
                df["Coverage"] = cov_split.str[0].astype(float)
                
                # Get second value (fraction or percentage)
                # pandas str accessor will return NaN for missing values
                fraction_val = cov_split.str[1].astype(float).fillna(0.0)
                
                # Determine if second value is fraction (0-1) or percentage (0-100)
                # If any value is > 1, assume it's percentage, otherwise assume fraction
                is_percentage = (fraction_val > 1.0).any() if len(fraction_val) > 0 else False
                
                if is_percentage:
                    df["Modified_Fraction"] = fraction_val
                    df["Fraction_Modified"] = df["Modified_Fraction"] / 100.0  # Convert percentage to fraction
                else:
                    df["Fraction_Modified"] = fraction_val
                    df["Modified_Fraction"] = df["Fraction_Modified"] * 100.0  # Convert fraction to percentage
                
                # Convert to new format columns for consistency
                df["Nvalid_cov"] = df["Coverage"]
                df["Nmod"] = df["Coverage"] * df["Fraction_Modified"]
                
                # Ensure Start is integer for proper comparison
                df["Start"] = df["Start"].astype(int)
            else:
                return []
            cpg_pairs = [
                (129467255, 129467256),
                (129467258, 129467259),
                (129467262, 129467263),
                (129467272, 129467273),
            ]
            rows: List[Dict[str, Any]] = []
            label_map = {
                "129467255/129467256": "1",
                "129467258/129467259": "2",
                "129467262/129467263": "3",
                "129467272/129467273": "4",
            }
            
            for p1, p2 in cpg_pairs:
                pos_key = f"{p1}/{p2}"
                site_label = label_map.get(pos_key, "Unknown")
                
                # For a CpG pair (p1, p2) where p1 and p2 are consecutive:
                # The CpG site consists of two cytosines:
                # - Forward strand: C at p1, G at p1+1 (p2)
                # - Reverse strand: C at p2 (reverse complement of G at p1+1), G at p1 (reverse complement of C at p1)
                #
                # In bedmethyl format, methylation is reported at the reference position of the C.
                # For this CpG site, we need to check:
                # - Forward strand reads: methylation at position p1 (the C on forward strand)
                # - Reverse strand reads: methylation at position p2 (the C on reverse strand)
                #
                # IMPORTANT: We must ensure we're checking the correct positions for this specific CpG site
                # and not accidentally assigning methylation from adjacent sites.
                
                # Check forward strand reads at position p1 (the C on forward strand for this CpG)
                fwd_p1 = df[
                    (df["Chromosome"] == "chr10")
                    & (df["Start"] == p1 - 1)
                    & (df["Strand"] == "+")
                ]
                
                # Check reverse strand reads at position p2 (the C on reverse strand for this CpG)
                # This is correct because reverse strand reads see the C at p2 for this CpG site
                rev_p2 = df[
                    (df["Chromosome"] == "chr10")
                    & (df["Start"] == p2 - 1)
                    & (df["Strand"] == "-")
                ]
                
                # IMPORTANT: We should NOT check reverse strand at p1 or forward strand at p2,
                # as those would represent methylation from adjacent CpG sites or the wrong cytosine.
                # For example, reverse strand at p1 would be the G position (not a C), and
                # forward strand at p2 would be the G position (not a C) for this CpG site.
                
                # Get forward strand data from p1
                if not fwd_p1.empty:
                    cov_f = float(fwd_p1["Nvalid_cov"].iloc[0])
                    mf = float(fwd_p1["Fraction_Modified"].iloc[0])  # Fraction (0-1), not percentage
                    nmod_f = float(fwd_p1["Nmod"].iloc[0])  # Direct count from bedmethyl file
                    meth_fwd_count = int(round(nmod_f))  # Use Nmod directly from bedmethyl
                else:
                    cov_f = 0.0
                    mf = 0.0
                    meth_fwd_count = 0
                
                # Get reverse strand data from p2 (the C position on reverse strand for this CpG)
                if not rev_p2.empty:
                    cov_r = float(rev_p2["Nvalid_cov"].iloc[0])
                    mr = float(rev_p2["Fraction_Modified"].iloc[0])  # Fraction (0-1), not percentage
                    nmod_r = float(rev_p2["Nmod"].iloc[0])  # Direct count from bedmethyl file
                    meth_rev_count = int(round(nmod_r))  # Use Nmod directly from bedmethyl
                else:
                    cov_r = 0.0
                    mr = 0.0
                    meth_rev_count = 0
                
                # Only add row if we have data
                if cov_f > 0 or cov_r > 0:
                    tot = cov_f + cov_r
                    # Calculate weighted average methylation fraction (0-1), then convert to percentage
                    # weighted = ((cov_f * mf) + (cov_r * mr)) / tot if tot > 0 else 0.0
                    # Or equivalently: weighted = (nmod_f + nmod_r) / tot
                    weighted = ((cov_f * mf) + (cov_r * mr)) / tot if tot > 0 else 0.0
                    weighted_pct = weighted * 100.0  # Convert to percentage for display
                    
                    rows.append(
                        {
                            "site": f"{site_label} (CpG {pos_key})",
                            "chr": "chr10",
                            "pos": pos_key,
                            "cov_fwd": int(cov_f),
                            "cov_rev": int(cov_r),
                            "cov_total": int(tot),
                            "meth": round(weighted_pct, 2),  # Store as percentage for display
                            "meth_fwd": int(meth_fwd_count),  # Direct from Nmod column
                            "meth_rev": int(meth_rev_count),  # Direct from Nmod column
                            "notes": "Combined methylation from both strands of CpG pair",
                        }
                    )
            
            return rows
        except Exception:
            return []

    def _refresh_mgmt() -> None:
        """Refresh MGMT data."""
        try:
            # Simple directory check
            if not sample_dir or not sample_dir.exists():
                logging.warning(f"[MGMT] Sample directory not found: {sample_dir}")
                return
            
            # Run MGMT refresh directly (already in background)
            _refresh_mgmt_sync(sample_dir, launcher)
                
        except Exception as e:
            logging.exception(f"[MGMT] Refresh failed: {e}")

    def _refresh_mgmt_sync(sample_dir: Path, launcher: Any) -> None:
        """Synchronous MGMT refresh."""
        try:
            # Define the helper function first
            def _count_from_name(p: Path) -> int:
                try:
                    return int(p.name.split("_")[0])
                except Exception:
                    return -1

            # First check for final_mgmt.csv files (highest priority)
            final_files = list(sample_dir.glob("final_mgmt.csv"))
            is_final_file = False
            if final_files:
                latest_csv = final_files[0]
                is_final_file = True
            else:
                # Fallback to numeric-prefixed files
                csv_files = list(sample_dir.glob("*_mgmt.csv"))
                if not csv_files:
                    return

                latest_csv = max(csv_files, key=_count_from_name)
            key = str(sample_dir)
            state = launcher._mgmt_state.get(key, {})
            
            # Check if this is a fresh page visit - force updates on fresh visits
            is_fresh_visit = "last_visit_time" not in state
            if is_fresh_visit:
                state["last_visit_time"] = time.time()
            
            current_count = _count_from_name(latest_csv)
            csv_mtime = latest_csv.stat().st_mtime
            # Force update on fresh page visit or if file changed
            csv_needs_update = (
                is_fresh_visit
                or state.get("csv_path") != str(latest_csv)
                or state.get("csv_mtime") != csv_mtime
            )
            
            # Always read CSV data (needed for plot even if CSV hasn't changed)
            try:
                df = pd.read_csv(latest_csv)
            except Exception as e:
                logging.error(f"[MGMT] Failed to read CSV file {latest_csv}: {e}")
                return  # Cannot proceed without CSV data
            
            # Only update CSV-related UI elements if CSV has changed
            if csv_needs_update:
                try:
                    status = str(df.get("status", pd.Series(["Unknown"])).iloc[0])
                    average = float(df.get("average", pd.Series([0.0])).iloc[0])
                    pred = float(df.get("pred", pd.Series([0.0])).iloc[0])
                    mgmt_status.set_text(f"MGMT status: {status}")
                    mgmt_avg.set_text(f"Average: {average:.2f}%")
                    mgmt_pred.set_text(f"Prediction: {pred:.2f}%")
                    try:
                        mgmt_results_table.rows = [
                            {
                                "average": f"{average:.2f}",
                                "pred": f"{pred:.2f}",
                                "status": status,
                            }
                        ]
                        mgmt_results_table.update()
                    except Exception:
                        pass
                except Exception as e:
                    logging.error(f"[MGMT] Failed to update CSV UI elements: {e}")
                    pass
            
            # Gather site_rows for table and plot annotations
            site_rows: List[Dict[str, Any]] = []
            # Handle bed file lookup based on CSV file type
            if is_final_file:
                # For final_mgmt.csv, look for final_mgmt.bed
                bed_path = sample_dir / "final_mgmt.bed"
                if not bed_path.exists():
                    alt_bed = sample_dir / "final_mgmt_mgmt.bed"
                    bed_path = alt_bed if alt_bed.exists() else bed_path
            else:
                # For numeric-prefixed CSV files, use the count
                bed_path = sample_dir / f"{current_count}_mgmt.bed"
                if not bed_path.exists():
                    alt_bed = sample_dir / f"{current_count}_mgmt_mgmt.bed"
                    bed_path = alt_bed if alt_bed.exists() else bed_path
            
            if bed_path.exists():
                bed_mtime = bed_path.stat().st_mtime
                bed_needs_update = (
                    is_fresh_visit
                    or state.get("bed_path") != str(bed_path)
                    or state.get("bed_mtime") != bed_mtime
                )
                
                # Always extract site_rows for plotting (needed even if bed hasn't changed)
                site_rows = _extract_mgmt_specific_sites(bed_path)
                
                # Only update table if bed file has changed
                if bed_needs_update:
                    try:
                        mgmt_sites_table.rows = site_rows
                        mgmt_sites_table.update()
                    except Exception:
                        pass
                else:
                    # Still update table even if unchanged (for consistency)
                    try:
                        mgmt_sites_table.rows = site_rows
                        mgmt_sites_table.update()
                    except Exception:
                        pass

            bam_path = sample_dir / "mgmt_sorted.bam"
            if bam_path.exists():
                try:
                    from robin.analysis.methylation_wrapper import (
                        locus_figure,
                    )

                    logging.debug(f"[MGMT] Creating locus figure for {bam_path}")
                    fig = locus_figure(
                        interval="chr10:129466536-129467536",
                        bam_path=str(bam_path),
                        motif="CG",
                        mods="m",
                        extra_cli=[
                            # Set figure size to match our UI element (width 2x height)
                            "--width", "18",
                            "--height", "8",
                            # "--minqual","20", "--reads","2000"
                        ],
                        site_rows=site_rows,
                    )
                    logging.debug(f"[MGMT] Locus figure created successfully, figure number: {fig.number}")
                    
                    # Update matplotlib element with the figure
                    # Suppress GridSpec warnings when updating figure
                    import warnings
                    import matplotlib.pyplot as plt
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", UserWarning)
                        # Close previous figure if it exists before assigning new one
                        if hasattr(mgmt_mpl, 'figure') and mgmt_mpl.figure is not None:
                            try:
                                plt.close(mgmt_mpl.figure)
                            except Exception:
                                pass
                        logging.debug(f"[MGMT] Assigning figure to matplotlib element")
                        mgmt_mpl.figure = fig
                        mgmt_mpl.update()
                        logging.debug(f"[MGMT] Plot updated successfully")
                except Exception as e:
                    # If methylartist fails, create a simple placeholder plot
                    logging.exception(f"[MGMT] Failed to create methylation plot: {e}")
                    import matplotlib.pyplot as plt
                    import matplotlib.patches as patches
                    
                    # Close previous figure if it exists
                    if hasattr(mgmt_mpl, 'figure') and mgmt_mpl.figure is not None:
                        try:
                            plt.close(mgmt_mpl.figure)
                        except Exception:
                            pass
                    
                    fig, ax = plt.subplots(figsize=(24, 12))
                    ax.text(0.5, 0.5, f"Methylation plot unavailable\n({str(e)})", 
                           ha='center', va='center', transform=ax.transAxes,
                           fontsize=10, color='red')
                    ax.set_xlim(0, 1)
                    ax.set_ylim(0, 1)
                    ax.axis('off')
                    ax.set_title("MGMT Methylation Plot")
                    
                    mgmt_mpl.figure = fig
                    mgmt_mpl.update()
            else:
                logging.warning(f"[MGMT] BAM file not found: {bam_path}")
                # Create a placeholder plot indicating BAM file is missing
                import matplotlib.pyplot as plt
                if hasattr(mgmt_mpl, 'figure') and mgmt_mpl.figure is not None:
                    try:
                        plt.close(mgmt_mpl.figure)
                    except Exception:
                        pass
                fig, ax = plt.subplots(figsize=(24, 12))
                ax.text(0.5, 0.5, "MGMT BAM file not found\n(mgmt_sorted.bam)", 
                       ha='center', va='center', transform=ax.transAxes,
                       fontsize=10, color='orange')
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.axis('off')
                ax.set_title("MGMT Methylation Plot")
                mgmt_mpl.figure = fig
                mgmt_mpl.update()

            launcher._mgmt_state[key] = {
                "last_count": current_count,
                "csv_path": str(latest_csv),
                "csv_mtime": csv_mtime,
                #'png_path': str(img_path) if img_path.exists() else '', 'png_mtime': img_path.stat().st_mtime if img_path.exists() else 0,
                "bed_path": str(bed_path) if bed_path.exists() else "",
                "bed_mtime": bed_path.stat().st_mtime if bed_path.exists() else 0,
            }
        except Exception as e:
            raise Exception(f"Failed to refresh MGMT section: {e}")

    # Start the refresh timer (every 30 seconds)
    ui.timer(30.0, _refresh_mgmt, active=True, immediate=True)
