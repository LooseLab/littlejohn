from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List
import logging

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
        # Persistent Matplotlib element (NiceGUI integration)
        mgmt_mpl = ui.matplotlib(figsize=(4, 3)).classes("mx-auto flex justify-center")
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
                {"name": "meth_fwd", "label": "Forward %", "field": "meth_fwd"},
                {"name": "meth_rev", "label": "Reverse %", "field": "meth_rev"},
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
            if df.shape[1] < 10:
                return []
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
            df = df.iloc[:, : len(cols)]
            df.columns = cols
            cov_split = df["Coverage_Info"].astype(str).str.split()
            df["Coverage"] = cov_split.str[0].astype(float)
            df["Modified_Fraction"] = cov_split.str[1].astype(float)
            cpg_pairs = [
                (129467255, 129467256),
                (129467258, 129467259),
                (129467262, 129467263),
                (129467272, 129467273),
            ]
            rows: List[Dict[str, Any]] = []
            for p1, p2 in cpg_pairs:
                fwd = df[
                    (df["Chromosome"] == "chr10")
                    & (df["Start"] == p1 - 1)
                    & (df["Strand"] == "+")
                ]
                rev = df[
                    (df["Chromosome"] == "chr10")
                    & (df["Start"] == p2 - 1)
                    & (df["Strand"] == "-")
                ]
                if not fwd.empty and not rev.empty:
                    cov_f = float(fwd["Coverage"].iloc[0])
                    cov_r = float(rev["Coverage"].iloc[0])
                    tot = cov_f + cov_r
                    mf = float(fwd["Modified_Fraction"].iloc[0])
                    mr = float(rev["Modified_Fraction"].iloc[0])
                    weighted = ((cov_f * mf) + (cov_r * mr)) / tot if tot > 0 else 0.0
                    label_map = {
                        "129467255/129467256": "Site 1",
                        "129467258/129467259": "Site 2",
                        "129467262/129467263": "Site 3",
                        "129467272/129467273": "Site 4",
                    }
                    pos_key = f"{p1}/{p2}"
                    rows.append(
                        {
                            "site": f"{label_map.get(pos_key,'Unknown')} (CpG {pos_key})",
                            "chr": "chr10",
                            "pos": pos_key,
                            "cov_fwd": int(cov_f),
                            "cov_rev": int(cov_r),
                            "cov_total": int(tot),
                            "meth": round(weighted * 100.0, 2),
                            "meth_fwd": round(mf * 100.0, 2),
                            "meth_rev": round(mr * 100.0, 2),
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
            if final_files:
                latest_csv = final_files[0]
            else:
                # Fallback to numeric-prefixed files
                csv_files = list(sample_dir.glob("*_mgmt.csv"))
                if not csv_files:
                    return

                latest_csv = max(csv_files, key=_count_from_name)
            key = str(sample_dir)
            state = launcher._mgmt_state.get(key, {})
            current_count = _count_from_name(latest_csv)
            csv_mtime = latest_csv.stat().st_mtime
            if (
                state.get("csv_path") == str(latest_csv)
                and state.get("csv_mtime") == csv_mtime
            ):
                return
            df = pd.read_csv(latest_csv)
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
            
            """
            img_path = sample_dir / f"{current_count}_mgmt.png"
            if img_path.exists():
                img_mtime = img_path.stat().st_mtime
                if state.get('png_path') != str(img_path) or state.get('png_mtime') != img_mtime:
                    try:
                        mgmt_img.set_source(str(img_path))
                    except Exception:
                        mgmt_img.source = str(img_path)
            """
            # Gather site_rows for table and plot annotations
            site_rows: List[Dict[str, Any]] = []
            bed_path = sample_dir / f"{current_count}_mgmt.bed"
            if not bed_path.exists():
                alt_bed = sample_dir / f"{current_count}_mgmt_mgmt.bed"
                bed_path = alt_bed if alt_bed.exists() else bed_path
            if bed_path.exists():
                bed_mtime = bed_path.stat().st_mtime
                if (
                    state.get("bed_path") != str(bed_path)
                    or state.get("bed_mtime") != bed_mtime
                ):
                    site_rows = _extract_mgmt_specific_sites(bed_path)
                    try:
                        mgmt_sites_table.rows = site_rows
                        mgmt_sites_table.update()
                    except Exception:
                        pass
                else:
                    # still populate for plotting if unchanged
                    site_rows = _extract_mgmt_specific_sites(bed_path)

            bam_path = sample_dir / "mgmt_sorted.bam"
            if bam_path.exists():
                try:
                    from robin.analysis.methylation_wrapper import (
                        locus_figure,
                    )

                    fig = locus_figure(
                        interval="chr10:129466536-129467536",
                        bam_path=str(bam_path),
                        motif="CG",
                        mods="m",
                        extra_cli=[
                            # "--minqual","20", "--reads","2000", "--height","6", "--width","10"
                        ],
                        site_rows=site_rows,
                    )
                    # Update NiceGUI Matplotlib element
                    mgmt_mpl.figure = fig
                    mgmt_mpl.update()
                except Exception as e:
                    # If methylartist fails, create a simple placeholder plot
                    import matplotlib.pyplot as plt
                    import matplotlib.patches as patches
                    
                    fig, ax = plt.subplots(figsize=(4, 3))
                    ax.text(0.5, 0.5, f"Methylation plot unavailable\n({str(e)})", 
                           ha='center', va='center', transform=ax.transAxes,
                           fontsize=10, color='red')
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

    # Start the refresh timer (every 15 seconds)
    ui.timer(15.0, _refresh_mgmt, active=True, immediate=True)
