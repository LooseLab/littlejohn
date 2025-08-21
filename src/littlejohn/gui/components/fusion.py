from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path
import os
import logging
import pickle

import pandas as pd

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


def _load_processed_pickle(file_path: Path) -> Optional[Dict[str, Any]]:
    """Load the preprocessed fusion data structure written by preprocess_fusion_data_standalone.

    Note: Although the extension is .csv, the file is a pickle per current pipeline.
    """
    try:
        if not file_path.exists() or file_path.stat().st_size == 0:
            return None
        with open(file_path, "rb") as f:
            data = pickle.load(f)
        # Expected keys: annotated_data (DataFrame), goodpairs (Series), gene_groups (list), candidate_count (int)
        return data if isinstance(data, dict) else None
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
    """Simple plot for selected gene group using matplotlib via ui.pyplot.

    Shows read mapping start (mS) along x-axis and gene labels on y-axis, colored by the 'Color' assigned during preprocessing.
    """
    try:
        subset = annotated_data[goodpairs]
        subset = subset[subset["col4"].isin(gene_group)]
        if subset.empty:
            with container:
                ui.label("No reads for selected gene group").classes("text-gray-600")
            return

        # Use matplotlib for a quick scatter style plot
        with container:
            with ui.pyplot(figsize=(10, 4)).classes("w-full"):
                import matplotlib.pyplot as plt  # local import to avoid global dependency during module import

                plt.rcParams["figure.constrained_layout.use"] = True
                fig, ax = plt.subplots(1, 1, figsize=(10, 4))

                # Map gene names to y positions
                genes = list(sorted(subset["col4"].unique()))
                y_map = {g: i for i, g in enumerate(genes)}
                x = subset["reference_start"].astype(int)
                y = subset["col4"].map(y_map).astype(int)
                colors = subset.get("Color", pd.Series(["#1f77b4"] * len(subset)))
                ax.scatter(x, y, c=colors, s=15, alpha=0.8)
                ax.set_yticks(list(y_map.values()))
                ax.set_yticklabels(list(y_map.keys()))
                ax.set_xlabel("Mapping Start (reference_start)")
                ax.set_title(
                    f"Reads supporting fusion group: {', '.join(gene_group)}",
                    fontsize=12,
                    fontweight="bold",
                )
                ax.grid(True, axis="x", linestyle=":", alpha=0.4)
    except Exception as e:
        logging.exception(f"[Fusion] Failed to plot gene group {gene_group}: {e}")


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
            "table_container": None,
            "table": None,
            "plot_container": None,
            "card": None,
        },
        "genome": {
            "data": None,
            "mtime": None,
            "table_container": None,
            "table": None,
            "plot_container": None,
            "card": None,
        },
    }

    def refresh() -> None:
        try:
            if sample_dir is None or not sample_dir.exists():
                return
            # Load target panel processed
            target_file = sample_dir / "fusion_candidates_master_processed.csv"
            genome_file = sample_dir / "fusion_candidates_all_processed.csv"

            t = _load_processed_pickle(target_file)
            g = _load_processed_pickle(genome_file)

            # Update target panel UI
            target_mtime = target_file.stat().st_mtime if target_file.exists() else None
            if t is not None and target_mtime != state["target"].get("mtime"):
                state["target"]["data"] = t
                state["target"]["mtime"] = target_mtime
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
                            ).classes("w-40")
                        with ui.row().classes("w-full"):
                            state["target"]["card"] = ui.card()
                            with state["target"]["card"]:
                                ui.label("Select gene pair to see results.").tailwind(
                                    "drop-shadow", "font-bold"
                                )

            # Update genome-wide UI
            genome_mtime = genome_file.stat().st_mtime if genome_file.exists() else None
            if g is not None and genome_mtime != state["genome"].get("mtime"):
                state["genome"]["data"] = g
                state["genome"]["mtime"] = genome_mtime
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
                # plot area
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
                            ).classes("w-40")
                        with ui.row().classes("w-full"):
                            state["genome"]["card"] = ui.card()
                            with state["genome"]["card"]:
                                ui.label("Select gene pair to see results.").tailwind(
                                    "drop-shadow", "font-bold"
                                )
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
        state["target"]["table_container"] = ui.column().classes("w-full")
        state["target"]["plot_container"] = ui.column().classes("w-full mt-2")

        # Genome-wide
        ui.separator()
        ui.label("Genome-wide").classes("text-base font-medium mt-2 mb-1")
        state["genome"]["table_container"] = ui.column().classes("w-full")
        state["genome"]["plot_container"] = ui.column().classes("w-full mt-2")

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
        refresh()
        ui.timer(10.0, refresh)
    except Exception:
        pass
