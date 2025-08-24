from __future__ import annotations

from pathlib import Path
from typing import Any

import natsort
import numpy as np
import pandas as pd
import logging

try:
    from nicegui import ui, app
except ImportError:  # pragma: no cover
    ui = None
    app = None

from littlejohn.gui.theme import styled_table


def add_coverage_section(launcher: Any, sample_dir: Path) -> None:
    """Build the Coverage UI section and attach refresh timers.

    This mirrors the existing inline implementation but lives in a reusable module.
    """
    with ui.card().classes("w-full"):
        ui.label("🧫 Coverage").classes("text-lg font-semibold mb-2")
        # Summary row (quality + metrics)
        with ui.row().classes("w-full items-center justify-between mb-2"):
            with ui.column().classes("gap-1"):
                ui.label("Coverage Analysis").classes("text-sm font-medium")
                with ui.row().classes("items-center gap-2"):
                    cov_quality_label = ui.label("Quality: --").classes(
                        "text-gray-600 font-medium"
                    )
                    cov_quality_badge = ui.label("--x").classes(
                        "px-2 py-1 rounded bg-gray-100 text-gray-600"
                    )
            with ui.column().classes("gap-1 items-end"):
                cov_global_lbl = ui.label("Global Estimated Coverage: --x").classes(
                    "text-sm text-gray-600"
                )
                cov_target_lbl = ui.label("Targets Estimated Coverage: --x").classes(
                    "text-sm text-gray-600"
                )
                cov_enrich_lbl = ui.label("Estimated enrichment: --x").classes(
                    "text-sm text-gray-600"
                )
        with ui.card().classes("w-full"):
            with ui.grid(columns=2).classes("w-full gap-4"):
                # Per Chromosome Coverage (bar)
                echart_chr_cov = ui.echart(
                    {
                        "backgroundColor": "transparent",
                        "title": {
                            "text": "Per Chromosome Coverage",
                            "left": "center",
                            "top": 10,
                            "textStyle": {"fontSize": 16, "color": "#000"},
                        },
                        "tooltip": {"trigger": "axis"},
                        "grid": {
                            "left": "5%",
                            "right": "5%",
                            "bottom": "10%",
                            "top": "20%",
                            "containLabel": True,
                        },
                        "xAxis": {
                            "type": "category",
                            "data": [],
                            "axisLabel": {"rotate": 45},
                        },
                        "yAxis": {"type": "value", "name": "Coverage (x)"},
                        "series": [],
                    }
                ).classes("w-full h-64")
                # Per Chromosome Target Coverage (scatter)
                echart_target_cov = ui.echart(
                    {
                        "backgroundColor": "transparent",
                        "title": {
                            "text": "Per Chromosome Target Coverage",
                            "left": "center",
                            "top": 10,
                            "textStyle": {"fontSize": 16, "color": "#000"},
                        },
                        "legend": {"data": ["Off Target", "On Target"], "top": 45},
                        "tooltip": {"trigger": "axis"},
                        "grid": {
                            "left": "5%",
                            "right": "5%",
                            "bottom": "10%",
                            "top": "25%",
                            "containLabel": True,
                        },
                        "xAxis": {
                            "type": "category",
                            "data": [],
                            "axisLabel": {"rotate": 45},
                        },
                        "yAxis": {"type": "value", "name": "Coverage (x)"},
                        "series": [],
                    }
                ).classes("w-full h-64")

        with ui.card().classes("w-full"):
            ui.label("📈 Coverage Over Time").classes("text-lg font-semibold mb-2")
            echart_time = ui.echart(
                {
                    "backgroundColor": "transparent",
                    "title": {
                        "text": "Coverage Over Time",
                        "left": "center",
                        "top": 10,
                        "textStyle": {"fontSize": 16, "color": "#000"},
                    },
                    "tooltip": {"trigger": "axis"},
                    "grid": {
                        "left": "5%",
                        "right": "5%",
                        "bottom": "10%",
                        "top": "20%",
                        "containLabel": True,
                    },
                    "xAxis": {"type": "time"},
                    "yAxis": {"type": "value", "name": "Coverage (x)"},
                    "series": [{"type": "line", "smooth": True, "data": []}],
                }
            ).classes("w-full h-64")

        with ui.card().classes("w-full"):
            ui.label("🎯 Target Coverage").classes("text-lg font-semibold mb-2")
            target_boxplot = ui.echart(
                {
                    "backgroundColor": "transparent",
                    "title": {
                        "text": "Target Coverage",
                        "left": "center",
                        "top": 10,
                        "textStyle": {"fontSize": 16, "color": "#000"},
                    },
                    "grid": {
                        "left": "5%",
                        "right": "5%",
                        "bottom": "10%",
                        "top": "20%",
                        "containLabel": True,
                    },
                    "dataset": [
                        {
                            "id": "raw",
                            "dimensions": [
                                "chrom",
                                "min",
                                "Q1",
                                "median",
                                "Q3",
                                "max",
                                "chrom_index",
                            ],
                            "source": [],
                        },
                        {
                            "id": "rawdata",
                            "dimensions": ["chrom", "coverage", "name"],
                            "source": [],
                        },
                        {
                            "id": "outliers",
                            "dimensions": ["chrom", "coverage", "name"],
                            "source": [],
                        },
                        {
                            "id": "globaloutliers",
                            "dimensions": ["chrom", "coverage", "name"],
                            "source": [],
                        },
                    ],
                    "xAxis": {
                        "type": "category",
                        "name": "Chromosome",
                        "nameGap": 30,
                        "axisLabel": {"interval": 0, "rotate": 45},
                    },
                    "yAxis": {"type": "value", "name": "Coverage (x)"},
                    "legend": {
                        "top": 50,
                        "selected": {
                            "box plot": True,
                            "outliers": True,
                            "global outliers": True,
                            "raw data": False,
                        },
                    },
                    "series": [
                        {
                            "name": "box plot",
                            "type": "boxplot",
                            "datasetId": "raw",
                            "encode": {
                                "x": "chrom",
                                "y": ["min", "Q1", "median", "Q3", "max"],
                                "itemName": ["chrom"],
                                "tooltip": ["min", "Q1", "median", "Q3", "max"],
                            },
                        },
                        {
                            "name": "outliers",
                            "type": "scatter",
                            "datasetId": "outliers",
                            "symbolSize": 6,
                            "label": {
                                "show": True,
                                "position": "right",
                                "formatter": "{@name}",
                            },
                            "encode": {
                                "x": "chrom",
                                "y": "coverage",
                                "label": ["name"],
                                "tooltip": ["name", "coverage"],
                            },
                        },
                        {
                            "name": "global outliers",
                            "type": "scatter",
                            "datasetId": "globaloutliers",
                            "symbolSize": 6,
                            "label": {
                                "show": True,
                                "position": "right",
                                "formatter": "{@name}",
                            },
                            "encode": {
                                "x": "chrom",
                                "y": "coverage",
                                "label": ["name"],
                                "tooltip": ["name", "coverage"],
                            },
                        },
                    ],
                }
            ).classes("w-full h-80")
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-3"):
                target_search = ui.input("Search targets…").props(
                    "borderless dense clearable"
                )
            _, target_cov_table = styled_table(
                columns=[
                    {
                        "name": "chrom",
                        "label": "Chrom",
                        "field": "chrom",
                        "sortable": True,
                    },
                    {
                        "name": "startpos",
                        "label": "Start",
                        "field": "startpos",
                        "sortable": True,
                    },
                    {
                        "name": "endpos",
                        "label": "End",
                        "field": "endpos",
                        "sortable": True,
                    },
                    {
                        "name": "name",
                        "label": "Name",
                        "field": "name",
                        "sortable": True,
                    },
                    {
                        "name": "coverage",
                        "label": "Coverage (x)",
                        "field": "coverage",
                        "sortable": True,
                    },
                ],
                rows=[],
                pagination=20,
                class_size="table-xs",
            )
            try:
                target_cov_table.props(
                    'multi-sort rows-per-page-options="[10,20,50,0]"'
                )
                target_search.bind_value(target_cov_table, "filter")
            except Exception:
                pass
            target_cov_table.add_slot(
                "body-cell-coverage",
                """
    <q-td key="coverage" :props="props">
    <q-badge :color="props.value >= 30 ? 'green' : props.value >= 20 ? 'blue' : props.value >= 10 ? 'orange' : 'red'">
        {{ Number(props.value).toFixed(2) }}
    </q-badge>
    </q-td>
    """,
            )

    # IGV viewer section
    with ui.card().classes("w-full"):
        ui.label("🧬 IGV").classes("text-lg font-semibold mb-2")
        igv_div = ui.element("div").classes("w-full h-96 border")
        igv_status = ui.label("IGV will load once an IGV-ready BAM is available.").classes(
            "text-sm text-gray-600"
        )
        


        async def _trigger_build_sorted_bam() -> None:
            try:
                # Debug: Check what's in the launcher
                debug_info = f"Launcher type: {type(launcher).__name__}"
                if hasattr(launcher, 'workflow_runner'):
                    debug_info += f", Has workflow_runner: {launcher.workflow_runner is not None}"
                    if launcher.workflow_runner is not None:
                        debug_info += f", Runner type: {type(launcher.workflow_runner).__name__}"
                else:
                    debug_info += ", No workflow_runner attribute"
                
                # Check if we have a workflow runner
                if hasattr(launcher, 'workflow_runner') and launcher.workflow_runner:
                    # Check if the target.bam file exists
                    target_bam = sample_dir / "target.bam"
                    if not target_bam.exists():
                        ui.notify("No target.bam file found. Please run target analysis first.", type="warning")
                        return
                    
                    # Check if IGV BAM already exists
                    igv_bam = sample_dir / "igv" / "igv_ready.bam"
                    if igv_bam.exists() and (sample_dir / "igv" / "igv_ready.bam.bai").exists():
                        ui.notify("IGV BAM already exists and is ready.", type="positive")
                        return
                    
                    # Submit the IGV BAM generation job to the workflow system
                    try:
                        runner = launcher.workflow_runner
                        sample_id = sample_dir.name
                        
                        # Check if it's a Ray workflow or Simple workflow
                        if hasattr(runner, 'submit_sample_job'):
                            # Simple workflow
                            success = runner.submit_sample_job(str(sample_dir), 'igv_bam', sample_id)
                        elif hasattr(runner, 'manager') and hasattr(runner.manager, 'submit_sample_job'):
                            # Ray workflow - need to get the coordinator
                            # This is more complex for Ray, so we'll provide a fallback
                            ui.notify(
                                "Ray workflow detected. IGV BAM generation should happen automatically "
                                "after target analysis completes. If you need to regenerate it, "
                                "please restart the workflow or check the logs.", 
                                type="info"
                            )
                            return
                        else:
                            ui.notify("Unknown workflow type. Cannot submit IGV BAM job.", type="warning")
                            return
                        
                        if success:
                            ui.notify("IGV BAM generation job submitted to workflow queue!", type="positive")
                        else:
                            ui.notify("Failed to submit IGV BAM generation job. Check logs for details.", type="warning")
                            
                    except Exception as e:
                        ui.notify(f"Error submitting IGV BAM job: {e}", type="negative")
                        
                else:
                    ui.notify(f"No workflow runner available. Debug: {debug_info}. GUI is running but workflow integration is not available. This may be a configuration issue.", type="warning")
            except Exception as e:
                try:
                    ui.notify(f"Error checking IGV BAM status: {e}", type="negative")
                except Exception:
                    # Client may have disconnected, ignore UI errors
                    pass
        with ui.row().classes("items-center gap-2 mt-2"):
            ui.button("Generate IGV BAM", on_click=_trigger_build_sorted_bam)

    # Helpers
    def _log_notify(message: str, level: str = "warning", notify: bool = False) -> None:
        try:
            # Python logging
            log_func = getattr(logging, level, logging.warning)
            log_func(f"[coverage] {message}")
            # Append to GUI log buffer if available
            try:
                launcher._log_buffer.append(f"[coverage] {message}\n")  # type: ignore[attr-defined]
                if hasattr(launcher, "log_area") and launcher.log_area is not None:  # type: ignore[attr-defined]
                    launcher.log_area.set_value("".join(launcher._log_buffer))  # type: ignore[attr-defined]
            except Exception:
                pass
            # User notification (optional)
            if notify and ui is not None:
                n_type = "warning" if level in {"warning", "debug", "info"} else "error"
                try:
                    ui.notify(message, type=n_type)
                except Exception:
                    pass
        except Exception:
            # Never let logging itself break the UI
            pass

    def _update_chr_cov(cov_df: pd.DataFrame) -> None:
        try:

            def chr_key(label: str) -> int:
                s = str(label)
                if s.startswith("chr"):
                    s = s[3:]
                mapping = {"X": 23, "Y": 24}
                try:
                    return int(s)
                except Exception:
                    return mapping.get(s, 1000)

            pattern = r"^chr([0-9]+|X|Y)$"
            name_col = "#rname" if "#rname" in cov_df.columns else "rname"
            temp_df = cov_df[cov_df[name_col].astype(str).str.match(pattern)]
            temp_df = temp_df[temp_df[name_col] != "chrM"]
            names = sorted(temp_df[name_col].astype(str).unique(), key=chr_key)
            temp_df = temp_df.set_index(name_col).loc[names].reset_index()
            echart_chr_cov.options["xAxis"]["data"] = names
            echart_chr_cov.options["series"] = [
                {
                    "type": "bar",
                    "name": "Chromosome",
                    "barWidth": "60%",
                    "data": [float(v) for v in temp_df["meandepth"].tolist()],
                }
            ]
            echart_chr_cov.update()
        except Exception as e:
            _log_notify(
                f"Chromosome coverage update failed: {e}", level="warning", notify=False
            )

    def _update_target_cov(cov_df: pd.DataFrame, bed_df: pd.DataFrame) -> None:
        try:

            def chr_key(label: str) -> int:
                s = str(label)
                if s.startswith("chr"):
                    s = s[3:]
                mapping = {"X": 23, "Y": 24}
                try:
                    return int(s)
                except Exception:
                    return mapping.get(s, 1000)

            bed_df = bed_df.copy()
            bed_df["length"] = (bed_df["endpos"] - bed_df["startpos"] + 1).astype(float)
            grouped = (
                bed_df.groupby("chrom")
                .agg({"bases": "sum", "length": "sum"})
                .reset_index()
            )
            grouped["meandepth"] = grouped["bases"] / grouped["length"]
            pattern = r"^chr([0-9]+|X|Y)$"
            name_col = "#rname" if "#rname" in cov_df.columns else "rname"
            temp_df = cov_df[cov_df[name_col].astype(str).str.match(pattern)]
            temp_df = temp_df[temp_df[name_col] != "chrM"]
            names = sorted(temp_df[name_col].astype(str).unique(), key=chr_key)
            echart_target_cov.options["xAxis"]["data"] = names
            grouped = grouped.set_index("chrom").reindex(names).reset_index()
            echart_target_cov.options["series"] = [
                {
                    "type": "scatter",
                    "name": "Off Target",
                    "symbolSize": 8,
                    "data": [
                        float(v)
                        for v in temp_df.set_index(name_col)
                        .loc[names]["meandepth"]
                        .tolist()
                    ],
                },
                {
                    "type": "scatter",
                    "name": "On Target",
                    "symbolSize": 8,
                    "data": [float(v) for v in grouped["meandepth"].fillna(0).tolist()],
                },
            ]
            echart_target_cov.update()
        except Exception as e:
            _log_notify(
                f"Target vs Off-target coverage update failed: {e}",
                level="warning",
                notify=False,
            )

    def _update_boxplot(bed_df: pd.DataFrame) -> None:
        try:
            df = bed_df.copy()
            if "coverage" not in df.columns:
                df["length"] = (df["endpos"] - df["startpos"] + 1).astype(float)
                df["coverage"] = df["bases"] / df["length"]
            chroms = natsort.natsorted(df["chrom"].astype(str).unique())
            chrom_lookup = {chrom: idx for idx, chrom in enumerate(chroms)}
            df["chrom_index"] = df["chrom"].map(chrom_lookup)
            agg = (
                df.groupby("chrom")
                .agg(
                    min=("coverage", "min"),
                    Q1=("coverage", lambda x: float(np.percentile(x, 25))),
                    median=("coverage", "median"),
                    Q3=("coverage", lambda x: float(np.percentile(x, 75))),
                    max=("coverage", "max"),
                    chrom_index=("chrom_index", "first"),
                )
                .reset_index()
            )
            agg["chrom"] = pd.Categorical(agg["chrom"], categories=chroms, ordered=True)
            agg = agg.sort_values("chrom").reset_index(drop=True)
            result = [
                ["chrom", "min", "Q1", "median", "Q3", "max", "chrom_index"]
            ] + agg.values.tolist()

            def iqr_bounds(sub):
                q1 = np.percentile(sub["coverage"], 25)
                q3 = np.percentile(sub["coverage"], 75)
                iqr = q3 - q1
                return q1 - 1.5 * iqr, q3 + 1.5 * iqr

            out_rows = []
            for c in chroms:
                sub = df[df["chrom"].astype(str) == c]
                lb, ub = iqr_bounds(sub)
                out = sub[(sub["coverage"] < lb) | (sub["coverage"] > ub)]
                out_rows += out[["chrom", "coverage", "name"]].values.tolist()
            lb_g, ub_g = iqr_bounds(df)
            glob = df[(df["coverage"] < lb_g) | (df["coverage"] > ub_g)]
            glob_rows = glob[["chrom", "coverage", "name"]].values.tolist()
            target_boxplot.options["dataset"][0]["source"] = result
            target_boxplot.options["dataset"][1]["source"] = df[
                ["chrom", "coverage", "name"]
            ].values.tolist()
            target_boxplot.options["dataset"][2]["source"] = out_rows
            target_boxplot.options["dataset"][3]["source"] = glob_rows
            target_boxplot.options["xAxis"]["data"] = chroms
            target_boxplot.update()
        except Exception as e:
            _log_notify(
                f"Target boxplot update failed: {e}", level="warning", notify=False
            )

    def _update_time(npy_path: Path) -> None:
        try:
            arr = np.load(npy_path)
            echart_time.options["series"][0]["data"] = arr.tolist()
            echart_time.update()
        except Exception as e:
            _log_notify(
                f"Coverage-over-time update failed: {e}", level="warning", notify=False
            )

    def _update_target_table(df: pd.DataFrame) -> None:
        try:
            dfx = df.copy()
            if "coverage" not in dfx.columns:
                dfx["length"] = (dfx["endpos"] - dfx["startpos"] + 1).astype(float)
                dfx["coverage"] = dfx["bases"] / dfx["length"]
            dfx["coverage"] = dfx["coverage"].astype(float).round(2)
            rows = dfx[["chrom", "startpos", "endpos", "name", "coverage"]].to_dict(
                orient="records"
            )
            target_cov_table.rows = rows
            target_cov_table.update()
        except Exception as e:
            _log_notify(
                f"Target table update failed: {e}", level="warning", notify=False
            )

    def _refresh_coverage() -> None:
        try:
            if not sample_dir or not sample_dir.exists():
                _log_notify(
                    f"Sample directory not found: {sample_dir}",
                    level="warning",
                    notify=True,
                )
                return
            key = str(sample_dir)
            state = launcher._coverage_state.get(key, {})
            cov_main = sample_dir / "coverage_main.csv"
            bed_cov = sample_dir / "bed_coverage_main.csv"
            target_cov = sample_dir / "target_coverage.csv"
            cov_time = sample_dir / "coverage_time_chart.npy"
            if cov_main.exists():
                m = cov_main.stat().st_mtime
                if state.get("cov_main_mtime") != m or "cov_df" not in state:
                    try:
                        cov_df = pd.read_csv(cov_main)
                        # Backfill meandepth if missing
                        if "meandepth" not in cov_df.columns and {
                            "covbases",
                            "endpos",
                        }.issubset(set(cov_df.columns)):
                            cov_df = cov_df.copy()
                            with np.errstate(divide="ignore", invalid="ignore"):
                                cov_df["meandepth"] = (
                                    cov_df["covbases"]
                                    / cov_df["endpos"].replace(0, np.nan)
                                ).fillna(0)
                        _update_chr_cov(cov_df)
                        state["cov_df"] = cov_df
                        state["cov_main_mtime"] = m  # only set on success
                    except Exception as e:
                        _log_notify(
                            f"Failed to read coverage_main.csv: {e}",
                            level="error",
                            notify=True,
                        )
            else:
                _log_notify("coverage_main.csv not found", level="debug", notify=False)
            if bed_cov.exists():
                m = bed_cov.stat().st_mtime
                if state.get("bed_cov_mtime") != m or "bed_df" not in state:
                    try:
                        bed_df = pd.read_csv(bed_cov)
                        state["bed_df"] = bed_df
                        _update_boxplot(bed_df)
                        if state.get("cov_df") is not None:
                            _update_target_cov(state["cov_df"], bed_df)
                        state["bed_cov_mtime"] = m  # only set on success
                    except Exception as e:
                        _log_notify(
                            f"Failed to read bed_coverage_main.csv: {e}",
                            level="error",
                            notify=True,
                        )
            else:
                _log_notify(
                    "bed_coverage_main.csv not found", level="debug", notify=False
                )
            if target_cov.exists():
                m = target_cov.stat().st_mtime
                if state.get("target_cov_mtime") != m:
                    try:
                        tdf = pd.read_csv(target_cov)
                        _update_target_table(tdf)
                        state["target_cov_mtime"] = m  # only set on success
                    except Exception as e:
                        _log_notify(
                            f"Failed to read target_coverage.csv: {e}",
                            level="error",
                            notify=True,
                        )
            else:
                _log_notify(
                    "target_coverage.csv not found", level="debug", notify=False
                )
            if cov_time.exists():
                m = cov_time.stat().st_mtime
                if state.get("cov_time_mtime") != m:
                    try:
                        _update_time(cov_time)
                        state["cov_time_mtime"] = m
                    except Exception as e:
                        _log_notify(
                            f"Failed to load coverage_time_chart.npy: {e}",
                            level="warning",
                            notify=False,
                        )
            # Summary
            try:
                global_cov = None
                target_cov_v = None
                enrich_v = None
                if state.get("cov_df") is not None:
                    cdf = state["cov_df"]
                    if (
                        "covbases" in cdf.columns
                        and "endpos" in cdf.columns
                        and cdf["endpos"].sum() > 0
                    ):
                        global_cov = float(cdf["covbases"].sum()) / float(
                            cdf["endpos"].sum()
                        )
                if state.get("bed_df") is not None:
                    bdf = state["bed_df"].copy()
                    if "length" not in bdf.columns:
                        bdf["length"] = (bdf["endpos"] - bdf["startpos"] + 1).astype(
                            float
                        )
                    if bdf["length"].sum() > 0:
                        target_cov_v = float(bdf["bases"].sum()) / float(
                            bdf["length"].sum()
                        )
                if (
                    global_cov is not None
                    and target_cov_v is not None
                    and global_cov > 0
                ):
                    enrich_v = target_cov_v / global_cov
                if target_cov_v is not None:
                    if target_cov_v >= 30:
                        q_text, q_cls, q_bg = (
                            "Excellent",
                            "text-green-600",
                            "bg-green-100",
                        )
                    elif target_cov_v >= 20:
                        q_text, q_cls, q_bg = "Good", "text-blue-600", "bg-blue-100"
                    elif target_cov_v >= 10:
                        q_text, q_cls, q_bg = (
                            "Moderate",
                            "text-yellow-600",
                            "bg-yellow-100",
                        )
                    else:
                        q_text, q_cls, q_bg = (
                            "Insufficient",
                            "text-red-600",
                            "bg-red-100",
                        )
                    cov_quality_label.set_text(f"Quality: {q_text}")
                    try:
                        cov_quality_label.classes(
                            replace=f"text-sm font-medium {q_cls}"
                        )
                        cov_quality_badge.classes(
                            replace=f"px-2 py-1 rounded {q_bg} {q_cls}"
                        )
                    except Exception:
                        pass
                    cov_quality_badge.set_text(f"{target_cov_v:.2f}x")
                if global_cov is not None:
                    cov_global_lbl.set_text(
                        f"Global Estimated Coverage: {global_cov:.2f}x"
                    )
                if target_cov_v is not None:
                    cov_target_lbl.set_text(
                        f"Targets Estimated Coverage: {target_cov_v:.2f}x"
                    )
                if enrich_v is not None:
                    cov_enrich_lbl.set_text(f"Estimated enrichment: {enrich_v:.2f}x")
            except Exception as e:
                _log_notify(
                    f"Coverage summary update failed: {e}",
                    level="warning",
                    notify=False,
                )
            launcher._coverage_state[key] = state
            # IGV auto-init and auto-load when BAM becomes available (UI-only, light work)
            try:
                # Determine candidate BAMs (prefer standardized IGV path)
                igv_dir = sample_dir / "igv"
                clair_dir = sample_dir / "clair3"
                candidates = [
                    igv_dir / "igv_ready.bam",
                    clair_dir / "sorted_targets_exceeding_rerun.bam",
                    clair_dir / "sorted_targets_exceeding.bam",
                    sample_dir / "target.bam",
                ]
                bam_path = None
                for p in candidates:
                    bai = Path(f"{p}.bai")
                    if p.exists() and bai.exists():
                        bam_path = p
                        break

                # Mount directory once
                if bam_path is not None and app is not None:
                    if bam_path.parent == igv_dir:
                        mount = f"/samples/{sample_dir.name}/igv"
                    elif bam_path.parent == clair_dir:
                        mount = f"/samples/{sample_dir.name}/clair3"
                    else:
                        mount = f"/samples/{sample_dir.name}"
                    try:
                        app.add_static_files(mount, str(bam_path.parent))
                    except Exception:
                        # Already mounted or in tests without app
                        pass

                    # Define URLs first
                    bam_url = f"{mount}/{bam_path.name}"
                    bai_url = f"{mount}/{bam_path.name}.bai"
                    
                    # Initialize IGV once per page (include track immediately to avoid race)
                    if not state.get("igv_initialized"):
                        js_create = f"""
                            const el = getElement('{igv_div.id}');
                            const track = {{ name: '{sample_dir.name}', url: '{bam_url}', indexURL: '{bai_url}', format: 'bam', visibilityWindow: 1000000 }};
                            const options = {{ genome: 'hg38', locus: 'chr1:1-1', tracks: [track] }};
                            igv.createBrowser(el, options).then(b => {{ window.lj_igv = b; }});
                        """
                        try:
                            ui.run_javascript(js_create, timeout=30.0)
                            state["igv_initialized"] = True
                            igv_status.set_text("IGV initialised.")
                        except Exception:
                            pass

                    # Load/refresh BAM track if changed
                    if state.get("igv_initialized") and state.get("igv_loaded_bam") != bam_url:
                        js_track = f"""
                            if (window.lj_igv) {{
                                const name = '{sample_dir.name}';
                                const track = {{ name, url: '{bam_url}', indexURL: '{bai_url}', format: 'bam', visibilityWindow: 1000000 }};
                                window.lj_igv.loadTrack(track);
                            }}
                        """
                        try:
                            ui.run_javascript(js_track, timeout=60.0)
                            state["igv_loaded_bam"] = bam_url
                            igv_status.set_text("BAM loaded in IGV.")
                        except Exception:
                            pass
                else:
                    try:
                        igv_status.set_text("No IGV-viewable BAM found yet (need coordinate-sorted BAM with .bai).")
                    except Exception:
                        pass
                

            except Exception:
                # Never let IGV logic break coverage refresh
                pass
        except Exception as e:
            _log_notify(
                f"Unexpected coverage refresh error: {e}", level="error", notify=True
            )

    # Trigger an immediate refresh on page load, then continue with periodic refreshes
    try:
        ui.timer(0.5, _refresh_coverage, once=True)
    except Exception:
        pass
    ui.timer(30.0, _refresh_coverage, active=True)
