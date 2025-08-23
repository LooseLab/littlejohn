from __future__ import annotations

from typing import Any, Dict, List, Tuple
from pathlib import Path

import natsort
import numpy as np
from functools import lru_cache
import logging
import importlib.resources as importlib_resources
import pandas as pd

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None


def add_cnv_section(launcher: Any, sample_dir: Path) -> None:
    """Build the CNV UI section and attach refresh timers.

    Uses `launcher._cnv_state` for per-sample cache/state.
    Expects CNV.npy, CNV3.npy, CNV_dict.npy, XYestimate.pkl, cnv_data_array.npy in sample folder.
    
    Controls now trigger immediate refresh instead of waiting for timer updates.
    """
    with ui.card().classes("w-full"):
        ui.label("🧬 Copy Number Variation (CNV)").classes("text-lg font-semibold mb-2")
        with ui.row().classes("w-full items-center justify-between mb-2"):
            with ui.column().classes("gap-1"):
                cnv_status = ui.label("Status: Awaiting Data").classes("text-gray-600")
                cnv_xy = ui.label("Genetic Sex: --").classes("text-gray-600")
            with ui.column().classes("gap-1 items-end"):
                cnv_bin = ui.label("Bin Width: --").classes("text-sm text-gray-600")
                cnv_var = ui.label("Variance: --").classes("text-sm text-gray-600")
        # Controls similar to original module
        with ui.row().classes("gap-4 items-center mb-2"):
            ui.label("Select Chromosome").classes("text-sm")
            cnv_chrom_select = ui.select(options={"All": "All"}, value="All").style(
                "width: 160px"
            )
            ui.label("Select Gene").classes("text-sm ml-4")
            cnv_gene_select = ui.select(options={"All": "All"}, value="All").style(
                "width: 200px"
            )
            ui.label("Color by").classes("text-sm ml-4")
            cnv_color = ui.toggle(
                options={"chromosome": "Chromosome", "value": "Up/Down"},
                value="chromosome",
            ).classes("mt-1")
            ui.label("Y-axis scale").classes("text-sm ml-4")
            cnv_scale = ui.toggle(
                options={"linear": "Linear", "log": "Log"}, value="linear"
            ).classes("mt-1")
            ui.label("Breakpoints").classes("text-sm ml-4")
            cnv_bp = ui.toggle(
                options={"hide": "Hide", "show": "Show"}, value="show"
            ).classes("mt-1")
        cnv_abs = ui.echart(
            {
                "backgroundColor": "transparent",
                "title": {"text": "CNV Scatter Plot", "left": "center", "top": 10},
                "grid": {
                    "left": "5%",
                    "right": "5%",
                    "bottom": "10%",
                    "top": "20%",
                    "containLabel": True,
                },
                "tooltip": {"trigger": "axis"},
                "xAxis": {"type": "value", "max": "dataMax"},
                "yAxis": [
                    {"type": "value", "name": "Ploidy"},
                    {
                        "type": "value",
                        "name": "Breakpoint density",
                        "position": "right",
                    },
                ],
                "dataZoom": [
                    {"type": "slider", "xAxisIndex": [0]},
                    {"type": "slider", "yAxisIndex": [0, 1], "right": 20},
                ],
                "series": [
                    {"type": "scatter", "name": "CNV", "symbolSize": 3, "data": []},
                    {
                        "type": "scatter",
                        "name": "centromeres_highlight",
                        "data": [],
                        "symbolSize": 3,
                        "markArea": {
                            "itemStyle": {"color": "rgba(135, 206, 250, 0.4)"},
                            "data": [],
                        },
                    },
                    {
                        "type": "scatter",
                        "name": "cytobands_highlight",
                        "data": [],
                        "symbolSize": 3,
                        "markArea": {
                            "itemStyle": {"color": "rgba(200, 200, 200, 0.4)"},
                            "data": [],
                        },
                        "markLine": {"symbol": "none", "data": []},
                    },
                ],
            }
        ).classes("w-full h-72")
        cnv_diff = ui.echart(
            {
                "backgroundColor": "transparent",
                "title": {"text": "Difference Plot", "left": "center", "top": 10},
                "grid": {
                    "left": "5%",
                    "right": "5%",
                    "bottom": "10%",
                    "top": "20%",
                    "containLabel": True,
                },
                "tooltip": {"trigger": "axis"},
                "xAxis": {"type": "value", "max": "dataMax"},
                "yAxis": [
                    {"type": "value", "name": "Relative"},
                    {
                        "type": "value",
                        "name": "Breakpoint density",
                        "position": "right",
                    },
                ],
                "dataZoom": [
                    {"type": "slider", "xAxisIndex": [0]},
                    {"type": "slider", "yAxisIndex": [0, 1], "right": 20},
                ],
                "series": [
                    {"type": "scatter", "name": "CNV Δ", "symbolSize": 3, "data": []},
                    {
                        "type": "scatter",
                        "name": "centromeres_highlight",
                        "data": [],
                        "symbolSize": 3,
                        "markArea": {
                            "itemStyle": {"color": "rgba(135, 206, 250, 0.4)"},
                            "data": [],
                        },
                    },
                    {
                        "type": "scatter",
                        "name": "cytobands_highlight",
                        "data": [],
                        "symbolSize": 3,
                        "markArea": {
                            "itemStyle": {"color": "rgba(200, 200, 200, 0.4)"},
                            "data": [],
                        },
                        "markLine": {"symbol": "none", "data": []},
                    },
                ],
            }
        ).classes("w-full h-72")

        # Cytoband analysis table (updated when a single chromosome is selected)
        ui.label("CNV Cytoband Analysis").classes("text-md font-semibold mt-4")
        cyto_summary = ui.label("Select a chromosome to view cytoband events").classes(
            "text-gray-600 mb-1"
        )
        cyto_columns = [
            {"name": "chrom", "label": "Chr", "field": "chrom", "sortable": True},
            {"name": "region", "label": "Region", "field": "region", "sortable": True},
            {
                "name": "start_mb",
                "label": "Start (Mb)",
                "field": "start_mb",
                "sortable": True,
                "align": "right",
            },
            {
                "name": "end_mb",
                "label": "End (Mb)",
                "field": "end_mb",
                "sortable": True,
                "align": "right",
            },
            {
                "name": "length_mb",
                "label": "Length (Mb)",
                "field": "length_mb",
                "sortable": True,
                "align": "right",
            },
            {
                "name": "mean_cnv",
                "label": "Mean CNV",
                "field": "mean_cnv",
                "sortable": True,
                "align": "right",
            },
            {
                "name": "state",
                "label": "State",
                "field": "state",
                "sortable": True,
                "align": "center",
            },
            {"name": "genes", "label": "Genes", "field": "genes"},
        ]
        cyto_table = ui.table(columns=cyto_columns, rows=[], pagination=20).classes(
            "w-full"
        )
        try:
            cyto_table.props('multi-sort rows-per-page-options="[10,20,50,0]"')
        except Exception:
            pass

    # Adaptive thinning helpers
    MAX_POINTS_PER_CHART = 2000

    def _get_visible_range(chart, series_list):
        try:
            # Compute overall data span
            min_x = None
            max_x = None
            for s in series_list:
                data = s.get("data") or []
                if not data:
                    continue
                sx = data[0][0] if isinstance(data[0], (list, tuple)) else None
                ex = data[-1][0] if isinstance(data[-1], (list, tuple)) else None
                if sx is None or ex is None:
                    # Fallback compute min/max
                    for p in data:
                        try:
                            x = float(p[0])
                        except Exception:
                            continue
                        min_x = x if min_x is None else min(min_x, x)
                        max_x = x if max_x is None else max(max_x, x)
                else:
                    vmin = min(float(sx), float(ex))
                    vmax = max(float(sx), float(ex))
                    min_x = vmin if min_x is None else min(min_x, vmin)
                    max_x = vmax if max_x is None else max(max_x, vmax)
            dz = None
            try:
                dzo = chart.options.get("dataZoom")
                if isinstance(dzo, list) and dzo:
                    dz = dzo[0]
            except Exception:
                dz = None
            if not dz:
                return (min_x, max_x)
            # Prefer explicit values
            sv = dz.get("startValue") if isinstance(dz, dict) else None
            ev = dz.get("endValue") if isinstance(dz, dict) else None
            if sv is not None or ev is not None:
                left = float(sv) if sv is not None else min_x
                right = float(ev) if ev is not None else max_x
                return (left, right)
            # Fallback to percentage range
            sp = dz.get("start") if isinstance(dz, dict) else None
            ep = dz.get("end") if isinstance(dz, dict) else None
            if (
                (sp is not None or ep is not None)
                and min_x is not None
                and max_x is not None
            ):
                width = (
                    max_x - min_x if max_x is not None and min_x is not None else None
                )
                if width and width > 0:
                    left = min_x + (float(sp or 0) / 100.0) * width
                    right = min_x + (float(ep or 100) / 100.0) * width
                    return (left, right)
            return (min_x, max_x)
        except Exception:
            return (None, None)

    def _evenly_sample(seq, k):
        try:
            n = len(seq)
            if k <= 0 or n <= k:
                return list(seq)
            if k == 1:
                return [seq[n // 2]]
            # Choose k indices evenly across [0, n-1]
            return [seq[int(round(i * (n - 1) / (k - 1)))] for i in range(k)]
        except Exception:
            return list(seq)[:k]

    def _thin_chart_series(chart, max_points: int = MAX_POINTS_PER_CHART) -> None:
        try:
            series = chart.options.get("series", [])
            # Identify data series to thin (exclude overlays)
            data_idx = []
            data_series = []
            for idx, s in enumerate(series):
                name = s.get("name", "")
                if s.get("type") == "scatter" and name not in (
                    "centromeres_highlight",
                    "cytobands_highlight",
                ):
                    data = s.get("data") or []
                    if isinstance(data, list) and data:
                        data_idx.append(idx)
                        data_series.append(s)
            if not data_series:
                return
            x_range = _get_visible_range(chart, data_series)
            # Gather visible counts and data within range
            vis_data = []
            total = 0
            left, right = x_range
            for s in data_series:
                pts = s.get("data") or []
                if left is not None and right is not None:
                    sub = [
                        p
                        for p in pts
                        if isinstance(p, (list, tuple)) and left <= float(p[0]) <= right
                    ]
                else:
                    sub = list(pts)
                vis_data.append(sub)
                total += len(sub)
            if total <= max_points:
                return
            # Allocate budgets proportional to visible counts with a small floor
            budgets = []
            remaining = max_points
            for sub in vis_data:
                share = int(max(1, round((len(sub) / total) * max_points)))
                budgets.append(share)
            # Normalize budgets to exactly max_points
            adj = sum(budgets) - max_points
            i = 0
            while adj != 0 and budgets:
                if adj > 0 and budgets[i] > 1:
                    budgets[i] -= 1
                    adj -= 1
                elif adj < 0:
                    budgets[i] += 1
                    adj += 1
                i = (i + 1) % len(budgets)
            # Apply sampling and replace data (preserve points outside range sparsely)
            for (idx, s), sub, k in zip(zip(data_idx, data_series), vis_data, budgets):
                original = s.get("data") or []
                # Keep outside-range points sparsely so context remains when zoomed out/in
                if left is not None and right is not None:
                    outside = [
                        p
                        for p in original
                        if isinstance(p, (list, tuple))
                        and not (left <= float(p[0]) <= right)
                    ]
                    outside_keep = _evenly_sample(
                        outside, max(0, k // 10)
                    )  # at most 10% of budget
                else:
                    outside_keep = []
                inside_keep = _evenly_sample(sub, max(1, k - len(outside_keep)))
                new_data = inside_keep + outside_keep
                chart.options["series"][idx]["data"] = new_data
        except Exception:
            pass

    @lru_cache(maxsize=1)
    def _load_centromere_regions() -> Dict[str, List[Tuple[int, int, str]]]:
        """Load centromere/satellite regions from packaged resources.
        Returns mapping: chrom -> list of (start_bp, end_bp, name).
        """
        regions: Dict[str, List[Tuple[int, int, str]]] = {}
        try:
            res_path = (
                importlib_resources.files("littlejohn.resources") / "cenSatRegions.bed"
            )
            with res_path.open("r") as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) < 4:
                        continue
                    chrom, start, end, name = (
                        parts[0],
                        int(parts[1]),
                        int(parts[2]),
                        parts[3],
                    )
                    regions.setdefault(chrom, []).append((start, end, name))
        except Exception:
            pass
        return regions

    @lru_cache(maxsize=1)
    def _load_cytobands_df() -> pd.DataFrame:
        try:
            res_path = (
                importlib_resources.files("littlejohn.resources") / "cytoBand.txt"
            )
            df = pd.read_csv(
                res_path,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "name", "stain"],
            )
            return df
        except Exception:
            return pd.DataFrame(columns=["chrom", "start", "end", "name", "stain"])

    @lru_cache(maxsize=1)
    def _load_gene_bed() -> pd.DataFrame:
        # Use a compact gene list; fall back to all_genes.bed
        for fname in ["unique_genes.bed", "all_genes.bed"]:
            try:
                res_path = importlib_resources.files("littlejohn.resources") / fname
                if res_path.exists():
                    return pd.read_csv(
                        res_path,
                        sep="\t",
                        header=None,
                        names=["chrom", "start", "end", "gene"],
                    )
            except Exception:
                continue
        return pd.DataFrame(columns=["chrom", "start", "end", "gene"])

    def _sex_label(xy_val: Any) -> str:
        try:
            s = str(xy_val).strip().upper()
            if s in ("MALE", "XY"):
                return "Male"
            if s in ("FEMALE", "XX"):
                return "Female"
        except Exception:
            pass
        return "Unknown"

    def _analyze_cytoband_cnv(
        cnv_data: Dict[str, np.ndarray],
        chromosome: str,
        bin_width: int,
        sex_estimate: str,
    ) -> pd.DataFrame:
        try:
            import numpy as _np
            import pandas as _pd

            logger = logging.getLogger(__name__)
            if not cnv_data or chromosome not in cnv_data:
                return _pd.DataFrame()
            if bin_width > 10_000_000:
                return _pd.DataFrame()

            cyto_df = _load_cytobands_df()
            chromosome_cytobands = cyto_df[cyto_df["chrom"] == chromosome].copy()
            if chromosome_cytobands.empty:
                return _pd.DataFrame()

            # Centromere mask
            centro = _load_centromere_regions()
            mask = _np.ones(len(cnv_data[chromosome]), dtype=bool)
            cent_regions = centro.get(chromosome, [])
            if cent_regions:
                # Combine all annotated centromere/satellite spans
                for s_bp, e_bp, _nm in cent_regions:
                    s_bin = max(0, int(s_bp // bin_width))
                    e_bin = min(len(mask), int(e_bp // bin_width))
                    if e_bin > s_bin:
                        mask[s_bin:e_bin] = False
            chr_cnv = _np.asarray(cnv_data[chromosome])
            chr_cnv = chr_cnv[mask] if mask.any() else _np.asarray(cnv_data[chromosome])

            # Chromosome-wide stats
            chr_mean = float(_np.mean(chr_cnv)) if chr_cnv.size else 0.0
            chr_std = float(_np.std(chr_cnv)) if chr_cnv.size else 1.0

            # Autosomal distribution for whole-chrom thresholds
            chrom_means: List[float] = []
            for ck, arr in cnv_data.items():
                if ck.startswith("chr") and ck[3:].isdigit():
                    a = _np.asarray(arr)
                    if a.size:
                        chrom_means.append(float(_np.mean(a)))
            means_mean = float(_np.mean(chrom_means)) if chrom_means else 0.0
            means_std = float(_np.std(chrom_means)) if chrom_means else 1.0

            # Thresholds
            sex_lbl = sex_estimate
            if chromosome.startswith("chr") and chromosome[3:].isdigit():
                gain_threshold = means_mean + (1.0 * means_std)
                loss_threshold = means_mean - (1.0 * means_std)
                cyto_gain_th = chr_mean + (1.0 * chr_std)
                cyto_loss_th = chr_mean - (1.0 * chr_std)
            elif chromosome == "chrX":
                gain_threshold = means_mean + (1.0 * means_std)
                loss_threshold = means_mean - (1.0 * means_std)
                cyto_gain_th = chr_mean + (1.0 * chr_std)
                cyto_loss_th = chr_mean - (1.0 * chr_std)
            else:  # chrY
                if sex_lbl in ("Male", "XY"):
                    gain_threshold = means_mean + (1.0 * means_std)
                    loss_threshold = means_mean - (1.0 * means_std)
                    cyto_gain_th = chr_mean + (1.0 * chr_std)
                    cyto_loss_th = chr_mean - (1.0 * chr_std)
                else:
                    gain_threshold = means_mean + (1.2 * means_std)
                    loss_threshold = means_mean - (1.2 * means_std)
                    cyto_gain_th = chr_mean + (1.2 * chr_std)
                    cyto_loss_th = chr_mean - (1.2 * chr_std)

            # Whole chromosome event detection
            bins_above_gain = float(
                (_np.asarray(cnv_data[chromosome]) > gain_threshold).sum()
            ) / max(1, len(cnv_data[chromosome]))
            bins_below_loss = float(
                (_np.asarray(cnv_data[chromosome]) < loss_threshold).sum()
            ) / max(1, len(cnv_data[chromosome]))
            min_prop = 0.7
            whole_chr_event = False
            whole_chr_state = "NORMAL"
            if bins_above_gain > min_prop:
                whole_chr_event = True
                whole_chr_state = "GAIN"
            elif bins_below_loss > min_prop:
                whole_chr_event = True
                whole_chr_state = "LOSS"

            merged_rows: List[dict] = []
            if whole_chr_event:
                gene_df = _load_gene_bed()
                genes_in_chr = (
                    gene_df[gene_df["chrom"] == chromosome]["gene"].astype(str).tolist()
                )
                merged_rows.append(
                    {
                        "chrom": chromosome,
                        "start_pos": int(chromosome_cytobands["start"].min()),
                        "end_pos": int(chromosome_cytobands["end"].max()),
                        "name": f"{chromosome} WHOLE CHROMOSOME {whole_chr_state}",
                        "mean_cnv": chr_mean,
                        "cnv_state": whole_chr_state,
                        "length": int(chromosome_cytobands["end"].max())
                        - int(chromosome_cytobands["start"].min()),
                        "genes": genes_in_chr,
                    }
                )

            # Group contiguous cytobands by state
            current_group = None
            vals = _np.asarray(cnv_data[chromosome])
            for _, band in chromosome_cytobands.iterrows():
                s_bp = int(band["start"])
                e_bp = int(band["end"])
                s_bin = max(0, s_bp // bin_width)
                e_bin = min(len(vals) - 1, max(0, e_bp // bin_width))
                region = (
                    vals[s_bin : e_bin + 1]
                    if len(vals) and e_bin >= s_bin
                    else _np.array([])
                )
                mean_val = float(_np.mean(region)) if region.size else 0.0
                if whole_chr_event:
                    if whole_chr_state == "GAIN":
                        state = (
                            "LOSS"
                            if mean_val < (chr_mean - 2.0 * chr_std)
                            else "NORMAL"
                        )
                    else:
                        state = (
                            "GAIN"
                            if mean_val > (chr_mean + 2.0 * chr_std)
                            else "NORMAL"
                        )
                else:
                    state = (
                        "GAIN"
                        if mean_val > cyto_gain_th
                        else ("LOSS" if mean_val < cyto_loss_th else "NORMAL")
                    )

                if current_group is None:
                    current_group = {
                        "chrom": chromosome,
                        "start_pos": s_bp,
                        "end_pos": e_bp,
                        "bands": [str(band["name"])],
                        "mean_vals": [mean_val],
                        "cnv_state": state,
                    }
                elif state == current_group["cnv_state"]:
                    current_group["end_pos"] = e_bp
                    current_group["bands"].append(str(band["name"]))
                    current_group["mean_vals"].append(mean_val)
                else:
                    # finalize
                    mean_cnv = (
                        float(_np.mean(current_group["mean_vals"]))
                        if current_group["mean_vals"]
                        else 0.0
                    )
                    row = {
                        "chrom": current_group["chrom"],
                        "start_pos": int(current_group["start_pos"]),
                        "end_pos": int(current_group["end_pos"]),
                        "name": f"{current_group['chrom']} {current_group['bands'][0]}-{current_group['bands'][-1]}",
                        "mean_cnv": mean_cnv,
                        "cnv_state": current_group["cnv_state"],
                        "length": int(current_group["end_pos"])
                        - int(current_group["start_pos"]),
                    }
                    if row["cnv_state"] in ("GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"):
                        gene_df = _load_gene_bed()
                        genes = (
                            gene_df[
                                (gene_df["chrom"] == row["chrom"])
                                & (gene_df["start"] <= row["end_pos"])
                                & (gene_df["end"] >= row["start_pos"])
                            ]["gene"]
                            .astype(str)
                            .tolist()
                        )
                        row["genes"] = genes
                    merged_rows.append(row)
                    # start new group
                    current_group = {
                        "chrom": chromosome,
                        "start_pos": s_bp,
                        "end_pos": e_bp,
                        "bands": [str(band["name"])],
                        "mean_vals": [mean_val],
                        "cnv_state": state,
                    }
            # finalize last
            if current_group is not None:
                mean_cnv = (
                    float(_np.mean(current_group["mean_vals"]))
                    if current_group["mean_vals"]
                    else 0.0
                )
                row = {
                    "chrom": current_group["chrom"],
                    "start_pos": int(current_group["start_pos"]),
                    "end_pos": int(current_group["end_pos"]),
                    "name": f"{current_group['chrom']} {current_group['bands'][0]}-{current_group['bands'][-1]}",
                    "mean_cnv": mean_cnv,
                    "cnv_state": current_group["cnv_state"],
                    "length": int(current_group["end_pos"])
                    - int(current_group["start_pos"]),
                }
                if row["cnv_state"] in ("GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"):
                    gene_df = _load_gene_bed()
                    genes = (
                        gene_df[
                            (gene_df["chrom"] == row["chrom"])
                            & (gene_df["start"] <= row["end_pos"])
                            & (gene_df["end"] >= row["start_pos"])
                        ]["gene"]
                        .astype(str)
                        .tolist()
                    )
                    row["genes"] = genes
                merged_rows.append(row)

            df = _pd.DataFrame(merged_rows)
            if not df.empty:
                df = df.sort_values("start_pos")
            return df
        except Exception:
            return pd.DataFrame()

    def _get_cytoband_cnv_summary(
        cnv_data: Dict[str, np.ndarray],
        chromosome: str,
        bin_width: int,
        sex_estimate: str,
    ) -> str:
        try:
            df = _analyze_cytoband_cnv(cnv_data, chromosome, bin_width, sex_estimate)
            if df.empty:
                return "No significant CNV changes detected"
            gains = df[df["cnv_state"] == "GAIN"]
            losses = df[df["cnv_state"] == "LOSS"]
            parts: List[str] = []
            if not gains.empty:
                parts.append(
                    "Gains: "
                    + ", ".join(
                        f"{r['name']} ({r['mean_cnv']:.2f})"
                        for _, r in gains.iterrows()
                    )
                )
            if not losses.empty:
                parts.append(
                    "Losses: "
                    + ", ".join(
                        f"{r['name']} ({r['mean_cnv']:.2f})"
                        for _, r in losses.iterrows()
                    )
                )
            return "\n".join(parts) if parts else "No significant CNV changes detected"
        except Exception:
            return "No CNV data available"

    def _compute_all_cytoband_df(
        cnv_data: Dict[str, np.ndarray], bin_width: int, sex_estimate: str
    ) -> pd.DataFrame:
        try:
            frames: List[pd.DataFrame] = []
            for chrom in natsort.natsorted(cnv_data.keys()):
                if (
                    chrom == "chrM"
                    or not isinstance(chrom, str)
                    or not chrom.startswith("chr")
                ):
                    continue
                # restrict to autosomes and sex chromosomes
                tail = chrom[3:]
                if not (tail.isdigit() or tail in ("X", "Y")):
                    continue
                df = _analyze_cytoband_cnv(cnv_data, chrom, bin_width, sex_estimate)
                if not df.empty:
                    frames.append(df)
            if frames:
                out = pd.concat(frames, ignore_index=True)
                if not out.empty:
                    # Natural chromosome sort: chr1..chr22, chrX, chrY
                    def _rank(label: Any) -> int:
                        try:
                            s = str(label)
                            if s.startswith("chr"):
                                s = s[3:]
                            mapping = {"X": 23, "Y": 24, "M": 25}
                            return int(s) if s.isdigit() else mapping.get(s, 1000)
                        except Exception:
                            return 1000

                    out["_chrom_rank"] = out["chrom"].map(_rank)
                    out = out.sort_values(["_chrom_rank", "start_pos"]).drop(
                        columns=["_chrom_rank"]
                    )
                return out
            return pd.DataFrame()
        except Exception:
            return pd.DataFrame()

    def _build_cyto_rows(df: pd.DataFrame) -> List[Dict[str, Any]]:
        rows: List[Dict[str, Any]] = []
        try:
            if df is None or df.empty:
                return rows
            for _, r in df.iterrows():
                if r.get("cnv_state") in ("GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"):
                    rows.append(
                        {
                            "chrom": str(r["chrom"]).replace("chr", ""),
                            "region": str(r["name"]).replace(f"{r['chrom']} ", ""),
                            "start_mb": f"{float(r['start_pos'])/1e6:.2f}",
                            "end_mb": f"{float(r['end_pos'])/1e6:.2f}",
                            "length_mb": f"{float(r['length'])/1e6:.2f}",
                            "mean_cnv": f"{float(r['mean_cnv']):.3f}",
                            "state": str(r["cnv_state"]),
                            "genes": ", ".join(
                                r.get("genes", [])
                                if isinstance(r.get("genes"), list)
                                else []
                            ),
                        }
                    )
        except Exception:
            return rows
        return rows

    def _render_cnv_from_state(state: Dict[str, Any]) -> None:
        try:
            cnv_map = state.get("cnv")
            cnv3_map = state.get("cnv3")
            if isinstance(cnv_map, dict) and "cnv" in cnv_map:
                cnv_map = cnv_map["cnv"]
            if isinstance(cnv3_map, dict) and "cnv" in cnv3_map:
                cnv3_map = cnv3_map["cnv"]
            if not cnv_map:
                return
            binw = state.get("cnv_dict", {}).get("bin_width", 1000000)
            selected = state.get("selected_chrom", "All")
            use_log = state.get("y_scale", "linear") == "log"
            raw_color_mode = state.get("color_mode", "chromosome")
            # normalize color mode to expected keys
            lval = str(raw_color_mode).strip().lower()
            if lval in ("chromosome", "chromosomes"):
                color_mode = "chromosome"
            elif lval in (
                "value",
                "up/down",
                "updown",
                "up_down",
                "updown ",
                "up down",
            ):
                color_mode = "value"
            else:
                color_mode = "chromosome"
            cnv_abs.options["yAxis"][0]["type"] = "log" if use_log else "value"
            cnv_abs.options["yAxis"][0]["logBase"] = 10 if use_log else None
            # ensure xAxis sane when switching modes
            if selected == "All":
                # Ensure full-range x-axes and clear any previous zoom constraints
                cnv_abs.options["xAxis"]["min"] = 0
                cnv_abs.options["xAxis"]["max"] = "dataMax"
                cnv_diff.options["xAxis"]["min"] = 0
                cnv_diff.options["xAxis"]["max"] = "dataMax"
                # also clear any leftover x-zoom from prior gene view and align diff zoom to full range
                try:
                    if (
                        isinstance(cnv_abs.options.get("dataZoom"), list)
                        and cnv_abs.options["dataZoom"]
                    ):
                        dz = cnv_abs.options["dataZoom"][0]
                        dz.pop("startValue", None)
                        dz.pop("endValue", None)
                        dz.update({"start": 0, "end": 100})
                    if (
                        isinstance(cnv_diff.options.get("dataZoom"), list)
                        and cnv_diff.options["dataZoom"]
                    ):
                        dz2 = cnv_diff.options["dataZoom"][0]
                        dz2.pop("startValue", None)
                        dz2.pop("endValue", None)
                        dz2.update({"start": 0, "end": 100})
                        # also clear any markLine/Area induced constraints by ensuring xAxis remains full dataMax
                        cnv_diff.options["xAxis"]["max"] = "dataMax"
                except Exception:
                    pass
            else:
                cnv_abs.options["xAxis"]["min"] = 0
                cnv_abs.options["xAxis"]["max"] = "dataMax"
                cnv_diff.options["xAxis"]["min"] = 0
                cnv_diff.options["xAxis"]["max"] = "dataMax"
            logging.debug(
                f"CNV render: selected={selected}, y_scale={state.get('y_scale')}, color_mode={state.get('color_mode')}"
            )
            # Absolute plot
            series_abs = []
            # Prepare chromosome partitions for labels/areas when viewing All
            chrom_bounds = []  # list of (name, start_bp, end_bp)
            chrom_offsets: Dict[str, float] = {}
            if selected == "All":
                total = 0
                for contig, cnv in natsort.natsorted(cnv_map.items()):
                    if contig == "chrM":
                        continue
                    x = (np.arange(len(cnv)) + total) * binw
                    pts = list(zip(x.tolist(), [float(v) for v in cnv]))
                    start_bp = total * binw
                    total += len(cnv)
                    end_bp = total * binw
                    chrom_offsets[contig] = start_bp
                    chrom_bounds.append((contig, start_bp, end_bp))
                    if color_mode == "chromosome":
                        series_abs.append(
                            {
                                "type": "scatter",
                                "name": contig,
                                "symbolSize": 3,
                                "data": pts,
                            }
                        )
                    else:
                        try:
                            autosome_vals = [
                                v
                                for k, arr in cnv_map.items()
                                if k.startswith("chr") and k[3:].isdigit()
                                for v in arr
                            ]
                            mean_val = (
                                float(np.mean(autosome_vals)) if autosome_vals else 2.0
                            )
                            std_val = (
                                float(np.std(autosome_vals)) if autosome_vals else 1.0
                            )
                        except Exception:
                            mean_val, std_val = 2.0, 1.0
                        high, low, norm = [], [], []
                        for xi, vi in pts:
                            z = (vi - mean_val) / std_val if std_val > 0 else 0.0
                            (high if z > 0.5 else low if z < -0.5 else norm).append(
                                [xi, vi]
                            )
                        if high:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"High {contig}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": "#007AFF"},
                                    "data": high,
                                }
                            )
                        if low:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Low {contig}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": "#FF3B30"},
                                    "data": low,
                                }
                            )
                        if norm:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Normal {contig}",
                                    "symbolSize": 2,
                                    "itemStyle": {"color": "#8E8E93"},
                                    "data": norm,
                                }
                            )
            else:
                cnv = cnv_map.get(selected)
                if cnv is not None:
                    x = (np.arange(len(cnv))) * binw
                    pts = list(zip(x.tolist(), [float(v) for v in cnv]))
                    if color_mode == "chromosome":
                        series_abs.append(
                            {
                                "type": "scatter",
                                "name": selected,
                                "symbolSize": 3,
                                "data": pts,
                            }
                        )
                    else:
                        expected = 2.0
                        try:
                            if selected in ("chrX", "chrY") and str(
                                state.get("xy", "")
                            ).upper().startswith("MALE"):
                                expected = 1.0
                        except Exception:
                            pass
                        vals = [v for _, v in pts]
                        std_val = float(np.std(vals)) if vals else 1.0
                        high, low, norm = [], [], []
                        for xi, vi in pts:
                            z = (vi - expected) / std_val if std_val > 0 else 0
                            (high if z > 0.5 else low if z < -0.5 else norm).append(
                                [xi, vi]
                            )
                        if high:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"High {selected}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": "#007AFF"},
                                    "data": high,
                                }
                            )
                        if low:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Low {selected}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": "#FF3B30"},
                                    "data": low,
                                }
                            )
                        if norm:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Normal {selected}",
                                    "symbolSize": 2,
                                    "itemStyle": {"color": "#8E8E93"},
                                    "data": norm,
                                }
                            )
            # Preserve highlight series (centromeres, cytobands) and replace data series only
            keep = [
                s
                for s in cnv_abs.options["series"]
                if s.get("name") in ("centromeres_highlight", "cytobands_highlight")
            ]
            cnv_abs.options["series"] = series_abs + keep
            # Build background chromosome areas and vertical labels when showing All
            try:
                if selected == "All" and chrom_bounds:
                    # Alternating shaded bands per chromosome for readability
                    areas_data = []
                    lines_data = []
                    for contig, start_bp, end_bp in chrom_bounds:
                        areas_data.append(
                            [{"xAxis": float(start_bp)}, {"xAxis": float(end_bp)}]
                        )
                        # Label at the boundary (end of chromosome region)
                        lines_data.append(
                            {
                                "xAxis": float(end_bp),
                                "lineStyle": {"type": "dashed", "color": "#A0A0A0"},
                                "label": {"show": True, "formatter": contig},
                            }
                        )

                    # Helper to set overlays by series name regardless of index
                    def _apply_overlays(
                        chart, band_areas, band_lines, centro_areas_list
                    ):
                        try:
                            idx_cyto = next(
                                (
                                    i
                                    for i, s in enumerate(chart.options["series"])
                                    if s.get("name") == "cytobands_highlight"
                                ),
                                None,
                            )
                            idx_centro = next(
                                (
                                    i
                                    for i, s in enumerate(chart.options["series"])
                                    if s.get("name") == "centromeres_highlight"
                                ),
                                None,
                            )
                            if idx_cyto is not None:
                                # Remove background shading per request; keep dashed vertical lines only
                                chart.options["series"][idx_cyto]["markArea"][
                                    "data"
                                ] = []
                                chart.options["series"][idx_cyto].setdefault(
                                    "markLine", {}
                                )
                                chart.options["series"][idx_cyto]["markLine"][
                                    "data"
                                ] = band_lines
                            # Do not show centromeres in All view
                            if idx_centro is not None:
                                chart.options["series"][idx_centro]["markArea"][
                                    "data"
                                ] = []
                        except Exception:
                            pass

                    _apply_overlays(cnv_abs, areas_data, lines_data, [])
                    _apply_overlays(cnv_diff, areas_data, lines_data, [])
                else:
                    # Clear overlays when focusing on a single chromosome
                    def _clear_overlays(chart):
                        try:
                            for s in chart.options["series"]:
                                if s.get("name") in (
                                    "cytobands_highlight",
                                    "centromeres_highlight",
                                ):
                                    if "markArea" in s and "data" in s["markArea"]:
                                        s["markArea"]["data"] = []
                                    if "markLine" in s and "data" in s["markLine"]:
                                        s["markLine"]["data"] = []
                        except Exception:
                            pass

                    _clear_overlays(cnv_abs)
                    _clear_overlays(cnv_diff)
                    # In single-chromosome view, draw centromeres, cytobands (by state), genes, and breakpoint candidates
                    if selected != "All":
                        try:
                            centro = _load_centromere_regions()
                            idx_centro_abs = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_abs.options["series"])
                                    if s.get("name") == "centromeres_highlight"
                                ),
                                None,
                            )
                            idx_centro_diff = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_diff.options["series"])
                                    if s.get("name") == "centromeres_highlight"
                                ),
                                None,
                            )
                            areas = []
                            for s, e, _n in centro.get(selected, []):
                                areas.append([{"xAxis": float(s)}, {"xAxis": float(e)}])
                            if idx_centro_abs is not None:
                                cnv_abs.options["series"][idx_centro_abs]["markArea"][
                                    "data"
                                ] = areas
                            if idx_centro_diff is not None:
                                cnv_diff.options["series"][idx_centro_diff]["markArea"][
                                    "data"
                                ] = areas
                        except Exception:
                            pass
                        # Cytobands colored by CNV state (from CNV3 values)
                        try:
                            idx_cyto_abs = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_abs.options["series"])
                                    if s.get("name") == "cytobands_highlight"
                                ),
                                None,
                            )
                            if (
                                idx_cyto_abs is not None
                                and cnv3_map
                                and selected in cnv3_map
                            ):
                                cyto_df = _load_cytobands_df()
                                bands = cyto_df[cyto_df["chrom"] == selected]
                                vals = (
                                    np.array(cnv3_map[selected])
                                    if isinstance(
                                        cnv3_map[selected], (list, np.ndarray)
                                    )
                                    else np.array([])
                                )
                                band_areas = []
                                for _, row in bands.iterrows():
                                    s_bp, e_bp = int(row["start"]), int(row["end"])
                                    s_bin = max(0, s_bp // binw)
                                    e_bin = min(len(vals) - 1, max(0, e_bp // binw))
                                    if len(vals) > 0 and e_bin >= s_bin:
                                        mean_val = float(
                                            np.mean(vals[s_bin : e_bin + 1])
                                        )
                                    else:
                                        mean_val = 0.0
                                    if mean_val > 0.5:
                                        color = "rgba(52, 199, 89, 0.12)"
                                    elif mean_val < -0.5:
                                        color = "rgba(255, 45, 85, 0.12)"
                                    else:
                                        color = "rgba(0, 0, 0, 0.03)"
                                    band_areas.append(
                                        [
                                            {
                                                "name": str(row["name"]),
                                                "xAxis": float(s_bp),
                                                "itemStyle": {"color": color},
                                                "label": {
                                                    "show": True,
                                                    "position": "insideTop",
                                                    "color": "#555",
                                                    "fontSize": 11,
                                                },
                                            },
                                            {"xAxis": float(e_bp)},
                                        ]
                                    )
                                cnv_abs.options["series"][idx_cyto_abs]["markArea"][
                                    "data"
                                ] = band_areas
                        except Exception:
                            pass
                        # Genes of interest and gene selector options
                        try:
                            gene_df = _load_gene_bed()
                            gchr = gene_df[gene_df["chrom"] == selected]
                            # limit to reduce clutter; still add many labels
                            gene_opts = {"All": "All"}
                            for _, gr in gchr.iterrows():
                                gene_opts[str(gr["gene"])] = str(gr["gene"])
                            try:
                                cnv_gene_select.set_options(gene_opts)
                            except Exception:
                                pass
                            # annotate genes on main series
                            if series_abs:
                                main = series_abs[0]
                                # attach gene regions via markArea on main series after replacement
                                mark = []
                                for _, gr in gchr.iterrows():
                                    mark.append(
                                        [
                                            {
                                                "name": str(gr["gene"]),
                                                "xAxis": float(gr["start"]),
                                                "label": {
                                                    "position": "insideTop",
                                                    "color": "#000",
                                                    "fontSize": 11,
                                                },
                                            },
                                            {"xAxis": float(gr["end"])},
                                        ]
                                    )
                                main.setdefault("markArea", {"data": []})
                                main["markArea"]["data"] = (
                                    main["markArea"]["data"] or []
                                ) + mark
                                # ensure series_abs[0] updated
                                series_abs[0] = main
                            # zoom to selected gene if any
                            sel_gene = launcher._cnv_state.setdefault(
                                str(sample_dir), {}
                            ).get("selected_gene", "All")
                            if sel_gene and sel_gene != "All":
                                row = gchr[gchr["gene"] == sel_gene]
                                if not row.empty:
                                    s_bp = int(row.iloc[0]["start"])
                                    e_bp = int(row.iloc[0]["end"])
                                    pad = int(0.1 * (e_bp - s_bp + 1))
                                    for chart in (cnv_abs, cnv_diff):
                                        try:
                                            chart.options["dataZoom"][0].update(
                                                {
                                                    "startValue": max(0, s_bp - pad),
                                                    "endValue": e_bp + pad,
                                                }
                                            )
                                        except Exception:
                                            pass
                        except Exception:
                            pass
                        # Breakpoint candidates as dashed vertical lines
                        try:
                            idx_cyto_abs = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_abs.options["series"])
                                    if s.get("name") == "cytobands_highlight"
                                ),
                                None,
                            )
                            if (
                                idx_cyto_abs is not None
                                and state.get("bp_array") is not None
                            ):
                                arr = state["bp_array"]
                                pos = [
                                    int(r["end"]) for r in arr if r["name"] == selected
                                ]
                                lines = [
                                    {
                                        "xAxis": float(p),
                                        "lineStyle": {
                                            "type": "dashed",
                                            "color": "#E0162B",
                                        },
                                    }
                                    for p in pos
                                ]
                                cnv_abs.options["series"][idx_cyto_abs].setdefault(
                                    "markLine", {"symbol": "none", "data": []}
                                )
                                cnv_abs.options["series"][idx_cyto_abs]["markLine"][
                                    "data"
                                ] = lines
                        except Exception:
                            pass
                # debug label removed
            except Exception:
                pass
            # Adaptive thinning based on current zoom and cap total points
            _thin_chart_series(cnv_abs, MAX_POINTS_PER_CHART)
            cnv_abs.update()
            # Difference plot
            if cnv3_map:
                series_diff = []
                if selected == "All":
                    total = 0
                    for contig, cnv in natsort.natsorted(cnv3_map.items()):
                        if contig == "chrM":
                            continue
                        x = (np.arange(len(cnv)) + total) * binw
                        pts = list(zip(x.tolist(), [float(v) for v in cnv]))
                        total += len(cnv)
                        series_diff.append(
                            {
                                "type": "scatter",
                                "name": contig,
                                "symbolSize": 3,
                                "data": pts,
                            }
                        )
                else:
                    cnv = cnv3_map.get(selected)
                    if cnv is not None:
                        x = (np.arange(len(cnv))) * binw
                        pts = list(zip(x.tolist(), [float(v) for v in cnv]))
                        series_diff.append(
                            {
                                "type": "scatter",
                                "name": selected,
                                "symbolSize": 3,
                                "data": pts,
                            }
                        )
                # Preserve highlight series by name
                try:
                    base_series = cnv_diff.options["series"]
                    keep = [
                        s
                        for s in base_series
                        if s.get("name")
                        in ("centromeres_highlight", "cytobands_highlight")
                    ]
                except Exception:
                    keep = []
                cnv_diff.options["series"] = series_diff + keep
                _thin_chart_series(cnv_diff, MAX_POINTS_PER_CHART)
                cnv_diff.update()

            # Cytoband CNV table update (whole-genome table with per-chromosome subsetting)
            try:
                selected = state.get("selected_chrom", "All")
                binw = state.get("cnv_dict", {}).get("bin_width", 1_000_000)
                sex_lbl = _sex_label(state.get("xy"))
                # Prefer difference map (CNV3) for calling; fall back to absolute
                if isinstance(cnv3_map, dict):
                    data = cnv3_map
                    source = "cnv3"
                else:
                    data = cnv_map if isinstance(cnv_map, dict) else None
                    source = "cnv"
                if data and binw:
                    cache_key = (
                        f"{source}:{state.get(source+'_m')}:{int(binw)}:{sex_lbl}"
                    )
                    if state.get("cyto_cache_key") != cache_key:
                        df_all = _compute_all_cytoband_df(data, int(binw), sex_lbl)
                        state["cyto_df_all"] = df_all
                        state["cyto_cache_key"] = cache_key
                    df_all = state.get("cyto_df_all")
                    if isinstance(df_all, pd.DataFrame) and not df_all.empty:
                        if selected and selected != "All":
                            df_show = df_all[df_all["chrom"] == selected]
                        else:
                            df_show = df_all
                        cyto_table.rows = _build_cyto_rows(df_show)
                        try:
                            cyto_table.update()
                        except Exception:
                            pass
                        if selected and selected != "All":
                            cyto_summary.set_text(
                                _get_cytoband_cnv_summary(
                                    data, selected, int(binw), sex_lbl
                                )
                            )
                        else:
                            cyto_summary.set_text(
                                f"Whole genome cytoband events: {len(cyto_table.rows)}"
                            )
                    else:
                        cyto_table.rows = []
                        try:
                            cyto_table.update()
                        except Exception:
                            pass
                        cyto_summary.set_text("No significant CNV changes detected")
                else:
                    cyto_summary.set_text("CNV data not available")
            except Exception:
                pass
        except Exception:
            pass

    def _refresh_cnv() -> None:
        try:
            if not sample_dir or not sample_dir.exists():
                return
            key = str(sample_dir)
            state = launcher._cnv_state.get(key, {})
            
            # Simple debouncing to prevent rapid successive calls
            import time
            current_time = time.time()
            last_refresh = state.get("_last_refresh", 0)
            if current_time - last_refresh < 0.1:  # 100ms debounce
                return
            state["_last_refresh"] = current_time
            # Sync state from UI controls in case events are not firing in this environment
            try:
                ui_changed = False
                ui_sel = getattr(cnv_chrom_select, "value", None)
                if ui_sel and ui_sel != state.get("selected_chrom"):
                    state["selected_chrom"] = ui_sel
                    ui_changed = True
                ui_scale = getattr(cnv_scale, "value", None)
                if ui_scale and ui_scale != state.get("y_scale"):
                    state["y_scale"] = ui_scale
                    ui_changed = True
                ui_bp = getattr(cnv_bp, "value", None)
                if ui_bp is not None:
                    desired = ui_bp == "show"
                    if desired != state.get("show_bp", True):
                        state["show_bp"] = desired
                        ui_changed = True
                ui_color = getattr(cnv_color, "value", None)
                if ui_color and ui_color != state.get("color_mode"):
                    state["color_mode"] = ui_color
                    ui_changed = True
            except Exception:
                pass
            cnv_npy = sample_dir / "CNV.npy"
            cnv3_npy = sample_dir / "CNV3.npy"
            cnv_dict_npy = sample_dir / "CNV_dict.npy"
            data_array_npy = sample_dir / "cnv_data_array.npy"
            xy_pkl = sample_dir / "XYestimate.pkl"
            changed = False
            if cnv_dict_npy.exists():
                m = cnv_dict_npy.stat().st_mtime
                if state.get("dict_m") != m:
                    state["cnv_dict"] = np.load(cnv_dict_npy, allow_pickle=True).item()
                    cnv_bin.set_text(
                        f"Bin Width: {state['cnv_dict'].get('bin_width', '--'):,}"
                    )
                    cnv_var.set_text(
                        f"Variance: {state['cnv_dict'].get('variance','--'):.3f}"
                        if isinstance(state["cnv_dict"].get("variance"), (int, float))
                        else "Variance: --"
                    )
                    state["dict_m"] = m
            if xy_pkl.exists():
                m = xy_pkl.stat().st_mtime
                if state.get("xy_m") != m:
                    try:
                        import pickle

                        with xy_pkl.open("rb") as f:
                            xy = pickle.load(f)
                        cnv_xy.set_text(f"Genetic Sex: {xy}")
                        state["xy"] = xy
                    except Exception:
                        pass
                    state["xy_m"] = m

            def _load_npy(path_key, file_path):
                m = file_path.stat().st_mtime
                if state.get(path_key + "_m") != m:
                    try:
                        state[path_key] = np.load(file_path, allow_pickle=True).item()
                    except Exception:
                        state[path_key] = None
                    state[path_key + "_m"] = m
                    return True
                return False

            if cnv_npy.exists() and _load_npy("cnv", cnv_npy):
                changed = True
            if cnv3_npy.exists() and _load_npy("cnv3", cnv3_npy):
                changed = True
            if state.get("cnv"):
                if changed:
                    cnv_status.set_text("Status: CNV data loaded")
                # Populate selector at least once, or when data changed
                if changed or not state.get("_chrom_opts_set"):
                    try:
                        cnv_map = state["cnv"].get("cnv", state["cnv"])
                        chrom_opts = {"All": "All"}
                        for contig in natsort.natsorted(cnv_map.keys()):
                            if contig != "chrM":
                                chrom_opts[contig] = contig
                        cnv_chrom_select.set_options(chrom_opts)
                        # Keep selected option if still valid
                        chrom_opts = {"All": "All"}
                        for contig in natsort.natsorted(cnv_map.keys()):
                            if contig != "chrM":
                                chrom_opts[contig] = contig
                        cnv_chrom_select.set_options(chrom_opts)
                        # Keep selected option if still valid
                        sel = launcher._cnv_state.setdefault(str(sample_dir), {}).get(
                            "selected_chrom", "All"
                        )
                        if sel not in chrom_opts:
                            sel = "All"
                        cnv_chrom_select.value = sel
                        try:
                            cnv_chrom_select.update()
                        except Exception:
                            pass
                        state["_chrom_opts_set"] = True
                    except Exception:
                        pass
                # Only re-render when data or UI state changed, or on first render
                if changed or ui_changed or not state.get("_rendered_once"):
                    _render_cnv_from_state(state)
                    state["_rendered_once"] = True
            # Breakpoint density overlay
            if data_array_npy.exists():
                try:
                    arr = np.load(data_array_npy, allow_pickle=True)
                    if hasattr(arr, "dtype") and "name" in arr.dtype.names:
                        binw = state.get("cnv_dict", {}).get("bin_width", 1000000)
                        state["bp_array"] = arr
                        dens: Dict[int, int] = {}
                        for r in arr:
                            end = int(r["end"])
                            bin_idx = end // binw
                            dens[bin_idx] = dens.get(bin_idx, 0) + 1
                        dens_pts = [[k * binw, v] for k, v in sorted(dens.items())]
                        # Only add density overlay in single-chromosome view; preserve existing diff series
                        selected = launcher._cnv_state.setdefault(
                            str(sample_dir), {}
                        ).get("selected_chrom", "All")
                        current_series = [
                            s
                            for s in cnv_diff.options["series"]
                            if s.get("name") != "Breakpoint Density"
                        ]
                        if selected != "All" and state.get("show_bp", True):
                            current_series.append(
                                {
                                    "type": "scatter",
                                    "name": "Breakpoint Density",
                                    "yAxisIndex": 1,
                                    "symbolSize": 8,
                                    "data": dens_pts,
                                }
                            )
                        cnv_diff.options["series"] = current_series
                        cnv_diff.update()
                except Exception:
                    pass
            launcher._cnv_state[key] = state
        except Exception:
            pass

    # Bind control events
    try:

        def _val(ev, default=None):
            return (
                getattr(ev, "value", None)
                if hasattr(ev, "value")
                else (getattr(ev, "args", None) or default)
            )

        def _on_chrom(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["selected_chrom"] = _val(ev, "All") or "All"
            logging.debug(f"CNV select changed -> {st['selected_chrom']}")
            # reset x zoom when switching scope
            try:
                for chart in (cnv_abs, cnv_diff):
                    if (
                        isinstance(chart.options.get("dataZoom"), list)
                        and chart.options["dataZoom"]
                    ):
                        chart.options["dataZoom"][0].pop("startValue", None)
                        chart.options["dataZoom"][0].pop("endValue", None)
            except Exception:
                pass
            # Trigger immediate refresh to update all UI elements
            _refresh_cnv()

        def _on_scale(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["y_scale"] = _val(ev, "linear") or "linear"
            # Trigger immediate refresh to update all UI elements
            _refresh_cnv()

        def _on_bp(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["show_bp"] = _val(ev, "show") == "show"
            # Trigger immediate refresh to update all UI elements
            _refresh_cnv()

        def _on_color(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            val = _val(ev, "chromosome") or "chromosome"
            # Accept either keys or labels from the toggle
            vlow = str(val).strip().lower()
            if vlow in ("chromosome", "chromosomes"):
                st["color_mode"] = "chromosome"
            elif vlow in ("value", "up/down", "updown", "up_down", "up down"):
                st["color_mode"] = "value"
            else:
                st["color_mode"] = "chromosome"
            logging.debug(f"CNV color mode -> {st['color_mode']} (raw={val})")
            # Trigger immediate refresh to update all UI elements
            _refresh_cnv()

        # Bind both native change and model-value updates for robustness
        cnv_chrom_select.on("change", _on_chrom)
        cnv_chrom_select.on("update:model-value", _on_chrom)

        # Gene selection zoom
        def _on_gene(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["selected_gene"] = _val(ev, "All") or "All"
            # Trigger immediate refresh to update all UI elements
            _refresh_cnv()

        cnv_gene_select.on("change", _on_gene)
        cnv_gene_select.on("update:model-value", _on_gene)
        cnv_scale.on("change", _on_scale)
        cnv_scale.on("update:model-value", _on_scale)
        cnv_bp.on("change", _on_bp)
        cnv_bp.on("update:model-value", _on_bp)
        cnv_color.on("change", _on_color)
        cnv_color.on("update:model-value", _on_color)
    except Exception:
        pass

    ui.timer(30.0, _refresh_cnv, active=True)
