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
    """
    with ui.card().classes('w-full'):
        ui.label('🧬 Copy Number Variation (CNV)').classes('text-lg font-semibold mb-2')
        with ui.row().classes('w-full items-center justify-between mb-2'):
            with ui.column().classes('gap-1'):
                cnv_status = ui.label('Status: Awaiting Data').classes('text-gray-600')
                cnv_xy = ui.label('Genetic Sex: --').classes('text-gray-600')
            with ui.column().classes('gap-1 items-end'):
                cnv_bin = ui.label('Bin Width: --').classes('text-sm text-gray-600')
                cnv_var = ui.label('Variance: --').classes('text-sm text-gray-600')
        # Controls similar to original module
        with ui.row().classes('gap-4 items-center mb-2'):
            ui.label('Select Chromosome').classes('text-sm')
            cnv_chrom_select = ui.select(options={'All': 'All'}, value='All').style('width: 160px')
            ui.label('Select Gene').classes('text-sm ml-4')
            cnv_gene_select = ui.select(options={'All': 'All'}, value='All').style('width: 200px')
            ui.label('Color by').classes('text-sm ml-4')
            cnv_color = ui.toggle(options={'chromosome': 'Chromosome', 'value': 'Up/Down'}, value='chromosome').classes('mt-1')
            ui.label('Y-axis scale').classes('text-sm ml-4')
            cnv_scale = ui.toggle(options={'linear': 'Linear', 'log': 'Log'}, value='linear').classes('mt-1')
            ui.label('Breakpoints').classes('text-sm ml-4')
            cnv_bp = ui.toggle(options={'hide': 'Hide', 'show': 'Show'}, value='show').classes('mt-1')
        cnv_abs = ui.echart({
            'backgroundColor': 'transparent',
            'title': {'text': 'CNV Scatter Plot', 'left': 'center', 'top': 10},
            'grid': {'left': '5%', 'right': '5%', 'bottom': '10%', 'top': '20%', 'containLabel': True},
            'tooltip': {'trigger': 'axis'},
            'xAxis': {'type': 'value', 'max': 'dataMax'},
            'yAxis': [{'type': 'value', 'name': 'Ploidy'}, {'type': 'value', 'name': 'Breakpoint density', 'position': 'right'}],
            'dataZoom': [{'type': 'slider', 'xAxisIndex': [0]}, {'type': 'slider', 'yAxisIndex': [0,1], 'right': 20}],
            'series': [
                {'type': 'scatter', 'name': 'CNV', 'symbolSize': 3, 'data': []},
                {'type': 'scatter', 'name': 'centromeres_highlight', 'data': [], 'symbolSize': 3, 'markArea': {'itemStyle': {'color': 'rgba(135, 206, 250, 0.4)'}, 'data': []}},
                {'type': 'scatter', 'name': 'cytobands_highlight', 'data': [], 'symbolSize': 3, 'markArea': {'itemStyle': {'color': 'rgba(200, 200, 200, 0.4)'}, 'data': []}, 'markLine': {'symbol': 'none', 'data': []}},
            ],
        }).classes('w-full h-72')
        cnv_diff = ui.echart({
            'backgroundColor': 'transparent',
            'title': {'text': 'Difference Plot', 'left': 'center', 'top': 10},
            'grid': {'left': '5%', 'right': '5%', 'bottom': '10%', 'top': '20%', 'containLabel': True},
            'tooltip': {'trigger': 'axis'},
            'xAxis': {'type': 'value', 'max': 'dataMax'},
            'yAxis': [{'type': 'value', 'name': 'Relative'}, {'type': 'value', 'name': 'Breakpoint density', 'position': 'right'}],
            'dataZoom': [{'type': 'slider', 'xAxisIndex': [0]}, {'type': 'slider', 'yAxisIndex': [0,1], 'right': 20}],
            'series': [
                {'type': 'scatter', 'name': 'CNV Δ', 'symbolSize': 3, 'data': []},
                {'type': 'scatter', 'name': 'centromeres_highlight', 'data': [], 'symbolSize': 3, 'markArea': {'itemStyle': {'color': 'rgba(135, 206, 250, 0.4)'}, 'data': []}},
                {'type': 'scatter', 'name': 'cytobands_highlight', 'data': [], 'symbolSize': 3, 'markArea': {'itemStyle': {'color': 'rgba(200, 200, 200, 0.4)'}, 'data': []}, 'markLine': {'symbol': 'none', 'data': []}},
            ],
        }).classes('w-full h-72')

    @lru_cache(maxsize=1)
    def _load_centromere_regions() -> Dict[str, List[Tuple[int, int, str]]]:
        """Load centromere/satellite regions from packaged resources.
        Returns mapping: chrom -> list of (start_bp, end_bp, name).
        """
        regions: Dict[str, List[Tuple[int, int, str]]] = {}
        try:
            res_path = importlib_resources.files('littlejohn.resources') / 'cenSatRegions.bed'
            with res_path.open('r') as fh:
                for line in fh:
                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue
                    chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
                    regions.setdefault(chrom, []).append((start, end, name))
        except Exception:
            pass
        return regions

    @lru_cache(maxsize=1)
    def _load_cytobands_df() -> pd.DataFrame:
        try:
            res_path = importlib_resources.files('littlejohn.resources') / 'cytoBand.txt'
            df = pd.read_csv(res_path, sep='\t', header=None,
                             names=['chrom','start','end','name','stain'])
            return df
        except Exception:
            return pd.DataFrame(columns=['chrom','start','end','name','stain'])

    @lru_cache(maxsize=1)
    def _load_gene_bed() -> pd.DataFrame:
        # Use a compact gene list; fall back to all_genes.bed
        for fname in ['unique_genes.bed', 'all_genes.bed']:
            try:
                res_path = importlib_resources.files('littlejohn.resources') / fname
                if res_path.exists():
                    return pd.read_csv(res_path, sep='\t', header=None,
                                       names=['chrom','start','end','gene'])
            except Exception:
                continue
        return pd.DataFrame(columns=['chrom','start','end','gene'])

    def _render_cnv_from_state(state: Dict[str, Any]) -> None:
        try:
            cnv_map = state.get('cnv')
            cnv3_map = state.get('cnv3')
            if isinstance(cnv_map, dict) and 'cnv' in cnv_map:
                cnv_map = cnv_map['cnv']
            if isinstance(cnv3_map, dict) and 'cnv' in cnv3_map:
                cnv3_map = cnv3_map['cnv']
            if not cnv_map:
                return
            binw = state.get('cnv_dict', {}).get('bin_width', 1000000)
            selected = state.get('selected_chrom', 'All')
            use_log = state.get('y_scale', 'linear') == 'log'
            raw_color_mode = state.get('color_mode', 'chromosome')
            # normalize color mode to expected keys
            lval = str(raw_color_mode).strip().lower()
            if lval in ('chromosome', 'chromosomes'):
                color_mode = 'chromosome'
            elif lval in ('value', 'up/down', 'updown', 'up_down', 'updown ', 'up down'):
                color_mode = 'value'
            else:
                color_mode = 'chromosome'
            cnv_abs.options['yAxis'][0]['type'] = 'log' if use_log else 'value'
            cnv_abs.options['yAxis'][0]['logBase'] = 10 if use_log else None
            # ensure xAxis sane when switching modes
            if selected == 'All':
                # Ensure full-range x-axes and clear any previous zoom constraints
                cnv_abs.options['xAxis']['min'] = 0
                cnv_abs.options['xAxis']['max'] = 'dataMax'
                cnv_diff.options['xAxis']['min'] = 0
                cnv_diff.options['xAxis']['max'] = 'dataMax'
                # also clear any leftover x-zoom from prior gene view and align diff zoom to full range
                try:
                    if isinstance(cnv_abs.options.get('dataZoom'), list) and cnv_abs.options['dataZoom']:
                        dz = cnv_abs.options['dataZoom'][0]
                        dz.pop('startValue', None)
                        dz.pop('endValue', None)
                        dz.update({'start': 0, 'end': 100})
                    if isinstance(cnv_diff.options.get('dataZoom'), list) and cnv_diff.options['dataZoom']:
                        dz2 = cnv_diff.options['dataZoom'][0]
                        dz2.pop('startValue', None)
                        dz2.pop('endValue', None)
                        dz2.update({'start': 0, 'end': 100})
                        # also clear any markLine/Area induced constraints by ensuring xAxis remains full dataMax
                        cnv_diff.options['xAxis']['max'] = 'dataMax'
                except Exception:
                    pass
            else:
                cnv_abs.options['xAxis']['min'] = 0
                cnv_abs.options['xAxis']['max'] = 'dataMax'
                cnv_diff.options['xAxis']['min'] = 0
                cnv_diff.options['xAxis']['max'] = 'dataMax'
            logging.debug(f"CNV render: selected={selected}, y_scale={state.get('y_scale')}, color_mode={state.get('color_mode')}")
            # Absolute plot
            series_abs = []
            # Prepare chromosome partitions for labels/areas when viewing All
            chrom_bounds = []  # list of (name, start_bp, end_bp)
            chrom_offsets: Dict[str, float] = {}
            if selected == 'All':
                total = 0
                for contig, cnv in natsort.natsorted(cnv_map.items()):
                    if contig == 'chrM':
                        continue
                    x = (np.arange(len(cnv)) + total) * binw
                    pts = list(zip(x.tolist(), [float(v) for v in cnv]))
                    start_bp = total * binw
                    total += len(cnv)
                    end_bp = total * binw
                    chrom_offsets[contig] = start_bp
                    chrom_bounds.append((contig, start_bp, end_bp))
                    if color_mode == 'chromosome':
                        series_abs.append({
                            'type': 'scatter', 'name': contig, 'symbolSize': 3,
                            'data': pts[::max(1, len(pts)//10000 or 1)]
                        })
                    else:
                        try:
                            autosome_vals = [v for k, arr in cnv_map.items() if k.startswith('chr') and k[3:].isdigit() for v in arr]
                            mean_val = float(np.mean(autosome_vals)) if autosome_vals else 2.0
                            std_val = float(np.std(autosome_vals)) if autosome_vals else 1.0
                        except Exception:
                            mean_val, std_val = 2.0, 1.0
                        high, low, norm = [], [], []
                        for xi, vi in pts:
                            z = (vi - mean_val) / std_val if std_val > 0 else 0.0
                            (high if z > 0.5 else low if z < -0.5 else norm).append([xi, vi])
                        if high:
                            series_abs.append({'type': 'scatter', 'name': f'High {contig}', 'symbolSize': 4, 'itemStyle': {'color': '#007AFF'}, 'data': high[::max(1, len(high)//10000 or 1)]})
                        if low:
                            series_abs.append({'type': 'scatter', 'name': f'Low {contig}', 'symbolSize': 4, 'itemStyle': {'color': '#FF3B30'}, 'data': low[::max(1, len(low)//10000 or 1)]})
                        if norm:
                            series_abs.append({'type': 'scatter', 'name': f'Normal {contig}', 'symbolSize': 2, 'itemStyle': {'color': '#8E8E93'}, 'data': norm[::max(1, len(norm)//10000 or 1)]})
            else:
                cnv = cnv_map.get(selected)
                if cnv is not None:
                    x = (np.arange(len(cnv))) * binw
                    pts = list(zip(x.tolist(), [float(v) for v in cnv]))
                    if color_mode == 'chromosome':
                        series_abs.append({'type': 'scatter', 'name': selected, 'symbolSize': 3, 'data': pts[::max(1, len(pts)//10000 or 1)]})
                    else:
                        expected = 2.0
                        try:
                            if selected in ('chrX', 'chrY') and str(state.get('xy','')).upper().startswith('MALE'):
                                expected = 1.0
                        except Exception:
                            pass
                        vals = [v for _, v in pts]
                        std_val = float(np.std(vals)) if vals else 1.0
                        high, low, norm = [], [], []
                        for xi, vi in pts:
                            z = (vi - expected) / std_val if std_val > 0 else 0
                            (high if z > 0.5 else low if z < -0.5 else norm).append([xi, vi])
                        if high:
                            series_abs.append({'type': 'scatter', 'name': f'High {selected}', 'symbolSize': 4, 'itemStyle': {'color': '#007AFF'}, 'data': high})
                        if low:
                            series_abs.append({'type': 'scatter', 'name': f'Low {selected}', 'symbolSize': 4, 'itemStyle': {'color': '#FF3B30'}, 'data': low})
                        if norm:
                            series_abs.append({'type': 'scatter', 'name': f'Normal {selected}', 'symbolSize': 2, 'itemStyle': {'color': '#8E8E93'}, 'data': norm})
            # Preserve highlight series (centromeres, cytobands) and replace data series only
            keep = [s for s in cnv_abs.options['series'] if s.get('name') in ('centromeres_highlight', 'cytobands_highlight')]
            cnv_abs.options['series'] = series_abs + keep
            # Build background chromosome areas and vertical labels when showing All
            try:
                if selected == 'All' and chrom_bounds:
                    # Alternating shaded bands per chromosome for readability
                    areas_data = []
                    lines_data = []
                    for contig, start_bp, end_bp in chrom_bounds:
                        areas_data.append([
                            {'xAxis': float(start_bp)},
                            {'xAxis': float(end_bp)}
                        ])
                        # Label at the boundary (end of chromosome region)
                        lines_data.append({
                            'xAxis': float(end_bp),
                            'lineStyle': {'type': 'dashed', 'color': '#A0A0A0'},
                            'label': {'show': True, 'formatter': contig}
                        })
                    # Helper to set overlays by series name regardless of index
                    def _apply_overlays(chart, band_areas, band_lines, centro_areas_list):
                        try:
                            idx_cyto = next((i for i, s in enumerate(chart.options['series']) if s.get('name') == 'cytobands_highlight'), None)
                            idx_centro = next((i for i, s in enumerate(chart.options['series']) if s.get('name') == 'centromeres_highlight'), None)
                            if idx_cyto is not None:
                                # Remove background shading per request; keep dashed vertical lines only
                                chart.options['series'][idx_cyto]['markArea']['data'] = []
                                chart.options['series'][idx_cyto].setdefault('markLine', {})
                                chart.options['series'][idx_cyto]['markLine']['data'] = band_lines
                            # Do not show centromeres in All view
                            if idx_centro is not None:
                                chart.options['series'][idx_centro]['markArea']['data'] = []
                        except Exception:
                            pass
                    _apply_overlays(cnv_abs, areas_data, lines_data, [])
                    _apply_overlays(cnv_diff, areas_data, lines_data, [])
                else:
                    # Clear overlays when focusing on a single chromosome
                    def _clear_overlays(chart):
                        try:
                            for s in chart.options['series']:
                                if s.get('name') in ('cytobands_highlight', 'centromeres_highlight'):
                                    if 'markArea' in s and 'data' in s['markArea']:
                                        s['markArea']['data'] = []
                                    if 'markLine' in s and 'data' in s['markLine']:
                                        s['markLine']['data'] = []
                        except Exception:
                            pass
                    _clear_overlays(cnv_abs)
                    _clear_overlays(cnv_diff)
                    # In single-chromosome view, draw centromeres, cytobands (by state), genes, and breakpoint candidates
                    if selected != 'All':
                        try:
                            centro = _load_centromere_regions()
                            idx_centro_abs = next((i for i, s in enumerate(cnv_abs.options['series']) if s.get('name') == 'centromeres_highlight'), None)
                            idx_centro_diff = next((i for i, s in enumerate(cnv_diff.options['series']) if s.get('name') == 'centromeres_highlight'), None)
                            areas = []
                            for s, e, _n in centro.get(selected, []):
                                areas.append([{'xAxis': float(s)}, {'xAxis': float(e)}])
                            if idx_centro_abs is not None:
                                cnv_abs.options['series'][idx_centro_abs]['markArea']['data'] = areas
                            if idx_centro_diff is not None:
                                cnv_diff.options['series'][idx_centro_diff]['markArea']['data'] = areas
                        except Exception:
                            pass
                        # Cytobands colored by CNV state (from CNV3 values)
                        try:
                            idx_cyto_abs = next((i for i, s in enumerate(cnv_abs.options['series']) if s.get('name') == 'cytobands_highlight'), None)
                            if idx_cyto_abs is not None and cnv3_map and selected in cnv3_map:
                                cyto_df = _load_cytobands_df()
                                bands = cyto_df[cyto_df['chrom'] == selected]
                                vals = np.array(cnv3_map[selected]) if isinstance(cnv3_map[selected], (list, np.ndarray)) else np.array([])
                                band_areas = []
                                for _, row in bands.iterrows():
                                    s_bp, e_bp = int(row['start']), int(row['end'])
                                    s_bin = max(0, s_bp // binw)
                                    e_bin = min(len(vals)-1, max(0, e_bp // binw))
                                    if len(vals) > 0 and e_bin >= s_bin:
                                        mean_val = float(np.mean(vals[s_bin:e_bin+1]))
                                    else:
                                        mean_val = 0.0
                                    if mean_val > 0.5:
                                        color = 'rgba(52, 199, 89, 0.12)'
                                    elif mean_val < -0.5:
                                        color = 'rgba(255, 45, 85, 0.12)'
                                    else:
                                        color = 'rgba(0, 0, 0, 0.03)'
                                    band_areas.append([{'name': str(row['name']), 'xAxis': float(s_bp), 'itemStyle': {'color': color}, 'label': {'show': True, 'position': 'insideTop', 'color': '#555', 'fontSize': 11}},
                                                      {'xAxis': float(e_bp)}])
                                cnv_abs.options['series'][idx_cyto_abs]['markArea']['data'] = band_areas
                        except Exception:
                            pass
                        # Genes of interest and gene selector options
                        try:
                            gene_df = _load_gene_bed()
                            gchr = gene_df[gene_df['chrom'] == selected]
                            # limit to reduce clutter; still add many labels
                            gene_opts = {'All': 'All'}
                            for _, gr in gchr.iterrows():
                                gene_opts[str(gr['gene'])] = str(gr['gene'])
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
                                    mark.append([{'name': str(gr['gene']), 'xAxis': float(gr['start']), 'label': {'position': 'insideTop', 'color': '#000', 'fontSize': 11}}, {'xAxis': float(gr['end'])}])
                                main.setdefault('markArea', {'data': []})
                                main['markArea']['data'] = (main['markArea']['data'] or []) + mark
                                # ensure series_abs[0] updated
                                series_abs[0] = main
                            # zoom to selected gene if any
                            sel_gene = launcher._cnv_state.setdefault(str(sample_dir), {}).get('selected_gene', 'All')
                            if sel_gene and sel_gene != 'All':
                                row = gchr[gchr['gene'] == sel_gene]
                                if not row.empty:
                                    s_bp = int(row.iloc[0]['start'])
                                    e_bp = int(row.iloc[0]['end'])
                                    pad = int(0.1 * (e_bp - s_bp + 1))
                                    for chart in (cnv_abs, cnv_diff):
                                        try:
                                            chart.options['dataZoom'][0].update({'startValue': max(0, s_bp - pad), 'endValue': e_bp + pad})
                                        except Exception:
                                            pass
                        except Exception:
                            pass
                        # Breakpoint candidates as dashed vertical lines
                        try:
                            idx_cyto_abs = next((i for i, s in enumerate(cnv_abs.options['series']) if s.get('name') == 'cytobands_highlight'), None)
                            if idx_cyto_abs is not None and state.get('bp_array') is not None:
                                arr = state['bp_array']
                                pos = [int(r['end']) for r in arr if r['name'] == selected]
                                lines = [{'xAxis': float(p), 'lineStyle': {'type': 'dashed', 'color': '#E0162B'}} for p in pos]
                                cnv_abs.options['series'][idx_cyto_abs].setdefault('markLine', {'symbol': 'none', 'data': []})
                                cnv_abs.options['series'][idx_cyto_abs]['markLine']['data'] = lines
                        except Exception:
                            pass
                # debug label removed
            except Exception:
                pass
            cnv_abs.update()
            # Difference plot
            if cnv3_map:
                series_diff = []
                if selected == 'All':
                    total = 0
                    for contig, cnv in natsort.natsorted(cnv3_map.items()):
                        if contig == 'chrM':
                            continue
                        x = (np.arange(len(cnv)) + total) * binw
                        pts = list(zip(x.tolist(), [float(v) for v in cnv]))
                        total += len(cnv)
                        series_diff.append({'type': 'scatter', 'name': contig, 'symbolSize': 3, 'data': pts[::max(1, len(pts)//10000 or 1)]})
                else:
                    cnv = cnv3_map.get(selected)
                    if cnv is not None:
                        x = (np.arange(len(cnv))) * binw
                        pts = list(zip(x.tolist(), [float(v) for v in cnv]))
                        series_diff.append({'type': 'scatter', 'name': selected, 'symbolSize': 3, 'data': pts})
                # Preserve highlight series by name
                try:
                    base_series = cnv_diff.options['series']
                    keep = [s for s in base_series if s.get('name') in ('centromeres_highlight', 'cytobands_highlight')]
                except Exception:
                    keep = []
                cnv_diff.options['series'] = series_diff + keep
                cnv_diff.update()
        except Exception:
            pass

    def _refresh_cnv() -> None:
        try:
            if not sample_dir or not sample_dir.exists():
                return
            key = str(sample_dir)
            state = launcher._cnv_state.get(key, {})
            # Sync state from UI controls in case events are not firing in this environment
            try:
                ui_changed = False
                ui_sel = getattr(cnv_chrom_select, 'value', None)
                if ui_sel and ui_sel != state.get('selected_chrom'):
                    state['selected_chrom'] = ui_sel
                    ui_changed = True
                ui_scale = getattr(cnv_scale, 'value', None)
                if ui_scale and ui_scale != state.get('y_scale'):
                    state['y_scale'] = ui_scale
                    ui_changed = True
                ui_bp = getattr(cnv_bp, 'value', None)
                if ui_bp is not None:
                    desired = (ui_bp == 'show')
                    if desired != state.get('show_bp', True):
                        state['show_bp'] = desired
                        ui_changed = True
                ui_color = getattr(cnv_color, 'value', None)
                if ui_color and ui_color != state.get('color_mode'):
                    state['color_mode'] = ui_color
                    ui_changed = True
            except Exception:
                pass
            cnv_npy = sample_dir / 'CNV.npy'
            cnv3_npy = sample_dir / 'CNV3.npy'
            cnv_dict_npy = sample_dir / 'CNV_dict.npy'
            data_array_npy = sample_dir / 'cnv_data_array.npy'
            xy_pkl = sample_dir / 'XYestimate.pkl'
            changed = False
            if cnv_dict_npy.exists():
                m = cnv_dict_npy.stat().st_mtime
                if state.get('dict_m') != m:
                    state['cnv_dict'] = np.load(cnv_dict_npy, allow_pickle=True).item()
                    cnv_bin.set_text(f"Bin Width: {state['cnv_dict'].get('bin_width', '--'):,}")
                    cnv_var.set_text(f"Variance: {state['cnv_dict'].get('variance','--'):.3f}" if isinstance(state['cnv_dict'].get('variance'), (int,float)) else 'Variance: --')
                    state['dict_m'] = m
            if xy_pkl.exists():
                m = xy_pkl.stat().st_mtime
                if state.get('xy_m') != m:
                    try:
                        import pickle
                        with xy_pkl.open('rb') as f:
                            xy = pickle.load(f)
                        cnv_xy.set_text(f'Genetic Sex: {xy}')
                        state['xy'] = xy
                    except Exception:
                        pass
                    state['xy_m'] = m
            def _load_npy(path_key, file_path):
                m = file_path.stat().st_mtime
                if state.get(path_key+'_m') != m:
                    try:
                        state[path_key] = np.load(file_path, allow_pickle=True).item()
                    except Exception:
                        state[path_key] = None
                    state[path_key+'_m'] = m
                    return True
                return False
            if cnv_npy.exists() and _load_npy('cnv', cnv_npy):
                changed = True
            if cnv3_npy.exists() and _load_npy('cnv3', cnv3_npy):
                changed = True
            if state.get('cnv'):
                if changed:
                    cnv_status.set_text('Status: CNV data loaded')
                # Populate selector at least once, or when data changed
                if changed or not state.get('_chrom_opts_set'):
                    try:
                        cnv_map = state['cnv'].get('cnv', state['cnv'])
                        chrom_opts = {'All': 'All'}
                        for contig in natsort.natsorted(cnv_map.keys()):
                            if contig != 'chrM':
                                chrom_opts[contig] = contig
                        cnv_chrom_select.set_options(chrom_opts)
                        # Keep selected option if still valid
                        chrom_opts = {'All': 'All'}
                        for contig in natsort.natsorted(cnv_map.keys()):
                            if contig != 'chrM':
                                chrom_opts[contig] = contig
                        cnv_chrom_select.set_options(chrom_opts)
                        # Keep selected option if still valid
                        sel = launcher._cnv_state.setdefault(str(sample_dir), {}).get('selected_chrom', 'All')
                        if sel not in chrom_opts:
                            sel = 'All'
                        cnv_chrom_select.value = sel
                        try:
                            cnv_chrom_select.update()
                        except Exception:
                            pass
                        state['_chrom_opts_set'] = True
                    except Exception:
                        pass
                # Only re-render when data or UI state changed, or on first render
                if changed or ui_changed or not state.get('_rendered_once'):
                    _render_cnv_from_state(state)
                    state['_rendered_once'] = True
            # Breakpoint density overlay
            if data_array_npy.exists():
                try:
                    arr = np.load(data_array_npy, allow_pickle=True)
                    if hasattr(arr, 'dtype') and 'name' in arr.dtype.names:
                        binw = state.get('cnv_dict',{}).get('bin_width', 1000000)
                        state['bp_array'] = arr
                        dens: Dict[int,int] = {}
                        for r in arr:
                            end = int(r['end'])
                            bin_idx = end // binw
                            dens[bin_idx] = dens.get(bin_idx, 0) + 1
                        dens_pts = [[k*binw, v] for k, v in sorted(dens.items())]
                        # Only add density overlay in single-chromosome view; preserve existing diff series
                        selected = launcher._cnv_state.setdefault(str(sample_dir), {}).get('selected_chrom', 'All')
                        current_series = [s for s in cnv_diff.options['series'] if s.get('name') != 'Breakpoint Density']
                        if selected != 'All' and state.get('show_bp', True):
                            current_series.append({'type': 'scatter', 'name': 'Breakpoint Density', 'yAxisIndex': 1, 'symbolSize': 8, 'data': dens_pts})
                        cnv_diff.options['series'] = current_series
                        cnv_diff.update()
                except Exception:
                    pass
            launcher._cnv_state[key] = state
        except Exception:
            pass

    # Bind control events
    try:
        def _val(ev, default=None):
            return getattr(ev, 'value', None) if hasattr(ev, 'value') else (getattr(ev, 'args', None) or default)

        def _on_chrom(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st['selected_chrom'] = _val(ev, 'All') or 'All'
            logging.debug(f"CNV select changed -> {st['selected_chrom']}")
            # reset x zoom when switching scope
            try:
                for chart in (cnv_abs, cnv_diff):
                    if isinstance(chart.options.get('dataZoom'), list) and chart.options['dataZoom']:
                        chart.options['dataZoom'][0].pop('startValue', None)
                        chart.options['dataZoom'][0].pop('endValue', None)
            except Exception:
                pass
            _render_cnv_from_state(st)

        def _on_scale(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st['y_scale'] = _val(ev, 'linear') or 'linear'
            _render_cnv_from_state(st)

        def _on_bp(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st['show_bp'] = (_val(ev, 'show') == 'show')
            _render_cnv_from_state(st)

        def _on_color(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            val = _val(ev, 'chromosome') or 'chromosome'
            # Accept either keys or labels from the toggle
            vlow = str(val).strip().lower()
            if vlow in ('chromosome', 'chromosomes'):
                st['color_mode'] = 'chromosome'
            elif vlow in ('value', 'up/down', 'updown', 'up_down', 'up down'):
                st['color_mode'] = 'value'
            else:
                st['color_mode'] = 'chromosome'
            logging.debug(f"CNV color mode -> {st['color_mode']} (raw={val})")
            _render_cnv_from_state(st)

        # Bind both native change and model-value updates for robustness
        cnv_chrom_select.on('change', _on_chrom)
        cnv_chrom_select.on('update:model-value', _on_chrom)
        # Gene selection zoom
        def _on_gene(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st['selected_gene'] = _val(ev, 'All') or 'All'
            _render_cnv_from_state(st)
        cnv_gene_select.on('change', _on_gene)
        cnv_gene_select.on('update:model-value', _on_gene)
        cnv_scale.on('change', _on_scale)
        cnv_scale.on('update:model-value', _on_scale)
        cnv_bp.on('change', _on_bp)
        cnv_bp.on('update:model-value', _on_bp)
        cnv_color.on('change', _on_color)
        cnv_color.on('update:model-value', _on_color)
    except Exception:
        pass

    ui.timer(30.0, _refresh_cnv, active=True)


