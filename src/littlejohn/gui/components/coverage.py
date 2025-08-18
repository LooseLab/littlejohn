from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict, List

import natsort
import numpy as np
import pandas as pd

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None


def add_coverage_section(launcher: Any, sample_dir: Path) -> None:
    """Build the Coverage UI section and attach refresh timers.

    This mirrors the existing inline implementation but lives in a reusable module.
    """
    # Summary row (quality + metrics)
    with ui.card().classes('w-full'):
        ui.label('🧫 Coverage').classes('text-lg font-semibold mb-2')
        with ui.row().classes('w-full items-center justify-between mb-2'):
            with ui.column().classes('gap-1'):
                ui.label('Coverage Analysis').classes('text-sm font-medium')
                with ui.row().classes('items-center gap-2'):
                    cov_quality_label = ui.label('Quality: --').classes('text-gray-600 font-medium')
                    cov_quality_badge = ui.label('--x').classes('px-2 py-1 rounded bg-gray-100 text-gray-600')
            with ui.column().classes('gap-1 items-end'):
                cov_global_lbl = ui.label('Global Estimated Coverage: --x').classes('text-sm text-gray-600')
                cov_target_lbl = ui.label('Targets Estimated Coverage: --x').classes('text-sm text-gray-600')
                cov_enrich_lbl = ui.label('Estimated enrichment: --x').classes('text-sm text-gray-600')
        with ui.grid(columns=2).classes('w-full gap-4'):
            # Per Chromosome Coverage (bar)
            echart_chr_cov = ui.echart({
                'backgroundColor': 'transparent',
                'title': {'text': 'Per Chromosome Coverage', 'left': 'center', 'top': 10, 'textStyle': {'fontSize': 16, 'color': '#000'}},
                'tooltip': {'trigger': 'axis'},
                'grid': {'left': '5%', 'right': '5%', 'bottom': '10%', 'top': '20%', 'containLabel': True},
                'xAxis': {'type': 'category', 'data': [], 'axisLabel': {'rotate': 45}},
                'yAxis': {'type': 'value', 'name': 'Coverage (x)'},
                'series': []
            }).classes('w-full h-64')
            # Per Chromosome Target Coverage (scatter)
            echart_target_cov = ui.echart({
                'backgroundColor': 'transparent',
                'title': {'text': 'Per Chromosome Target Coverage', 'left': 'center', 'top': 10, 'textStyle': {'fontSize': 16, 'color': '#000'}},
                'legend': {'data': ['Off Target', 'On Target'], 'top': 45},
                'tooltip': {'trigger': 'axis'},
                'grid': {'left': '5%', 'right': '5%', 'bottom': '10%', 'top': '25%', 'containLabel': True},
                'xAxis': {'type': 'category', 'data': [], 'axisLabel': {'rotate': 45}},
                'yAxis': {'type': 'value', 'name': 'Coverage (x)'},
                'series': []
            }).classes('w-full h-64')

    with ui.card().classes('w-full'):
        ui.label('📈 Coverage Over Time').classes('text-lg font-semibold mb-2')
        echart_time = ui.echart({
            'backgroundColor': 'transparent',
            'title': {'text': 'Coverage Over Time', 'left': 'center', 'top': 10, 'textStyle': {'fontSize': 16, 'color': '#000'}},
            'tooltip': {'trigger': 'axis'},
            'grid': {'left': '5%', 'right': '5%', 'bottom': '10%', 'top': '20%', 'containLabel': True},
            'xAxis': {'type': 'time'},
            'yAxis': {'type': 'value', 'name': 'Coverage (x)'},
            'series': [{'type': 'line', 'smooth': True, 'data': []}]
        }).classes('w-full h-64')

    with ui.card().classes('w-full'):
        ui.label('🎯 Target Coverage').classes('text-lg font-semibold mb-2')
        with ui.grid(columns=2).classes('w-full gap-4'):
            target_boxplot = ui.echart({
                'backgroundColor': 'transparent',
                'title': {'text': 'Target Coverage', 'left': 'center', 'top': 10, 'textStyle': {'fontSize': 16, 'color': '#000'}},
                'grid': {'left': '5%', 'right': '5%', 'bottom': '10%', 'top': '20%', 'containLabel': True},
                'dataset': [
                    {'id': 'raw', 'dimensions': ['chrom','min','Q1','median','Q3','max','chrom_index'], 'source': []},
                    {'id': 'rawdata', 'dimensions': ['chrom','coverage','name'], 'source': []},
                    {'id': 'outliers', 'dimensions': ['chrom','coverage','name'], 'source': []},
                    {'id': 'globaloutliers', 'dimensions': ['chrom','coverage','name'], 'source': []},
                ],
                'xAxis': {'type': 'category', 'name': 'Chromosome', 'nameGap': 30, 'axisLabel': {'interval': 0, 'rotate': 45}},
                'yAxis': {'type': 'value', 'name': 'Coverage (x)'},
                'legend': {'top': 50, 'selected': {'box plot': True, 'outliers': True, 'global outliers': True, 'raw data': False}},
                'series': [
                    {'name': 'box plot', 'type': 'boxplot', 'datasetId': 'raw', 'encode': {'x': 'chrom', 'y': ['min','Q1','median','Q3','max'], 'itemName': ['chrom'], 'tooltip': ['min','Q1','median','Q3','max']}},
                    {'name': 'outliers', 'type': 'scatter', 'datasetId': 'outliers', 'symbolSize': 6, 'label': {'show': True, 'position': 'right', 'formatter': '{@name}'}, 'encode': {'x': 'chrom', 'y': 'coverage', 'label': ['name'], 'tooltip': ['name','coverage']}},
                    {'name': 'global outliers', 'type': 'scatter', 'datasetId': 'globaloutliers', 'symbolSize': 6, 'label': {'show': True, 'position': 'right', 'formatter': '{@name}'}, 'encode': {'x': 'chrom', 'y': 'coverage', 'label': ['name'], 'tooltip': ['name','coverage']}},
                ],
            }).classes('w-full h-80')
            target_cov_table = ui.table(
                columns=[
                    {'name': 'chrom', 'label': 'Chrom', 'field': 'chrom'},
                    {'name': 'startpos', 'label': 'Start', 'field': 'startpos'},
                    {'name': 'endpos', 'label': 'End', 'field': 'endpos'},
                    {'name': 'name', 'label': 'Name', 'field': 'name'},
                    {'name': 'coverage', 'label': 'Coverage (x)', 'field': 'coverage'},
                ],
                rows=[],
                pagination=20
            ).classes('w-full h-80')
            target_cov_table.add_slot('body-cell-coverage', """
<q-td key="coverage" :props="props">
  <q-badge :color="props.value >= 30 ? 'green' : props.value >= 20 ? 'blue' : props.value >= 10 ? 'orange' : 'red'">
    {{ Number(props.value).toFixed(2) }}
  </q-badge>
</q-td>
""")

    # Helpers
    def _update_chr_cov(cov_df: pd.DataFrame) -> None:
        try:
            def chr_key(label: str) -> int:
                s = str(label)
                if s.startswith('chr'):
                    s = s[3:]
                mapping = {'X': 23, 'Y': 24}
                try:
                    return int(s)
                except Exception:
                    return mapping.get(s, 1000)

            pattern = r'^chr([0-9]+|X|Y)$'
            name_col = '#rname' if '#rname' in cov_df.columns else 'rname'
            temp_df = cov_df[cov_df[name_col].astype(str).str.match(pattern)]
            temp_df = temp_df[temp_df[name_col] != 'chrM']
            names = sorted(temp_df[name_col].astype(str).unique(), key=chr_key)
            temp_df = temp_df.set_index(name_col).loc[names].reset_index()
            echart_chr_cov.options['xAxis']['data'] = names
            echart_chr_cov.options['series'] = [{
                'type': 'bar', 'name': 'Chromosome', 'barWidth': '60%',
                'data': [float(v) for v in temp_df['meandepth'].tolist()]
            }]
            echart_chr_cov.update()
        except Exception:
            pass

    def _update_target_cov(cov_df: pd.DataFrame, bed_df: pd.DataFrame) -> None:
        try:
            def chr_key(label: str) -> int:
                s = str(label)
                if s.startswith('chr'):
                    s = s[3:]
                mapping = {'X': 23, 'Y': 24}
                try:
                    return int(s)
                except Exception:
                    return mapping.get(s, 1000)
            bed_df = bed_df.copy()
            bed_df['length'] = (bed_df['endpos'] - bed_df['startpos'] + 1).astype(float)
            grouped = bed_df.groupby('chrom').agg({'bases': 'sum', 'length': 'sum'}).reset_index()
            grouped['meandepth'] = grouped['bases'] / grouped['length']
            pattern = r'^chr([0-9]+|X|Y)$'
            name_col = '#rname' if '#rname' in cov_df.columns else 'rname'
            temp_df = cov_df[cov_df[name_col].astype(str).str.match(pattern)]
            temp_df = temp_df[temp_df[name_col] != 'chrM']
            names = sorted(temp_df[name_col].astype(str).unique(), key=chr_key)
            echart_target_cov.options['xAxis']['data'] = names
            grouped = grouped.set_index('chrom').reindex(names).reset_index()
            echart_target_cov.options['series'] = [
                {'type': 'scatter', 'name': 'Off Target', 'symbolSize': 8, 'data': [float(v) for v in temp_df.set_index(name_col).loc[names]['meandepth'].tolist()]},
                {'type': 'scatter', 'name': 'On Target', 'symbolSize': 8, 'data': [float(v) for v in grouped['meandepth'].fillna(0).tolist()]},
            ]
            echart_target_cov.update()
        except Exception:
            pass

    def _update_boxplot(bed_df: pd.DataFrame) -> None:
        try:
            df = bed_df.copy()
            if 'coverage' not in df.columns:
                df['length'] = (df['endpos'] - df['startpos'] + 1).astype(float)
                df['coverage'] = df['bases'] / df['length']
            chroms = natsort.natsorted(df['chrom'].astype(str).unique())
            chrom_lookup = {chrom: idx for idx, chrom in enumerate(chroms)}
            df['chrom_index'] = df['chrom'].map(chrom_lookup)
            agg = (
                df.groupby('chrom')
                .agg(
                    min=('coverage','min'),
                    Q1=('coverage', lambda x: float(np.percentile(x, 25))),
                    median=('coverage','median'),
                    Q3=('coverage', lambda x: float(np.percentile(x, 75))),
                    max=('coverage','max'),
                    chrom_index=('chrom_index','first')
                )
                .reset_index()
            )
            agg['chrom'] = pd.Categorical(agg['chrom'], categories=chroms, ordered=True)
            agg = agg.sort_values('chrom').reset_index(drop=True)
            result = [['chrom','min','Q1','median','Q3','max','chrom_index']] + agg.values.tolist()
            def iqr_bounds(sub):
                q1 = np.percentile(sub['coverage'], 25)
                q3 = np.percentile(sub['coverage'], 75)
                iqr = q3 - q1
                return q1 - 1.5*iqr, q3 + 1.5*iqr
            out_rows = []
            for c in chroms:
                sub = df[df['chrom'].astype(str) == c]
                lb, ub = iqr_bounds(sub)
                out = sub[(sub['coverage'] < lb) | (sub['coverage'] > ub)]
                out_rows += out[['chrom','coverage','name']].values.tolist()
            lb_g, ub_g = iqr_bounds(df)
            glob = df[(df['coverage'] < lb_g) | (df['coverage'] > ub_g)]
            glob_rows = glob[['chrom','coverage','name']].values.tolist()
            target_boxplot.options['dataset'][0]['source'] = result
            target_boxplot.options['dataset'][1]['source'] = df[['chrom','coverage','name']].values.tolist()
            target_boxplot.options['dataset'][2]['source'] = out_rows
            target_boxplot.options['dataset'][3]['source'] = glob_rows
            target_boxplot.options['xAxis']['data'] = chroms
            target_boxplot.update()
        except Exception:
            pass

    def _update_time(npy_path: Path) -> None:
        try:
            arr = np.load(npy_path)
            echart_time.options['series'][0]['data'] = arr.tolist()
            echart_time.update()
        except Exception:
            pass

    def _update_target_table(df: pd.DataFrame) -> None:
        try:
            dfx = df.copy()
            if 'coverage' not in dfx.columns:
                dfx['length'] = (dfx['endpos'] - dfx['startpos'] + 1).astype(float)
                dfx['coverage'] = dfx['bases'] / dfx['length']
            dfx['coverage'] = dfx['coverage'].astype(float).round(2)
            rows = dfx[['chrom','startpos','endpos','name','coverage']].to_dict(orient='records')
            target_cov_table.rows = rows
            target_cov_table.update()
        except Exception:
            pass

    def _refresh_coverage() -> None:
        try:
            if not sample_dir or not sample_dir.exists():
                return
            key = str(sample_dir)
            state = launcher._coverage_state.get(key, {})
            cov_main = sample_dir / 'coverage_main.csv'
            bed_cov = sample_dir / 'bed_coverage_main.csv'
            target_cov = sample_dir / 'target_coverage.csv'
            cov_time = sample_dir / 'coverage_time_chart.npy'
            if cov_main.exists():
                m = cov_main.stat().st_mtime
                if state.get('cov_main_mtime') != m:
                    try:
                        cov_df = pd.read_csv(cov_main)
                        _update_chr_cov(cov_df)
                        state['cov_df'] = cov_df
                    except Exception:
                        pass
                    state['cov_main_mtime'] = m
            if bed_cov.exists():
                m = bed_cov.stat().st_mtime
                if state.get('bed_cov_mtime') != m:
                    try:
                        bed_df = pd.read_csv(bed_cov)
                        state['bed_df'] = bed_df
                    except Exception:
                        bed_df = None
                    if bed_df is not None:
                        _update_boxplot(bed_df)
                    if state.get('cov_df') is not None and bed_df is not None:
                        _update_target_cov(state['cov_df'], bed_df)
                    state['bed_cov_mtime'] = m
            if target_cov.exists():
                m = target_cov.stat().st_mtime
                if state.get('target_cov_mtime') != m:
                    try:
                        tdf = pd.read_csv(target_cov)
                        _update_target_table(tdf)
                    except Exception:
                        pass
                    state['target_cov_mtime'] = m
            if cov_time.exists():
                m = cov_time.stat().st_mtime
                if state.get('cov_time_mtime') != m:
                    _update_time(cov_time)
                    state['cov_time_mtime'] = m
            # Summary
            try:
                global_cov = None
                target_cov_v = None
                enrich_v = None
                if state.get('cov_df') is not None:
                    cdf = state['cov_df']
                    if 'covbases' in cdf.columns and 'endpos' in cdf.columns and cdf['endpos'].sum() > 0:
                        global_cov = float(cdf['covbases'].sum()) / float(cdf['endpos'].sum())
                if state.get('bed_df') is not None:
                    bdf = state['bed_df'].copy()
                    if 'length' not in bdf.columns:
                        bdf['length'] = (bdf['endpos'] - bdf['startpos'] + 1).astype(float)
                    if bdf['length'].sum() > 0:
                        target_cov_v = float(bdf['bases'].sum()) / float(bdf['length'].sum())
                if global_cov is not None and target_cov_v is not None and global_cov > 0:
                    enrich_v = target_cov_v / global_cov
                if target_cov_v is not None:
                    if target_cov_v >= 30:
                        q_text, q_cls, q_bg = 'Excellent', 'text-green-600', 'bg-green-100'
                    elif target_cov_v >= 20:
                        q_text, q_cls, q_bg = 'Good', 'text-blue-600', 'bg-blue-100'
                    elif target_cov_v >= 10:
                        q_text, q_cls, q_bg = 'Moderate', 'text-yellow-600', 'bg-yellow-100'
                    else:
                        q_text, q_cls, q_bg = 'Insufficient', 'text-red-600', 'bg-red-100'
                    cov_quality_label.set_text(f'Quality: {q_text}')
                    try:
                        cov_quality_label.classes(replace=f'text-sm font-medium {q_cls}')
                        cov_quality_badge.classes(replace=f'px-2 py-1 rounded {q_bg} {q_cls}')
                    except Exception:
                        pass
                    cov_quality_badge.set_text(f'{target_cov_v:.2f}x')
                if global_cov is not None:
                    cov_global_lbl.set_text(f'Global Estimated Coverage: {global_cov:.2f}x')
                if target_cov_v is not None:
                    cov_target_lbl.set_text(f'Targets Estimated Coverage: {target_cov_v:.2f}x')
                if enrich_v is not None:
                    cov_enrich_lbl.set_text(f'Estimated enrichment: {enrich_v:.2f}x')
            except Exception:
                pass
            launcher._coverage_state[key] = state
        except Exception:
            pass

    ui.timer(30.0, _refresh_coverage, active=True)


