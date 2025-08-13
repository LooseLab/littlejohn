from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List
import csv

import pandas as pd

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None


def add_classification_section(sample_dir: Path) -> None:
    """Build the Classification section (Sturgeon, NanoDX, PanNanoDX, RF)."""
    ui.label('🧪 Classification').classes('text-lg font-semibold mb-2')
    tool_to_file = {
        'Sturgeon': {'file': 'sturgeon_scores.csv', 'mode': 'fraction'},
        'NanoDX': {'file': 'NanoDX_scores.csv', 'mode': 'fraction'},
        'PanNanoDX': {'file': 'PanNanoDX_scores.csv', 'mode': 'fraction'},
        'Random Forest': {'file': 'random_forest_scores.csv', 'mode': 'percent'},
    }
    charts: Dict[str, Dict[str, Any]] = {}
    for tool_name, cfg in tool_to_file.items():
        exp = ui.expansion(tool_name, icon='analytics').classes('w-full')
        with exp:
            summary_labels = None
            if tool_name == 'Sturgeon':
                with ui.card().classes('w-full bg-gradient-to-r from-blue-50 to-indigo-50 mb-2 p-2'):
                    with ui.row().classes('gap-6 items-center'):
                        st_class = ui.label('Sturgeon classification: Unknown').classes('text-sm font-semibold text-blue-800')
                        st_conf = ui.label('Confidence: --%').classes('text-sm text-gray-700')
                        st_probes = ui.label('Probes: --').classes('text-sm text-gray-700')
                    summary_labels = {'class': st_class, 'conf': st_conf, 'probes': st_probes}
            elif tool_name in ('NanoDX', 'PanNanoDX'):
                with ui.card().classes('w-full bg-gradient-to-r from-blue-50 to-indigo-50 mb-2 p-2'):
                    with ui.row().classes('gap-6 items-center'):
                        ndx_class = ui.label(f'{tool_name} classification: Unknown').classes('text-sm font-semibold text-blue-800')
                        ndx_conf = ui.label('Confidence: --%').classes('text-sm text-gray-700')
                        ndx_feats = ui.label('Probes: --').classes('text-sm text-gray-700')
                    summary_labels = {'class': ndx_class, 'conf': ndx_conf, 'probes': ndx_feats}
            elif tool_name == 'Random Forest':
                with ui.card().classes('w-full bg-gradient-to-r from-blue-50 to-indigo-50 mb-2 p-2'):
                    with ui.row().classes('gap-6 items-center'):
                        rf_class = ui.label('Forest classification: Unknown').classes('text-sm font-semibold text-blue-800')
                        rf_conf = ui.label('Confidence: --%').classes('text-sm text-gray-700')
                        rf_feats = ui.label('Features: --').classes('text-sm text-gray-700')
                    summary_labels = {'class': rf_class, 'conf': rf_conf, 'probes': rf_feats}
            ui.label(f'{tool_name} current classification').classes('text-sm text-gray-700')
            bar = ui.echart({
                'backgroundColor': 'transparent',
                'title': {'text': f'{tool_name} (Top classes)', 'left': 'center', 'top': 10, 'textStyle': {'fontSize': 16, 'color': '#000'}},
                'tooltip': {'trigger': 'axis', 'axisPointer': {'type': 'shadow'}, 'formatter': '{b}: {c}%'},
                'grid': {'left': '5%', 'right': '5%', 'bottom': '5%', 'top': '25%', 'containLabel': True},
                'xAxis': {'type': 'value', 'min': 0, 'max': 100, 'interval': 20, 'axisLabel': {'formatter': '{value}%'}},
                'yAxis': {'type': 'category', 'inverse': True, 'data': []},
                'series': [{'type': 'bar', 'data': [], 'barMaxWidth': '60%', 'itemStyle': {'color': '#007AFF', 'borderRadius': [0, 4, 4, 0]}, 'label': {'show': True, 'position': 'right', 'formatter': '{c}%'}}],
            }).classes('w-full h-60')
            ui.label(f'{tool_name} confidence over time').classes('text-sm text-gray-700 mt-2')
            ts = ui.echart({
                'backgroundColor': 'transparent',
                'title': {'text': f'{tool_name} (time series)', 'left': 'center', 'top': 5, 'textStyle': {'fontSize': 16, 'color': '#000'}},
                'tooltip': {'trigger': 'axis'},
                'legend': {'type': 'scroll', 'top': 45},
                'grid': {'left': '5%', 'right': '5%', 'bottom': '5%', 'top': '30%', 'containLabel': True},
                'xAxis': {'type': 'time'},
                'yAxis': {'type': 'value', 'min': 0, 'max': 100, 'axisLabel': {'formatter': '{value}%'}, 'splitLine': {'show': True, 'lineStyle': {'type': 'dashed', 'color': '#E0E0E0'}}},
                'series': [],
            }).classes('w-full h-64')
            charts[tool_name] = {'bar': bar, 'ts': ts, 'file': cfg['file'], 'mode': cfg['mode'], 'summary': summary_labels, 'expansion': exp, 'last_mtime': None}

    def _read_scores_csv(csv_path: Path, mode: str):
        try:
            with csv_path.open('r', newline='') as fh:
                reader = csv.DictReader(fh)
                rows = list(reader)
                if not rows:
                    return None
                numeric_keys = []
                for k in reader.fieldnames or []:
                    if k.lower() in {'timestamp', 'number_probes'}:
                        continue
                    try:
                        float(rows[-1].get(k, ''))
                        numeric_keys.append(k)
                    except Exception:
                        continue
                if not numeric_keys:
                    return None
                last_scores = {}
                for k in numeric_keys:
                    try:
                        v = float(rows[-1].get(k, 0))
                        if mode == 'fraction':
                            v *= 100.0
                        last_scores[k] = round(v, 2)
                    except Exception:
                        pass
                time_key = 'timestamp' if 'timestamp' in (reader.fieldnames or []) else None
                x_labels = []
                series_map: Dict[str, List[List[Any]]] = {k: [] for k in numeric_keys}
                for idx, r in enumerate(rows):
                    if time_key:
                        raw_x = r.get(time_key)
                        try:
                            x_val = float(raw_x)
                        except Exception:
                            x_val = raw_x
                    else:
                        x_val = idx + 1
                    x_labels.append(x_val)
                    for k in numeric_keys:
                        try:
                            vv = float(r.get(k, 0))
                            if mode == 'fraction':
                                vv *= 100.0
                            series_map[k].append([x_val, round(vv, 2)])
                        except Exception:
                            series_map[k].append([x_val, 0])
                if time_key:
                    def _sort_key(p):
                        try:
                            return float(p[0])
                        except Exception:
                            return str(p[0])
                    for k in series_map:
                        series_map[k] = sorted(series_map[k], key=_sort_key)
                return {'last': last_scores, 'x': x_labels, 'series': series_map, 'has_time': bool(time_key)}
        except Exception:
            return None

    def _update_charts_from_file(tool_name: str, file_name: str):
        try:
            file_path = sample_dir / file_name if sample_dir else None
            if not file_path or not file_path.exists():
                return
            mode = charts[tool_name].get('mode', 'percent')
            data = _read_scores_csv(file_path, mode)
            if not data:
                return
            bar = charts[tool_name]['bar']
            ts = charts[tool_name]['ts']
            last = data['last']
            top = sorted(last.items(), key=lambda kv: kv[1], reverse=True)[:10]
            labels = [k for k, _ in top][::-1]
            values = [round(v, 2) for _, v in top][::-1]
            bar.options['yAxis']['data'] = labels
            bar.options['series'][0]['data'] = values
            bar.update()
            ts.options['series'] = []
            if data.get('has_time'):
                ts.options['xAxis'] = {'type': 'time'}
                for k, _ in top[:5]:
                    ts.options['series'].append({'name': k, 'type': 'line', 'smooth': True, 'animation': False, 'data': data['series'][k]})
            else:
                ts.options['xAxis'] = {'type': 'category', 'data': data['x']}
                for k, _ in top[:5]:
                    ts.options['series'].append({'name': k, 'type': 'line', 'smooth': True, 'animation': False, 'data': [y for _, y in data['series'][k]]})
            # Inline thresholds where required
            if tool_name == 'Sturgeon':
                ts.options.setdefault('series', [])
                ts.options['series'].append({'name': 'thresholds', 'type': 'line', 'data': [], 'markLine': {'silent': True, 'lineStyle': {'type': 'dashed', 'color': '#999'}, 'data': [{'yAxis': 60}, {'yAxis': 80}]}})
            if tool_name == 'Random Forest':
                ts.options.setdefault('series', [])
                ts.options['series'].append({'name': 'thresholds', 'type': 'line', 'data': [], 'markLine': {'silent': True, 'lineStyle': {'type': 'dashed', 'color': '#999'}, 'data': [{'yAxis': 65}, {'yAxis': 85}]}})
            ts.update()
            # Summaries
            if charts[tool_name].get('summary'):
                labels_map = charts[tool_name]['summary']
                if labels_map and top:
                    best_label, best_value = top[0]
                    if tool_name == 'Random Forest':
                        labels_map['class'].set_text(f'Forest classification: {best_label}')
                    else:
                        labels_map['class'].set_text(f'{tool_name} classification: {best_label}')
                    labels_map['conf'].set_text(f'Confidence: {best_value:.2f}%')
                    try:
                        with file_path.open('r', newline='') as fh:
                            reader = csv.DictReader(fh)
                            rows = list(reader)
                        if rows:
                            npv = rows[-1].get('number_probes') or rows[-1].get('number_probes'.lower())
                            if npv is not None:
                                label = 'Features' if tool_name == 'Random Forest' else 'Probes'
                                labels_map['probes'].set_text(f'{label}: {int(float(npv))}')
                    except Exception:
                        pass
                try:
                    exp = charts[tool_name].get('expansion')
                    if exp and top:
                        best_label, best_value = top[0]
                        exp.props(f'label="{tool_name} — {best_label} ({round(best_value, 2)}%)"')
                except Exception:
                    pass
        except Exception:
            pass

    def _refresh_classification() -> None:
        for tool_name, cfg in tool_to_file.items():
            file_path = sample_dir / cfg['file'] if sample_dir else None
            try:
                if file_path and file_path.exists():
                    mtime = file_path.stat().st_mtime
                    if charts[tool_name].get('last_mtime') is None or mtime > charts[tool_name].get('last_mtime', 0):
                        _update_charts_from_file(tool_name, cfg['file'])
                        charts[tool_name]['last_mtime'] = mtime
            except Exception:
                pass

    ui.timer(2.0, _refresh_classification, active=True)


