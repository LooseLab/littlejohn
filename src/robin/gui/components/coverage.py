from __future__ import annotations

from pathlib import Path
from typing import Any
import time

import natsort
import numpy as np
import pandas as pd
import logging

try:
    from nicegui import ui, app
except ImportError:  # pragma: no cover
    ui = None
    app = None

from robin.gui.theme import styled_table


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
        igv_div = ui.element("div").classes("w-full h-[600px] border")
        igv_status = ui.label("Checking for IGV-ready BAM files...").classes(
            "text-sm text-gray-600"
        )
        
        # IGV state management - prevent unnecessary redraws
        def _get_igv_state():
            key = str(sample_dir)
            if key not in launcher._coverage_state:
                launcher._coverage_state[key] = {}
            return launcher._coverage_state[key]
        
        def _is_igv_ready():
            """Check if IGV browser is properly initialized and ready."""
            state = _get_igv_state()
            return (
                state.get("igv_initialized", False) and 
                state.get("igv_browser_ready", False) and
                state.get("igv_loaded_bam") is not None
            )
        
        def _set_igv_ready(bam_url: str):
            """Mark IGV as ready with the current BAM."""
            state = _get_igv_state()
            state["igv_initialized"] = True
            state["igv_browser_ready"] = True
            state["igv_loaded_bam"] = bam_url
        
        def _clear_igv_state():
            """Reset IGV state to force reinitialization."""
            state = _get_igv_state()
            state["igv_initialized"] = False
            state["igv_browser_ready"] = False
            state["igv_loaded_bam"] = None
        
        # Immediate check for existing IGV BAM files
        def _check_existing_igv_bam():
            try:
                # Only proceed if IGV isn't already ready
                if _is_igv_ready():
                    igv_status.set_text("IGV is already loaded and ready.")
                    return
                
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
                
                if bam_path is not None:
                    igv_status.set_text(f"Found existing BAM: {bam_path.name}")
                    # Trigger immediate IGV loading
                    ui.timer(0.1, lambda: _load_igv_bam(bam_path), once=True)
                else:
                    igv_status.set_text("No IGV-ready BAM found yet (need coordinate-sorted BAM with .bai).")
                    
            except Exception as e:
                igv_status.set_text(f"Error checking for BAM files: {e}")
        
        # Function to load BAM into IGV
        def _load_igv_bam(bam_path: Path):
            try:
                # Check if we already have this BAM loaded
                state = _get_igv_state()
                bam_url = f"/samples/{sample_dir.name}/{bam_path.name}"
                if state.get("igv_loaded_bam") == bam_url and _is_igv_ready():
                    igv_status.set_text(f"BAM {bam_path.name} already loaded in IGV.")
                    return
                
                # Mount directory once
                if bam_path.parent == sample_dir / "igv":
                    mount = f"/samples/{sample_dir.name}/igv"
                elif bam_path.parent == sample_dir / "clair3":
                    mount = f"/samples/{sample_dir.name}/clair3"
                else:
                    mount = f"/samples/{sample_dir.name}"
                
                try:
                    if app is not None:
                        app.add_static_files(mount, str(bam_path.parent))
                except Exception:
                    # Already mounted or in tests without app
                    pass

                # Define URLs
                bam_url = f"{mount}/{bam_path.name}"
                bai_url = f"{mount}/{bam_path.name}.bai"
                
                # Check if we need to create a new browser or just add tracks
                if not state.get("igv_initialized") or not state.get("igv_browser_ready"):
                    # Create new IGV browser
                    js_create = f"""
                        const el = getElement('{igv_div.id}');
                        const track = {{ name: '{sample_dir.name}', url: '{bam_url}', indexURL: '{bai_url}', format: 'bam', visibilityWindow: 1000000 }};
                        const options = {{ genome: 'hg38', locus: 'chr1:1-1', tracks: [track] }};
                        igv.createBrowser(el, options).then(b => {{ 
                            window.lj_igv = b; 
                            window.lj_igv_browser_ready = true;
                            console.log('IGV browser created successfully');
                        }});
                    """
                    
                    try:
                        ui.run_javascript(js_create, timeout=30.0)
                        _set_igv_ready(bam_url)
                        igv_status.set_text(f"IGV browser created with {bam_path.name}")
                    except Exception as e:
                        igv_status.set_text(f"Failed to create IGV browser: {e}")
                        _clear_igv_state()
                else:
                    # Browser exists, just add/update the track
                    js_add_track = f"""
                        if (window.lj_igv && window.lj_igv_browser_ready) {{
                            // Remove existing track with same name if it exists
                            const existingTrack = window.lj_igv.findTracksByName('{sample_dir.name}')[0];
                            if (existingTrack) {{
                                window.lj_igv.removeTrack(existingTrack);
                            }}
                            
                            // Add new track
                            const track = {{ name: '{sample_dir.name}', url: '{bam_url}', indexURL: '{bai_url}', format: 'bam', visibilityWindow: 1000000 }};
                            window.lj_igv.loadTrack(track);
                            console.log('Track updated in existing IGV browser');
                        }}
                    """
                    
                    try:
                        ui.run_javascript(js_add_track, timeout=30.0)
                        _set_igv_ready(bam_url)
                        igv_status.set_text(f"Track updated in IGV: {bam_path.name}")
                    except Exception as e:
                        igv_status.set_text(f"Failed to update track: {e}")
                        
            except Exception as e:
                igv_status.set_text(f"Error loading IGV: {e}")
                _clear_igv_state()
        
        # Check for existing files immediately
        ui.timer(0.1, _check_existing_igv_bam, once=True)
        
        # Function to refresh IGV BAM check (useful after regeneration)
        def _refresh_igv_check():
            try:
                # Only refresh if IGV isn't ready or if we need to check for new files
                if _is_igv_ready():
                    # Check if the current BAM still exists
                    current_bam = _get_igv_state().get("igv_loaded_bam")
                    if current_bam:
                        # Extract filename from URL and check if it still exists
                        bam_name = current_bam.split("/")[-1]
                        igv_dir = sample_dir / "igv"
                        clair_dir = sample_dir / "clair3"
                        candidates = [
                            igv_dir / bam_name,
                            clair_dir / bam_name,
                            sample_dir / bam_name,
                        ]
                        
                        if any(p.exists() and Path(f"{p}.bai").exists() for p in candidates):
                            igv_status.set_text("IGV is ready and BAM file is current.")
                            return
                
                # If we get here, we need to refresh
                _clear_igv_state()
                _check_existing_igv_bam()
            except Exception as e:
                print(f"Error refreshing IGV check: {e}")
        
        # Function to clear data tracks from IGV (keeps reference tracks)
        def _clear_igv_tracks():
            try:
                # First check if IGV is actually ready
                if not _is_igv_ready():
                    igv_status.set_text("IGV is not ready - cannot clear tracks.")
                    return
                
                # More robust JavaScript to clear tracks
                js_clear = """
                    try {
                        console.log('Attempting to clear IGV tracks...');
                        console.log('window.lj_igv:', window.lj_igv);
                        console.log('window.lj_igv_browser_ready:', window.lj_igv_browser_ready);
                        
                        if (window.lj_igv) {
                            // Try different methods to get tracks
                            let tracks = [];
                            
                            // Method 1: Try getTracks() (standard IGV.js API)
                            if (typeof window.lj_igv.getTracks === 'function') {
                                tracks = window.lj_igv.getTracks();
                                console.log('Using getTracks() method, found:', tracks.length);
                            }
                            // Method 2: Try findTracksByName (alternative API)
                            else if (typeof window.lj_igv.findTracksByName === 'function') {
                                tracks = window.lj_igv.findTracksByName('*');
                                console.log('Using findTracksByName method, found:', tracks.length);
                            }
                            // Method 3: Try accessing tracks property directly
                            else if (window.lj_igv.tracks && Array.isArray(window.lj_igv.tracks)) {
                                tracks = window.lj_igv.tracks;
                                console.log('Using tracks property, found:', tracks.length);
                            }
                            // Method 4: Try to enumerate all properties to find tracks
                            else {
                                console.log('Available methods on IGV browser:');
                                for (let prop in window.lj_igv) {
                                    if (typeof window.lj_igv[prop] === 'function') {
                                        console.log('Method:', prop);
                                    }
                                }
                                console.log('Available properties on IGV browser:');
                                for (let prop in window.lj_igv) {
                                    if (typeof window.lj_igv[prop] !== 'function') {
                                        console.log('Property:', prop, '=', window.lj_igv[prop]);
                                    }
                                }
                                throw new Error('Could not find tracks using any known method');
                            }
                            
                            if (tracks.length > 0) {
                                console.log('Found tracks to remove:', tracks.length);
                                
                                // Only remove data tracks, keep reference tracks
                                const tracksToRemove = [];
                                const tracksToKeep = [];
                                
                                tracks.forEach((track, index) => {
                                    const trackName = track.name || track.id || 'unnamed';
                                    const trackType = track.type || 'unknown';
                                    
                                    // Keep essential IGV reference tracks
                                    if (trackName === 'ideogram' || 
                                        trackName === 'ruler' || 
                                        trackName === 'sequence' ||
                                        trackName === 'Refseq Select' ||
                                        trackType === 'ideogram' ||
                                        trackType === 'ruler' ||
                                        trackType === 'sequence') {
                                        tracksToKeep.push(track);
                                        console.log('Keeping track', index, ':', trackName, '(', trackType, ')');
                                    } else {
                                        // This is a data track (BAM, etc.) - remove it
                                        tracksToRemove.push(track);
                                        console.log('Will remove track', index, ':', trackName, '(', trackType, ')');
                                    }
                                });
                                
                                console.log('Tracks to keep:', tracksToKeep.length);
                                console.log('Tracks to remove:', tracksToRemove.length);
                                
                                // Remove only the data tracks
                                tracksToRemove.forEach((track, index) => {
                                    try {
                                        if (typeof window.lj_igv.removeTrack === 'function') {
                                            window.lj_igv.removeTrack(track);
                                            console.log('Removed data track:', track.name || track.id || 'unnamed');
                                        } else {
                                            console.log('removeTrack method not available');
                                        }
                                    } catch (e) {
                                        console.error('Error removing track:', e);
                                    }
                                });
                                
                                console.log('Data track removal completed');
                                
                                // Force a redraw
                                if (typeof window.lj_igv.redraw === 'function') {
                                    window.lj_igv.redraw();
                                    console.log('Forced IGV redraw');
                                } else if (typeof window.lj_igv.update === 'function') {
                                    window.lj_igv.update();
                                    console.log('Forced IGV update');
                                } else {
                                    console.log('No redraw/update method available');
                                }
                            } else {
                                console.log('No tracks found to remove');
                            }
                        } else {
                            console.log('No IGV browser found');
                        }
                    } catch (e) {
                        console.error('Error in track clearing:', e);
                    }
                """
                
                ui.run_javascript(js_clear, timeout=15.0)
                
                # Update status and state
                igv_status.set_text("IGV data tracks cleared. Reference tracks preserved.")
                
                # Don't clear the IGV state - just mark that data tracks are cleared
                state = _get_igv_state()
                state["data_tracks_cleared"] = True
                
            except Exception as e:
                igv_status.set_text(f"Error clearing tracks: {e}")
                print(f"Error in _clear_igv_tracks: {e}")
        
        # Function to reload BAM data into existing IGV browser
        def _reload_bam_track():
            try:
                # First check if IGV is actually ready
                if not _is_igv_ready():
                    igv_status.set_text("IGV is not ready - cannot reload BAM.")
                    return
                
                # Find the current BAM file
                state = _get_igv_state()
                current_bam_url = state.get("igv_loaded_bam")
                
                if not current_bam_url:
                    igv_status.set_text("No BAM file currently loaded - cannot reload.")
                    return
                
                # Extract BAM name for display
                bam_name = current_bam_url.split("/")[-1]
                
                # JavaScript to reload the BAM track
                js_reload = f"""
                    try {{
                        console.log('Attempting to reload BAM track...');
                        console.log('BAM URL:', '{current_bam_url}');
                        
                        if (window.lj_igv) {{
                            // Create the track configuration
                            const trackConfig = {{
                                name: '{bam_name}',
                                url: '{current_bam_url}',
                                indexURL: '{current_bam_url}.bai',
                                format: 'bam',
                                visibilityWindow: 1000000,
                                type: 'alignment'
                            }};
                            
                            console.log('Track config:', trackConfig);
                            
                            // Add the track to the browser
                            if (typeof window.lj_igv.loadTrack === 'function') {{
                                window.lj_igv.loadTrack(trackConfig);
                                console.log('BAM track reloaded successfully');
                            }} else {{
                                console.log('loadTrack method not available');
                                // Try alternative method
                                if (typeof window.lj_igv.addTrack === 'function') {{
                                    window.lj_igv.addTrack(trackConfig);
                                    console.log('BAM track added using addTrack method');
                                }} else {{
                                    throw new Error('No track loading method available');
                                }}
                            }}
                            
                            // Force a redraw
                            if (typeof window.lj_igv.redraw === 'function') {{
                                window.lj_igv.redraw();
                                console.log('Forced IGV redraw after reload');
                            }} else if (typeof window.lj_igv.update === 'function') {{
                                window.lj_igv.update();
                                console.log('Forced IGV update after reload');
                            }}
                        }} else {{
                            console.log('No IGV browser found');
                        }}
                    }} catch (e) {{
                        console.error('Error reloading BAM track:', e);
                    }}
                """
                
                ui.run_javascript(js_reload, timeout=20.0)
                
                # Update status
                igv_status.set_text(f"Reloading BAM track: {bam_name}")
                
                # Mark that data tracks are no longer cleared
                state["data_tracks_cleared"] = False
                
            except Exception as e:
                igv_status.set_text(f"Error reloading BAM: {e}")
                print(f"Error in _reload_bam_track: {e}")
        
        # Function to debug IGV state
        def _debug_igv_state():
            try:
                js_debug = """
                    console.log('=== IGV Debug Info ===');
                    console.log('window.lj_igv:', window.lj_igv);
                    console.log('window.lj_igv_browser_ready:', window.lj_igv_browser_ready);
                    
                    if (window.lj_igv) {
                        console.log('IGV browser exists');
                        console.log('Browser type:', typeof window.lj_igv);
                        console.log('Available methods:', Object.getOwnPropertyNames(window.lj_igv));
                        
                        // Try different methods to get tracks
                        let tracks = [];
                        let methodUsed = 'none';
                        
                        if (typeof window.lj_igv.getTracks === 'function') {
                            tracks = window.lj_igv.getTracks();
                            methodUsed = 'getTracks()';
                        } else if (typeof window.lj_igv.findTracksByName === 'function') {
                            tracks = window.lj_igv.findTracksByName('*');
                            methodUsed = 'findTracksByName()';
                        } else if (window.lj_igv.tracks && Array.isArray(window.lj_igv.tracks)) {
                            tracks = window.lj_igv.tracks;
                            methodUsed = 'tracks property';
                        }
                        
                        console.log('Method used to get tracks:', methodUsed);
                        console.log('Current tracks:', tracks.length);
                        
                        if (tracks.length > 0) {
                            tracks.forEach((track, i) => {
                                console.log(`Track ${i}:`, {
                                    name: track.name || 'unnamed',
                                    id: track.id || 'no-id',
                                    type: track.type || 'unknown-type',
                                    visible: track.visible !== undefined ? track.visible : 'unknown'
                                });
                            });
                        }
                        
                        // Check for common IGV methods
                        const commonMethods = ['getTracks', 'findTracksByName', 'removeTrack', 'addTrack', 'redraw', 'update'];
                        console.log('Available common methods:');
                        commonMethods.forEach(method => {
                            console.log(`  ${method}:`, typeof window.lj_igv[method] === 'function' ? '✓' : '✗');
                        });
                        
                    } else {
                        console.log('No IGV browser found');
                    }
                    
                    // Check Python state
                    console.log('Python state check requested');
                """
                
                ui.run_javascript(js_debug, timeout=10.0)
                
                # Also show Python state
                state = _get_igv_state()
                debug_info = f"Python IGV State: initialized={state.get('igv_initialized')}, browser_ready={state.get('igv_browser_ready')}, bam={state.get('igv_loaded_bam')}, tracks_cleared={state.get('tracks_cleared')}"
                igv_status.set_text(debug_info)
                print(debug_info)
                
            except Exception as e:
                igv_status.set_text(f"Debug error: {e}")
                print(f"Error in debug: {e}")
        
        # IGV control buttons
        with ui.row().classes("items-center gap-2 mt-2"):
            ui.button("🔄 Refresh IGV Check", on_click=_refresh_igv_check)
            ui.button("Clear Data Tracks", on_click=_clear_igv_tracks)
            ui.button("📁 Reload BAM", on_click=_reload_bam_track)
            ui.button("🐛 Debug IGV", on_click=_debug_igv_state)
        
        # BAM generation buttons
        async def _trigger_build_sorted_bam(force_regenerate: bool = False) -> None:
            try:
                # Debug: Check what's in the launcher
                debug_info = f"launcher type: {type(launcher)}, workflow_runner: {getattr(launcher, 'workflow_runner', 'None')}"
                
                if force_regenerate:
                    ui.notify(f"Force regenerate mode: will recreate IGV BAM even if it exists", type="info")
                
                # Check if target.bam exists (required for IGV BAM generation)
                target_bam = sample_dir / "target.bam"
                if not target_bam.exists():
                    ui.notify("No target.bam file found. Please run target analysis first.", type="warning")
                    return
                
                # Check if IGV BAM already exists (only for normal generation, not force regenerate)
                igv_bam = sample_dir / "igv" / "igv_ready.bam"
                if not force_regenerate and igv_bam.exists() and (sample_dir / "igv" / "igv_ready.bam.bai").exists():
                    ui.notify("IGV BAM already exists and is ready.", type="positive")
                    return
                
                if force_regenerate:
                    ui.notify(f"Proceeding with force regenerate - existing files will be overwritten", type="info")
                
                # Check if we have a workflow runner
                if not hasattr(launcher, 'workflow_runner') or not launcher.workflow_runner:
                    ui.notify(f"No workflow runner available. Debug: {debug_info}. GUI is running but workflow integration is not available. This may be a configuration issue.", type="warning")
                    return
                
                # Submit the IGV BAM generation job to the workflow system
                try:
                    runner = launcher.workflow_runner
                    sample_id = sample_dir.name
                    
                    # Check if it's a Ray workflow or Simple workflow
                    if hasattr(runner, 'submit_sample_job'):
                        # Simple workflow
                        success = runner.submit_sample_job(str(sample_dir), 'igv_bam', sample_id, force_regenerate)
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
                        action = "regenerated" if force_regenerate else "generated"
                        ui.notify(f"IGV BAM {action} job submitted to workflow queue!", type="positive")
                        # Refresh the IGV check after successful job submission
                        _refresh_igv_check()
                    else:
                        ui.notify("Failed to submit IGV BAM generation job. Check logs for details.", type="warning")
                        
                except Exception as e:
                    ui.notify(f"Error submitting IGV BAM job: {e}", type="negative")
                    
            except Exception as e:
                try:
                    ui.notify(f"Error checking IGV BAM status: {e}", type="negative")
                except Exception:
                    # Client may have disconnected, ignore UI errors
                    pass
        
        with ui.row().classes("items-center gap-2 mt-2"):
            ui.button("Generate IGV BAM (if missing)", on_click=_trigger_build_sorted_bam)
            ui.button("Force Regenerate IGV BAM", on_click=lambda: _trigger_build_sorted_bam(force_regenerate=True))

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
            
            # IGV management - only update if necessary
            try:
                # Check if IGV needs attention (only if not already ready)
                if not _is_igv_ready():
                    # Look for available BAM files
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
                    
                    if bam_path is not None:
                        # Found a BAM file, trigger IGV loading
                        ui.timer(0.1, lambda: _load_igv_bam(bam_path), once=True)
                else:
                    # IGV is ready, just update status occasionally
                    if state.get("last_igv_status_update", 0) < time.time() - 60:  # Update status every minute
                        try:
                            igv_status.set_text("IGV is ready and loaded.")
                            state["last_igv_status_update"] = time.time()
                        except Exception:
                            pass
            except Exception:
                # Never let IGV logic break coverage refresh
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
