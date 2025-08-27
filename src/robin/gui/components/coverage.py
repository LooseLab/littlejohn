from __future__ import annotations

from pathlib import Path
from typing import Any
import time
import os

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
        
        # Add IGV library status indicator
        igv_lib_status = ui.label("IGV library: Checking...").classes(
            "text-xs text-gray-500"
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
        
        # Function to check if IGV library is available
        def _check_igv_library():
            """Check if the IGV JavaScript library is available."""
            js_check = """
                try {
                    if (typeof igv === 'undefined') {
                        console.error('IGV library not loaded');
                        return false;
                    }
                    
                    // Check if IGV is fully initialized
                    if (typeof igv.createBrowser !== 'function') {
                        console.error('IGV library not fully initialized');
                        return false;
                    }
                    
                    console.log('IGV library is available and ready');
                    return true;
                } catch (e) {
                    console.error('Error checking IGV library:', e);
                    return false;
                }
            """
            try:
                result = ui.run_javascript(js_check, timeout=10.0)
                # Update the status indicator
                if result:
                    igv_lib_status.set_text("IGV library: ✓ Loaded and ready")
                    igv_lib_status.classes("text-xs text-green-600")
                else:
                    igv_lib_status.set_text("IGV library: ✗ Not ready")
                    igv_lib_status.classes("text-xs text-red-600")
                return result
            except Exception as e:
                igv_lib_status.set_text("IGV library: ✗ Error checking")
                igv_lib_status.classes("text-xs text-red-600")
                print(f"Error checking IGV library: {e}")
                return False
        
        # Function to retry IGV creation
        def _retry_igv_creation(bam_path: Path):
            """Retry IGV browser creation after a failure."""
            try:
                igv_status.set_text("Retrying IGV browser creation...")
                _load_igv_bam(bam_path)
            except Exception as e:
                igv_status.set_text(f"Retry failed: {e}")
                print(f"IGV retry error: {e}")
        
        # Function to wait for element to be ready
        def _wait_for_element_ready():
            """Wait for the IGV div element to be properly rendered."""
            js_wait = f"""
                return new Promise((resolve, reject) => {{
                    const el = document.getElementById('{igv_div.id}');
                    if (!el) {{
                        reject(new Error('Element not found'));
                        return;
                    }}
                    
                    if (el.offsetWidth > 0 && el.offsetHeight > 0) {{
                        resolve(true);
                        return;
                    }}
                    
                    // Wait for element to be ready
                    const checkReady = () => {{
                        if (el.offsetWidth > 0 && el.offsetHeight > 0) {{
                            resolve(true);
                        }} else {{
                            setTimeout(checkReady, 100);
                        }}
                    }};
                    checkReady();
                }});
            """
            try:
                return ui.run_javascript(js_wait, timeout=10.0)
            except Exception:
                return False
        
        # Function to load BAM into IGV
        def _load_igv_bam(bam_path: Path):
            try:
                # First check if IGV library is available
                if not _check_igv_library():
                    igv_status.set_text("IGV library not yet loaded. Please wait and try again.")
                    return
                
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
                    # Wait for element to be ready before creating IGV browser
                    if not _wait_for_element_ready():
                        igv_status.set_text("Waiting for IGV element to be ready...")
                        ui.timer(0.5, lambda: _load_igv_bam(bam_path), once=True)
                        return
                    
                    # Create new IGV browser
                    js_create = f"""
                        try {{
                            console.log('Starting IGV browser creation...');
                            console.log('IGV library available:', typeof igv !== 'undefined');
                            console.log('IGV library version:', igv.version || 'unknown');
                            
                            // Check for minimum version compatibility
                            if (igv.version && igv.version < '2.0.0') {{
                                console.warn('IGV version may be too old for this functionality');
                            }}
                            console.log('Target element:', document.getElementById('{igv_div.id}'));
                            
                            const el = document.getElementById('{igv_div.id}');
                            if (!el) {{
                                throw new Error('Target element not found');
                            }}
                            
                            // Check if element is properly rendered
                            if (el.offsetWidth === 0 || el.offsetHeight === 0) {{
                                throw new Error('Target element has no dimensions - not ready for IGV');
                            }}
                            
                            // Check if IGV library is fully initialized
                            if (typeof igv.createBrowser !== 'function') {{
                                throw new Error('IGV library not fully initialized');
                            }}
                            
                            const track = {{ name: '{sample_dir.name}', url: '{bam_url}', indexURL: '{bai_url}', format: 'bam', visibilityWindow: 1000000 }};
                            const options = {{ genome: 'hg38', locus: 'chr1:1-1', tracks: [track] }};
                            
                            console.log('Creating IGV browser with options:', options);
                            
                                                    igv.createBrowser(el, options).then(b => {{
                            window.lj_igv = b; 
                            window.lj_igv_browser_ready = true;
                            console.log('IGV browser created successfully');
                        }}).catch(error => {{
                            console.error('IGV browser creation failed:', error);
                            window.lj_igv_browser_ready = false;
                            
                            // Try fallback method
                            try {{
                                console.log('Trying fallback IGV creation...');
                                const fallbackOptions = {{ genome: 'hg38' }};
                                igv.createBrowser(el, fallbackOptions).then(b => {{
                                    window.lj_igv = b;
                                    window.lj_igv_browser_ready = true;
                                    console.log('IGV browser created with fallback method');
                                }});
                            }} catch (fallbackError) {{
                                console.error('Fallback IGV creation also failed:', fallbackError);
                                throw error;
                            }}
                        }});
                        }} catch (error) {{
                            console.error('Error in IGV browser creation:', error);
                            throw error;
                        }}
                    """
                    
                    try:
                        ui.run_javascript(js_create, timeout=30.0)
                        _set_igv_ready(bam_url)
                        igv_status.set_text(f"IGV browser created with {bam_path.name}")
                    except Exception as e:
                        igv_status.set_text(f"Failed to create IGV browser: {e}")
                        print(f"IGV browser creation error: {e}")
                        
                        # Try to retry after a delay
                        ui.timer(2.0, lambda: _retry_igv_creation(bam_path), once=True)
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
        
        # Add a delay to ensure IGV library is fully loaded
        ui.timer(1.0, _check_existing_igv_bam, once=True)
        
        # Check IGV library status periodically
        ui.timer(2.0, _check_igv_library, once=True)
        ui.timer(5.0, _check_igv_library, once=True)
        
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
        
        # Function to check console for errors
        def _check_console_errors():
            """Check the browser console for any JavaScript errors."""
            js_check_console = """
                try {
                    // Check if there are any console errors
                    const errors = [];
                    
                    // Override console.error to capture errors
                    const originalError = console.error;
                    console.error = function(...args) {
                        errors.push(args.join(' '));
                        originalError.apply(console, args);
                    };
                    
                    // Check for common IGV-related issues
                    if (typeof igv === 'undefined') {
                        errors.push('IGV library not loaded');
                    } else if (typeof igv.createBrowser !== 'function') {
                        errors.push('IGV createBrowser method not available');
                    } else if (igv.version) {
                        errors.push(`IGV version: ${igv.version}`);
                    }
                    
                    // Check for element issues
                    const el = document.getElementById('igv_div');
                    if (!el) {
                        errors.push('IGV div element not found');
                    } else if (el.offsetWidth === 0 || el.offsetHeight === 0) {
                        errors.push('IGV div element has no dimensions');
                    }
                    
                    // Check for IGV library issues
                    if (typeof igv === 'undefined') {
                        errors.push('IGV library not loaded');
                    } else if (typeof igv.createBrowser !== 'function') {
                        errors.push('IGV createBrowser method not available');
                    } else if (igv.version) {
                        errors.push(`IGV version: ${igv.version}`);
                    }
                    
                    // Restore original console.error
                    console.error = originalError;
                    
                    return errors;
                } catch (e) {
                    return ['Error checking console: ' + e.message];
                }
            """
            
            try:
                errors = ui.run_javascript(js_check_console, timeout=10.0)
                if errors and len(errors) > 0:
                    error_text = "Console errors found:\n" + "\n".join(errors)
                    ui.notify(error_text, type="warning")
                    print(error_text)
                else:
                    ui.notify("No console errors found", type="positive")
            except Exception as e:
                ui.notify(f"Error checking console: {e}", type="negative")
        
        # Function to test IGV browser creation
        def _test_igv_browser():
            """Test IGV browser creation with minimal configuration."""
            try:
                if not _check_igv_library():
                    ui.notify("IGV library not available", type="warning")
                    return
                
                js_test = f"""
                    try {{
                        console.log('Testing IGV browser creation...');
                        const el = document.getElementById('{igv_div.id}');
                        if (!el) {{
                            throw new Error('Test element not found');
                        }}
                        
                        const options = {{ 
                            genome: 'hg38', 
                            locus: 'chr1:1-1' 
                        }};
                        
                        console.log('Creating test IGV browser...');
                        igv.createBrowser(el, options).then(b => {{
                            window.lj_igv = b;
                            window.lj_igv_browser_ready = true;
                            console.log('Test IGV browser created successfully');
                        }}).catch(error => {{
                            console.error('Test IGV browser creation failed:', error);
                            window.lj_igv_browser_ready = false;
                        }});
                    }} catch (error) {{
                        console.error('Error in test IGV creation:', error);
                    }}
                """
                
                ui.run_javascript(js_test, timeout=30.0)
                ui.notify("Test IGV browser creation initiated", type="info")
                
            except Exception as e:
                ui.notify(f"Error testing IGV: {e}", type="negative")
                print(f"Error in test IGV: {e}")
        
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
            ui.button("🔄 Manual IGV Load", on_click=lambda: _check_existing_igv_bam())
            ui.button("🧪 Test IGV", on_click=lambda: _test_igv_browser())
            ui.button("📋 Check Console", on_click=lambda: _check_console_errors())
        
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

    # SNP Analysis section
    with ui.card().classes("w-full"):
        ui.label("🧬 SNP Analysis").classes("text-lg font-semibold mb-2")
        
        # SNP Analysis controls and status
        with ui.row().classes("w-full items-center justify-between mb-4"):
            with ui.column().classes("gap-2"):
                ui.label("Variant Calling").classes("text-sm font-medium")
                snp_status_label = ui.label("Ready to run SNP analysis").classes(
                    "text-sm text-gray-600"
                )
                
                # Add file availability status
                snp_files_status = ui.label("Checking required files...").classes(
                    "text-xs text-gray-500"
                )
            
            with ui.column().classes("gap-2 items-end"):
                # SNP Analysis button
                snp_analysis_button = ui.button(
                    "Run SNP Analysis",
                    icon="play_arrow",
                    on_click=lambda: _trigger_snp_analysis()
                ).props("color=primary data-snp-analysis-button")
                
                # Force regenerate checkbox
                force_regenerate_checkbox = ui.checkbox(
                    "Force regenerate",
                    value=False
                ).props("dense")
        
        # SNP Analysis requirements check
        with ui.expansion().classes("w-full").props("icon=info"):
            ui.label("Requirements").classes("text-sm font-medium mb-2")
            with ui.column().classes("gap-1 text-sm"):
                ui.label("• target.bam - Target regions BAM file (generated by target analysis)")
                ui.label("• targets_exceeding_threshold.bed - BED file defining regions (auto-generated by target analysis)")
                ui.label("• Reference genome (hg38 recommended) - provided via CLI --reference option")
                ui.label("• Docker running with hkubal/clairs-to:latest image")
                ui.label("• snpEff and SnpSift installed")
            
            # Add helpful note about file generation
            with ui.row().classes("w-full mt-2 p-2 bg-blue-50 rounded"):
                ui.icon("info").classes("text-blue-600")
                ui.label("Note: Both target.bam and targets_exceeding_threshold.bed are automatically generated when you run target analysis. You don't need to create them manually.").classes("text-sm text-blue-800")
        
        # SNP Analysis results display
        with ui.expansion().classes("w-full").props("icon=assessment"):
            ui.label("Results").classes("text-sm font-medium mb-2")
            
            # Results status
            snp_results_status = ui.label("No SNP analysis results yet").classes(
                "text-sm text-gray-600"
            ).props("data-snp-results-status")
            
            # Results table placeholder
            snp_results_container = ui.column().classes("w-full")
            
                            # Function to check and display SNP results
        def _check_snp_results():
            try:
                clair_dir = sample_dir / "clair3"
                snp_vcf = clair_dir / "snpsift_output.vcf"
                indel_vcf = clair_dir / "snpsift_indel_output.vcf"
                snp_csv = clair_dir / "snpsift_output.vcf.csv"
                indel_csv = clair_dir / "snpsift_indel_output.vcf.csv"
                
                if snp_vcf.exists() and indel_vcf.exists():
                    snp_results_status.set_text("SNP analysis completed successfully!")
                    
                    # Clear previous results
                    snp_results_container.clear()
                    
                    # Display results summary
                    with snp_results_container:
                        with ui.row().classes("w-full gap-4"):
                            with ui.card().classes("flex-1"):
                                ui.label("SNPs").classes("text-sm font-medium")
                                if snp_csv.exists():
                                    try:
                                        snp_df = pd.read_csv(snp_csv)
                                        ui.label(f"Total SNPs: {len(snp_df)}").classes("text-xs text-gray-600")
                                    except Exception:
                                        ui.label("SNP data available").classes("text-xs text-green-600")
                                else:
                                    ui.label("SNP data available").classes("text-xs text-green-600")
                            
                            with ui.card().classes("flex-1"):
                                ui.label("INDELs").classes("text-sm font-medium")
                                if indel_csv.exists():
                                    try:
                                        indel_df = pd.read_csv(indel_csv)
                                        ui.label(f"Total INDELs: {len(indel_df)}").classes("text-xs text-gray-600")
                                    except Exception:
                                        ui.label("INDEL data available").classes("text-xs text-green-600")
                                else:
                                    ui.label("INDEL data available").classes("text-xs text-green-600")
                    
                    # Update button state
                    snp_analysis_button.set_text("Rerun SNP Analysis")
                    snp_analysis_button.props("color=secondary")
                    
                    # Add detailed results viewer
                    with ui.expansion().classes("w-full").props("icon=table_chart"):
                        ui.label("Detailed Results").classes("text-sm font-medium mb-2")
                        
                        # Tabs for SNPs and INDELs
                        with ui.tabs().classes('w-full') as tabs:
                            with ui.tab('SNPs', icon='dna'):
                                _display_variant_table(snp_csv, "SNP", clair_dir)
                            
                            with ui.tab('INDELs', icon='straighten'):
                                _display_variant_table(indel_csv, "INDEL", clair_dir)
                    
                else:
                    snp_results_status.set_text("No SNP analysis results found")
                    
            except Exception as e:
                snp_results_status.set_text(f"Error checking results: {e}")
        
        # Helper function to detect gzipped files
        def _is_gzipped(file_path):
            """Check if a file is gzipped by reading first few bytes"""
            try:
                with open(file_path, 'rb') as f:
                    return f.read(2) == b'\x1f\x8b'
            except:
                return False
        
        # Function to display variant table with filtering and details
        def _display_variant_table(csv_file, variant_type, clair_dir):
            """Display variant data in a comprehensive table with filtering"""
            try:
                # Try to read the VCF file directly instead of the malformed CSV
                vcf_file = None
                if variant_type == "SNP":
                    vcf_file = clair_dir / "snpsift_output.vcf"
                else:
                    vcf_file = clair_dir / "snpsift_indel_output.vcf"
                
                if not vcf_file.exists():
                    ui.label(f"No {variant_type} VCF data available").classes("text-sm text-gray-500")
                    return
                
                # Read the VCF data directly
                df = None
                try:
                    # Check if file is compressed
                    if str(vcf_file).endswith('.gz') or _is_gzipped(vcf_file):
                        import gzip
                        with gzip.open(vcf_file, 'rt') as f:
                            lines = [line.strip() for line in f if not line.startswith('#')]
                    else:
                        with open(vcf_file, 'r') as f:
                            lines = [line.strip() for line in f if not line.startswith('#')]
                    
                    if not lines:
                        ui.label(f"No {variant_type} variants found").classes("text-sm text-gray-500")
                        return
                    
                    # Parse VCF data properly
                    data = []
                    for line in lines:
                        fields = line.split('\t')
                        if len(fields) >= 8:  # Minimum VCF format
                            data.append(fields[:8])  # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
                    
                    # Create DataFrame with proper column names
                    df = pd.DataFrame(data, columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
                    
                except Exception as e:
                    ui.label(f"Error reading {variant_type} VCF: {e}").classes("text-sm text-red-500")
                    return
                
                if df.empty:
                    ui.label(f"No {variant_type} variants found").classes("text-sm text-gray-500")
                    return
                
                # Display summary statistics
                with ui.row().classes("w-full gap-4 mb-4"):
                    with ui.card().classes("flex-1"):
                        ui.label(f"Total {variant_type}s").classes("text-sm font-medium")
                        ui.label(f"{len(df)}").classes("text-2xl font-bold text-blue-600")
                    
                    with ui.card().classes("flex-1"):
                        ui.label("Passing Filters").classes("text-sm font-medium")
                        # Check if we have filter data in the FILTER column
                        if 'FILTER' in df.columns:
                            try:
                                # Count PASS variants
                                passing = sum(1 for val in df['FILTER'] if str(val) == 'PASS')
                                ui.label(f"{passing}").classes("text-2xl font-bold text-green-600")
                            except:
                                ui.label("N/A").classes("text-2xl font-bold text-gray-400")
                        else:
                            ui.label("N/A").classes("text-2xl font-bold text-gray-400")
                    
                    with ui.card().classes("flex-1"):
                        ui.label("Avg Quality").classes("text-sm font-medium")
                        # Check if we have quality data in the QUAL column
                        if 'QUAL' in df.columns:
                            try:
                                # Calculate average quality from QUAL column
                                qual_values = pd.to_numeric(df['QUAL'], errors='coerce')
                                avg_qual = qual_values.mean()
                                if pd.notna(avg_qual):
                                    ui.label(f"{avg_qual:.1f}").classes("text-2xl font-bold text-purple-600")
                                else:
                                    ui.label("N/A").classes("text-2xl font-bold text-gray-400")
                            except:
                                ui.label("N/A").classes("text-2xl font-bold text-gray-400")
                        else:
                            ui.label("N/A").classes("text-2xl font-bold text-gray-400")
                
                # Add filtering controls
                with ui.row().classes("w-full gap-4 mb-4"):
                    # Quality filter
                    quality_threshold = ui.input(
                        value="0",
                        label="Min Quality"
                    ).classes("w-32")
                    
                    # Filter status
                    filter_status = ui.select(
                        options=["All", "PASS", "NonSomatic", "LowQual", "LowAltBQ", "LowAltMQ"],
                        value="All",
                        label="Filter Status"
                    ).classes("w-40")
                    
                    # Apply filters button
                    filter_button = ui.button("Apply Filters", icon="filter_list").classes("w-32")
                
                # Create filtered dataframe
                filtered_df = df.copy()
                
                def apply_filters():
                    nonlocal filtered_df
                    filtered_df = df.copy()
                    
                    # Apply quality filter
                    try:
                        quality_val = float(quality_threshold.value)
                        if quality_val > 0 and 'QUAL' in filtered_df.columns:
                            # Filter by quality in QUAL column
                            qual_values = pd.to_numeric(filtered_df['QUAL'], errors='coerce')
                            quality_mask = qual_values >= quality_val
                            filtered_df = filtered_df[quality_mask]
                    except (ValueError, TypeError):
                        pass  # Invalid quality value, skip filtering
                    
                    # Apply filter status
                    if filter_status.value != "All" and 'FILTER' in filtered_df.columns:
                        # Filter by status in FILTER column
                        filter_mask = filtered_df['FILTER'] == filter_status.value
                        filtered_df = filtered_df[filter_mask]
                    
                    # Update table
                    if variant_table:
                        variant_table.clear()
                        new_table = _create_variant_table(filtered_df, variant_type)
                        if new_table:
                            variant_table.replace(new_table)
                    
                    # Update summary
                    ui.notify(f"Showing {len(filtered_df)} of {len(df)} {variant_type}s")
                
                filter_button.on_click(apply_filters)
                
                # Create conventional table using ui.table.from_pandas
                def _create_variant_table(data_df, v_type):
                    """Create a conventional table from pandas DataFrame"""
                    if data_df.empty:
                        ui.label(f"No {v_type}s match the current filters").classes("text-sm text-gray-500")
                        return None
                    
                    # Prepare DataFrame for display
                    display_df = data_df.copy()
                    
                    # The DataFrame already has proper VCF column names
                    # Select relevant columns for display
                    display_columns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
                    if all(col in display_df.columns for col in display_columns):
                        display_df = display_df[display_columns]
                    else:
                        ui.label(f"Invalid VCF format for {v_type}s").classes("text-sm text-red-500")
                        return None
                        
                    # Clean up data
                    display_df = display_df.replace({np.nan: None})
                    
                    # Add quality color coding
                    def color_quality(qual):
                        try:
                            qual_val = float(qual)
                            if qual_val >= 30:
                                return f"<span style='color: #059669; font-weight: bold;'>{qual_val:.1f}</span>"
                            elif qual_val >= 20:
                                return f"<span style='color: #d97706; font-weight: bold;'>{qual_val:.1f}</span>"
                            else:
                                return f"<span style='color: #dc2626; font-weight: bold;'>{qual_val:.1f}</span>"
                        except:
                            return str(qual)
                    
                    # Apply quality formatting
                    display_df['QUAL'] = display_df['QUAL'].apply(color_quality)
                    
                    # Add filter color coding
                    def color_filter(filter_val):
                        if filter_val == 'PASS':
                            return f"<span style='color: #059669; font-weight: bold;'>{filter_val}</span>"
                        else:
                            return f"<span style='color: #dc2626; font-weight: bold;'>{filter_val}</span>"
                    
                    display_df['FILTER'] = display_df['FILTER'].apply(color_filter)
                    
                    # Truncate INFO field for display
                    def truncate_info(info_val):
                        if info_val and len(str(info_val)) > 50:
                            return str(info_val)[:50] + "..."
                        return str(info_val)
                    
                    display_df['INFO'] = display_df['INFO'].apply(truncate_info)
                    
                    # Create table
                    variant_table = (
                        ui.table.from_pandas(display_df, pagination=25)
                        .props("dense")
                        .style("height: 600px")
                        .style("font-size: 90%; font-weight: 300")
                    )
                    
                    # Make all columns sortable
                    for col in variant_table.columns:
                        col["sortable"] = True
                    
                    # Add table controls
                    with variant_table.add_slot("top-left"):
                        def toggle_fs():
                            variant_table.toggle_fullscreen()
                            fs_button.props(
                                "icon=fullscreen_exit"
                                if variant_table.is_fullscreen
                                else "icon=fullscreen"
                            )
                        
                        fs_button = ui.button(
                            "Toggle fullscreen",
                            icon="fullscreen",
                            on_click=toggle_fs,
                        ).props("flat")
                        
                        # Column visibility toggle
                        with ui.button(icon="menu").props("flat"):
                            with ui.menu(), ui.column().classes("gap-0 p-2"):
                                for column in variant_table.columns:
                                    ui.switch(
                                        column["label"],
                                        value=True,
                                        on_change=lambda e, column=column: toggle_column(column, e.value)
                                    )
                    
                    # Add search functionality
                    with variant_table.add_slot("top-right"):
                        with ui.input(placeholder=f"Search {v_type}s...").props("type=search").bind_value(variant_table, "filter").add_slot("append"):
                            ui.icon("search")
                    
                    # Column toggle function
                    def toggle_column(column: dict, visible: bool) -> None:
                        column["classes"] = "" if visible else "hidden"
                        column["headerClasses"] = "" if visible else "hidden"
                        variant_table.update()
                    
                    return variant_table
                    
                
                # Create and display the table
                variant_table = _create_variant_table(filtered_df, variant_type)
                
                # Add export functionality
                with ui.row().classes("w-full mt-4"):
                    ui.button(
                        f"Export {variant_type}s to CSV",
                        icon="download",
                        on_click=lambda: _export_variants(filtered_df, variant_type)
                    ).classes("w-full")
                    
                    # Add row count display
                    ui.label(f"Showing {len(filtered_df)} of {len(df)} {variant_type}s").classes("text-sm text-gray-600 ml-auto")
                
            except Exception as e:
                ui.label(f"Error displaying {variant_type} table: {e}").classes("text-sm text-red-500")
        
        # Function to show variant details
        def _show_variant_details(variant_row, variant_type, clair_dir):
            """Show detailed information for a specific variant"""
            try:
                with ui.dialog() as dialog, ui.card():
                    ui.label(f"{variant_type} Details").classes("text-lg font-semibold mb-4")
                    
                    # Display variant information
                    with ui.column().classes("gap-2"):
                        for col, val in variant_row.items():
                            if pd.notna(val) and str(val).strip():
                                with ui.row().classes("w-full"):
                                    ui.label(f"{col}:").classes("font-medium w-32")
                                    ui.label(str(val)).classes("flex-1")
                    
                    # Add close button
                    ui.button("Close", on_click=dialog.close).classes("w-full mt-4")
                    
            except Exception as e:
                ui.notify(f"Error showing variant details: {e}", type="error")
        
        # Function to export variants
        def _export_variants(data_df, variant_type):
            """Export filtered variants to CSV"""
            try:
                import tempfile
                import os
                
                # Create temporary file
                with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
                    data_df.to_csv(f.name, index=False)
                    temp_path = f.name
                
                # Download the file
                ui.download(temp_path, filename=f"{variant_type.lower()}_variants.csv")
                
                # Clean up
                os.unlink(temp_path)
                
                ui.notify(f"{variant_type} variants exported successfully!")
                
            except Exception as e:
                ui.notify(f"Export failed: {e}", type="error")
        
        # Check for existing results (after all functions are defined)
        _check_snp_results()
        
        # Function to check file availability
        def _check_snp_file_availability():
            try:
                target_bam = sample_dir / "target.bam"
                targets_bed = sample_dir / "targets_exceeding_threshold.bed"
                
                if not target_bam.exists():
                    snp_files_status.set_text("❌ Missing target.bam - run target analysis first")
                    snp_files_status.classes(replace="text-xs text-red-500")
                    snp_analysis_button.disable()
                    return False
                
                if not targets_bed.exists():
                    snp_files_status.set_text("❌ Missing targets_exceeding_threshold.bed - run target analysis first")
                    snp_files_status.classes(replace="text-xs text-red-500")
                    snp_analysis_button.disable()
                    return False
                
                # Check if BED file has content
                try:
                    with open(targets_bed, 'r') as f:
                        bed_content = f.read().strip()
                    
                    if not bed_content:
                        snp_files_status.set_text("⚠️ BED file is empty - no regions exceed threshold")
                        snp_files_status.classes(replace="text-xs text-yellow-500")
                        snp_analysis_button.disable()
                        return False
                    else:
                        snp_files_status.set_text("✅ All required files available")
                        snp_files_status.classes(replace="text-xs text-green-500")
                        snp_analysis_button.enable()
                        return True
                        
                except Exception as e:
                    snp_files_status.set_text(f"❌ Error reading BED file: {e}")
                    snp_files_status.classes(replace="text-xs text-red-500")
                    snp_analysis_button.disable()
                    return False
                    
            except Exception as e:
                snp_files_status.set_text(f"❌ Error checking files: {e}")
                snp_files_status.classes(replace="text-xs text-red-500")
                snp_analysis_button.disable()
                return False
        
        # Check file availability initially
        _check_snp_file_availability()
        
        # Function to trigger SNP analysis
        def _trigger_snp_analysis():
            try:
                # Check prerequisites
                target_bam = sample_dir / "target.bam"
                targets_bed = sample_dir / "targets_exceeding_threshold.bed"
                
                if not target_bam.exists():
                    ui.notify("Missing target.bam file. Please run target analysis first.", type="warning")
                    return
                
                if not targets_bed.exists():
                    ui.notify("Missing targets_exceeding_threshold.bed file. Please run target analysis first.", type="warning")
                    return
                
                # Check if BED file has content (not empty)
                try:
                    with open(targets_bed, 'r') as f:
                        bed_content = f.read().strip()
                    
                    if not bed_content:
                        ui.notify("targets_exceeding_threshold.bed file is empty. No regions exceed coverage threshold.", type="warning")
                        return
                        
                except Exception as e:
                    ui.notify(f"Error reading targets_exceeding_threshold.bed file: {e}", type="error")
                    return
                
                # Check for reference genome
                # Remove the unnecessary fallback reference check since we get it from CLI
                # The reference genome is passed via CLI and stored in workflow_runner.reference
                
                # Update UI state
                snp_analysis_button.disable()
                snp_status_label.set_text("Starting SNP analysis...")
                snp_status_label.classes(replace="text-sm text-blue-600")
                
                # Create and submit SNP analysis job
                try:
                    # Import required modules
                    from littlejohn.analysis.target_analysis import snp_analysis_handler
                    from littlejohn.workflow_simple import Job, WorkflowContext
                    
                    # Create job metadata with reference genome
                    # Get reference genome from workflow runner (set by CLI)
                    reference_genome = None
                    
                    # Try to get reference from workflow runner
                    print(f"=== REFERENCE GENOME DEBUGGING ===")
                    print(f"launcher type: {type(launcher)}")
                    print(f"launcher dir: {[attr for attr in dir(launcher) if not attr.startswith('_')]}")
                    print(f"hasattr(launcher, 'workflow_runner'): {hasattr(launcher, 'workflow_runner')}")
                    
                    if hasattr(launcher, 'workflow_runner'):
                        print(f"launcher.workflow_runner: {launcher.workflow_runner}")
                        if launcher.workflow_runner:
                            print(f"workflow_runner type: {type(launcher.workflow_runner)}")
                            print(f"workflow_runner dir: {[attr for attr in dir(launcher.workflow_runner) if not attr.startswith('_')]}")
                            print(f"hasattr(workflow_runner, 'reference'): {hasattr(launcher.workflow_runner, 'reference')}")
                            
                            if hasattr(launcher.workflow_runner, 'reference'):
                                reference_genome = launcher.workflow_runner.reference
                                print(f"workflow_runner.reference: {reference_genome}")
                                if reference_genome:
                                    print(f"✅ SUCCESS: Using reference genome from workflow runner: {reference_genome}")
                                else:
                                    print("❌ FAILED: No reference genome in workflow runner")
                            else:
                                print("❌ Workflow runner has no reference attribute")
                        else:
                            print("❌ Workflow runner is None")
                    else:
                        print("❌ No workflow runner attribute found")
                    
                    print(f"=== END REFERENCE GENOME DEBUGGING ===")
                    
                    # Fallback to environment variable
                    if not reference_genome:
                        env_reference = os.environ.get("LITTLEJOHN_REFERENCE")
                        if env_reference and os.path.exists(env_reference):
                            reference_genome = env_reference
                            print(f"Using reference genome from environment: {reference_genome}")
                    
                    if not reference_genome:
                        ui.notify("No reference genome found. SNP calling may fail. Please ensure --reference is provided via CLI.", type="warning")
                    else:
                        ui.notify(f"✅ Reference genome found: {os.path.basename(reference_genome)}", type="positive")
                    
                    # Use monitored_directory as work_dir, or fall back to sample_dir if not available
                    work_dir = launcher.monitored_directory if launcher.monitored_directory else str(sample_dir)
                    
                    metadata = {
                        "work_dir": work_dir,
                        "threads": 4,
                        "force_regenerate": force_regenerate_checkbox.value,
                        "reference": reference_genome
                    }
                    
                    print(f"=== JOB METADATA DEBUGGING ===")
                    print(f"Created metadata: {metadata}")
                    print(f"reference_genome value: {reference_genome}")
                    print(f"=== END JOB METADATA DEBUGGING ===")
                    
                    # Create workflow context
                    context = WorkflowContext(
                        filepath=str(sample_dir),
                        metadata=metadata
                    )
                    
                    # Add sample_id method
                    def get_sample_id():
                        return sample_dir.name
                    
                    context.get_sample_id = get_sample_id
                    context.add_result = lambda key, value: None  # Mock for now
                    context.add_error = lambda key, value: None   # Mock for now
                    
                    # Create job
                    job = Job(
                        job_id=hash(f"snp_analysis_{sample_dir.name}") % 1000000,  # Simple hash-based ID
                        job_type="snp_analysis",
                        context=context,
                        origin="fast",
                        workflow=["fast:snp_analysis"]
                    )
                    
                    # Run SNP analysis in background
                    def run_snp_analysis():
                        try:
                            # Use the work_dir variable we calculated earlier
                            snp_analysis_handler(job, work_dir=work_dir)
                            
                            # Update UI on completion using the main thread
                            # We'll use a different approach to avoid threading issues
                            print("SNP analysis completed successfully!")
                            
                        except Exception as e:
                            # Log error instead of trying to update UI from background thread
                            ui.run_javascript(f"""
                                // Update status
                                var statusElement = document.querySelector('{snp_status_label.id}');
                                if (statusElement) {{
                                    statusElement.textContent = "SNP analysis failed: " + "{str(e)}";
                                    statusElement.className = "text-sm text-red-600";
                                }}
                                
                                // Re-enable button
                                var buttonElement = document.querySelector('{snp_analysis_button.id}');
                                if (buttonElement) {{
                                    buttonElement.disabled = false;
                                }}
                            """)
                    
                    # Run in background thread
                    import threading
                    thread = threading.Thread(target=run_snp_analysis, daemon=True)
                    thread.start()
                    
                    ui.notify("SNP analysis started in background", type="info")
                    
                except Exception as e:
                    ui.notify(f"Failed to start SNP analysis: {e}", type="error")
                    snp_analysis_button.enable()
                    snp_status_label.set_text("Failed to start SNP analysis")
                    snp_status_label.classes(replace="text-sm text-red-600")
                    
            except Exception as e:
                ui.notify(f"Error triggering SNP analysis: {e}", type="error")
                snp_analysis_button.enable()
                snp_status_label.set_text("Error triggering SNP analysis")
                snp_status_label.classes(replace="text-sm text-red-600")

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
            
            # SNP Analysis results check - only update if necessary
            try:
                if state.get("last_snp_check", 0) < time.time() - 30:  # Check every 30 seconds
                    # Check for SNP analysis results
                    clair_dir = sample_dir / "clair3"
                    snp_vcf = clair_dir / "snpsift_output.vcf"
                    indel_vcf = clair_dir / "snpsift_indel_output.vcf"
                    
                    if snp_vcf.exists() and indel_vcf.exists():
                        # Update SNP results status if it exists
                        try:
                            # Find the SNP results status label and update it
                            snp_results_elements = document.querySelectorAll('[data-snp-results-status]')
                            if snp_results_elements.length > 0:
                                for element in snp_results_elements:
                                    element.textContent = "SNP analysis completed successfully!"
                                    element.className = "text-sm text-green-600"
                        except Exception:
                            pass
                        
                        # Update button state if it exists
                        try:
                            snp_button_elements = document.querySelectorAll('[data-snp-analysis-button]')
                            if snp_button_elements.length > 0:
                                for element in snp_button_elements:
                                    element.textContent = "Rerun SNP Analysis"
                                    element.disabled = false
                                    element.className = "q-btn q-btn--standard q-btn--rectangle q-btn--secondary"
                        except Exception:
                            pass
                    
                    # Also check file availability for SNP analysis
                    try:
                        target_bam = sample_dir / "target.bam"
                        targets_bed = sample_dir / "targets_exceeding_threshold.bed"
                        
                        if target_bam.exists() and targets_bed.exists():
                            # Check if BED file has content
                            try:
                                with open(targets_bed, 'r') as f:
                                    bed_content = f.read().strip()
                                
                                if bed_content:
                                    # Files are available, enable button if not already enabled
                                    try:
                                        snp_button_elements = document.querySelectorAll('[data-snp-analysis-button]')
                                        if snp_button_elements.length > 0:
                                            for element in snp_button_elements:
                                                if element.disabled:
                                                    element.disabled = false;
                                                    element.className = "q-btn q-btn--standard q-btn--rectangle q-btn--primary";
                                    except Exception:
                                        pass
                            except Exception:
                                pass
                    except Exception:
                        pass
                    
                    state["last_snp_check"] = time.time()
            except Exception:
                # Never let SNP logic break coverage refresh
                pass
            
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
