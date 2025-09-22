from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional
import json
import csv
from datetime import datetime
import hashlib
import logging
import asyncio
import concurrent.futures

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

from robin.gui.config import get_confidence_level as get_classifier_confidence_level


def add_summary_section(sample_dir: Path, sample_id: str) -> None:
    """Build the Summary section at the top of the sample detail page."""

    # State management for the summary section
    state = {
        "sample_dir": sample_dir,
        "sample_id": sample_id,
        "last_data_hash": None,
        "run_info_labels": {},
        "classification_labels": {},
        "coverage_labels": {},
        "cnv_labels": {},
        "mgmt_labels": {},
        "fusion_labels": {},
        "refresh_timer": None,
    }

    with ui.card().classes("w-full"):
        ui.label("Run Details").classes("text-lg font-semibold text-blue-800")
        ui.separator().classes().style("border: 1px solid var(--md-primary)")# Run Information Section
        with ui.row().classes("w-full gap-4"):
            state["run_info_labels"]["run"] = _create_dashboard_card(
                "Run Time", "Loading...", "schedule", "Sequencing run timestamp"
            )
            state["run_info_labels"]["model"] = _create_dashboard_card(
                "Basecall Model", "Loading...", "settings", "AI model used for basecalling"
            )
            state["run_info_labels"]["device"] = _create_dashboard_card(
                "Device", "Loading...", "smartphone", "Sequencing device identifier"
            )
            state["run_info_labels"]["flow_cell"] = _create_dashboard_card(
                "Flow Cell", "Loading...", "tag", "Flow cell identification"
            )
            state["run_info_labels"]["panel"] = _create_dashboard_card(
                "Analysis Panel", "Loading...", "science", "Target panel used for analysis"
            )
            state["run_info_labels"]["sample"] = _create_dashboard_card(
                "Sample ID", sample_id, "play_arrow", "Unique sample identifier"
            )

    with ui.card().classes("w-full"):
        ui.label("Classification Details").classes("text-lg font-semibold text-blue-800")
        ui.separator().classes().style("border: 1px solid var(--md-primary)")
        # Classification Summary Section
        with ui.row().classes("w-full gap-4"):
            state["classification_labels"]["sturgeon"] = _create_classification_dashboard_card(
                "Sturgeon", "Loading...", "psychology", "AI brain tumor classification"
            )
            state["classification_labels"]["nanodx"] = _create_classification_dashboard_card(
                "NanoDX", "Loading...", "biotech", "Molecular cancer diagnostics"
            )
            state["classification_labels"]["pannanodx"] = _create_classification_dashboard_card(
                "PanNanoDX", "Loading...", "science", "Pan-cancer analysis"
            )
            state["classification_labels"]["random_forest"] = _create_classification_dashboard_card(
                "Random Forest", "Loading...", "forest", "Machine learning classification"
            )

    # Analysis Details Section
    with ui.card().classes("w-full"):
        ui.label("Analysis Details").classes("text-lg font-semibold text-blue-800")
        ui.separator().classes().style("border: 1px solid var(--md-primary)")
        with ui.row().classes("w-full gap-4"):
            state["coverage_labels"] = _create_coverage_dashboard_card()
            state["cnv_labels"] = _create_cnv_dashboard_card()
            state["mgmt_labels"] = _create_mgmt_dashboard_card()
            state["fusion_labels"] = _create_fusion_dashboard_card()

    # Async data loading functions
    async def _initial_load_async():
        """Load initial data asynchronously."""
        try:
            await _refresh_summary_data_async()
        except Exception as e:
            logging.exception(f"[Summary] Initial load failed: {e}")

    async def _refresh_summary_data_async():
        """Refresh all summary data asynchronously if source files have changed."""
        try:
            # Run data hash calculation in background thread
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future = executor.submit(_calculate_data_hash, state["sample_dir"])
                current_hash = await asyncio.wrap_future(future)
            
            if current_hash == state["last_data_hash"]:
                return  # No changes, skip update

            # Update the data (these functions are already synchronous and safe to call)
            _update_run_information(state)
            _update_classification_data(state)
            _update_analysis_data(state)

            # Update the hash
            state["last_data_hash"] = current_hash

            logging.debug(f"[Summary] Refreshed data for {state['sample_id']}")

        except Exception as e:
            logging.exception(f"[Summary] Refresh failed: {e}")

    # Start the refresh timer (every 30 seconds) with async function
    state["refresh_timer"] = ui.timer(30.0, lambda: ui.timer(0.1, _refresh_summary_data_async, once=True), active=True)

    # Load initial data asynchronously
    ui.timer(0.1, _initial_load_async, once=True)


def _create_dashboard_card(title: str, value: str, icon: str, description: str) -> Dict[str, Any]:
    """Create a compact dashboard-style card matching the design pattern. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label(title).classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-pink-100 rounded-full flex items-center justify-center"):
                ui.icon(icon).classes("w-3.5 h-3.5 text-pink-600")
        
        # Main value
        value_label = ui.label(value).classes("text-xl font-bold text-gray-900 mb-1")
        
        # Description
        ui.label(description).classes("text-xs text-gray-500")
        
        return {"value": value_label}


def _create_classification_dashboard_card(title: str, classification: str, icon: str, description: str) -> Dict[str, Any]:
    """Create a compact classification dashboard card with detailed information. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label(title).classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-blue-100 rounded-full flex items-center justify-center"):
                ui.icon(icon).classes("w-3.5 h-3.5 text-blue-600")
        
        # Main classification result with confidence badge
        with ui.row().classes("flex items-center justify-between mb-1"):
            classification_label = ui.label(classification).classes("text-xl font-bold text-gray-900")
            # Confidence level badge
            confidence_badge = ui.label("Loading...").classes("px-2 py-1 text-xs font-medium rounded-full bg-gray-100 text-gray-600")
        
        # Compact details in a single row
        with ui.row().classes("flex items-center justify-between mb-1"):
            confidence_label = ui.label("Confidence: Loading...").classes("text-xs font-medium text-gray-700")
            features_label = ui.label("Features: Loading...").classes("text-xs text-gray-500")
        
        # Description
        ui.label(description).classes("text-xs text-gray-500")
        
        return {
            "classification": classification_label,
            "confidence": confidence_label,
            "confidence_badge": confidence_badge,
            "features": features_label
        }


def _create_classification_card(
    title: str,
    classification: str,
    confidence: str,
    confidence_level: str,
    features: str = "",
    features_label: str = "",
    model: str = "",
) -> Dict[str, Any]:
    """Create a classification summary card. Returns labels for updating."""
    labels = {}
    with ui.card().classes("flex-1 elevation-4 rounded-xl bg-gradient-to-br from-blue-50 to-indigo-50 border-l-4 border-blue-500"):
        ui.label(f"{title}").classes("font-bold text-blue-800 mb-2")
        ui.separator().classes().style("border: 1px solid var(--md-primary)")
        labels["classification"] = ui.label(f"Class: {classification}").classes("font-bold text-medium text-blue-600")
        confidence_color = _get_confidence_color(confidence_level)
        labels["confidence"] = ui.label(f"Confidence: {confidence}%").classes(f"text-sm text-{confidence_color}-600")
        labels["confidence_level"] = ui.label(confidence_level).classes(f"text-sm text-{confidence_color}-600")
        if model:
            ui.label(f"Model: {model}").classes("text-sm text-blue-600")
        if features and features_label:
            labels["features"] = ui.label(f"Features: {features_label}").classes("text-sm text-blue-600")

    return labels


def _create_coverage_dashboard_card() -> Dict[str, Any]:
    """Create a compact coverage analysis dashboard card. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label("Coverage Analysis").classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-blue-100 rounded-full flex items-center justify-center"):
                ui.icon("analytics").classes("w-3.5 h-3.5 text-blue-600")
        
        # Main quality result with badge
        with ui.row().classes("flex items-center justify-between mb-1"):
            quality_label = ui.label("Loading...").classes("text-xl font-bold text-gray-900")
            # Coverage badge
            coverage_badge = ui.label("--x").classes("px-2 py-1 text-xs font-medium rounded-full bg-gray-100 text-gray-600")
        
        # Coverage details in compact layout
        with ui.column().classes("mb-2"):
            global_coverage_label = ui.label("Global: Loading...").classes("text-xs text-gray-600")
            target_coverage_label = ui.label("Targets: Loading...").classes("text-xs text-gray-600")
            enrichment_label = ui.label("Enrichment: Loading...").classes("text-xs text-gray-600")
        
        # Coverage thresholds as small badges
        with ui.row().classes("gap-1 flex-wrap"):
            ui.label("≥30x").classes("px-1 py-0.5 text-xs bg-green-100 text-green-800 rounded")
            ui.label("≥20x").classes("px-1 py-0.5 text-xs bg-blue-100 text-blue-800 rounded")
            ui.label("≥10x").classes("px-1 py-0.5 text-xs bg-yellow-100 text-yellow-800 rounded")
            ui.label("<10x").classes("px-1 py-0.5 text-xs bg-red-100 text-red-800 rounded")
        
        return {
            "quality": quality_label,
            "coverage_badge": coverage_badge,
            "global_coverage": global_coverage_label,
            "target_coverage": target_coverage_label,
            "enrichment": enrichment_label
        }


def _create_cnv_dashboard_card() -> Dict[str, Any]:
    """Create a compact CNV analysis dashboard card. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label("Copy Number Analysis").classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-purple-100 rounded-full flex items-center justify-center"):
                ui.icon("person").classes("w-3.5 h-3.5 text-purple-600")
        
        # Main genetic sex result
        genetic_sex_label = ui.label("Loading...").classes("text-xl font-bold text-gray-900 mb-1")
        
        # Analysis details in compact layout
        with ui.column().classes("mb-2"):
            bin_width_label = ui.label("Bin Width: Loading...").classes("text-xs text-gray-600")
            variance_label = ui.label("Variance: Loading...").classes("text-xs text-gray-600")
        
        # CNV counts as badges
        with ui.row().classes("gap-2 mb-1"):
            gained_badge = ui.label("Gained: --").classes("px-2 py-1 text-xs font-medium rounded-full bg-green-100 text-green-800")
            lost_badge = ui.label("Lost: --").classes("px-2 py-1 text-xs font-medium rounded-full bg-red-100 text-red-800")
        
        # Description
        ui.label("Copy number analysis across genome with breakpoint detection").classes("text-xs text-gray-500")
        
        return {
            "genetic_sex": genetic_sex_label,
            "bin_width": bin_width_label,
            "variance": variance_label,
            "gained": gained_badge,
            "lost": lost_badge
        }


def _create_mgmt_dashboard_card() -> Dict[str, Any]:
    """Create a compact MGMT analysis dashboard card. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label("MGMT Analysis").classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-orange-100 rounded-full flex items-center justify-center"):
                ui.icon("science").classes("w-3.5 h-3.5 text-orange-600")
        
        # Main status result with badge
        with ui.row().classes("flex items-center justify-between mb-1"):
            status_label = ui.label("Loading...").classes("text-xl font-bold text-gray-900")
            # Methylation badge
            methylation_badge = ui.label("--%").classes("px-2 py-1 text-xs font-medium rounded-full bg-gray-100 text-gray-600")
        
        # Analysis details in compact layout
        with ui.column().classes("mb-2"):
            average_methylation_label = ui.label("Average: Loading...").classes("text-xs text-gray-600")
            prediction_score_label = ui.label("Score: Loading...").classes("text-xs text-gray-600")
        
        # Description
        cpg_sites_label = ui.label("MGMT status determined from methylation analysis of -- CpG sites").classes("text-xs text-gray-500")
        
        return {
            "status": status_label,
            "methylation_badge": methylation_badge,
            "average_methylation": average_methylation_label,
            "prediction_score": prediction_score_label,
            "cpg_sites": cpg_sites_label
        }


def _create_fusion_dashboard_card() -> Dict[str, Any]:
    """Create a compact fusion analysis dashboard card. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label("Fusion Analysis").classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-green-100 rounded-full flex items-center justify-center"):
                ui.icon("merge").classes("w-3.5 h-3.5 text-green-600")
        
        # Panel info and main result
        with ui.row().classes("flex items-center justify-between mb-1"):
            panel_label = ui.label("Panel: --").classes("text-sm font-medium text-gray-700")
            target_fusions_badge = ui.label("-- target fusions").classes("px-2 py-1 text-xs font-medium rounded-full bg-blue-100 text-blue-800")
        
        # Analysis details in compact layout
        with ui.column().classes("mb-2"):
            genome_fusions_label = ui.label("-- genome wide fusions").classes("text-xs text-gray-600")
        
        # Description
        ui.label("Fusion candidates identified from reads with supplementary alignments").classes("text-xs text-gray-500")
        
        return {
            "panel": panel_label,
            "target_fusions": target_fusions_badge,
            "genome_fusions": genome_fusions_label
        }




def _calculate_data_hash(sample_dir: Path) -> str:
    """Calculate a hash of key data files to detect changes."""
    try:
        hash_data = []

        # Key files to monitor for changes
        key_files = [
            "sturgeon_scores.csv",
            "NanoDX_scores.csv",
            "PanNanoDX_scores.csv",
            "random_forest_scores.csv",
            "coverage_summary.csv",
            "bed_coverage_main.csv",
            "cnv_analysis_counter.txt",
            "1_mgmt.csv",
            "sv_count.txt",
        ]

        for filename in key_files:
            file_path = sample_dir / filename
            if file_path.exists():
                try:
                    # Get file modification time and size
                    stat = file_path.stat()
                    hash_data.append(f"{filename}:{stat.st_mtime}:{stat.st_size}")
                except Exception:
                    pass

        # Create hash from all data (excluding directory modification time for stability)
        hash_string = "|".join(sorted(hash_data))
        return hashlib.md5(hash_string.encode()).hexdigest()

    except Exception:
        return "unknown"


def _update_run_information(state: Dict[str, Any]) -> None:
    """Update run information labels."""
    try:
        run_info = _extract_run_information(state["sample_dir"], state["sample_id"])

        if "run" in state["run_info_labels"]:
            state["run_info_labels"]["run"]["value"].text = run_info.get(
                "run_time", "Unknown"
            )
        if "model" in state["run_info_labels"]:
            state["run_info_labels"]["model"]["value"].text = run_info.get(
                "model", "Unknown"
            )
        if "device" in state["run_info_labels"]:
            state["run_info_labels"]["device"]["value"].text = run_info.get(
                "device", "Unknown"
            )
        if "flow_cell" in state["run_info_labels"]:
            state["run_info_labels"]["flow_cell"]["value"].text = run_info.get(
                "flow_cell", "Unknown"
            )
        if "panel" in state["run_info_labels"]:
            state["run_info_labels"]["panel"]["value"].text = run_info.get(
                "panel", "Unknown"
            )

    except Exception as e:
        logging.exception(f"[Summary] Failed to update run information: {e}")


def _update_classification_data(state: Dict[str, Any]) -> None:
    """Update classification data labels."""
    try:
        classification_data = _extract_classification_data(state["sample_dir"])

        # Update Sturgeon
        if "sturgeon" in state["classification_labels"]:
            labels = state["classification_labels"]["sturgeon"]
            data = classification_data.get("sturgeon", {})
            if "classification" in labels:
                classification = data.get('classification', 'Unknown')
                labels["classification"].text = classification
            if "confidence" in labels:
                confidence = data.get('confidence', 0)
                labels["confidence"].text = f"Confidence: {confidence:.1f}%"
            if "confidence_badge" in labels:
                confidence = data.get('confidence', 0)
                confidence_level = _get_confidence_level(confidence, "sturgeon")
                labels["confidence_badge"].text = confidence_level
                # Update badge color using style method
                if confidence >= 80:
                    labels["confidence_badge"].style("background-color: #dcfce7; color: #166534;")  # green
                elif confidence >= 50:
                    labels["confidence_badge"].style("background-color: #fef3c7; color: #92400e;")  # yellow
                else:
                    labels["confidence_badge"].style("background-color: #fee2e2; color: #991b1b;")  # red
            if "features" in labels:
                features = data.get('features', 0)
                labels["features"].text = f"Features: {features:,}"

        # Update NanoDX
        if "nanodx" in state["classification_labels"]:
            labels = state["classification_labels"]["nanodx"]
            data = classification_data.get("nanodx", {})
            if "classification" in labels:
                classification = data.get('classification', 'Unknown')
                labels["classification"].text = classification
            if "confidence" in labels:
                confidence = data.get('confidence', 0)
                labels["confidence"].text = f"Confidence: {confidence:.1f}%"
            if "confidence_badge" in labels:
                confidence = data.get('confidence', 0)
                confidence_level = _get_confidence_level(confidence, "nanodx")
                labels["confidence_badge"].text = confidence_level
                # Update badge color using style method
                if confidence >= 80:
                    labels["confidence_badge"].style("background-color: #dcfce7; color: #166534;")  # green
                elif confidence >= 50:
                    labels["confidence_badge"].style("background-color: #fef3c7; color: #92400e;")  # yellow
                else:
                    labels["confidence_badge"].style("background-color: #fee2e2; color: #991b1b;")  # red
            if "features" in labels:
                features = data.get('features', 0)
                labels["features"].text = f"Features: {features:,}"

        # Update PanNanoDX
        if "pannanodx" in state["classification_labels"]:
            labels = state["classification_labels"]["pannanodx"]
            data = classification_data.get("pannanodx", {})
            if "classification" in labels:
                classification = data.get('classification', 'Unknown')
                labels["classification"].text = classification
            if "confidence" in labels:
                confidence = data.get('confidence', 0)
                labels["confidence"].text = f"Confidence: {confidence:.1f}%"
            if "confidence_badge" in labels:
                confidence = data.get('confidence', 0)
                confidence_level = _get_confidence_level(confidence, "pannanodx")
                labels["confidence_badge"].text = confidence_level
                # Update badge color using style method
                if confidence >= 80:
                    labels["confidence_badge"].style("background-color: #dcfce7; color: #166534;")  # green
                elif confidence >= 50:
                    labels["confidence_badge"].style("background-color: #fef3c7; color: #92400e;")  # yellow
                else:
                    labels["confidence_badge"].style("background-color: #fee2e2; color: #991b1b;")  # red
            if "features" in labels:
                features = data.get('features', 0)
                labels["features"].text = f"Features: {features:,}"

        # Update Random Forest
        if "random_forest" in state["classification_labels"]:
            labels = state["classification_labels"]["random_forest"]
            data = classification_data.get("random_forest", {})
            if "classification" in labels:
                classification = data.get('classification', 'Unknown')
                labels["classification"].text = classification
            if "confidence" in labels:
                confidence = data.get('confidence', 0)
                labels["confidence"].text = f"Confidence: {confidence:.1f}%"
            if "confidence_badge" in labels:
                confidence = data.get('confidence', 0)
                confidence_level = _get_confidence_level(confidence, "random_forest")
                labels["confidence_badge"].text = confidence_level
                # Update badge color using style method
                if confidence >= 80:
                    labels["confidence_badge"].style("background-color: #dcfce7; color: #166534;")  # green
                elif confidence >= 50:
                    labels["confidence_badge"].style("background-color: #fef3c7; color: #92400e;")  # yellow
                else:
                    labels["confidence_badge"].style("background-color: #fee2e2; color: #991b1b;")  # red
            if "features" in labels:
                features = data.get('features', 0)
                labels["features"].text = f"Features: {features:,}"

    except Exception as e:
        logging.exception(f"[Summary] Failed to update classification data: {e}")


def _update_analysis_data(state: Dict[str, Any]) -> None:
    """Update analysis data labels."""
    try:
        # Update coverage data
        coverage_data = _extract_coverage_data(state["sample_dir"])
        if "coverage_labels" in state:
            labels = state["coverage_labels"]
            if "quality" in labels:
                quality_status = coverage_data.get("quality", "Unknown")
                labels["quality"].text = quality_status
            if "coverage_badge" in labels:
                coverage_value = coverage_data.get("target_coverage", "0.0x")
                labels["coverage_badge"].text = coverage_value
                # Update badge color using style method
                try:
                    coverage_num = float(coverage_value.replace("x", ""))
                    if coverage_num >= 30:
                        labels["coverage_badge"].style("background-color: #dcfce7; color: #166534;")  # green
                    elif coverage_num >= 20:
                        labels["coverage_badge"].style("background-color: #dbeafe; color: #1e40af;")  # blue
                    elif coverage_num >= 10:
                        labels["coverage_badge"].style("background-color: #fef3c7; color: #92400e;")  # yellow
                    else:
                        labels["coverage_badge"].style("background-color: #fee2e2; color: #991b1b;")  # red
                except (ValueError, AttributeError):
                    # Handle "Not available" or other non-numeric values
                    labels["coverage_badge"].style("background-color: #f3f4f6; color: #6b7280;")  # gray
            if "global_coverage" in labels:
                labels["global_coverage"].text = f"Global: {coverage_data.get('global_coverage', '0.0x')}"
            if "target_coverage" in labels:
                labels["target_coverage"].text = f"Targets: {coverage_data.get('target_coverage', '0.0x')}"
            if "enrichment" in labels:
                labels["enrichment"].text = f"Enrichment: {coverage_data.get('enrichment', '0.0x')}"

        # Update CNV data
        cnv_data = _extract_cnv_data(state["sample_dir"])
        if "cnv_labels" in state:
            labels = state["cnv_labels"]
            if "genetic_sex" in labels:
                labels["genetic_sex"].text = cnv_data.get('genetic_sex', 'Unknown')
            if "bin_width" in labels:
                labels["bin_width"].text = f"Bin Width: {cnv_data.get('bin_width', 'Unknown')}"
            if "variance" in labels:
                labels["variance"].text = f"Variance: {cnv_data.get('variance', 'Unknown')}"
            if "gained" in labels:
                labels["gained"].text = f"Gained: {cnv_data.get('gained', 0)}"
            if "lost" in labels:
                labels["lost"].text = f"Lost: {cnv_data.get('lost', 0)}"

        # Update MGMT data
        mgmt_data = _extract_mgmt_data(state["sample_dir"])
        if "mgmt_labels" in state:
            labels = state["mgmt_labels"]
            if "status" in labels:
                methylation_status = mgmt_data.get("status", "Unknown")
                labels["status"].text = methylation_status
            if "methylation_badge" in labels:
                methylation_value = mgmt_data.get("methylation_percent", "0.0%")
                labels["methylation_badge"].text = methylation_value
                # Update badge color using style method
                try:
                    meth_num = float(methylation_value.replace("%", ""))
                    if meth_num > 10:
                        labels["methylation_badge"].style("background-color: #dcfce7; color: #166534;")  # green
                    elif meth_num > 5:
                        labels["methylation_badge"].style("background-color: #fef3c7; color: #92400e;")  # yellow
                    else:
                        labels["methylation_badge"].style("background-color: #fee2e2; color: #991b1b;")  # red
                except (ValueError, AttributeError):
                    # Handle "Not available" or other non-numeric values
                    labels["methylation_badge"].style("background-color: #f3f4f6; color: #6b7280;")  # gray
            if "average_methylation" in labels:
                labels["average_methylation"].text = f"Average: {mgmt_data.get('average_methylation', '0.0%')}"
            if "prediction_score" in labels:
                labels["prediction_score"].text = f"Score: {mgmt_data.get('prediction_score', '0.0%')}"
            if "cpg_sites" in labels:
                labels["cpg_sites"].text = f"MGMT status determined from methylation analysis of {mgmt_data.get('cpg_sites', 0)} CpG sites"

        # Update fusion data
        fusion_data = _extract_fusion_data(state["sample_dir"])
        if "fusion_labels" in state:
            labels = state["fusion_labels"]
            if "target_fusions" in labels:
                labels["target_fusions"].text = f"{fusion_data.get('target_fusions', 0)} target fusions"
            if "genome_fusions" in labels:
                labels["genome_fusions"].text = f"{fusion_data.get('genome_fusions', 0)} genome wide fusions"
            if "panel" in labels:
                panel = _get_analysis_panel(state["sample_dir"])
                labels["panel"].text = f"Panel: {panel}"

    except Exception as e:
        logging.exception(f"[Summary] Failed to update analysis data: {e}")


def _get_analysis_panel(sample_dir: Path) -> str:
    """Get the analysis panel from master.csv"""
    try:
        master_csv_path = sample_dir / "master.csv"
        if master_csv_path.exists():
            import pandas as pd
            df = pd.read_csv(master_csv_path)
            if not df.empty and "analysis_panel" in df.columns:
                panel = df.iloc[0]["analysis_panel"]
                if panel and str(panel).strip() != "":
                    return str(panel).strip()
        return "rCNS2"  # Default fallback
    except Exception:
        return "rCNS2"  # Default fallback


def _extract_run_information(sample_dir: Path, sample_id: str) -> Dict[str, str]:
    """Extract run information from sample directory."""
    run_info = {
        "run_time": "Not available",
        "model": "Not available",
        "device": "Not available",
        "flow_cell": "Not available",
        "panel": "Not available",
    }

    try:
        # Try to get run time from file modification times
        if sample_dir.exists():
            latest_time = 0
            for file_path in sample_dir.rglob("*"):
                if file_path.is_file():
                    mtime = file_path.stat().st_mtime
                    if mtime > latest_time:
                        latest_time = mtime
            if latest_time > 0:
                run_info["run_time"] = datetime.fromtimestamp(latest_time).strftime(
                    "%Y-%m-%d %H:%M"
                )

        # Try to extract from various metadata files (search recursively)
        metadata_files = [
            "sample_metadata.json",
            "run_info.json",
            "sequencing_summary.txt",
            "master.csv",
        ]

        for metadata_file in metadata_files:
            candidate_paths: List[Path] = []
            direct_path = sample_dir / metadata_file
            if direct_path.exists():
                candidate_paths = [direct_path]
            else:
                candidate_paths = list(sample_dir.rglob(metadata_file))

            for file_path in candidate_paths:
                if metadata_file.endswith(".json"):
                    try:
                        with open(file_path, "r") as f:
                            data = json.load(f)
                            if "model" in data and data["model"]:
                                run_info["model"] = str(data["model"]).strip()
                            if "device" in data and data["device"]:
                                run_info["device"] = str(data["device"]).strip()
                            if "flow_cell" in data and data["flow_cell"]:
                                run_info["flow_cell"] = str(data["flow_cell"]).strip()
                    except Exception:
                        pass
                elif metadata_file == "master.csv":
                    try:
                        with open(file_path, "r") as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                # Helper to fetch by case-insensitive key
                                def get_ci(
                                    r: Dict[str, str], key: str
                                ) -> Optional[str]:
                                    for k, v in r.items():
                                        if k and k.strip().lower() == key:
                                            return v
                                    return None

                                model_val = get_ci(row, "run_info_model") or get_ci(
                                    row, "basecall_models"
                                )
                                if model_val and model_val.strip():
                                    run_info["model"] = model_val.strip()

                                device_val = get_ci(row, "run_info_device") or get_ci(
                                    row, "devices"
                                )
                                if device_val and device_val.strip():
                                    run_info["device"] = device_val.strip()

                                flowcell_val = get_ci(
                                    row, "run_info_flow_cell"
                                ) or get_ci(row, "flowcell_ids")
                                if flowcell_val and flowcell_val.strip():
                                    run_info["flow_cell"] = flowcell_val.strip()

                                panel_val = get_ci(row, "analysis_panel")
                                if panel_val and panel_val.strip():
                                    run_info["panel"] = panel_val.strip()

                                run_time_val = get_ci(row, "run_info_run_time")
                                if run_time_val and run_time_val.strip():
                                    try:
                                        dt = datetime.fromisoformat(
                                            run_time_val.replace("Z", "+00:00")
                                        )
                                        run_info["run_time"] = dt.strftime(
                                            "%Y-%m-%d %H:%M"
                                        )
                                    except Exception:
                                        pass
                    except Exception:
                        pass
                # Stop early if we have all fields
                if all(
                    run_info[k] != "Not available"
                    for k in ["run_time", "model", "device", "flow_cell"]
                ):
                    break
            # Stop early if complete
            if all(
                run_info[k] != "Not available"
                for k in ["run_time", "model", "device", "flow_cell"]
            ):
                break
    except Exception:
        pass

    return run_info


def _extract_coverage_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract coverage analysis data using the same logic as the coverage component."""
    coverage_data = {
        "quality": "Not available",
        "global_coverage": "Not available",
        "target_coverage": "Not available",
        "enrichment": "Not available",
    }

    try:
        # Use the same file reading logic as the coverage component
        cov_main = sample_dir / "coverage_main.csv"
        bed_cov = sample_dir / "bed_coverage_main.csv"

        global_cov = None
        target_cov_v = None
        enrich_v = None

        # Read coverage_main.csv for global coverage (same as coverage component)
        if cov_main.exists():
            try:
                import pandas as pd

                cov_df = pd.read_csv(cov_main)

                # Calculate global coverage using the same formula as coverage component
                if (
                    "covbases" in cov_df.columns
                    and "endpos" in cov_df.columns
                    and cov_df["endpos"].sum() > 0
                ):
                    global_cov = float(cov_df["covbases"].sum()) / float(
                        cov_df["endpos"].sum()
                    )
                    coverage_data["global_coverage"] = f"{global_cov:.2f}x"
            except Exception as e:
                logging.debug(f"   MGMT: <access denied>: {e}")
                pass

        # Read bed_coverage_main.csv for target coverage (same as coverage component)
        if bed_cov.exists():
            try:
                import pandas as pd

                bed_df = pd.read_csv(bed_cov)

                # Calculate target coverage using the same formula as coverage component
                if "bases" in bed_df.columns:
                    # Calculate length if not present (same as coverage component)
                    if "length" not in bed_df.columns:
                        bed_df["length"] = (
                            bed_df["endpos"] - bed_df["startpos"] + 1
                        ).astype(float)

                    if bed_df["length"].sum() > 0:
                        target_cov_v = float(bed_df["bases"].sum()) / float(
                            bed_df["length"].sum()
                        )
                        coverage_data["target_coverage"] = f"{target_cov_v:.2f}x"

                        # Determine quality based on target coverage (same logic as coverage component)
                        if target_cov_v >= 30:
                            coverage_data["quality"] = "Excellent"
                        elif target_cov_v >= 20:
                            coverage_data["quality"] = "Good"
                        elif target_cov_v >= 10:
                            coverage_data["quality"] = "Moderate"
                        else:
                            coverage_data["quality"] = "Insufficient"
            except Exception as e:
                logging.debug(f"   MGMT: <access denied>: {e}")
                pass

        # Calculate enrichment if we have both values (same logic as coverage component)
        if global_cov is not None and target_cov_v is not None and global_cov > 0:
            enrich_v = target_cov_v / global_cov
            coverage_data["enrichment"] = f"{enrich_v:.2f}x"

    except Exception as e:
        logging.debug(f"   MGMT: <access denied>: {e}")
        pass

    return coverage_data


def _extract_cnv_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract CNV analysis data using the same methods as the CNV component."""
    cnv_data = {
        "genetic_sex": "Not available",
        "bin_width": "Not available",
        "variance": "Not available",
        "gained": 0,
        "lost": 0,
    }

    try:
        # Load CNV_dict.npy for bin width and variance (same as CNV component)
        cnv_dict_npy = sample_dir / "CNV_dict.npy"
        if cnv_dict_npy.exists():
            try:
                import numpy as np

                cnv_dict = np.load(cnv_dict_npy, allow_pickle=True).item()
                if isinstance(cnv_dict, dict):
                    if "bin_width" in cnv_dict:
                        bin_width = cnv_dict["bin_width"]
                        if isinstance(bin_width, (int, float)):
                            cnv_data["bin_width"] = f"{bin_width:,}"
                    if "variance" in cnv_dict:
                        variance = cnv_dict["variance"]
                        if isinstance(variance, (int, float)):
                            cnv_data["variance"] = f"{variance:.3f}"
            except Exception as e:
                logging.debug(f"   MGMT: <access denied>: {e}")
                pass

        # Load XYestimate.pkl for genetic sex (same as CNV component)
        xy_pkl = sample_dir / "XYestimate.pkl"
        if xy_pkl.exists():
            try:
                import pickle

                with xy_pkl.open("rb") as f:
                    xy = pickle.load(f)
                if xy:
                    cnv_data["genetic_sex"] = str(xy)
            except Exception as e:
                logging.debug(f"   MGMT: <access denied>: {e}")
                pass

        # Look for CNV result files for gained/lost counts
        cnv_files = [
            "cnv_analysis_counter.txt",
            "cnv_results.csv",
            "cnv_summary.txt",
            "cnv_analysis_results.pkl",
        ]

        for cnv_file in cnv_files:
            file_path = sample_dir / cnv_file
            if file_path.exists():
                if cnv_file.endswith(".txt"):
                    try:
                        with open(file_path, "r") as f:
                            content = f.read().strip()
                            # Parse counter file for basic CNV info
                            if content.isdigit():
                                cnv_data["gained"] = int(
                                    content
                                )  # Assume this is gained regions
                                cnv_data["lost"] = 0  # Default value
                    except Exception as e:
                        logging.debug(f"   MGMT: <access denied>: {e}")
                        pass
                elif cnv_file.endswith(".csv"):
                    try:
                        with open(file_path, "r") as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                if "genetic_sex" in row:
                                    cnv_data["genetic_sex"] = row["genetic_sex"]
                                if "bin_width" in row:
                                    cnv_data["bin_width"] = row["bin_width"]
                                if "variance" in row:
                                    cnv_data["variance"] = row["variance"]
                                if "gained" in row:
                                    cnv_data["gained"] = int(row["gained"])
                                if "lost" in row:
                                    cnv_data["lost"] = int(row["lost"])
                    except Exception as e:
                        logging.debug(f"   MGMT: <access denied>: {e}")
                        pass
                break

    except Exception as e:
        logging.debug(f"   MGMT: <access denied>: {e}")
        pass

    return cnv_data


def _extract_mgmt_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract MGMT analysis data using the same logic as the MGMT component."""
    mgmt_data = {
        "status": "Not available",
        "methylation_percent": "Not available",
        "average_methylation": "Not available",
        "prediction_score": "Not available",
        "cpg_sites": "Not available",
    }

    try:
        # Use the same file discovery logic as the MGMT component
        csv_files = list(sample_dir.glob("*_mgmt.csv"))
        if not csv_files:
            return mgmt_data

        # Find the latest MGMT CSV file (same logic as MGMT component)
        def _count_from_name(p: Path) -> int:
            try:
                return int(p.name.split("_")[0])
            except Exception as e:
                logging.debug(f"   MGMT: <access denied>: {e}")
                return -1

        latest_csv = max(csv_files, key=_count_from_name)

        # Read the latest MGMT CSV file
        try:
            import pandas as pd

            df = pd.read_csv(latest_csv)

            # Extract data using the same column names as the MGMT component
            if "status" in df.columns:
                status = str(df.get("status", pd.Series(["Unknown"])).iloc[0])
                mgmt_data["status"] = status
            elif "pred" in df.columns:
                # Convert prediction score to status (same logic as MGMT component)
                try:
                    pred_val = float(df.get("pred", pd.Series([0.0])).iloc[0])
                    if pred_val > 0.1:
                        mgmt_data["status"] = "Methylated"
                    else:
                        mgmt_data["status"] = "Unmethylated"
                except Exception as e:
                    logging.debug(f"   MGMT: <access denied>: {e}")
                    pass

            if "average" in df.columns:
                try:
                    avg_val = float(df.get("average", pd.Series([0.0])).iloc[0])
                    mgmt_data["average_methylation"] = f"{avg_val:.1f}%"
                    mgmt_data["methylation_percent"] = f"{avg_val:.1f}%"
                except Exception as e:
                    logging.debug(f"   MGMT: <access denied>: {e}")
                    pass

            if "pred" in df.columns:
                try:
                    pred_val = float(df.get("pred", pd.Series([0.0])).iloc[0])
                    mgmt_data["prediction_score"] = f"{pred_val:.1f}%"
                except Exception as e:
                    logging.debug(f"   MGMT: <access denied>: {e}")
                    pass

            # Count CPG sites from the BED file if available
            current_count = _count_from_name(latest_csv)
            bed_path = sample_dir / f"{current_count}_mgmt.bed"
            if not bed_path.exists():
                alt_bed = sample_dir / f"{current_count}_mgmt_mgmt.bed"
                bed_path = alt_bed if alt_bed.exists() else bed_path

            if bed_path.exists():
                try:
                    # Count lines in BED file to estimate CPG sites
                    with open(bed_path, "r") as f:
                        line_count = sum(1 for line in f if line.strip())
                    mgmt_data["cpg_sites"] = line_count
                except Exception as e:
                    logging.debug(f"   MGMT: <access denied>: {e}")
                    pass

        except Exception as e:
            logging.debug(f"   MGMT: <access denied>: {e}")
            pass

    except Exception as e:
        logging.debug(f"   MGMT: <access denied>: {e}")
        pass

    return mgmt_data


def _extract_classification_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract classification data from various score files."""
    classification_data = {
        "sturgeon": {
            "classification": "Not available",
            "confidence": 0.0,
            "confidence_level": "Not available",
            "features": 0,
        },
        "nanodx": {
            "classification": "Not available",
            "confidence": 0.0,
            "confidence_level": "Not available",
            "features": 0,
        },
        "pannanodx": {
            "classification": "Not available",
            "confidence": 0.0,
            "confidence_level": "Not available",
            "features": 0,
        },
        "random_forest": {
            "classification": "Not available",
            "confidence": 0.0,
            "confidence_level": "Not available",
            "features": 0,
        },
    }

    try:
        # Extract Sturgeon data
        sturgeon_file = sample_dir / "sturgeon_scores.csv"
        if sturgeon_file.exists():
            try:
                with open(sturgeon_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        # Get the latest row (last one)
                        features = int(row.get("number_probes", 0))
                        # Find the highest scoring classification
                        max_score = 0.0
                        best_class = "Unknown"
                        for col, value in row.items():
                            if col not in ["timestamp", "number_probes"]:
                                try:
                                    score = float(value)
                                    if score > max_score:
                                        max_score = score
                                        best_class = col
                                except Exception as e:
                                    logging.debug(f"   Sturgeon: <access denied>: {e}")
                                    pass

                        classification_data["sturgeon"] = {
                            "classification": best_class,
                            "confidence": max_score * 100,
                            "confidence_level": _get_confidence_level(max_score * 100, "sturgeon"),
                            "features": features,
                        }
            except Exception as e:
                logging.debug(f"   Sturgeon: <access denied>: {e}")
                pass

        # Extract NanoDX data
        nanodx_file = sample_dir / "NanoDX_scores.csv"
        if nanodx_file.exists():
            try:
                with open(nanodx_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        # Get the latest row (last one)
                        features = int(row.get("number_probes", 0))
                        # Find the highest scoring classification
                        max_score = 0.0
                        best_class = "Unknown"
                        for col, value in row.items():
                            if col not in ["timestamp", "number_probes"]:
                                try:
                                    score = float(value)
                                    if score > max_score:
                                        max_score = score
                                        best_class = col
                                except Exception as e:
                                    logging.debug(f"   NanoDX: <access denied>: {e}")
                                    pass

                        classification_data["nanodx"] = {
                            "classification": best_class,
                            "confidence": max_score * 100,
                            "confidence_level": _get_confidence_level(max_score * 100, "nanodx"),
                            "features": features,
                        }
            except Exception as e:
                logging.debug(f"   NanoDX: <access denied>: {e}")
                pass

        # Extract PanNanoDX data
        pannanodx_file = sample_dir / "PanNanoDX_scores.csv"
        if pannanodx_file.exists():
            try:
                with open(pannanodx_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        # Get the latest row (last one)
                        features = int(row.get("number_probes", 0))
                        # Find the highest scoring classification
                        max_score = 0.0
                        best_class = "Unknown"
                        for col, value in row.items():
                            if col not in ["timestamp", "number_probes"]:
                                try:
                                    score = float(value)
                                    if score > max_score:
                                        max_score = score
                                        best_class = col
                                except Exception as e:
                                    logging.debug(f"   PanNanoDX: <access denied>: {e}")
                                    pass

                        classification_data["pannanodx"] = {
                            "classification": best_class,
                            "confidence": max_score * 100,
                            "confidence_level": _get_confidence_level(max_score * 100, "pannanodx"),
                            "features": features,
                        }
            except Exception as e:
                logging.debug(f"   PanNanoDX: <access denied>: {e}")
                pass

        # Extract Random Forest data
        rf_file = sample_dir / "random_forest_scores.csv"
        if rf_file.exists():
            try:
                with open(rf_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        # Get the latest row (last one)
                        features = int(row.get("number_probes", 0))
                        # Find the highest scoring classification
                        max_score = 0.0
                        best_class = "Unknown"
                        for col, value in row.items():
                            if col not in ["timestamp", "number_probes"]:
                                try:
                                    score = float(value)
                                    if score > max_score:
                                        max_score = score
                                        best_class = col
                                except Exception as e:
                                    logging.debug(f"   Random Forest: <access denied>: {e}")
                                    pass

                        # Some random forest outputs are already in percent (0-100),
                        # while others are fractional (0-1). Normalize to percent.
                        confidence_percent = (
                            max_score * 100 if max_score <= 1 else max_score
                        )
                        classification_data["random_forest"] = {
                            "classification": best_class,
                            "confidence": confidence_percent,
                            "confidence_level": _get_confidence_level(
                                confidence_percent, "random_forest"
                            ),
                            "features": features,
                        }
            except Exception as e:
                logging.debug(f"   Random Forest: <access denied>: {e}")
                pass

    except Exception as e:
        logging.debug(f"   Random Forest: <access denied>: {e}")
        pass

    return classification_data


def _extract_fusion_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract fusion analysis data."""
    fusion_data = {"target_fusions": "Not available", "genome_fusions": "Not available"}

    try:
        # Look for fusion result files
        fusion_files = [
            "sv_count.txt",
            "fusion_results.csv",
            "fusion_summary.txt",
            "fusion_analysis.csv",
        ]

        for fusion_file in fusion_files:
            file_path = sample_dir / fusion_file
            if file_path.exists():
                if fusion_file.endswith(".txt"):
                    try:
                        with open(file_path, "r") as f:
                            content = f.read().strip()
                            # Parse counter file for basic fusion info
                            if content.isdigit():
                                fusion_data["genome_fusions"] = int(
                                    content
                                )  # Assume this is genome-wide fusions
                                fusion_data["target_fusions"] = 0  # Default value
                    except Exception as e:
                        logging.debug(f"   Fusion: <access denied>: {e}")
                        pass
                elif fusion_file.endswith(".csv"):
                    try:
                        with open(file_path, "r") as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                if "target_fusions" in row:
                                    try:
                                        fusion_data["target_fusions"] = int(
                                            row["target_fusions"]
                                        )
                                    except Exception as e:
                                        logging.debug(f"   Fusion: <access denied>: {e}")
                                        pass
                                if "genome_fusions" in row:
                                    try:
                                        fusion_data["genome_fusions"] = int(
                                            row["genome_fusions"]
                                        )
                                    except Exception as e:
                                        logging.debug(f"   Fusion: <access denied>: {e}")
                                        pass
                    except Exception as e:
                        logging.debug(f"   Fusion: <access denied>: {e}")
                        pass
                break

    except Exception:
        pass

    return fusion_data


def _get_confidence_color(confidence_level: str) -> str:
    """Get the appropriate color class for confidence level."""
    if "High" in confidence_level:
        return "text-green-600"
    elif "Medium" in confidence_level:
        return "text-blue-600"
    elif "Low" in confidence_level:
        return "text-orange-600"
    else:
        return "text-gray-600"


def _get_quality_color(quality: str) -> str:
    """Get the appropriate color class for quality status."""
    if quality == "Excellent":
        return "text-green-600"
    elif quality == "Good":
        return "text-blue-600"
    elif quality == "Moderate":
        return "text-orange-600"
    elif quality == "Insufficient":
        return "text-red-600"
    else:
        return "text-gray-600"


def _get_coverage_badge_color(coverage: str) -> str:
    """Get the appropriate badge color for coverage value."""
    try:
        cov_val = float(coverage.replace("x", ""))
        if cov_val >= 30:
            return "positive"
        elif cov_val >= 20:
            return "primary"
        elif cov_val >= 10:
            return "warning"
        else:
            return "negative"
    except Exception as e:
        logging.debug(f"   Coverage: <access denied>: {e}")
        return "grey"


def _get_methylation_color(status: str) -> str:
    """Get the appropriate color class for methylation status."""
    if "Methylated" in status:
        return "text-green-600"
    elif "Unmethylated" in status:
        return "text-orange-600"
    else:
        return "text-gray-600"


def _get_confidence_level(confidence: float, classifier: str = "default") -> str:
    """Get the confidence level based on classifier-specific thresholds."""
    return get_classifier_confidence_level(classifier, confidence)


def _get_confidence_badge_classes(confidence: float) -> str:
    """Get CSS classes for confidence level badge based on confidence score."""
    if confidence >= 80:
        return "px-2 py-1 text-xs font-medium rounded-full bg-green-100 text-green-800"
    elif confidence >= 50:
        return "px-2 py-1 text-xs font-medium rounded-full bg-yellow-100 text-yellow-800"
    else:
        return "px-2 py-1 text-xs font-medium rounded-full bg-red-100 text-red-800"


def _get_methylation_badge_color(methylation: str) -> str:
    """Get the appropriate badge color for methylation value."""
    try:
        meth_val = float(methylation.replace("%", ""))
        if meth_val > 10:
            return "positive"
        elif meth_val > 5:
            return "warning"
        else:
            return "negative"
    except Exception as e:
        logging.debug(f"   Methylation: <access denied>: {e}")
        return "grey"


