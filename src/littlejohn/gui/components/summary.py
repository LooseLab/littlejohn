from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional
import json
import csv
from datetime import datetime
import hashlib
import logging

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None


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
        "refresh_timer": None
    }
    
    # Run Information Section
    with ui.card().classes("w-full bg-gradient-to-r from-blue-50 to-indigo-50"):
        ui.label("Run Information").classes("text-lg font-semibold mb-3 text-blue-800")
        
        with ui.row().classes("w-full gap-8 items-center flex-wrap"):
            state["run_info_labels"]["run"] = _create_info_item("Run", "Loading...", "schedule")
            state["run_info_labels"]["model"] = _create_info_item("Model", "Loading...", "settings")
            state["run_info_labels"]["device"] = _create_info_item("Device", "Loading...", "smartphone")
            state["run_info_labels"]["flow_cell"] = _create_info_item("Flow Cell", "Loading...", "tag")
            state["run_info_labels"]["sample"] = _create_info_item("Sample", sample_id, "play_arrow")
    
    # Classification Summary Section
    with ui.card().classes("w-full"):
        ui.label("Classification Summary").classes("text-lg font-semibold mb-3")
        
        with ui.row().classes("w-full gap-4 flex-wrap"):
            state["classification_labels"]["sturgeon"] = _create_classification_card(
                "Sturgeon classification", 
                "Loading...",
                "0.0%",
                "Unknown",
                "0",
                "Features"
            )
            state["classification_labels"]["nanodx"] = _create_classification_card(
                "NanoDX classification", 
                "Loading...",
                "0.0%",
                "Unknown",
                "0",
                "Features",
                "Capper_et_al_NN.pkl"
            )
            state["classification_labels"]["pannanodx"] = _create_classification_card(
                "NanoDX classification", 
                "Loading...",
                "0.0%",
                "Unknown",
                "0",
                "Features",
                "pancan_devel_v5i_NN.pkl"
            )
            state["classification_labels"]["random_forest"] = _create_classification_card(
                "Forest classification", 
                "Loading...",
                "0.0%",
                "Unknown",
                "0",
                "Features"
            )
    
    # Analysis Details Section
    with ui.card().classes("w-full"):
        ui.label("Analysis Details").classes("text-lg font-semibold mb-3")
        
        with ui.row().classes("w-full gap-4 flex-wrap"):
            state["coverage_labels"] = _create_coverage_analysis_card()
            state["cnv_labels"] = _create_cnv_analysis_card()
            state["mgmt_labels"] = _create_mgmt_analysis_card()
            state["fusion_labels"] = _create_fusion_analysis_card()
    
    # Initial data load
    def _initial_load():
        """Load initial data immediately."""
        try:
            _refresh_summary_data()
        except Exception as e:
            logging.exception(f"[Summary] Initial load failed: {e}")
    
    # Refresh function
    def _refresh_summary_data():
        """Refresh all summary data if source files have changed."""
        try:
            # Check if data has actually changed
            current_hash = _calculate_data_hash(state["sample_dir"])
            if current_hash == state["last_data_hash"]:
                return  # No changes, skip update
            
            # Update the data
            _update_run_information(state)
            _update_classification_data(state)
            _update_analysis_data(state)
            
            # Update the hash
            state["last_data_hash"] = current_hash
            
            logging.debug(f"[Summary] Refreshed data for {state['sample_id']}")
            
        except Exception as e:
            logging.exception(f"[Summary] Refresh failed: {e}")
    
    # Start the refresh timer (every 30 seconds)
    state["refresh_timer"] = ui.timer(30.0, _refresh_summary_data, active=True)
    
    # Load initial data
    ui.timer(0.1, _initial_load, once=True)


def _create_info_item(label: str, value: str, icon: str) -> Dict[str, Any]:
    """Create a single information item with icon and value. Returns labels for updating."""
    with ui.column().classes("items-center min-w-24"):
        ui.icon(icon).classes("text-2xl text-blue-600 mb-1")
        ui.label(label).classes("text-xs text-gray-600")
        value_label = ui.label(value).classes("text-sm font-semibold text-gray-800")
        return {"value": value_label}


def _create_classification_card(
    title: str, 
    classification: str, 
    confidence: str, 
    confidence_level: str, 
    features: str = "", 
    features_label: str = "", 
    model: str = ""
) -> Dict[str, Any]:
    """Create a classification summary card. Returns labels for updating."""
    labels = {}
    
    with ui.card().classes("flex-1 min-w-64 bg-gray-50"):
        ui.label(title).classes("text-sm font-semibold text-gray-800 mb-1")
        # Display the current top classification for this classifier
        labels["classification"] = ui.label(classification).classes("text-base font-semibold text-gray-900 mb-2")
        
        if model:
            ui.label(f"Model: {model}").classes("text-xs text-gray-600 mb-1")
        
        if features and features_label:
            labels["features"] = ui.label(f"{features_label}: {features}").classes("text-xs text-gray-600 mb-2")
        
        # Confidence percentage with color coding
        confidence_color = _get_confidence_color(confidence_level)
        labels["confidence"] = ui.label(confidence).classes(f"text-2xl font-bold {confidence_color} mb-1")
        
        labels["confidence_level"] = ui.label(confidence_level).classes("text-xs text-gray-600")
    
    return labels


def _create_coverage_analysis_card() -> Dict[str, Any]:
    """Create the coverage analysis card. Returns labels for updating."""
    labels = {}
    
    with ui.card().classes("flex-1 min-w-64"):
        # Header with quality status
        with ui.row().classes("items-center justify-between mb-2"):
            labels["quality"] = ui.label("Quality: Loading...").classes("text-sm font-semibold text-gray-600")
            labels["coverage_badge"] = ui.badge("--x").props("color=grey")
        
        ui.label("Coverage Depths").classes("text-xs font-semibold text-gray-700 mb-2")
        
        with ui.column().classes("gap-1 mb-3"):
            labels["global_coverage"] = ui.label("Global Estimated Coverage: --x").classes("text-xs text-gray-600")
            labels["target_coverage"] = ui.label("Targets Estimated Coverage: --x").classes("text-xs text-gray-600")
            labels["enrichment"] = ui.label("Estimated enrichment: --x").classes("text-xs text-gray-600")
        
        # Coverage thresholds
        with ui.row().classes("gap-2 flex-wrap"):
            ui.badge("≥30x Excellent").props("color=positive size=sm")
            ui.badge("≥20x Good").props("color=primary size=sm")
            ui.badge("≥10x Moderate").props("color=warning size=sm")
            ui.badge("<10x Insufficient").props("color=negative size=sm")
    
    return labels


def _create_cnv_analysis_card() -> Dict[str, Any]:
    """Create the CNV analysis card. Returns labels for updating."""
    labels = {}
    
    with ui.card().classes("flex-1 min-w-64"):
        with ui.row().classes("items-center gap-2 mb-2"):
            ui.icon("person").classes("text-blue-600")
            labels["genetic_sex"] = ui.label("Genetic Sex: Loading...").classes("text-sm font-semibold text-blue-600")
        
        ui.label("Analysis Details").classes("text-xs font-semibold text-gray-700 mb-2")
        
        with ui.column().classes("gap-1 mb-3"):
            labels["bin_width"] = ui.label("Bin Width: Loading...").classes("text-xs text-gray-600")
            labels["variance"] = ui.label("Variance: Loading...").classes("text-xs text-gray-600")
        
        with ui.row().classes("gap-2 mb-2"):
            labels["gained"] = ui.badge("Gained: --").props("color=positive")
            labels["lost"] = ui.badge("Lost: --").props("color=negative")
        
        ui.label("Copy number analysis across genome with breakpoint detection").classes("text-xs text-gray-600")
    
    return labels


def _create_mgmt_analysis_card() -> Dict[str, Any]:
    """Create the MGMT analysis card. Returns labels for updating."""
    labels = {}
    
    with ui.card().classes("flex-1 min-w-64"):
        # Header with methylation status
        with ui.row().classes("items-center justify-between mb-2"):
            labels["status"] = ui.label("Status: Loading...").classes("text-sm font-semibold text-gray-600")
            labels["methylation_badge"] = ui.badge("--%").props("color=grey")
        
        ui.label("Analysis Details").classes("text-xs font-semibold text-gray-700 mb-2")
        
        with ui.column().classes("gap-1 mb-3"):
            labels["average_methylation"] = ui.label("Average Methylation: --%").classes("text-xs text-gray-600")
            labels["prediction_score"] = ui.label("Prediction Score: --%").classes("text-xs text-gray-600")
        
        labels["cpg_sites"] = ui.label("MGMT status determined from methylation analysis of -- CpG sites").classes("text-xs text-gray-600")
    
    return labels


def _create_fusion_analysis_card() -> Dict[str, Any]:
    """Create the fusion analysis card. Returns labels for updating."""
    labels = {}
    
    with ui.card().classes("flex-1 min-w-64"):
        with ui.row().classes("items-center justify-between mb-2"):
            ui.label("Panel: rCNS2").classes("text-sm font-semibold text-gray-800")
            labels["target_fusions"] = ui.badge("-- between target fusions").props("color=primary")
        
        ui.label("Analysis Details").classes("text-xs font-semibold text-gray-700 mb-2")
        
        with ui.column().classes("gap-1 mb-3"):
            labels["genome_fusions"] = ui.label("-- genome wide fusions").classes("text-xs text-gray-600")
        
        ui.label("Fusion candidates identified from reads with supplementary alignments").classes("text-xs text-gray-600")
    
    return labels


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
            "sv_count.txt"
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
            state["run_info_labels"]["run"]["value"].text = run_info.get("run_time", "Unknown")
        if "model" in state["run_info_labels"]:
            state["run_info_labels"]["model"]["value"].text = run_info.get("model", "Unknown")
        if "device" in state["run_info_labels"]:
            state["run_info_labels"]["device"]["value"].text = run_info.get("device", "Unknown")
        if "flow_cell" in state["run_info_labels"]:
            state["run_info_labels"]["flow_cell"]["value"].text = run_info.get("flow_cell", "Unknown")
            
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
                labels["classification"].text = data.get("classification", "Unknown")
            if "confidence" in labels:
                labels["confidence"].text = f"{data.get('confidence', 0):.1f}%"
                # Update confidence color
                confidence_level = data.get("confidence_level", "Unknown")
                confidence_color = _get_confidence_color(confidence_level)
                labels["confidence"].classes(f"text-2xl font-bold {confidence_color} mb-1")
            if "confidence_level" in labels:
                labels["confidence_level"].text = data.get("confidence_level", "Unknown")
            if "features" in labels:
                labels["features"].text = f"Features: {data.get('features', 0):,}"
        
        # Update NanoDX
        if "nanodx" in state["classification_labels"]:
            labels = state["classification_labels"]["nanodx"]
            data = classification_data.get("nanodx", {})
            if "classification" in labels:
                labels["classification"].text = data.get("classification", "Unknown")
            if "confidence" in labels:
                labels["confidence"].text = f"{data.get('confidence', 0):.1f}%"
                confidence_level = data.get("confidence_level", "Unknown")
                confidence_color = _get_confidence_color(confidence_level)
                labels["confidence"].classes(f"text-2xl font-bold {confidence_color} mb-1")
            if "confidence_level" in labels:
                labels["confidence_level"].text = data.get("confidence_level", "Unknown")
            if "features" in labels:
                labels["features"].text = f"Features: {data.get('features', 0):,}"
        
        # Update PanNanoDX
        if "pannanodx" in state["classification_labels"]:
            labels = state["classification_labels"]["pannanodx"]
            data = classification_data.get("pannanodx", {})
            if "classification" in labels:
                labels["classification"].text = data.get("classification", "Unknown")
            if "confidence" in labels:
                labels["confidence"].text = f"{data.get('confidence', 0):.1f}%"
                confidence_level = data.get("confidence_level", "Unknown")
                confidence_color = _get_confidence_color(confidence_level)
                labels["confidence"].classes(f"text-2xl font-bold {confidence_color} mb-1")
            if "confidence_level" in labels:
                labels["confidence_level"].text = data.get("confidence_level", "Unknown")
            if "features" in labels:
                labels["features"].text = f"Features: {data.get('features', 0):,}"
        
        # Update Random Forest
        if "random_forest" in state["classification_labels"]:
            labels = state["classification_labels"]["random_forest"]
            data = classification_data.get("random_forest", {})
            if "classification" in labels:
                labels["classification"].text = data.get("classification", "Unknown")
            if "confidence" in labels:
                labels["confidence"].text = f"{data.get('confidence', 0):.1f}%"
                confidence_level = data.get("confidence_level", "Unknown")
                confidence_color = _get_confidence_color(confidence_level)
                labels["confidence"].classes(f"text-2xl font-bold {confidence_color} mb-1")
            if "confidence_level" in labels:
                labels["confidence_level"].text = data.get("confidence_level", "Unknown")
            if "features" in labels:
                labels["features"].text = f"Features: {data.get('features', 0):,}"
                
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
                quality_color = _get_quality_color(quality_status)
                labels["quality"].text = f"Quality: {quality_status}"
                labels["quality"].classes(f"text-sm font-semibold {quality_color}")
            if "coverage_badge" in labels:
                coverage_value = coverage_data.get("target_coverage", "0.0x")
                labels["coverage_badge"].props(f"color={_get_coverage_badge_color(coverage_value)}")
                labels["coverage_badge"].text = coverage_value
            if "global_coverage" in labels:
                labels["global_coverage"].text = f"Global Estimated Coverage: {coverage_data.get('global_coverage', '0.0x')}"
            if "target_coverage" in labels:
                labels["target_coverage"].text = f"Targets Estimated Coverage: {coverage_data.get('target_coverage', '0.0x')}"
            if "enrichment" in labels:
                labels["enrichment"].text = f"Estimated enrichment: {coverage_data.get('enrichment', '0.0x')}"
        
        # Update CNV data
        cnv_data = _extract_cnv_data(state["sample_dir"])
        if "cnv_labels" in state:
            labels = state["cnv_labels"]
            if "genetic_sex" in labels:
                labels["genetic_sex"].text = f"Genetic Sex: {cnv_data.get('genetic_sex', 'Unknown')}"
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
                methylation_color = _get_methylation_color(methylation_status)
                labels["status"].text = f"Status: {methylation_status}"
                labels["status"].classes(f"text-sm font-semibold {methylation_color}")
            if "methylation_badge" in labels:
                methylation_value = mgmt_data.get("methylation_percent", "0.0%")
                labels["methylation_badge"].props(f"color={_get_methylation_badge_color(methylation_value)}")
                labels["methylation_badge"].text = methylation_value
            if "average_methylation" in labels:
                labels["average_methylation"].text = f"Average Methylation: {mgmt_data.get('average_methylation', '0.0%')}"
            if "prediction_score" in labels:
                labels["prediction_score"].text = f"Prediction Score: {mgmt_data.get('prediction_score', '0.0%')}"
            if "cpg_sites" in labels:
                labels["cpg_sites"].text = f"MGMT status determined from methylation analysis of {mgmt_data.get('cpg_sites', 0)} CpG sites"
        
        # Update fusion data
        fusion_data = _extract_fusion_data(state["sample_dir"])
        if "fusion_labels" in state:
            labels = state["fusion_labels"]
            if "target_fusions" in labels:
                labels["target_fusions"].text = f"{fusion_data.get('target_fusions', 0)} between target fusions"
            if "genome_fusions" in labels:
                labels["genome_fusions"].text = f"{fusion_data.get('genome_fusions', 0)} genome wide fusions"
                
    except Exception as e:
        logging.exception(f"[Summary] Failed to update analysis data: {e}")


def _extract_run_information(sample_dir: Path, sample_id: str) -> Dict[str, str]:
    """Extract run information from sample directory."""
    run_info = {
        "run_time": "Not available",
        "model": "Not available", 
        "device": "Not available",
        "flow_cell": "Not available"
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
                run_info["run_time"] = datetime.fromtimestamp(latest_time).strftime("%Y-%m-%d %H:%M")

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
                if metadata_file.endswith('.json'):
                    try:
                        with open(file_path, 'r') as f:
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
                        with open(file_path, 'r') as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                # Helper to fetch by case-insensitive key
                                def get_ci(r: Dict[str, str], key: str) -> Optional[str]:
                                    for k, v in r.items():
                                        if k and k.strip().lower() == key:
                                            return v
                                    return None

                                model_val = get_ci(row, "run_info_model") or get_ci(row, "basecall_models")
                                if model_val and model_val.strip():
                                    run_info["model"] = model_val.strip()

                                device_val = get_ci(row, "run_info_device") or get_ci(row, "devices")
                                if device_val and device_val.strip():
                                    run_info["device"] = device_val.strip()

                                flowcell_val = get_ci(row, "run_info_flow_cell") or get_ci(row, "flowcell_ids")
                                if flowcell_val and flowcell_val.strip():
                                    run_info["flow_cell"] = flowcell_val.strip()

                                run_time_val = get_ci(row, "run_info_run_time")
                                if run_time_val and run_time_val.strip():
                                    try:
                                        dt = datetime.fromisoformat(run_time_val.replace('Z', '+00:00'))
                                        run_info["run_time"] = dt.strftime("%Y-%m-%d %H:%M")
                                    except Exception:
                                        pass
                    except Exception:
                        pass
                # Stop early if we have all fields
                if all(run_info[k] != "Not available" for k in ["run_time", "model", "device", "flow_cell"]):
                    break
            # Stop early if complete
            if all(run_info[k] != "Not available" for k in ["run_time", "model", "device", "flow_cell"]):
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
        "enrichment": "Not available"
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
                if ("covbases" in cov_df.columns and 
                    "endpos" in cov_df.columns and 
                    cov_df["endpos"].sum() > 0):
                    global_cov = float(cov_df["covbases"].sum()) / float(cov_df["endpos"].sum())
                    coverage_data["global_coverage"] = f"{global_cov:.2f}x"
            except Exception:
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
                        bed_df["length"] = (bed_df["endpos"] - bed_df["startpos"] + 1).astype(float)
                    
                    if bed_df["length"].sum() > 0:
                        target_cov_v = float(bed_df["bases"].sum()) / float(bed_df["length"].sum())
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
            except Exception:
                pass
        
        # Calculate enrichment if we have both values (same logic as coverage component)
        if (global_cov is not None and 
            target_cov_v is not None and 
            global_cov > 0):
            enrich_v = target_cov_v / global_cov
            coverage_data["enrichment"] = f"{enrich_v:.2f}x"
                        
    except Exception:
        pass
    
    return coverage_data


def _extract_cnv_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract CNV analysis data using the same methods as the CNV component."""
    cnv_data = {
        "genetic_sex": "Not available",
        "bin_width": "Not available",
        "variance": "Not available",
        "gained": 0,
        "lost": 0
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
            except Exception:
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
            except Exception:
                pass
        
        # Look for CNV result files for gained/lost counts
        cnv_files = [
            "cnv_analysis_counter.txt",
            "cnv_results.csv",
            "cnv_summary.txt",
            "cnv_analysis_results.pkl"
        ]
        
        for cnv_file in cnv_files:
            file_path = sample_dir / cnv_file
            if file_path.exists():
                if cnv_file.endswith('.txt'):
                    try:
                        with open(file_path, 'r') as f:
                            content = f.read().strip()
                            # Parse counter file for basic CNV info
                            if content.isdigit():
                                cnv_data["gained"] = int(content)  # Assume this is gained regions
                                cnv_data["lost"] = 0  # Default value
                    except:
                        pass
                elif cnv_file.endswith('.csv'):
                    try:
                        with open(file_path, 'r') as f:
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
                    except:
                        pass
                break
                        
    except Exception:
        pass
    
    return cnv_data


def _extract_mgmt_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract MGMT analysis data using the same logic as the MGMT component."""
    mgmt_data = {
        "status": "Not available",
        "methylation_percent": "Not available",
        "average_methylation": "Not available",
        "prediction_score": "Not available",
        "cpg_sites": "Not available"
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
            except Exception:
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
                except:
                    pass
            
            if "average" in df.columns:
                try:
                    avg_val = float(df.get("average", pd.Series([0.0])).iloc[0])
                    mgmt_data["average_methylation"] = f"{avg_val:.1f}%"
                    mgmt_data["methylation_percent"] = f"{avg_val:.1f}%"
                except:
                    pass
            
            if "pred" in df.columns:
                try:
                    pred_val = float(df.get("pred", pd.Series([0.0])).iloc[0])
                    mgmt_data["prediction_score"] = f"{pred_val:.1f}%"
                except:
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
                    with open(bed_path, 'r') as f:
                        line_count = sum(1 for line in f if line.strip())
                    mgmt_data["cpg_sites"] = line_count
                except:
                    pass
                    
        except Exception:
            pass
                        
    except Exception:
        pass
    
    return mgmt_data


def _extract_classification_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract classification data from various score files."""
    classification_data = {
        "sturgeon": {"classification": "Not available", "confidence": 0.0, "confidence_level": "Not available", "features": 0},
        "nanodx": {"classification": "Not available", "confidence": 0.0, "confidence_level": "Not available", "features": 0},
        "pannanodx": {"classification": "Not available", "confidence": 0.0, "confidence_level": "Not available", "features": 0},
        "random_forest": {"classification": "Not available", "confidence": 0.0, "confidence_level": "Not available", "features": 0}
    }
    
    try:
        # Extract Sturgeon data
        sturgeon_file = sample_dir / "sturgeon_scores.csv"
        if sturgeon_file.exists():
            try:
                with open(sturgeon_file, 'r') as f:
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
                                except:
                                    pass
                        
                        classification_data["sturgeon"] = {
                            "classification": best_class,
                            "confidence": max_score * 100,
                            "confidence_level": _get_confidence_level(max_score * 100),
                            "features": features
                        }
            except:
                pass
        
        # Extract NanoDX data
        nanodx_file = sample_dir / "NanoDX_scores.csv"
        if nanodx_file.exists():
            try:
                with open(nanodx_file, 'r') as f:
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
                                except:
                                    pass
                        
                        classification_data["nanodx"] = {
                            "classification": best_class,
                            "confidence": max_score * 100,
                            "confidence_level": _get_confidence_level(max_score * 100),
                            "features": features
                        }
            except:
                pass
        
        # Extract PanNanoDX data
        pannanodx_file = sample_dir / "PanNanoDX_scores.csv"
        if pannanodx_file.exists():
            try:
                with open(pannanodx_file, 'r') as f:
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
                                except:
                                    pass
                        
                        classification_data["pannanodx"] = {
                            "classification": best_class,
                            "confidence": max_score * 100,
                            "confidence_level": _get_confidence_level(max_score * 100),
                            "features": features
                        }
            except:
                pass
        
        # Extract Random Forest data
        rf_file = sample_dir / "random_forest_scores.csv"
        if rf_file.exists():
            try:
                with open(rf_file, 'r') as f:
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
                                except:
                                    pass
                        
                        # Some random forest outputs are already in percent (0-100),
                        # while others are fractional (0-1). Normalize to percent.
                        confidence_percent = max_score * 100 if max_score <= 1 else max_score
                        classification_data["random_forest"] = {
                            "classification": best_class,
                            "confidence": confidence_percent,
                            "confidence_level": _get_confidence_level(confidence_percent),
                            "features": features
                        }
            except:
                pass
                        
    except Exception:
        pass
    
    return classification_data


def _extract_fusion_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract fusion analysis data."""
    fusion_data = {
        "target_fusions": "Not available",
        "genome_fusions": "Not available"
    }
    
    try:
        # Look for fusion result files
        fusion_files = [
            "sv_count.txt",
            "fusion_results.csv",
            "fusion_summary.txt",
            "fusion_analysis.csv"
        ]
        
        for fusion_file in fusion_files:
            file_path = sample_dir / fusion_file
            if file_path.exists():
                if fusion_file.endswith('.txt'):
                    try:
                        with open(file_path, 'r') as f:
                            content = f.read().strip()
                            # Parse counter file for basic fusion info
                            if content.isdigit():
                                fusion_data["genome_fusions"] = int(content)  # Assume this is genome-wide fusions
                                fusion_data["target_fusions"] = 0  # Default value
                    except:
                        pass
                elif fusion_file.endswith('.csv'):
                    try:
                        with open(file_path, 'r') as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                if "target_fusions" in row:
                                    try:
                                        fusion_data["target_fusions"] = int(row["target_fusions"])
                                    except:
                                        pass
                                if "genome_fusions" in row:
                                    try:
                                        fusion_data["genome_fusions"] = int(row["genome_fusions"])
                                    except:
                                        pass
                    except:
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
        cov_val = float(coverage.replace('x', ''))
        if cov_val >= 30:
            return "positive"
        elif cov_val >= 20:
            return "primary"
        elif cov_val >= 10:
            return "warning"
        else:
            return "negative"
    except:
        return "grey"


def _get_methylation_color(status: str) -> str:
    """Get the appropriate color class for methylation status."""
    if "Methylated" in status:
        return "text-green-600"
    elif "Unmethylated" in status:
        return "text-orange-600"
    else:
        return "text-gray-600"


def _get_confidence_level(confidence: float) -> str:
    """Get the confidence level based on confidence percentage."""
    if confidence >= 80:
        return "High confidence"
    elif confidence >= 50:
        return "Medium confidence"
    elif confidence >= 20:
        return "Low confidence"
    else:
        return "Very low confidence"


def _get_methylation_badge_color(methylation: str) -> str:
    """Get the appropriate badge color for methylation value."""
    try:
        meth_val = float(methylation.replace('%', ''))
        if meth_val > 10:
            return "positive"
        elif meth_val > 5:
            return "warning"
        else:
            return "negative"
    except:
        return "grey"
