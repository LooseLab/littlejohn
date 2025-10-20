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

from robin.gui.config import get_confidence_level as get_classifier_confidence_level


@ui.refreshable
def _run_info_section(sample_dir: Path, sample_id: str):
    """Create refreshable run information section."""
    run_info = _extract_run_information(sample_dir, sample_id)
    
    with ui.card().classes("w-full"):
        ui.label("Run Details").classes("text-lg font-semibold text-blue-800")
        ui.separator().classes().style("border: 1px solid var(--md-primary)")
        
        with ui.row().classes("w-full gap-1 flex-wrap"):
            _create_dashboard_card(
                "Run Time", run_info.get("run_time", "Not available"), "schedule", "Sequencing run timestamp"
            )
            _create_dashboard_card(
                "Basecall Model", run_info.get("model", "Not available"), "settings", "AI model used for basecalling"
            )
            _create_dashboard_card(
                "Device", run_info.get("device", "Not available"), "smartphone", "Sequencing device identifier"
            )
            _create_dashboard_card(
                "Flow Cell", run_info.get("flow_cell", "Not available"), "tag", "Flow cell identification"
            )
            _create_dashboard_card(
                "Analysis Panel", run_info.get("panel", "Not available"), "science", "Target panel used for analysis"
            )
            _create_dashboard_card(
                "Sample ID", sample_id, "play_arrow", "Unique sample identifier"
            )


@ui.refreshable
def _classification_section(sample_dir: Path):
    """Create refreshable classification section."""
    classification_data = _extract_classification_data(sample_dir)
    
    with ui.card().classes("w-full"):
        ui.label("Classification Details").classes("text-lg font-semibold text-blue-800")
        ui.separator().classes().style("border: 1px solid var(--md-primary)")
        
        with ui.row().classes("w-full gap-1 flex-wrap"):
            # Sturgeon
            sturgeon_data = classification_data.get("sturgeon", {})
            _create_classification_dashboard_card_with_data(
                "Sturgeon", 
                sturgeon_data.get("classification", "Not available"),
                sturgeon_data.get("confidence", 0.0),
                sturgeon_data.get("confidence_level", "Not available"),
                sturgeon_data.get("features", 0),
                "psychology", 
                "AI brain tumor classification"
            )
            
            # NanoDX
            nanodx_data = classification_data.get("nanodx", {})
            _create_classification_dashboard_card_with_data(
                "NanoDX", 
                nanodx_data.get("classification", "Not available"),
                nanodx_data.get("confidence", 0.0),
                nanodx_data.get("confidence_level", "Not available"),
                nanodx_data.get("features", 0),
                "biotech", 
                "Molecular cancer diagnostics"
            )
            
            # PanNanoDX
            pannanodx_data = classification_data.get("pannanodx", {})
            _create_classification_dashboard_card_with_data(
                "PanNanoDX", 
                pannanodx_data.get("classification", "Not available"),
                pannanodx_data.get("confidence", 0.0),
                pannanodx_data.get("confidence_level", "Not available"),
                pannanodx_data.get("features", 0),
                "science", 
                "Pan-cancer analysis"
            )
            
            # Random Forest
            rf_data = classification_data.get("random_forest", {})
            _create_classification_dashboard_card_with_data(
                "Random Forest", 
                rf_data.get("classification", "Not available"),
                rf_data.get("confidence", 0.0),
                rf_data.get("confidence_level", "Not available"),
                rf_data.get("features", 0),
                "forest", 
                "Machine learning classification"
            )


@ui.refreshable
def _analysis_section(sample_dir: Path):
    """Create refreshable analysis section."""
    coverage_data = _extract_coverage_data(sample_dir)
    cnv_data = _extract_cnv_data(sample_dir)
    mgmt_data = _extract_mgmt_data(sample_dir)
    fusion_data = _extract_fusion_data(sample_dir)
    
    with ui.card().classes("w-full"):
        ui.label("Analysis Details").classes("text-lg font-semibold text-blue-800")
        ui.separator().classes().style("border: 1px solid var(--md-primary)")
        
        with ui.row().classes("w-full gap-1 flex-wrap"):
            # Coverage Analysis
            _create_coverage_dashboard_card_with_data(
                coverage_data.get("quality", "Not available"),
                coverage_data.get("target_coverage", "Not available"),
                coverage_data.get("global_coverage", "Not available"),
                coverage_data.get("enrichment", "Not available")
            )
            
            # CNV Analysis
            _create_cnv_dashboard_card_with_data(
                cnv_data.get("genetic_sex", "Not available"),
                cnv_data.get("bin_width", "Not available"),
                cnv_data.get("variance", "Not available"),
                cnv_data.get("gained", 0),
                cnv_data.get("lost", 0)
            )
            
            # MGMT Analysis
            _create_mgmt_dashboard_card_with_data(
                mgmt_data.get("status", "Not available"),
                mgmt_data.get("methylation_percent", "Not available"),
                mgmt_data.get("average_methylation", "Not available"),
                mgmt_data.get("prediction_score", "Not available"),
                mgmt_data.get("cpg_sites", "Not available")
            )
            
            # Fusion Analysis
            panel = _get_analysis_panel(sample_dir)
            _create_fusion_dashboard_card_with_data(
                panel,
                fusion_data.get("target_fusions", 0),
                fusion_data.get("genome_fusions", 0)
            )


def add_summary_section(sample_dir: Path, sample_id: str) -> None:
    """Build the Summary section at the top of the sample detail page using ui.refreshable."""
    
    # Create refreshable sections synchronously for immediate rendering
    # This ensures they appear at the top of the page
    _run_info_section(sample_dir, sample_id)
    _classification_section(sample_dir)
    _analysis_section(sample_dir)
    
    # Start periodic refresh timer (every 30 seconds)
    ui.timer(30.0, lambda: [
        _run_info_section.refresh(),
        _classification_section.refresh(),
        _analysis_section.refresh()
    ], active=True, immediate=True)


def _create_dashboard_card(title: str, value: str, icon: str, description: str) -> None:
    """Create a compact dashboard-style card matching the Mosaic design pattern."""
    with ui.card().classes("mosaic-card flex-1 min-w-0"):
        with ui.row().classes("mosaic-card__header"):
            # Title
            ui.label(title).classes("mosaic-card__title truncate")
            
            # Icon in circular background
            with ui.row().classes("mosaic-card__icon"):
                ui.icon(icon).classes("w-4 h-4")
        
        # Main value with text wrapping
        ui.label(value).classes("mosaic-card__content text-sm font-bold mb-1 break-words")
        
        # Description
        ui.label(description).classes("mosaic-card__subtitle")


def _create_classification_dashboard_card_with_data(
    title: str, 
    classification: str, 
    confidence: float, 
    confidence_level: str, 
    features: int, 
    icon: str, 
    description: str
) -> None:
    """Create a compact classification dashboard card with Mosaic styling."""
    with ui.card().classes("mosaic-card flex-1 min-w-0"):
        with ui.row().classes("mosaic-card__header"):
            # Title
            ui.label(title).classes("mosaic-card__title truncate")
            
            # Icon in circular background
            with ui.row().classes("mosaic-card__icon"):
                ui.icon(icon).classes("w-4 h-4")
        
        # Main classification result with confidence badge
        with ui.row().classes("flex items-start justify-between mb-1 gap-1"):
            ui.label(classification).classes("mosaic-card__content text-sm font-bold break-words flex-1")
            # Confidence level badge with color coding
            confidence_badge_class = "status-badge"
            if confidence >= 80:
                confidence_badge_class += " status-badge--success"
            elif confidence >= 50:
                confidence_badge_class += " status-badge--warning"
            else:
                confidence_badge_class += " status-badge--error"
            
            ui.label(confidence_level).classes(confidence_badge_class)
        
        # Compact details in a single row
        with ui.row().classes("flex items-center justify-between mb-1 gap-1"):
            ui.label(f"Confidence: {confidence:.1f}%").classes("mosaic-card__subtitle")
            ui.label(f"Features: {features:,}").classes("mosaic-card__subtitle")
        
        # Description
        ui.label(description).classes("mosaic-card__subtitle")


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


def _create_coverage_dashboard_card_with_data(
    quality: str, 
    target_coverage: str, 
    global_coverage: str, 
    enrichment: str
) -> None:
    """Create a compact coverage analysis dashboard card with data."""
    with ui.card().classes("mosaic-card flex-1 min-w-0"):
        with ui.row().classes("mosaic-card__header"):
            # Title
            ui.label("Coverage Analysis").classes("mosaic-card__title truncate")
            
            # Icon in circular background
            with ui.row().classes("mosaic-card__icon"):
                ui.icon("analytics").classes("w-5 h-5")
        
        # Main quality result with badge
        with ui.row().classes("flex items-start justify-between mb-1 gap-1"):
            ui.label(quality).classes("mosaic-card__content text-sm font-bold break-words flex-1")
            # Coverage badge with color coding
            coverage_badge_class = "status-badge"
            try:
                coverage_num = float(target_coverage.replace("x", ""))
                if coverage_num >= 30:
                    coverage_badge_class += " status-badge--success"
                elif coverage_num >= 20:
                    coverage_badge_class += " status-badge--info"
                elif coverage_num >= 10:
                    coverage_badge_class += " status-badge--warning"
                else:
                    coverage_badge_class += " status-badge--error"
            except (ValueError, AttributeError):
                coverage_badge_class += " status-badge--info"
            
            ui.label(target_coverage).classes(coverage_badge_class)
        
        # Coverage details in compact layout
        with ui.column().classes("mb-1 gap-0.5"):
            ui.label(f"Global: {global_coverage}").classes("mosaic-card__subtitle")
            ui.label(f"Targets: {target_coverage}").classes("mosaic-card__subtitle")
            ui.label(f"Enrichment: {enrichment}").classes("mosaic-card__subtitle")
        
        # Coverage thresholds as small badges
        with ui.row().classes("gap-1 flex-wrap"):
            ui.label("≥30x").classes("px-1 py-0.5 text-xs bg-green-100 text-green-800 rounded")
            ui.label("≥20x").classes("px-1 py-0.5 text-xs bg-blue-100 text-blue-800 rounded")
            ui.label("≥10x").classes("px-1 py-0.5 text-xs bg-yellow-100 text-yellow-800 rounded")
            ui.label("<10x").classes("px-1 py-0.5 text-xs bg-red-100 text-red-800 rounded")


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


def _create_cnv_dashboard_card_with_data(
    genetic_sex: str, 
    bin_width: str, 
    variance: str, 
    gained: int, 
    lost: int
) -> None:
    """Create a compact CNV analysis dashboard card with data."""
    with ui.card().classes("mosaic-card flex-1 min-w-0"):
        with ui.row().classes("mosaic-card__header"):
            # Title
            ui.label("Copy Number Analysis").classes("mosaic-card__title truncate")
            
            # Icon in circular background
            with ui.row().classes("mosaic-card__icon"):
                ui.icon("person").classes("w-5 h-5")
        
        # Main genetic sex result
        ui.label(genetic_sex).classes("mosaic-card__content text-sm font-bold mb-1 break-words")
        
        # Analysis details in compact layout
        with ui.column().classes("mb-1 gap-0.5"):
            ui.label(f"Bin Width: {bin_width}").classes("mosaic-card__subtitle")
            ui.label(f"Variance: {variance}").classes("mosaic-card__subtitle")
        
        # CNV counts as badges
        with ui.row().classes("gap-1 mb-1"):
            ui.label(f"Gained: {gained}").classes("px-1 py-0.5 text-xs font-medium rounded-full bg-green-100 text-green-800")
            ui.label(f"Lost: {lost}").classes("px-1 py-0.5 text-xs font-medium rounded-full bg-red-100 text-red-800")
        
        # Description
        ui.label("Copy number analysis across genome with breakpoint detection").classes("mosaic-card__subtitle")


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


def _create_mgmt_dashboard_card_with_data(
    status: str, 
    methylation_percent: str, 
    average_methylation: str, 
    prediction_score: str, 
    cpg_sites: str
) -> None:
    """Create a compact MGMT analysis dashboard card with data."""
    with ui.card().classes("mosaic-card flex-1 min-w-0"):
        with ui.row().classes("mosaic-card__header"):
            # Title
            ui.label("MGMT Analysis").classes("mosaic-card__title truncate")
            
            # Icon in circular background
            with ui.row().classes("mosaic-card__icon"):
                ui.icon("science").classes("w-5 h-5")
        
        # Main status result with badge
        with ui.row().classes("flex items-start justify-between mb-1 gap-1"):
            ui.label(status).classes("mosaic-card__content text-sm font-bold break-words flex-1")
            # Methylation badge with color coding
            methylation_badge_class = "status-badge"
            try:
                meth_num = float(methylation_percent.replace("%", ""))
                if meth_num > 10:
                    methylation_badge_class += " status-badge--success"
                elif meth_num > 5:
                    methylation_badge_class += " status-badge--warning"
                else:
                    methylation_badge_class += " status-badge--error"
            except (ValueError, AttributeError):
                methylation_badge_class += " status-badge--info"
            
            ui.label(methylation_percent).classes(methylation_badge_class)
        
        # Analysis details in compact layout
        with ui.column().classes("mb-1 gap-0.5"):
            ui.label(f"Average: {average_methylation}").classes("mosaic-card__subtitle")
            ui.label(f"Score: {prediction_score}").classes("mosaic-card__subtitle")
        
        # Description
        ui.label(f"MGMT status determined from methylation analysis of {cpg_sites} CpG sites").classes("mosaic-card__subtitle")


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


def _create_fusion_dashboard_card_with_data(
    panel: str, 
    target_fusions: int, 
    genome_fusions: int
) -> None:
    """Create a compact fusion analysis dashboard card with data."""
    with ui.card().classes("mosaic-card flex-1 min-w-0"):
        with ui.row().classes("mosaic-card__header"):
            # Title
            ui.label("Fusion Analysis").classes("mosaic-card__title truncate")
            
            # Icon in circular background
            with ui.row().classes("mosaic-card__icon"):
                ui.icon("merge").classes("w-5 h-5")
        
        # Panel info and main result
        with ui.row().classes("flex items-start justify-between mb-1 gap-1"):
            ui.label(f"Panel: {panel}").classes("mosaic-card__content text-sm font-medium break-words flex-1")
            ui.label(f"{target_fusions} target fusions").classes("px-1 py-0.5 text-xs font-medium rounded-full bg-blue-100 text-blue-800")
        
        # Analysis details in compact layout
        with ui.column().classes("mb-1"):
            ui.label(f"{genome_fusions} genome wide fusions").classes("mosaic-card__subtitle")
        
        # Description
        ui.label("Fusion candidates identified from reads with supplementary alignments").classes("mosaic-card__subtitle")


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
        # No fallback to rCNS2 - return empty string if not found
        return ""
    except Exception:
        # No fallback to rCNS2 - return empty string if error
        return ""


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
        # First check for final_mgmt.csv files (highest priority)
        final_files = list(sample_dir.glob("final_mgmt.csv"))
        if final_files:
            latest_csv = final_files[0]
        else:
            # Fallback to numeric-prefixed files
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
    """Extract fusion analysis data from generated summary files."""
    fusion_data = {"target_fusions": 0, "genome_fusions": 0}
    
    logging.info(f"[Summary] _extract_fusion_data() called with sample_dir: {sample_dir}")

    try:
        # Debug: List all fusion-related files
        fusion_files = list(sample_dir.glob("*fusion*"))
        logging.info(f"[Summary] Found fusion files: {[f.name for f in fusion_files]}")
        
        # Debug: Check if genome-wide processed file exists
        genome_file = sample_dir / "fusion_candidates_all_processed.csv"
        logging.info(f"[Summary] Genome-wide processed file exists: {genome_file.exists()}")
        if genome_file.exists():
            logging.info(f"[Summary] Genome-wide processed file size: {genome_file.stat().st_size} bytes")
        # First try to read from the new fusion_summary.csv file
        summary_file = sample_dir / "fusion_summary.csv"
        if summary_file.exists():
            try:
                with open(summary_file, "r") as f:
                    content = f.read()
                    
                with open(summary_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        fusion_data["target_fusions"] = int(row.get("target_fusions", 0))
                        fusion_data["genome_fusions"] = int(row.get("genome_fusions", 0))
                        break
                
                # Check if the summary file has incorrect genome-wide count (0)
                if fusion_data["genome_fusions"] == 0:
                    # Try to regenerate the summary file with correct data
                    try:
                        from robin.gui.components.fusion import _generate_summary_files_from_pickle
                        if _generate_summary_files_from_pickle(sample_dir, force_regenerate=True):
                            # Re-read the regenerated summary file
                            with open(summary_file, "r") as f:
                                reader = csv.DictReader(f)
                                for row in reader:
                                    fusion_data["target_fusions"] = int(row.get("target_fusions", 0))
                                    fusion_data["genome_fusions"] = int(row.get("genome_fusions", 0))
                                    break
                            return fusion_data
                    except Exception as e:
                        pass  # Don't return early, continue to pickle file loading
                else:
                    logging.info(f"[Summary] Fusion data extracted from summary file - target: {fusion_data['target_fusions']}, genome: {fusion_data['genome_fusions']}")
                    logging.info(f"[Summary] Summary file path: {summary_file}")
                    return fusion_data
            except Exception as e:
                logging.debug(f"   Fusion: Failed to read summary file: {e}")
        
        # If summary files don't exist, try to generate them from pickle files
        try:
            from robin.gui.components.fusion import _generate_summary_files_from_pickle
            if _generate_summary_files_from_pickle(sample_dir, force_regenerate=False):
                # Try reading the newly generated summary file
                if summary_file.exists():
                    with open(summary_file, "r") as f:
                        reader = csv.DictReader(f)
                        for row in reader:
                            fusion_data["target_fusions"] = int(row.get("target_fusions", 0))
                            fusion_data["genome_fusions"] = int(row.get("genome_fusions", 0))
                            break
                    logging.info(f"[Summary] Fusion data extracted from generated summary file - target: {fusion_data['target_fusions']}, genome: {fusion_data['genome_fusions']}")
                    return fusion_data
        except Exception as e:
            logging.debug(f"   Fusion: Failed to generate summary files from pickle: {e}")
        
        # If still no data, try to load directly from pickle files and count gene_pairs
        try:
            from robin.gui.components.fusion import _load_processed_pickle
            target_file = sample_dir / "fusion_candidates_master_processed.csv"
            genome_file = sample_dir / "fusion_candidates_all_processed.csv"
            
            
            target_data = _load_processed_pickle(target_file)
            genome_data = _load_processed_pickle(genome_file)
            
            
            if target_data and isinstance(target_data, dict):
                fusion_data["target_fusions"] = target_data.get("candidate_count", 0)
            
            if genome_data and isinstance(genome_data, dict):
                # Use gene_pairs count if candidate_count is 0 (same logic as GUI)
                genome_count = genome_data.get("candidate_count", 0)
                if genome_count == 0 and genome_data.get("gene_pairs"):
                    genome_count = len(genome_data.get("gene_pairs", []))
                fusion_data["genome_fusions"] = genome_count
            
            logging.info(f"[Summary] Fusion data loaded directly from pickle files - target: {fusion_data['target_fusions']}, genome: {fusion_data['genome_fusions']}")
            return fusion_data
        except Exception as e:
            logging.debug(f"   Fusion: Failed to load directly from pickle files: {e}")
        
        # Fallback to individual fusion_results.csv file
        fusion_results_file = sample_dir / "fusion_results.csv"
        if fusion_results_file.exists():
            try:
                with open(fusion_results_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        fusion_data["target_fusions"] = int(row.get("target_fusions", 0))
                        fusion_data["genome_fusions"] = int(row.get("genome_fusions", 0))
                        break
                logging.debug(f"[Summary] Fusion data extracted from results file - target: {fusion_data['target_fusions']}, genome: {fusion_data['genome_fusions']}")
                return fusion_data
            except Exception as e:
                logging.debug(f"   Fusion: Failed to read results file: {e}")
        
        # Final fallback to legacy sv_count.txt file
        sv_count_file = sample_dir / "sv_count.txt"
        if sv_count_file.exists():
            try:
                with open(sv_count_file, "r") as f:
                    content = f.read().strip()
                    if content.isdigit():
                        # Legacy behavior - assume this is genome-wide count
                        fusion_data["genome_fusions"] = int(content)
                        fusion_data["target_fusions"] = 0
                logging.debug(f"[Summary] Fusion data extracted from legacy sv_count file - target: {fusion_data['target_fusions']}, genome: {fusion_data['genome_fusions']}")
            except Exception as e:
                logging.debug(f"   Fusion: Failed to read legacy sv_count file: {e}")

    except Exception as e:
        logging.debug(f"   Fusion: <access denied>: {e}")

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


