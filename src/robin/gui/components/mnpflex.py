from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional
import base64
import json
import logging
import os
import threading
import time

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

from robin.utils.mnpflex_client_standalone import MNPFlexClient
from robin.analysis.utilities.matkit import reconstruct_full_bedmethyl_for_mnpflex
from robin.analysis.utilities.mnp_flex import APIClient as MnpFlexApiClient
from robin import resources
from robin.gui.theme import styled_table


def add_mnpflex_section(launcher: Any, sample_dir: Path, sample_id: str) -> None:
    """Display and refresh MNP-Flex results for the sample."""
    if ui is None:
        return
    username = os.getenv("MNPFLEX_USERNAME") or os.getenv("EPIGNOSTIX_USERNAME")
    password = os.getenv("MNPFLEX_PASSWORD") or os.getenv("EPIGNOSTIX_PASSWORD")
    if not username or not password:
        return

    state: Dict[str, Any] = {
        "running": False,
        "last_error": "",
        "last_updated": None,
        "auto_fetch_attempted": False,
    }

    def _sample_is_complete() -> bool:
        master_csv = sample_dir / "master.csv"
        if not master_csv.exists():
            return False
        try:
            import csv

            with master_csv.open("r", newline="") as fh:
                reader = csv.DictReader(fh)
                first_row = next(reader, None)
            if not first_row:
                return False
            active_jobs = int(first_row.get("samples_overview_active_jobs", 0) or 0)
            pending_jobs = int(first_row.get("samples_overview_pending_jobs", 0) or 0)
            total_jobs = int(first_row.get("samples_overview_total_jobs", 0) or 0)
            completed_jobs = int(
                first_row.get("samples_overview_completed_jobs", 0) or 0
            )
            return total_jobs > 0 and completed_jobs >= total_jobs and active_jobs == 0 and pending_jobs == 0
        except Exception:
            return False

    def _find_results_dir() -> Optional[Path]:
        candidates = [
            sample_dir / f"mnpflex_results_{sample_id}",
            sample_dir / "mnpflex_results",
            sample_dir,
        ]
        for candidate in candidates:
            if (candidate / "bundle_summary.json").exists():
                return candidate
        for candidate in sorted(sample_dir.glob("mnpflex_results_*")):
            if (candidate / "bundle_summary.json").exists():
                return candidate
        return None

    def _find_parquet_file() -> Optional[Path]:
        preferred = sample_dir / f"{sample_id}.parquet"
        if preferred.exists():
            return preferred
        matches = list(sample_dir.glob("*.parquet"))
        return matches[0] if matches else None

    def _build_bed_file_from_parquet() -> Path:
        parquet_path = _find_parquet_file()
        if not parquet_path or not parquet_path.exists():
            raise RuntimeError("No parquet data found for this sample.")
        bed_path = sample_dir / f"{sample_id}.mnpflex.bed"
        bed_df = reconstruct_full_bedmethyl_for_mnpflex(str(parquet_path))
        bed_df.to_csv(bed_path, sep="\t", index=False, header=False)
        return bed_path

    def _build_subset_bed_from_parquet() -> Path:
        bed_path = _build_bed_file_from_parquet()
        subset_path = sample_dir / f"{sample_id}.MNPFlex.subset.bed"
        reference_bed = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "mnp_flex_sample_clean.bed",
        )
        api_client = MnpFlexApiClient(base_url="https://mnp-flex.org", verify_ssl=True)
        api_client.process_streaming(reference_bed, str(bed_path), str(subset_path))
        return subset_path

    def _read_image_base64(path: Path) -> Optional[str]:
        if not path.exists():
            return None
        try:
            encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            return f"data:image/png;base64,{encoded}"
        except Exception:
            logging.exception(f"[MNPFlex] Failed to load image: {path}")
            return None

    def _fetch_results_from_api(output_dir: Path) -> None:
        username = os.getenv("MNPFLEX_USERNAME") or os.getenv("EPIGNOSTIX_USERNAME")
        password = os.getenv("MNPFLEX_PASSWORD") or os.getenv("EPIGNOSTIX_PASSWORD")
        if not username or not password:
            raise RuntimeError(
                "Missing MNP-Flex credentials. Set MNPFLEX_USERNAME/MNPFLEX_PASSWORD."
            )
        workflow_id_env = os.getenv("MNPFLEX_WORKFLOW_ID", "18")
        try:
            workflow_id = int(workflow_id_env)
        except ValueError:
            raise RuntimeError(
                f"Invalid MNPFLEX_WORKFLOW_ID: {workflow_id_env}. Must be an integer."
            )

        # Use the same subset format as mnpflex_preprocess_modkit.sh
        bed_path = _build_subset_bed_from_parquet()

        base_url = os.getenv("MNPFLEX_BASE_URL", "https://app.epignostix.com")
        verify_ssl_env = os.getenv("MNPFLEX_VERIFY_SSL", "true").lower()
        verify_ssl = verify_ssl_env not in {"0", "false", "no"}

        client = MNPFlexClient(
            base_url=base_url,
            username=username,
            password=password,
            verify_ssl=verify_ssl,
        )
        client.upload_retrieve_cleanup(
            bed_file_path=str(bed_path),
            sample_identifier=sample_id,
            workflow_id=workflow_id,
            output_dir=str(output_dir),
        )

    with ui.card().classes("w-full"):
        with ui.row().classes("w-full items-center justify-between"):
            ui.label("MNP-Flex results").classes("text-sm font-semibold text-gray-800")
            status_label = ui.label("Idle").classes("text-xs text-gray-600")

        with ui.row().classes("w-full items-center gap-4"):
            last_updated_label = ui.label("Last updated: --").classes(
                "text-xs text-gray-500"
            )
            error_label = ui.label("").classes("text-xs text-red-600")

        with ui.row().classes("w-full items-center gap-2"):
            fetch_button = ui.button("Run MNP-Flex analysis").classes(
                "text-sm font-semibold px-3 py-1 rounded"
            )
            build_button = ui.button("Generate MNP-Flex subset BED").classes(
                "text-sm font-semibold px-3 py-1 rounded bg-gray-100"
            )

        with ui.row().classes("w-full gap-6"):
            with ui.column().classes("w-full md:w-1/2 gap-2"):
                qc_status = ui.label("QC status: --").classes("text-xs text-gray-700")
                qc_coverage = ui.label("Average coverage: --").classes(
                    "text-xs text-gray-700"
                )
                qc_missing = ui.label("Missing sites: --").classes(
                    "text-xs text-gray-700"
                )
            with ui.column().classes("w-full md:w-1/2 gap-2"):
                mgmt_status = ui.label("MGMT status: --").classes(
                    "text-xs text-gray-700"
                )
                mgmt_average = ui.label("MGMT average: --").classes(
                    "text-xs text-gray-700"
                )
                mgmt_sites = ui.label("MGMT sites: --").classes("text-xs text-gray-700")

        ui.separator()
        ui.label("Hierarchical summary").classes("text-sm font-semibold text-gray-800")
        with ui.row().classes("w-full items-center gap-2"):
            ui.label("Top path:").classes("text-xs text-gray-700")
            top_path_value = ui.label("--").classes("text-xs text-gray-700")
            top_path_badge = ui.badge("--").classes(
                "text-xs bg-gray-100 text-gray-700"
            )
        ui.label("Hierarchy aggregates (top entry)").classes(
            "text-xs text-gray-700 mt-1"
        )
        with ui.row().classes("w-full gap-6"):
            with ui.column().classes("w-full md:w-1/2 gap-2"):
                with ui.row().classes("w-full items-center gap-2"):
                    ui.label("Subclass:").classes("text-xs text-gray-700")
                    agg_subclass_name = ui.label("--").classes("text-xs text-gray-700")
                    agg_subclass_badge = ui.badge("--").classes(
                        "text-xs bg-gray-100 text-gray-700"
                    )
                with ui.row().classes("w-full items-center gap-2"):
                    ui.label("Class:").classes("text-xs text-gray-700")
                    agg_class_name = ui.label("--").classes("text-xs text-gray-700")
                    agg_class_badge = ui.badge("--").classes(
                        "text-xs bg-gray-100 text-gray-700"
                    )
            with ui.column().classes("w-full md:w-1/2 gap-2"):
                with ui.row().classes("w-full items-center gap-2"):
                    ui.label("Family:").classes("text-xs text-gray-700")
                    agg_family_name = ui.label("--").classes("text-xs text-gray-700")
                    agg_family_badge = ui.badge("--").classes(
                        "text-xs bg-gray-100 text-gray-700"
                    )
                with ui.row().classes("w-full items-center gap-2"):
                    ui.label("Superfamily:").classes("text-xs text-gray-700")
                    agg_superfamily_name = ui.label("--").classes(
                        "text-xs text-gray-700"
                    )
                    agg_superfamily_badge = ui.badge("--").classes(
                        "text-xs bg-gray-100 text-gray-700"
                    )
        _, hierarchy_table = styled_table(
            columns=[
                {"name": "group", "label": "Group", "field": "group"},
                {"name": "score", "label": "Score", "field": "score"},
                {"name": "description", "label": "Description", "field": "description"},
            ],
            rows=[],
            pagination=0,
            class_size="table-xs",
        )

        ui.separator()
        ui.label("Classifier summary").classes("text-sm font-semibold text-gray-800")
        with ui.row().classes("w-full gap-6"):
            with ui.column().classes("w-full md:w-1/2 gap-2"):
                classifier_name = ui.label("Classifier: --").classes(
                    "text-xs text-gray-700"
                )
                classifier_version = ui.label("Version: --").classes(
                    "text-xs text-gray-700"
                )
                classifier_type = ui.label("Type: --").classes(
                    "text-xs text-gray-700"
                )

        ui.label("Top 10 classifier scores").classes("text-xs text-gray-700 mt-1")
        _, classifier_scores_table = styled_table(
            columns=[
                {"name": "score", "label": "Score", "field": "score"},
                {"name": "subclass", "label": "Subclass", "field": "subclass"},
                {"name": "class", "label": "Class", "field": "class"},
                {"name": "family", "label": "Family", "field": "family"},
                {"name": "superfamily", "label": "Superfamily", "field": "superfamily"},
            ],
            rows=[],
            pagination=0,
            class_size="table-xs",
        )

        images_container = ui.row().classes("w-full")

        def _format_score(value: Optional[float]) -> str:
            try:
                return f"{float(value):.4f}"
            except Exception:
                return "--"

        def _score_badge_classes(value: Optional[float]) -> str:
            try:
                score = float(value)
            except Exception:
                return "text-xs bg-gray-100 text-gray-700"
            if score >= 0.8:
                return "text-xs bg-green-100 text-green-800"
            if score >= 0.6:
                return "text-xs bg-teal-100 text-teal-800"
            if score >= 0.4:
                return "text-xs bg-yellow-100 text-yellow-800"
            if score >= 0.2:
                return "text-xs bg-orange-100 text-orange-800"
            return "text-xs bg-red-100 text-red-800"

        def _set_badge_value(badge, value: Optional[float]) -> None:
            badge.set_text(_format_score(value))
            classes = _score_badge_classes(value)
            try:
                badge.classes(replace=classes)
            except Exception:
                badge.classes(classes)

        def _extract_classifier_rows(classifier_summary: Dict[str, Any], limit: int = 10):
            scores = classifier_summary.get("scores") or []
            rows = []
            for item in scores:
                group = item.get("reference_group") or {}
                rows.append(
                    {
                        "score": _format_score(item.get("score")),
                        "subclass": group.get("molecular_subclass")
                        or group.get("name")
                        or "",
                        "class": group.get("molecular_class") or "",
                        "family": group.get("molecular_family") or "",
                        "superfamily": group.get("molecular_superfamily") or "",
                    }
                )
            rows.sort(
                key=lambda r: float(r["score"]) if r["score"] != "--" else -1,
                reverse=True,
            )
            return rows[:limit]

        def _sum_scores_by_field(scores, field_name, match_value):
            if not match_value:
                return None
            total = 0.0
            for item in scores:
                ref = item.get("reference_group") or {}
                if ref.get(field_name) == match_value:
                    total += float(item.get("score", 0) or 0)
            return total

        def _flatten_hierarchy(nodes, path=None):
            if path is None:
                path = []
            flat = []
            for node in nodes or []:
                group = node.get("group", "Unknown")
                score = node.get("score")
                current = path + [group]
                members = node.get("members") or []
                if members:
                    flat.extend(_flatten_hierarchy(members, current))
                else:
                    flat.append((score, current))
            return flat

        def _collect_descriptions(node) -> str:
            descriptions = []
            desc = (node.get("description") or "").strip()
            if desc:
                descriptions.append(desc)
            for member in node.get("members") or []:
                member_desc = _collect_descriptions(member)
                if member_desc:
                    descriptions.append(member_desc)
            # Deduplicate while preserving order
            seen = set()
            combined = []
            for text in descriptions:
                if text not in seen:
                    seen.add(text)
                    combined.append(text)
            return " ".join(combined).strip()

        def _extract_hierarchy_rows(nodes, limit: int = 10):
            rows = []
            for node in nodes or []:
                rows.append(
                    {
                        "group": node.get("group", "Unknown"),
                        "score": _format_score(node.get("score")),
                        "description": _collect_descriptions(node),
                    }
                )
            rows.sort(key=lambda r: float(r["score"]) if r["score"] != "--" else -1, reverse=True)
            return rows[:limit]

        def _update_labels(summary: Optional[Dict[str, Any]], summary_path: Optional[Path]) -> None:
            if summary is None:
                qc_status.set_text("QC status: --")
                qc_coverage.set_text("Average coverage: --")
                qc_missing.set_text("Missing sites: --")
                mgmt_status.set_text("MGMT status: --")
                mgmt_average.set_text("MGMT average: --")
                mgmt_sites.set_text("MGMT sites: --")
                classifier_name.set_text("Classifier: --")
                classifier_version.set_text("Version: --")
                classifier_type.set_text("Type: --")

                classifier_scores_table.rows = []
                classifier_scores_table.update()
                hierarchy_table.rows = []
                hierarchy_table.update()
                top_path_value.set_text("--")
                _set_badge_value(top_path_badge, None)
                agg_subclass_name.set_text("--")
                agg_class_name.set_text("--")
                agg_family_name.set_text("--")
                agg_superfamily_name.set_text("--")
                _set_badge_value(agg_subclass_badge, None)
                _set_badge_value(agg_class_badge, None)
                _set_badge_value(agg_family_badge, None)
                _set_badge_value(agg_superfamily_badge, None)
                last_updated_label.set_text("Last updated: --")
                return

            qc = summary.get("qc", {}) or {}
            mgmt = summary.get("mgmt", {}) or {}
            classifier_summary = summary.get("classifier_summary", {}) or {}
            classifier = classifier_summary.get("classifier", {}) or {}
            qc_status.set_text(f"QC status: {qc.get('status', 'Unknown')}")
            qc_coverage.set_text(
                f"Average coverage: {qc.get('avg_coverage', 'Unknown')}"
            )
            qc_missing.set_text(
                f"Missing sites: {qc.get('missing_site_count', 'Unknown')}"
            )
            mgmt_status.set_text(f"MGMT status: {mgmt.get('status', 'Unknown')}")
            mgmt_average.set_text(
                f"MGMT average: {mgmt.get('average', 'Unknown')}"
            )
            mgmt_sites.set_text(f"MGMT sites: {mgmt.get('site_count', 'Unknown')}")
            classifier_name.set_text(
                f"Classifier: {classifier.get('name', 'Unknown')}"
            )
            classifier_version.set_text(
                f"Version: {classifier.get('version', 'Unknown')}"
            )
            classifier_type.set_text(
                f"Type: {classifier.get('classifier_type', 'Unknown')}"
            )
            scores = classifier_summary.get("scores") or []
            classifier_scores_table.rows = _extract_classifier_rows(classifier_summary)
            classifier_scores_table.update()
            if scores:
                top = sorted(
                    scores,
                    key=lambda item: float(item.get("score", 0) or 0),
                    reverse=True,
                )[:10]
                top_ref = top[0].get("reference_group") or {}
                top_subclass = top_ref.get("molecular_subclass") or top_ref.get("name")
                top_class = top_ref.get("molecular_class")
                top_family = top_ref.get("molecular_family")
                top_superfamily = top_ref.get("molecular_superfamily")

                subclass_sum = _sum_scores_by_field(
                    scores, "molecular_subclass", top_subclass
                )
                class_sum = _sum_scores_by_field(scores, "molecular_class", top_class)
                family_sum = _sum_scores_by_field(scores, "molecular_family", top_family)
                superfamily_sum = _sum_scores_by_field(
                    scores, "molecular_superfamily", top_superfamily
                )

                agg_subclass_name.set_text(top_subclass or "N/A")
                agg_class_name.set_text(top_class or "N/A")
                agg_family_name.set_text(top_family or "N/A")
                agg_superfamily_name.set_text(top_superfamily or "N/A")
                _set_badge_value(agg_subclass_badge, subclass_sum)
                _set_badge_value(agg_class_badge, class_sum)
                _set_badge_value(agg_family_badge, family_sum)
                _set_badge_value(agg_superfamily_badge, superfamily_sum)
            else:
                agg_subclass_name.set_text("--")
                agg_class_name.set_text("--")
                agg_family_name.set_text("--")
                agg_superfamily_name.set_text("--")
                _set_badge_value(agg_subclass_badge, None)
                _set_badge_value(agg_class_badge, None)
                _set_badge_value(agg_family_badge, None)
                _set_badge_value(agg_superfamily_badge, None)

            hierarchy = classifier_summary.get("summary_hierarchical", []) or []
            hierarchy_table.rows = _extract_hierarchy_rows(hierarchy)
            hierarchy_table.update()
            flat = _flatten_hierarchy(hierarchy)
            if flat:
                best_score, best_path = max(flat, key=lambda x: x[0] or 0)
                top_path_value.set_text(" > ".join(best_path))
                _set_badge_value(top_path_badge, best_score)
            else:
                top_path_value.set_text("--")
                _set_badge_value(top_path_badge, None)

            if summary_path and summary_path.exists():
                updated = time.strftime(
                    "%Y-%m-%d %H:%M:%S", time.localtime(summary_path.stat().st_mtime)
                )
                last_updated_label.set_text(f"Last updated: {updated}")

        def _refresh_from_disk() -> None:
            results_dir = _find_results_dir()
            summary_path = None
            summary = None
            if results_dir:
                summary_path = results_dir / "bundle_summary.json"
                if summary_path.exists():
                    try:
                        summary = json.loads(summary_path.read_text())
                    except Exception:
                        logging.exception(
                            f"[MNPFlex] Failed to parse summary: {summary_path}"
                        )
                        summary = None

            _update_labels(summary, summary_path)

            images_container.clear()
            if results_dir:
                image_specs = [
                    ("QC coverage plot", results_dir / "qc_coverage_plot.png"),
                    (
                        "QC methylation density plot",
                        results_dir / "qc_methylation_density_plot.png",
                    ),
                    ("MGMT region plot", results_dir / "mgmt_region_plot.png"),
                ]
                with images_container:
                    for title, img_path in image_specs:
                        data_url = _read_image_base64(img_path)
                        if data_url:
                            with ui.column().classes("w-1/4 gap-3"):
                                ui.label(title).classes("text-xs text-gray-700")
                                ui.image(data_url).classes("w-full")

        def _run_fetch(auto: bool = False) -> None:
            if state["running"]:
                return
            state["running"] = True
            state["last_error"] = ""
            status_label.set_text("Running MNP-Flex analysis...")
            error_label.set_text("")
            fetch_button.disable()
            build_button.disable()

            def _worker() -> None:
                try:
                    output_dir = sample_dir / f"mnpflex_results_{sample_id}"
                    output_dir.mkdir(parents=True, exist_ok=True)
                    _fetch_results_from_api(output_dir)
                    state["last_updated"] = time.time()
                except Exception as exc:
                    state["last_error"] = str(exc)
                    logging.error(f"[MNPFlex] Fetch failed: {exc}")
                finally:
                    state["running"] = False

            thread = threading.Thread(target=_worker, daemon=True)
            thread.start()

        def _handle_fetch_click() -> None:
            _run_fetch(auto=False)

        def _handle_build_click() -> None:
            if state["running"]:
                return
            state["last_error"] = ""
            status_label.set_text("Generating MNP-Flex subset BED...")
            error_label.set_text("")
            fetch_button.disable()
            build_button.disable()

            def _worker() -> None:
                try:
                    bed_path = _build_subset_bed_from_parquet()
                    state["last_updated"] = time.time()
                    logging.info(f"[MNPFlex] BED file generated at {bed_path}")
                except Exception as exc:
                    state["last_error"] = str(exc)
                    logging.error(f"[MNPFlex] BED generation failed: {exc}")
                finally:
                    status_label.set_text("Idle")
                    fetch_button.enable()
                    build_button.enable()

            thread = threading.Thread(target=_worker, daemon=True)
            thread.start()

        def _poll_status() -> None:
            if state["running"]:
                status_label.set_text("Running MNP-Flex analysis...")
            else:
                status_label.set_text("Idle")
                fetch_button.enable()
                build_button.enable()

            if state["last_error"]:
                error_label.set_text(state["last_error"])
            else:
                error_label.set_text("")

            if not state["auto_fetch_attempted"]:
                results_dir = _find_results_dir()
                if results_dir is None and _sample_is_complete():
                    state["auto_fetch_attempted"] = True
                    _run_fetch(auto=True)
            _refresh_from_disk()

        fetch_button.on_click(_handle_fetch_click)
        build_button.on_click(_handle_build_click)
        ui.timer(30.0, _poll_status, active=True, immediate=True)
