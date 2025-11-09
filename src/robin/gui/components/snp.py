from __future__ import annotations

from typing import Any, Dict, List
from pathlib import Path
import logging
import json

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

logger = logging.getLogger(__name__)


def navigate_igv_to_snp(chrom: str, pos: int, flank: int = 100) -> None:
    """
    Navigate IGV browser to a specific SNP location.
    
    Args:
        chrom: Chromosome name (e.g., "chr1")
        pos: Position on the chromosome
        flank: Number of bases to include on each side (default: 100 bp)
    """
    try:
        # Ensure chromosome name has 'chr' prefix if needed
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        
        # Calculate window around the SNP
        start = max(1, pos - flank)
        end = pos + flank
        
        region = f"{chrom}:{start}-{end}"
        
        # Escape region string for JavaScript
        escaped_region = region.replace('"', '\\"').replace("'", "\\'")
        
        js_navigate = f"""
            (function() {{
                try {{
                    if (window.lj_igv && window.lj_igv_browser_ready) {{
                        console.log('[IGV] Navigating to SNP region: {escaped_region}');
                        window.lj_igv.search('{escaped_region}');
                    }} else {{
                        console.warn('[IGV] Browser not ready yet, will navigate when ready');
                        setTimeout(function() {{
                            if (window.lj_igv && window.lj_igv_browser_ready) {{
                                window.lj_igv.search('{escaped_region}');
                            }}
                        }}, 1000);
                    }}
                }} catch (error) {{
                    console.error('[IGV] Error navigating to SNP: ' + error);
                }}
            }})();
        """
        
        ui.run_javascript(js_navigate, timeout=5.0)
        
    except Exception as e:
        logger.error(f"Error navigating IGV to SNP {chrom}:{pos}: {e}")


def add_snp_section(launcher: Any, sample_dir: Path) -> None:
    """
    Add SNP analysis section to the sample details page.
    
    Args:
        launcher: The GUI launcher instance
        sample_dir: Path to the sample directory
    """
    if not sample_dir or not sample_dir.exists():
        return
    
    # Look for VCF files in clair3 directory
    clair3_dir = sample_dir / "clair3"
    if not clair3_dir.exists():
        return
    
    # Look for snpsift_output.vcf (preferred) or other VCF files
    display_file = clair3_dir / "snpsift_output_display.json"

    if not display_file.exists():
        with ui.card().classes("w-full"):
            ui.label("SNP Analysis").classes("text-lg font-semibold text-blue-800")
            ui.separator().classes().style("border: 1px solid var(--md-primary)")
            ui.label(
                "Precomputed SNP display data was not found. Ensure SNP analysis has completed."
            ).classes("text-body-medium text-gray-600")
        return

    try:
        with display_file.open("r", encoding="utf-8") as f_in:
            snp_display = json.load(f_in)
    except Exception as exc:
        logger.error(f"Failed to load SNP display data: {exc}")
        with ui.card().classes("w-full"):
            ui.label("SNP Analysis").classes("text-lg font-semibold text-blue-800")
            ui.separator().classes().style("border: 1px solid var(--md-primary)")
            ui.label(
                "Could not load SNP variant data. Check logs for details."
            ).classes("text-body-medium text-red-600")
        return

    columns: List[Dict[str, Any]] = snp_display.get("columns", [])
    rows_all: List[Dict[str, Any]] = snp_display.get("rows_all", [])
    rows_pathogenic: List[Dict[str, Any]] = snp_display.get("rows_pathogenic", [])
    summary: Dict[str, Any] = snp_display.get("summary", {})
    snp_regions_map: Dict[str, Dict[str, int]] = snp_display.get("snp_regions_map", {})

    total_variants = summary.get("total_variants", len(rows_all))
    pathogenic_count = summary.get("pathogenic_variants", len(rows_pathogenic))

    js_snp_regions_json = json.dumps(snp_regions_map)

    def navigate_to_snp_region(snp_key: str) -> None:
        if snp_key in snp_regions_map:
            snp_data = snp_regions_map[snp_key]
            navigate_igv_to_snp(snp_data["chrom"], snp_data["pos"])

    js_init_snp_regions = f"""
        (function() {{
            window.snpRegionsMap = window.snpRegionsMap || {{}};
            Object.assign(window.snpRegionsMap, {js_snp_regions_json});
            console.log('[SNP] Initialized SNP regions map with', Object.keys(window.snpRegionsMap).length, 'pathogenic SNPs');
        }})();
    """
    ui.run_javascript(js_init_snp_regions, timeout=5.0)

    with ui.card().classes("w-full"):
        ui.label("SNP Analysis").classes("text-lg font-semibold text-blue-800")
        ui.separator().classes().style("border: 1px solid var(--md-primary)")

        with ui.row().classes("w-full gap-4 mb-4"):
            ui.label(f"Total Variants: {total_variants}").classes(
                "text-body-medium font-semibold"
            )
            if pathogenic_count > 0:
                ui.label(f"Pathogenic Variants: {pathogenic_count}").classes(
                    "text-body-medium font-semibold text-red-600"
                )

        from robin.gui.theme import styled_table

        table_container, snp_table = styled_table(
            columns=columns,
            rows=rows_all,
            pagination=25,
            class_size="table-xs",
        )

        with snp_table.add_slot("top-left"):
            pathogenic_filter = ui.switch(
                "Show pathogenic variants only",
                value=False,
            )

        with snp_table.add_slot("top-right"):
            with ui.input(placeholder="Search SNPs...").props("type=search").bind_value(
                snp_table, "filter"
            ).add_slot("append"):
                ui.icon("search")

        if any(col.get("field") == "action" for col in columns):
            try:
                snp_table.add_slot(
                    "body-cell-action",
                    """
<q-td key="action" :props="props">
  <q-btn 
    v-if="props.row.is_pathogenic === 'Yes'"
    icon="visibility" 
    size="sm" 
    dense 
    flat 
    color="primary"
    @click="$parent.$emit('snp-view-igv', props.row.CHROM + ':' + props.row.POS)"
    title="View in IGV"
  />
</q-td>
""",
                )

                def on_snp_view_igv(e):
                    try:
                        snp_key = (
                            e.args if isinstance(e.args, str) else getattr(e, "args", None)
                        )
                        if snp_key:
                            logging.debug(f"[SNP] Button clicked for: {snp_key}")
                            navigate_to_snp_region(snp_key)
                    except Exception as ex:
                        logger.debug(f"Error handling SNP IGV view: {ex}")

                snp_table.on("snp-view-igv", on_snp_view_igv)
            except Exception as ex:
                logger.warning(f"Could not add action button slot: {ex}")

        def apply_pathogenic_filter(value: bool) -> None:
            snp_table.rows = rows_pathogenic if value else rows_all
            snp_table.update()

        pathogenic_filter.on(
            "update:model-value", lambda e: apply_pathogenic_filter(bool(e.args))
        )

        for col in snp_table.columns:
            col["sortable"] = True
