from __future__ import annotations

from typing import Any, Dict, List, Optional
from pathlib import Path
import logging
import sys
import os

import pandas as pd
import numpy as np

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

logger = logging.getLogger(__name__)


def is_pathogenic(info_str: str) -> bool:
    """
    Check if a variant is pathogenic based on CLNSIG field.

    Args:
        info_str (str): The INFO field string from the VCF file

    Returns:
        bool: True if the variant is pathogenic, False otherwise
    """
    if "CLNSIG=" not in info_str:
        return False

    # Extract CLNSIG value
    for field in info_str.split(";"):
        if field.startswith("CLNSIG="):
            clnsig_value = field.split("=")[1]
            # Check for PATHOGENIC case-insensitive
            return "PATHOGENIC" in clnsig_value.upper()
    return False


def process_annotations(record: dict) -> tuple[dict, dict]:
    """
    This function takes a dictionary record from a vcf file and explodes the record into multiple records based on the contents of the INFO field.
    We expect some unit of 16 entries in the INFO field. Where there are multiples of 16 entries, we split them into a new record entry for that specific mutation.
    """
    if "INFO" not in record.keys():
        return {}, {}
    annotations = record["INFO"]
    # This dictionary holds the information for a single record
    rec_dict = {}

    # Check for pathogenic variants in CLNSIG field
    rec_dict["is_pathogenic"] = False
    if "CLNSIG=" in annotations:
        for field in annotations.split(";"):
            if field.startswith("CLNSIG="):
                clnsig_value = field.split("=")[1]
                if "PATHOGENIC" in clnsig_value.upper():
                    rec_dict["is_pathogenic"] = True
                break

    # This dictionary holds one or more records derived from the annotation field.
    ann_dict = {}

    for ann in annotations.split(";"):
        if "=" in ann:
            mykey = ann.split("=")[0]
            myvalue = ann.split("=")[1]
            if "|" in myvalue:
                if mykey == "ANN":
                    if len(myvalue.split("|")) == 16:
                        try:
                            myvalues = myvalue.split("|")
                            count = 0
                            ann_dict[count] = dict()
                            ann_dict[count]["Allele"] = myvalues[0]
                            ann_dict[count]["Annotation"] = myvalues[1]
                            ann_dict[count]["Annotation_Impact"] = myvalues[2]
                            ann_dict[count]["Gene_Name"] = myvalues[3]
                            ann_dict[count]["Gene_ID"] = myvalues[4]
                            ann_dict[count]["Feature_Type"] = myvalues[5]
                            ann_dict[count]["Feature_ID"] = myvalues[6]
                            ann_dict[count]["Transcript_BioType"] = myvalues[7]
                            ann_dict[count]["Rank"] = myvalues[8]
                            ann_dict[count]["HGVS.c"] = myvalues[9]
                            ann_dict[count]["HGVS.p"] = myvalues[10]
                            ann_dict[count]["cDNA.pos / cDNA.length"] = myvalues[11]
                            ann_dict[count]["CDS.pos / CDS.length"] = myvalues[12]
                            ann_dict[count]["AA.pos / AA.length"] = myvalues[13]
                            ann_dict[count]["Distance"] = myvalues[14]
                            ann_dict[count]["ERRORS / WARNINGS / INFO"] = myvalues[15]
                        except Exception as e:
                            logger.error(f"Error parsing annotation: {e}")
                    elif len(myvalue.split("|")) > 16:
                        count = 0
                        for chunk in myvalue.split(","):
                            try:
                                myvalues = chunk.split("|")
                                if len(myvalues) >= 16:
                                    ann_dict[count] = dict()
                                    ann_dict[count]["Allele"] = myvalues[0]
                                    ann_dict[count]["Annotation"] = myvalues[1]
                                    ann_dict[count]["Annotation_Impact"] = myvalues[2]
                                    ann_dict[count]["Gene_Name"] = myvalues[3]
                                    ann_dict[count]["Gene_ID"] = myvalues[4]
                                    ann_dict[count]["Feature_Type"] = myvalues[5]
                                    ann_dict[count]["Feature_ID"] = myvalues[6]
                                    ann_dict[count]["Transcript_BioType"] = myvalues[7]
                                    ann_dict[count]["Rank"] = myvalues[8]
                                    ann_dict[count]["HGVS.c"] = myvalues[9]
                                    ann_dict[count]["HGVS.p"] = myvalues[10]
                                    ann_dict[count]["cDNA.pos / cDNA.length"] = myvalues[11]
                                    ann_dict[count]["CDS.pos / CDS.length"] = myvalues[12]
                                    ann_dict[count]["AA.pos / AA.length"] = myvalues[13]
                                    ann_dict[count]["Distance"] = myvalues[14]
                                    ann_dict[count]["ERRORS / WARNINGS / INFO"] = myvalues[15]
                                    count += 1
                            except Exception as e:
                                logger.error(f"Error parsing annotation chunk: {e}")
            else:
                rec_dict[mykey] = myvalue

    return ann_dict, rec_dict


def parse_vcf(vcf_path: Path) -> Optional[pd.DataFrame]:
    """
    Parse a VCF file and return a processed DataFrame.
    
    Args:
        vcf_path: Path to the VCF file
        
    Returns:
        Processed DataFrame with annotations exploded, or None if parsing fails
    """
    try:
        if not vcf_path.exists():
            logger.warning(f"VCF file not found: {vcf_path}")
            return None
            
        header = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT GT".split()
        vcf = pd.read_csv(vcf_path, delimiter="\t", comment="#", names=header)
        
        if len(vcf) == 0:
            logger.info(f"VCF file is empty: {vcf_path}")
            return pd.DataFrame()
        
        explodedvcf = []
        for record in vcf.to_dict("records"):
            result, result2 = process_annotations(record)
            if len(result) > 0:
                for res in result:
                    dat = {**record, **result[res], **result2}
                    explodedvcf.append(dat)
            elif len(result2) > 0:
                # Include records even if no annotations
                dat = {**record, **result2}
                explodedvcf.append(dat)
            else:
                # Include record even if no INFO processing
                explodedvcf.append(record)

        vcf_df = pd.DataFrame.from_records(explodedvcf)
        if "INFO" in vcf_df.columns:
            vcf_df = vcf_df.drop(columns=["INFO"]).drop_duplicates()
        else:
            vcf_df = vcf_df.drop_duplicates()

        def set_unique_values(series):
            set_series = set(series.dropna())
            if len(set_series) > 0:
                return "{}".format(", ".join(map(str, set_series)))
            return None

        if len(vcf_df) > 0:
            shared_columns = [
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "FORMAT",
                "GT",
            ]
            
            # Add Allele if it exists
            if "Allele" in vcf_df.columns:
                shared_columns.append("Allele")

            # Define columns to be aggregated
            non_shared_columns = [
                col for col in vcf_df.columns if col not in shared_columns
            ]

            vcf_df = vcf_df.replace({np.nan: None})
            result = (
                vcf_df.groupby(shared_columns)[non_shared_columns]
                .agg(set_unique_values)
                .reset_index()
            )

            return result
        else:
            return pd.DataFrame()
            
    except Exception as e:
        logger.error(f"Error parsing VCF file {vcf_path}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None


def navigate_igv_to_snp(chrom: str, pos: int, window_size: int = 1000) -> None:
    """
    Navigate IGV browser to a specific SNP location.
    
    Args:
        chrom: Chromosome name (e.g., "chr1")
        pos: Position on the chromosome
        window_size: Size of the window to display (default: 1000 bp)
    """
    try:
        # Ensure chromosome name has 'chr' prefix if needed
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        
        # Calculate window around the SNP
        start = max(1, pos - window_size // 2)
        end = pos + window_size // 2
        
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
    vcf_files = [
        clair3_dir / "snpsift_output.vcf",
        clair3_dir / "snpsift_indel_output.vcf",
    ]
    
    vcf_file = None
    for vcf_path in vcf_files:
        if vcf_path.exists():
            vcf_file = vcf_path
            break
    
    if not vcf_file:
        # Show placeholder indicating SNP analysis is in progress or not available
        with ui.card().classes("w-full"):
            ui.label("SNP Analysis").classes("text-lg font-semibold text-blue-800")
            ui.separator().classes().style("border: 1px solid var(--md-primary)")
            ui.label(
                "SNP analysis results are not yet available. "
                "SNP analysis will run automatically after target.bam is generated."
            ).classes("text-body-medium text-gray-600")
        return
    
    # Parse VCF file
    vcf_df = parse_vcf(vcf_file)
    
    if vcf_df is None or vcf_df.empty:
        with ui.card().classes("w-full"):
            ui.label("SNP Analysis").classes("text-lg font-semibold text-blue-800")
            ui.separator().classes().style("border: 1px solid var(--md-primary)")
            ui.label("No variants found in SNP analysis.").classes("text-body-medium text-gray-600")
        return
    
    # Create SNP section
    with ui.card().classes("w-full"):
        ui.label("SNP Analysis").classes("text-lg font-semibold text-blue-800")
        ui.separator().classes().style("border: 1px solid var(--md-primary)")
        
        # Count pathogenic variants
        pathogenic_count = 0
        if "is_pathogenic" in vcf_df.columns:
            if vcf_df["is_pathogenic"].dtype == bool:
                pathogenic_count = vcf_df["is_pathogenic"].sum()
            else:
                # Handle string values
                pathogenic_count = vcf_df["is_pathogenic"].apply(
                    lambda x: str(x).upper() in ("YES", "TRUE", "1") if pd.notna(x) else False
                ).sum()
        
        # Summary statistics
        with ui.row().classes("w-full gap-4 mb-4"):
            ui.label(f"Total Variants: {len(vcf_df)}").classes("text-body-medium font-semibold")
            if pathogenic_count > 0:
                ui.label(f"Pathogenic Variants: {pathogenic_count}").classes(
                    "text-body-medium font-semibold text-red-600"
                )
        
        # Create table columns
        columns = []
        
        # Essential columns
        essential_cols = ["CHROM", "POS", "REF", "ALT", "QUAL", "Gene_Name"]
        for col in essential_cols:
            if col in vcf_df.columns:
                columns.append({
                    "name": col,
                    "label": col.replace("_", " "),
                    "field": col,
                    "sortable": True,
                })
        
        # Add annotation columns if available
        annotation_cols = [
            "Annotation",
            "Annotation_Impact",
            "HGVS.c",
            "HGVS.p",
            "is_pathogenic",
        ]
        for col in annotation_cols:
            if col in vcf_df.columns:
                columns.append({
                    "name": col,
                    "label": col.replace("_", " "),
                    "field": col,
                    "sortable": True,
                })
        
        # Add IGV link column for pathogenic variants
        if "is_pathogenic" in vcf_df.columns:
            columns.append({
                "name": "action",
                "label": "View in IGV",
                "field": "action",
                "sortable": False,
            })
        
        # Prepare rows for display
        rows = []
        for idx, row in vcf_df.iterrows():
            row_dict = {}
            for col in columns:
                if col["field"] != "action":
                    value = row.get(col["field"])
                    if pd.isna(value):
                        value = ""
                    elif isinstance(value, bool):
                        value = "Yes" if value else "No"
                    else:
                        value = str(value)
                    row_dict[col["field"]] = value
            
            # Add IGV link for pathogenic variants
            is_pathogenic_value = row.get("is_pathogenic", False)
            if isinstance(is_pathogenic_value, bool):
                is_pathogenic_value = is_pathogenic_value
            elif isinstance(is_pathogenic_value, str):
                is_pathogenic_value = is_pathogenic_value.upper() in ("YES", "TRUE", "1")
            else:
                is_pathogenic_value = False
                
            if "is_pathogenic" in vcf_df.columns and is_pathogenic_value:
                row_dict["action"] = "🔍"
            else:
                row_dict["action"] = ""
            
            rows.append(row_dict)
        
        # Create JavaScript map of SNP regions for IGV navigation
        import json
        snp_regions_map = {}
        for idx, row_data in enumerate(rows):
            # Check if pathogenic (handle both boolean and string)
            is_path = row_data.get("is_pathogenic", False)
            if isinstance(is_path, str):
                is_path = is_path.upper() in ("YES", "TRUE", "1")
            
            if is_path:  # Only add pathogenic SNPs
                chrom = row_data.get("CHROM", "")
                pos_str = row_data.get("POS", "")
                if chrom and pos_str:
                    try:
                        pos = int(pos_str)
                        # Create a unique key for this SNP
                        snp_key = f"{chrom}:{pos}"
                        snp_regions_map[snp_key] = {"chrom": chrom, "pos": pos}
                    except (ValueError, TypeError):
                        pass
        
        js_snp_regions_json = json.dumps(snp_regions_map)
        
        # Function to navigate IGV to a SNP region
        def navigate_to_snp_region(snp_key: str):
            """Navigate IGV browser to the specified SNP region."""
            if snp_key in snp_regions_map:
                snp_data = snp_regions_map[snp_key]
                chrom = snp_data["chrom"]
                pos = snp_data["pos"]
                navigate_igv_to_snp(chrom, pos)
        
        # Initialize SNP regions map in window BEFORE creating table
        js_init_snp_regions = f"""
            (function() {{
                window.snpRegionsMap = window.snpRegionsMap || {{}};
                Object.assign(window.snpRegionsMap, {js_snp_regions_json});
                console.log('[SNP] Initialized SNP regions map with', Object.keys(window.snpRegionsMap).length, 'pathogenic SNPs');
            }})();
        """
        ui.run_javascript(js_init_snp_regions, timeout=5.0)
        
        # Create table using styled_table
        from robin.gui.theme import styled_table
        
        table_container, snp_table = styled_table(
            columns=columns,
            rows=rows,
            pagination=25,
            class_size="table-xs"
        )
        
        # Add filter for pathogenic variants
        with snp_table.add_slot("top-left"):
            pathogenic_filter = ui.switch(
                "Show pathogenic variants only",
                value=False,
            )
        
        # Add search functionality
        with snp_table.add_slot("top-right"):
            with ui.input(placeholder="Search SNPs...").props("type=search").bind_value(
                snp_table, "filter"
            ).add_slot("append"):
                ui.icon("search")
        
        # Add clickable action button column using slot for pathogenic variants
        if "is_pathogenic" in vcf_df.columns:
            try:
                # Use Vue event emission which is more reliable than window functions
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
"""
                )
                
                # Handle the event from the slot
                def on_snp_view_igv(e):
                    """Handle SNP view IGV event from table button."""
                    try:
                        snp_key = e.args if isinstance(e.args, str) else getattr(e, 'args', None)
                        if snp_key:
                            logging.debug(f"[SNP] Button clicked for: {snp_key}")
                            navigate_to_snp_region(snp_key)
                    except Exception as ex:
                        logger.debug(f"Error handling SNP IGV view: {ex}")
                
                snp_table.on("snp-view-igv", on_snp_view_igv)
                
            except Exception as ex:
                logger.warning(f"Could not add action button slot: {ex}")
        
        # Apply filter
        def apply_pathogenic_filter(value: bool):
            if value:
                # Filter to show only pathogenic variants
                filtered_rows = []
                for r in rows:
                    is_path = r.get("is_pathogenic", False)
                    if isinstance(is_path, str):
                        is_path = is_path.upper() in ("YES", "TRUE", "1")
                    if is_path:
                        filtered_rows.append(r)
                snp_table.rows = filtered_rows
            else:
                # Show all variants
                snp_table.rows = rows
            snp_table.update()
        
        pathogenic_filter.on("update:model-value", lambda e: apply_pathogenic_filter(e.args))
        
        # Make columns sortable
        for col in snp_table.columns:
            col["sortable"] = True

