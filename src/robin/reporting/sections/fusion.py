"""
fusion.py

This module contains the fusion analysis section of the report.
"""

import os
import logging
import pickle
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.units import inch
from .base import ReportSection
import pandas as pd
import networkx as nx
from reportlab.lib.styles import ParagraphStyle

logger = logging.getLogger(__name__)


class FusionSection(ReportSection):
    """Section containing the fusion analysis results."""

    def _get_gene_pairs(self, data):
        """Get unique gene fusion pairs from the processed data."""
        if data is None or data.empty:
            return []

        # For processed data, we can use the 'tag' column which contains gene pairs
        # or fall back to grouping by read_id if needed
        if 'tag' in data.columns:
            # Extract gene pairs from the tag column (format: "gene1,gene2")
            gene_pairs = []
            for tag in data['tag'].unique():
                if ',' in tag:
                    genes = tag.split(',')
                    if len(genes) >= 2:
                        # Sort genes to ensure consistent ordering
                        gene_pair = tuple(sorted(genes))
                        gene_pairs.append(gene_pair)
            
            # Get unique pairs
            unique_pairs = list(set(gene_pairs))
            return unique_pairs
        else:
            # Fallback to original method if tag column not available
            return self._get_gene_pairs_fallback(data)

    def _get_gene_pairs_fallback(self, data):
        """Fallback method for getting gene pairs from raw data."""
        if data is None or data.empty:
            return []

        # Group reads by readID to find pairs
        read_groups = data.groupby("readID")
        gene_pairs = []

        # First collect all read IDs for each gene pair
        gene_pair_reads = {}
        for _, group in read_groups:
            genes = group["Gene"].unique()
            if len(genes) >= 2:
                # Sort genes to ensure consistent ordering
                genes = sorted(genes)
                for i in range(len(genes) - 1):
                    for j in range(i + 1, len(genes)):
                        pair = (genes[i], genes[j])
                        if pair not in gene_pair_reads:
                            gene_pair_reads[pair] = set()
                        gene_pair_reads[pair].add(group["readID"].iloc[0])

        # Only include pairs with 3 or more supporting reads
        for pair, reads in gene_pair_reads.items():
            if len(reads) >= 3:
                gene_pairs.append(pair)

        # Get unique pairs
        unique_pairs = list(set(gene_pairs))
        return unique_pairs

    def _get_gene_networks(self, gene_pairs):
        """Create networks of connected genes."""
        G = nx.Graph()
        for pair in gene_pairs:
            G.add_edge(pair[0], pair[1])
        connected_components = list(nx.connected_components(G))
        return [list(component) for component in connected_components]

    def _load_fusion_data(self):
        """Load fusion data from the processed pickle files."""
        fusion_data = {
            "master_candidates": None,
            "all_candidates": None,
            "master_path": os.path.join(
                self.report.output, "fusion_candidates_master_processed.csv"
            ),
            "all_path": os.path.join(self.report.output, "fusion_candidates_all_processed.csv"),
        }

        try:
            # Load master fusion candidates from processed pickle
            if os.path.exists(fusion_data["master_path"]):
                with open(fusion_data["master_path"], "rb") as f:
                    try:
                        processed_data = pickle.load(f)
                        # Filter to only good pairs to reduce memory usage and processing time
                        annotated_data = processed_data.get("annotated_data", pd.DataFrame())
                        goodpairs = processed_data.get("goodpairs", pd.Series())
                        
                        if not annotated_data.empty and not goodpairs.empty:
                            # Only keep the good pairs (same as GUI does)
                            fusion_data["master_candidates"] = annotated_data[goodpairs]
                        else:
                            fusion_data["master_candidates"] = annotated_data
                        
                        logger.debug(
                            f"Loaded master fusion candidates from {fusion_data['master_path']} "
                            f"({len(fusion_data['master_candidates'])} good pairs from {len(annotated_data)} total candidates)"
                        )
                    except (pickle.UnpicklingError, EOFError) as e:
                        logger.warning(f"Error loading master fusion pickle: {e}")
                        fusion_data["master_candidates"] = pd.DataFrame()

            # Load all fusion candidates from processed pickle
            if os.path.exists(fusion_data["all_path"]):
                with open(fusion_data["all_path"], "rb") as f:
                    try:
                        processed_data = pickle.load(f)
                        # Filter to only good pairs to reduce memory usage and processing time
                        annotated_data = processed_data.get("annotated_data", pd.DataFrame())
                        goodpairs = processed_data.get("goodpairs", pd.Series())
                        
                        if not annotated_data.empty and not goodpairs.empty:
                            # Only keep the good pairs (same as GUI does)
                            fusion_data["all_candidates"] = annotated_data[goodpairs]
                        else:
                            fusion_data["all_candidates"] = annotated_data
                        
                        logger.debug(
                            f"Loaded all fusion candidates from {fusion_data['all_path']} "
                            f"({len(fusion_data['all_candidates'])} good pairs from {len(annotated_data)} total candidates)"
                        )
                    except (pickle.UnpicklingError, EOFError) as e:
                        logger.warning(f"Error loading all fusion pickle: {e}")
                        fusion_data["all_candidates"] = pd.DataFrame()

        except Exception as e:
            logger.error(f"Error loading fusion data: {str(e)}")
            logger.debug("Exception details:", exc_info=True)

        return fusion_data

    def _format_fusion_table(self, data, title):
        """Create a formatted table for fusion data."""
        if data is None or data.empty:
            return None

        # Create paragraph styles for cells
        header_style = ParagraphStyle(
            "HeaderStyle",
            parent=self.styles.styles["Normal"],
            fontName="Helvetica-Bold",
            fontSize=8,
            textColor=self.styles.COLORS["primary"],
            leading=10,
            spaceBefore=0,
            spaceAfter=0,
        )

        cell_style = ParagraphStyle(
            "CellStyle",
            parent=self.styles.styles["Normal"],
            fontName="Helvetica",
            fontSize=8,
            textColor=colors.black,
            leading=10,
            spaceBefore=0,
            spaceAfter=0,
        )

        # Prepare table data with updated headers using Paragraph objects
        table_data = [
            [
                Paragraph("Fusion Pair", header_style),
                Paragraph("Chr 1", header_style),
                Paragraph("Chr 2", header_style),
                Paragraph("Gene 1 Position", header_style),
                Paragraph("Gene 2 Position", header_style),
                Paragraph("Reads", header_style),
            ]
        ]

        # Track processed pairs to avoid duplicates
        processed_pairs = set()
        fusion_details = {}

        if 'tag' in data.columns:
            # Use processed data structure
            for tag in data['tag'].unique():
                if ',' in tag and tag not in processed_pairs:
                    processed_pairs.add(tag)
                    genes = tag.split(',')
                    if len(genes) >= 2:
                        gene_pair = f"{genes[0]}-{genes[1]}"
                        
                        # Get supporting reads count
                        supporting_reads = len(data[data['tag'] == tag]['read_id'].unique())
                        
                        # Get gene information from the data
                        gene_info = {}
                        for gene in genes:
                            gene_data = data[data['col4'] == gene]
                            if not gene_data.empty:
                                first_row = gene_data.iloc[0]
                                gene_info[gene] = {
                                    "chrom": first_row.get('reference_id', 'Unknown'),
                                    "position": f"{first_row.get('reference_start', 'N/A')}-{first_row.get('reference_end', 'N/A')}",
                                }

                        # Only store fusion details if we have complete information
                        if len(gene_info) == 2 and supporting_reads >= 3:
                            fusion_details[gene_pair] = {
                                "chrom1": gene_info[genes[0]]["chrom"],
                                "chrom2": gene_info[genes[1]]["chrom"],
                                "pos1": gene_info[genes[0]]["position"],
                                "pos2": gene_info[genes[1]]["position"],
                                "supporting_reads": supporting_reads,
                            }
        else:
            # Fallback to original method for raw data
            return self._format_fusion_table_fallback(data, title, header_style, cell_style)

        # Add rows to table using Paragraph objects for wrappable text
        for gene_pair, details in fusion_details.items():
            table_data.append(
                [
                    Paragraph(gene_pair, cell_style),
                    Paragraph(details["chrom1"], cell_style),
                    Paragraph(details["chrom2"], cell_style),
                    Paragraph(details["pos1"], cell_style),
                    Paragraph(details["pos2"], cell_style),
                    Paragraph(str(details["supporting_reads"]), cell_style),
                ]
            )

        # Return None if no fusions meet the threshold
        if len(table_data) == 1:  # Only header row
            return None

        # Define column widths (in inches) - adjusted for better proportions
        col_widths = [1.8, 0.5, 0.5, 1.6, 1.6, 0.6]

        # Create and style the table
        table = Table(
            table_data, colWidths=[inch * width for width in col_widths], repeatRows=1
        )
        table.setStyle(
            TableStyle(
                [
                    # Inherit modern table style
                    *self.MODERN_TABLE_STYLE._cmds,
                    # Preserve specific alignments
                    ("ALIGN", (0, 0), (0, -1), "LEFT"),  # Left-align Fusion Pair
                    ("ALIGN", (1, 0), (2, -1), "CENTER"),  # Center-align Chromosomes
                    ("ALIGN", (3, 0), (4, -1), "LEFT"),  # Left-align positions
                    ("ALIGN", (5, 0), (5, -1), "RIGHT"),  # Right-align Supporting Reads
                ]
            )
        )

        return table

    def _format_fusion_table_fallback(self, data, title, header_style, cell_style):
        """Fallback method for formatting fusion table from raw data."""
        # Group by readID to get fusion pairs
        read_groups = data.groupby("readID")

        # Prepare table data with updated headers using Paragraph objects
        table_data = [
            [
                Paragraph("Fusion Pair", header_style),
                Paragraph("Chr 1", header_style),
                Paragraph("Chr 2", header_style),
                Paragraph("Gene 1 Position", header_style),
                Paragraph("Gene 2 Position", header_style),
                Paragraph("Reads", header_style),
            ]
        ]

        # Track processed pairs to avoid duplicates
        processed_pairs = set()
        fusion_details = {}

        # First pass: collect all read IDs for each gene pair
        gene_pair_reads = {}
        for read_id, group in read_groups:
            if len(group) >= 2:
                genes = sorted(
                    group["Gene"].unique()
                )  # Sort to ensure consistent ordering
                if len(genes) >= 2:
                    for i in range(len(genes) - 1):
                        for j in range(i + 1, len(genes)):
                            gene_pair = f"{genes[i]}-{genes[j]}"
                            if gene_pair not in gene_pair_reads:
                                gene_pair_reads[gene_pair] = set()
                            gene_pair_reads[gene_pair].add(read_id)

        # Second pass: process each gene pair with 3 or more supporting reads
        for gene_pair, reads in gene_pair_reads.items():
            if len(reads) >= 3 and gene_pair not in processed_pairs:
                processed_pairs.add(gene_pair)
                genes = gene_pair.split("-")

                # Get details for each gene in the pair
                gene_info = {}
                for gene in genes:
                    gene_filtered = data[data["Gene"] == gene]
                    if gene_filtered.empty:
                        logger.warning(f"No data found for gene: {gene}")
                        continue
                    gene_data = gene_filtered.iloc[0]
                    gene_info[gene] = {
                        "chrom": gene_data["chromBED"],
                        "position": f"{gene_data['BS']}-{gene_data['BE']}",
                    }

                # Only store fusion details if we have complete information for both genes
                if (
                    len(gene_info) == 2
                    and genes[0] in gene_info
                    and genes[1] in gene_info
                ):
                    fusion_details[gene_pair] = {
                        "chrom1": gene_info[genes[0]]["chrom"],
                        "chrom2": gene_info[genes[1]]["chrom"],
                        "pos1": gene_info[genes[0]]["position"],
                        "pos2": gene_info[genes[1]]["position"],
                        "supporting_reads": len(reads),
                    }
                else:
                    logger.warning(
                        f"Incomplete gene information for fusion pair: {gene_pair}"
                    )

        # Add rows to table using Paragraph objects for wrappable text
        for gene_pair, details in fusion_details.items():
            table_data.append(
                [
                    Paragraph(gene_pair, cell_style),
                    Paragraph(details["chrom1"], cell_style),
                    Paragraph(details["chrom2"], cell_style),
                    Paragraph(details["pos1"], cell_style),
                    Paragraph(details["pos2"], cell_style),
                    Paragraph(str(details["supporting_reads"]), cell_style),
                ]
            )

        # Return None if no fusions meet the threshold
        if len(table_data) == 1:  # Only header row
            return None

        # Define column widths (in inches) - adjusted for better proportions
        col_widths = [1.8, 0.5, 0.5, 1.6, 1.6, 0.6]

        # Create and style the table
        table = Table(
            table_data, colWidths=[inch * width for width in col_widths], repeatRows=1
        )
        table.setStyle(
            TableStyle(
                [
                    # Inherit modern table style
                    *self.MODERN_TABLE_STYLE._cmds,
                    # Preserve specific alignments
                    ("ALIGN", (0, 0), (0, -1), "LEFT"),  # Left-align Fusion Pair
                    ("ALIGN", (1, 0), (2, -1), "CENTER"),  # Center-align Chromosomes
                    ("ALIGN", (3, 0), (4, -1), "LEFT"),  # Left-align positions
                    ("ALIGN", (5, 0), (5, -1), "RIGHT"),  # Right-align Supporting Reads
                ]
            )
        )

        return table

    def add_content(self):
        """Add the fusion analysis content to the report."""
        logger.debug("Starting fusion section content generation")

        # Add section title
        self.elements.append(
            Paragraph("Gene Fusion Analysis", self.styles.styles["Heading1"])
        )
        self.elements.append(Spacer(1, 12))

        # Load fusion data
        fusion_data = self._load_fusion_data()

        # Get unique fusion pairs
        master_pairs = self._get_gene_pairs(fusion_data["master_candidates"])
        all_pairs = self._get_gene_pairs(fusion_data["all_candidates"])

        # Get gene networks
        master_networks = self._get_gene_networks(master_pairs)
        all_networks = self._get_gene_networks(all_pairs)

        # Create summary text
        summary_text = []
        if master_networks:
            summary_text.append(
                f"Found {len(master_networks)} fusion networks within the target panel:<br/>"
            )
            for network in master_networks:
                summary_text.append(f"• {' - '.join(network)}")
        else:
            summary_text.append(
                "No fusion networks identified within the target panel."
            )

        if all_networks and len(all_networks) > len(master_networks):
            summary_text.append(
                f"<br/>Additional {len(all_networks) - len(master_networks)} genome-wide fusion networks identified."
            )

        summary_paragraph = "\n".join(summary_text)

        # Add summary to summary section of report
        self.summary_elements.append(
            Paragraph("Gene Fusions", self.styles.styles["Heading3"])
        )
        self.summary_elements.append(
            Paragraph(summary_paragraph, self.styles.styles["Normal"])
        )

        # Add detailed analysis section if there are any candidates
        if master_pairs or all_pairs:
            self.elements.append(
                Paragraph("Detailed Analysis", self.styles.styles["Heading2"])
            )
            self.elements.append(Spacer(1, 12))

            # Add target panel fusions
            if master_pairs:
                self.elements.append(
                    Paragraph("Target Panel Fusions", self.styles.styles["Heading3"])
                )
                table = self._format_fusion_table(
                    fusion_data["master_candidates"], "Target Panel Fusions"
                )
                if table:
                    self.elements.append(table)
                    self.elements.append(Spacer(1, 12))

            # Add genome-wide fusions
            if all_pairs and len(all_pairs) > len(master_pairs):
                self.elements.append(
                    Paragraph("Genome-wide Fusions", self.styles.styles["Heading3"])
                )
                table = self._format_fusion_table(
                    fusion_data["all_candidates"], "Genome-wide Fusions"
                )
                if table:
                    self.elements.append(table)
                    self.elements.append(Spacer(1, 12))

            # Add source information
            source_text = "Data sources:\n"
            source_text += f"- Target panel fusions: {fusion_data['master_path']}\n"
            source_text += f"- Genome-wide fusions: {fusion_data['all_path']}\n"
            source_text += "\nNote: Fusion candidates are identified from reads with supplementary alignments and filtered based on mapping quality (>50) and alignment length (>200bp). Processed data includes pre-computed gene pairs and networks."

            self.elements.append(Paragraph(source_text, self.styles.styles["Normal"]))

            # Build export DataFrames
            try:
                import pandas as pd

                # Summaries
                networks_rows = []
                for scope, nets in (
                    ("target", master_networks),
                    ("genome_wide", all_networks),
                ):
                    for net in nets:
                        networks_rows.append(
                            {
                                "Scope": scope,
                                "Network": " - ".join(net),
                                "Genes": ", ".join(net),
                            }
                        )
                if networks_rows:
                    self.export_frames["fusions_networks"] = pd.DataFrame(networks_rows)

                # Detailed tables if available
                def _fusion_rows(df):
                    if df is None or df.empty:
                        return []
                    
                    if 'tag' in df.columns:
                        # Use processed data structure
                        rows = []
                        for tag in df['tag'].unique():
                            if ',' in tag:
                                genes = tag.split(',')
                                if len(genes) >= 2:
                                    supporting_reads = len(df[df['tag'] == tag]['read_id'].unique())
                                    if supporting_reads >= 3:
                                        gene_pair = f"{genes[0]}-{genes[1]}"
                                        
                                        # Get gene information
                                        gene_info = {}
                                        for gene in genes:
                                            gene_data = df[df['col4'] == gene]
                                            if not gene_data.empty:
                                                first_row = gene_data.iloc[0]
                                                gene_info[gene] = {
                                                    "chrom": first_row.get('reference_id', 'Unknown'),
                                                    "position": f"{first_row.get('reference_start', 'N/A')}-{first_row.get('reference_end', 'N/A')}",
                                                }
                                        
                                        if len(gene_info) == 2:
                                            rows.append({
                                                "FusionPair": gene_pair,
                                                "Chrom1": gene_info[genes[0]]["chrom"],
                                                "Chrom2": gene_info[genes[1]]["chrom"],
                                                "Gene1Pos": gene_info[genes[0]]["position"],
                                                "Gene2Pos": gene_info[genes[1]]["position"],
                                                "SupportingReads": supporting_reads,
                                            })
                        return rows
                    else:
                        # Fallback to original method for raw data
                        read_groups = df.groupby("readID")
                        gene_pair_reads = {}
                        for read_id, group in read_groups:
                            genes = sorted(group["Gene"].unique())
                            if len(genes) >= 2:
                                for i in range(len(genes) - 1):
                                    for j in range(i + 1, len(genes)):
                                        gene_pair = f"{genes[i]}-{genes[j]}"
                                        gene_pair_reads.setdefault(gene_pair, set()).add(
                                            read_id
                                        )
                        rows = []
                        for gene_pair, reads in gene_pair_reads.items():
                            if len(reads) >= 3:
                                g1, g2 = gene_pair.split("-")
                                g1row = df[df["Gene"] == g1].iloc[0]
                                g2row = df[df["Gene"] == g2].iloc[0]
                                rows.append(
                                    {
                                        "FusionPair": gene_pair,
                                        "Chrom1": g1row["chromBED"],
                                        "Chrom2": g2row["chromBED"],
                                        "Gene1Pos": f"{g1row['BS']}-{g1row['BE']}",
                                        "Gene2Pos": f"{g2row['BS']}-{g2row['BE']}",
                                        "SupportingReads": len(reads),
                                    }
                                )
                        return rows

                master_rows = _fusion_rows(fusion_data["master_candidates"])
                if master_rows:
                    self.export_frames["fusions_target_panel"] = pd.DataFrame(
                        master_rows
                    )

                all_rows = _fusion_rows(fusion_data["all_candidates"])
                if all_rows:
                    self.export_frames["fusions_genome_wide"] = pd.DataFrame(all_rows)
            except Exception:
                pass
        else:
            self.elements.append(
                Paragraph(
                    "No fusion candidates were identified.",
                    self.styles.styles["Normal"],
                )
            )
