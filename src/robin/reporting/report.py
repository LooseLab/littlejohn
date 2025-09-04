"""
report.py

This module contains the main report class that coordinates the generation of the PDF report.
"""

import os
import json
import logging
import pandas as pd
from datetime import datetime
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from .styling.styles import ReportStyles
from robin.gui import fonts

logger = logging.getLogger(__name__)


class RobinReport:
    """Main class for generating ROBIN PDF reports."""

    def __init__(self, filename, output, center: str):
        """Initialize the report generator.

        Args:
            filename: Output PDF filename
            output: Directory containing analysis output files
            center: Center ID running the analysis
        """
        self.filename = filename
        self.output = output
        self.center = center
        self.sample_id = os.path.basename(os.path.normpath(output))

        # Handle filename with None prefix
        if filename.startswith("None"):
            final_folder = os.path.basename(os.path.normpath(output))
            self.filename = filename.replace("None", final_folder, 1)
            self.sample_id = final_folder

        # Initialize styling
        self.fonts_dir = os.path.join(os.path.dirname(os.path.abspath(fonts.__file__)))
        self.styles = ReportStyles(self.fonts_dir)

        # Initialize document elements
        self.elements_summary = []
        self.elements = []
        self.end_of_report_elements = []

        # Load master data
        self.masterdf = self._load_master_data()
        self.centreID = self._get_centre_id()

        # Create document
        self.doc = self._create_document()

        # Initialize sections
        self.sections = []
        self._initialize_sections()

    def _load_master_data(self):
        """Load the master data file if it exists."""
        master_path = os.path.join(self.output, "master.csv")
        if os.path.exists(master_path):
            return pd.read_csv(master_path)
        return None

    def _get_centre_id(self):
        """Get the centre ID from the provided center parameter."""
        return self.center

    def _create_document(self):
        """Create the PDF document with enhanced M3 margins and settings."""
        return SimpleDocTemplate(
            self.filename,
            pagesize=A4,
            rightMargin=1.0 * inch,  # Enhanced from 0.75 inch
            leftMargin=1.0 * inch,  # Enhanced from 0.75 inch
            topMargin=1.35 * inch,  # Keep existing for header compatibility
            bottomMargin=1.0 * inch,  # Enhanced from 0.75 inch
        )

    def _initialize_sections(self):
        """Initialize all report sections."""
        # Import sections here to avoid circular imports
        from .sections.classification import ClassificationSection
        from .sections.cnv import CNVSection
        from .sections.fusion import FusionSection
        from .sections.coverage import CoverageSection
        from .sections.mgmt import MGMTSection
        from .sections.run_data import RunDataSection
        from .sections.disclaimer import DisclaimerSection
        from .sections.variants import VariantsSection

        # Add sections in order
        self.sections = [
            ClassificationSection(self),
            CNVSection(self),
            VariantsSection(self),
            FusionSection(self),
            CoverageSection(self),
            MGMTSection(self),
            RunDataSection(self),
            DisclaimerSection(self),
        ]

    def generate_report(
        self,
        report_type="detailed",
        export_csv_dir=None,
        export_xlsx=False,
        export_zip=False,
    ):
        """Generate the complete PDF report.

        Args:
            report_type: Type of report to generate ('summary' or 'detailed')

        Returns:
            Path to the generated PDF file
        """
        try:
            logger.info("Starting report generation")

            # Add summary section header with enhanced M3 typography
            self.elements_summary.insert(
                0,
                Paragraph(
                    f"Summary - {self.sample_id}", self.styles.styles["DisplayMedium"]
                ),
            )

            # Add enhanced summary card with M3 styling
            summary_card_content = f"""
            <b>ROBIN Analysis Report</b><br/>
            Sample ID: {self.sample_id}<br/>
            Centre ID: {self.centreID if self.centreID else 'Not specified'}<br/>
            Report Type: {report_type.title()}<br/>
            Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
            """
            self.elements_summary.insert(
                1, Paragraph(summary_card_content, self.styles.styles["InfoCard"])
            )
            self.elements_summary.insert(2, Spacer(1, 16))

            # Process each section
            for section in self.sections:
                try:
                    section.add_content()
                    summary_elements, main_elements = section.get_elements()

                    # Always include summary elements
                    self.elements_summary.extend(summary_elements)

                    # For detailed report or if it's the disclaimer section, include main elements
                    if (
                        report_type == "detailed"
                        or section.__class__.__name__ == "DisclaimerSection"
                    ):
                        self.elements.extend(main_elements)
                except Exception as e:
                    logger.error(
                        f"Error processing section {section.__class__.__name__}: {e}",
                        exc_info=True,
                    )
                    # Add error card to report for better user feedback
                    error_content = f"""
                    <b>⚠️ Section Processing Error</b><br/>
                    Section: {section.__class__.__name__}<br/>
                    Error: {str(e)[:100]}{'...' if len(str(e)) > 100 else ''}<br/>
                    <i>This section was skipped due to processing errors.</i>
                    """
                    self.elements_summary.append(
                        Paragraph(error_content, self.styles.styles["Error"])
                    )

            # Add detailed analysis header and elements only for detailed reports
            if report_type == "detailed":
                self.elements.insert(0, PageBreak())
                self.elements.insert(
                    1,
                    Paragraph("Detailed Analysis", self.styles.styles["HeadlineLarge"]),
                )
                self.elements.insert(2, Spacer(1, 16))  # Enhanced spacing

            # Add page break before end of report elements
            if self.end_of_report_elements:
                self.end_of_report_elements.insert(0, PageBreak())

            # Combine all elements
            logger.info("Combining elements for final PDF")
            if report_type == "detailed":
                final_elements = (
                    self.elements_summary + self.elements + self.end_of_report_elements
                )
            else:
                final_elements = (
                    self.elements_summary + self.elements + self.end_of_report_elements
                )

            # Build the PDF
            from .header_footer import header_footer_canvas_factory

            self.doc.multiBuild(
                final_elements,
                canvasmaker=header_footer_canvas_factory(
                    self.sample_id, self.centreID, self.styles, self.fonts_dir
                ),
            )

            logger.info(f"PDF created: {self.filename}")

            # Add success message to report
            success_content = f"""
            <b>✅ Report Generation Complete</b><br/>
            PDF file: {os.path.basename(self.filename)}<br/>
            Total sections processed: {len(self.sections)}<br/>
            Report type: {report_type.title()}
            """
            self.end_of_report_elements.append(
                Paragraph(success_content, self.styles.styles["Success"])
            )
            self.end_of_report_elements.append(Spacer(1, 12))

            # Optionally export CSV/XLSX/ZIP
            if export_csv_dir:
                try:
                    os.makedirs(export_csv_dir, exist_ok=True)
                    manifest = {
                        "sample_id": self.sample_id,
                        "centre_id": self.centreID,
                        "report_type": report_type,
                        "files": [],
                    }

                    # Collect frames from sections
                    for section in self.sections:
                        frames = getattr(section, "get_export_frames", lambda: {})()
                        for name, df in frames.items():
                            safe_name = name.replace(" ", "_")
                            csv_path = os.path.join(
                                export_csv_dir,
                                f"{self.sample_id}_{safe_name}.csv",
                            )
                            try:
                                df.to_csv(csv_path, index=False)
                                manifest["files"].append(
                                    {
                                        "name": name,
                                        "path": csv_path,
                                        "rows": int(df.shape[0]),
                                        "cols": int(df.shape[1]) if df.shape else 0,
                                    }
                                )
                            except Exception as ex:
                                logger.error(
                                    f"Error writing CSV for frame {name}: {str(ex)}",
                                    exc_info=True,
                                )

                    # Write manifest
                    manifest_path = os.path.join(
                        export_csv_dir, f"{self.sample_id}_manifest.json"
                    )
                    with open(manifest_path, "w", encoding="utf-8") as fh:
                        json.dump(manifest, fh, indent=2)

                    # Optional XLSX workbook
                    if export_xlsx:
                        try:
                            xlsx_path = os.path.join(
                                export_csv_dir, f"{self.sample_id}_report_data.xlsx"
                            )
                            with pd.ExcelWriter(xlsx_path) as writer:
                                for f in manifest["files"]:
                                    # Reload to avoid potential dtype issues
                                    try:
                                        df = (
                                            pd.read_csv(f["path"])
                                            if f["path"].endswith(".csv")
                                            else None
                                        )
                                    except Exception:
                                        df = None
                                    if df is not None:
                                        sheet_name = (
                                            os.path.basename(f["path"])
                                            .replace(f"{self.sample_id}_", "")
                                            .replace(".csv", "")[:31]
                                        )
                                        df.to_excel(
                                            writer, index=False, sheet_name=sheet_name
                                        )
                            logger.info(f"XLSX written: {xlsx_path}")
                        except Exception as ex:
                            logger.error(
                                "Error writing XLSX: %s", str(ex), exc_info=True
                            )

                    # Optional ZIP archive
                    if export_zip:
                        try:
                            import zipfile

                            zip_path = os.path.join(
                                export_csv_dir, f"{self.sample_id}_report_data.zip"
                            )
                            with zipfile.ZipFile(
                                zip_path, "w", zipfile.ZIP_DEFLATED
                            ) as zf:
                                for f in manifest["files"]:
                                    if os.path.exists(f["path"]):
                                        zf.write(
                                            f["path"],
                                            arcname=os.path.basename(f["path"]),
                                        )
                                if os.path.exists(manifest_path):
                                    zf.write(
                                        manifest_path,
                                        arcname=os.path.basename(manifest_path),
                                    )
                            logger.info(f"ZIP written: {zip_path}")
                        except Exception as ex:
                            logger.error(
                                "Error writing ZIP: %s", str(ex), exc_info=True
                            )
                except Exception as ex:
                    logger.error(
                        "Error exporting CSV/XLSX/ZIP artifacts: %s",
                        str(ex),
                        exc_info=True,
                    )

            return self.filename
        except Exception as e:
            logger.error(f"Error generating report: {e}", exc_info=True)
            raise


def create_pdf(
    filename,
    output,
    center: str,
    report_type="detailed",
    export_csv_dir=None,
    export_xlsx=False,
    export_zip=False,
):
    """Create a PDF report from ROBIN analysis results.

    Args:
        filename: Output PDF filename
        output: Directory containing analysis output files
        report_type: Type of report to generate ('summary' or 'detailed')

    Returns:
        Path to the generated PDF file
    """
    report = RobinReport(filename, output, center)
    return report.generate_report(
        report_type=report_type,
        export_csv_dir=export_csv_dir,
        export_xlsx=export_xlsx,
        export_zip=export_zip,
    )
