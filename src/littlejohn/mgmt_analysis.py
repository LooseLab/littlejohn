#!/usr/bin/env python3
"""
MGMT (O6-methylguanine-DNA methyltransferase) promoter methylation analysis.

This module provides automated processing of BAM files to analyze methylation
patterns in the MGMT promoter region (chr10:129466536-129467536). The analysis
integrates with LittleJohn's preprocessing pipeline and provides rich metadata
output.

Features:
- Automated BAM processing with LittleJohn's preprocessing pipeline
- MGMT region extraction using bedtools
- Methylation analysis using matkit from robin package
- R-based prediction model for methylation status
- Methylation visualization using methylartist
- Comprehensive metadata extraction and logging
"""

import os
import time
import subprocess
import tempfile
import shutil
import pysam
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional
from pathlib import Path
from littlejohn.logging_config import get_job_logger

# Try to find HV path from robin module, fallback to common locations
try:
    HVPATH = os.path.join(
        os.path.dirname(os.path.abspath(__import__('robin').__file__)), "submodules", "hv_rapidCNS2"
    )
except ImportError:
    # Fallback paths if robin module is not available
    HVPATH = None


def has_reads(bam_file, chrom, start, end):
    """
    Quickly checks if any reads span a specific genomic locus.

    Args:
        bam_file (str): Path to your BAM file.
        chrom (str): Chromosome name (e.g., 'chr10').
        start (int): Start coordinate (0-based, inclusive).
        end (int): End coordinate (0-based, exclusive).

    Returns:
        bool: True if at least one read spans the region, False otherwise.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for _ in bam.fetch(chrom, start, end):
            return True  # Found at least one read
    return False  # No reads found


@dataclass
class MGMTMetadata:
    """Container for MGMT analysis metadata and results"""
    sample_id: str
    bam_path: str
    analysis_timestamp: float
    processing_steps: List[str] = field(default_factory=list)
    error_message: Optional[str] = None
    results: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.results is None:
            self.results = {}


class MGMTAnalysis:
    """MGMT promoter methylation analysis worker"""
    
    def __init__(self, work_dir=None, threads=4):
        self.work_dir = work_dir or os.getcwd()
        self.threads = threads
        
        # Find required files and paths
        self.mgmt_bed_path = self._find_mgmt_bed()
        self.hv_path = self._find_hv_path()
        
        # Initialize counter for incremental file naming
        self.file_counter = 1
        
        # Use logger for initialization
        import logging
        logger = logging.getLogger("littlejohn.mgmt")
        logger.info(f"MGMT Analysis initialized")
        logger.debug(f"Work directory: {self.work_dir}")
        logger.debug(f"MGMT BED file: {self.mgmt_bed_path}")
        logger.debug(f"HV path: {self.hv_path}")
        logger.debug(f"Threads: {self.threads}")
        
    def _find_mgmt_bed(self) -> str:
        """Find or create the MGMT BED file"""
        possible_paths = [
            "mgmt_hg38.bed",
            "bin/mgmt_hg38.bed", 
            "data/mgmt_hg38.bed",
            "/usr/local/share/mgmt_hg38.bed"
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        
        # Create default MGMT BED file
        default_path = os.path.join(self.work_dir, "mgmt_hg38.bed")
        with open(default_path, 'w') as f:
            f.write("chr10\t129466536\t129467536\tMGMT_promoter\t1000\t+\n")
        
        return default_path

    def _find_hv_path(self) -> str:
        """Find the HV rapidCNS2 path from robin module"""
        # Use the already defined HVPATH constant if available
        if HVPATH and os.path.exists(HVPATH):
            return HVPATH
        
        # Fallback: create a basic structure if HVPATH doesn't exist
        hv_path = os.path.join(self.work_dir, "hv_rapidCNS2")
        os.makedirs(hv_path, exist_ok=True)
        os.makedirs(os.path.join(hv_path, "bin"), exist_ok=True)
        
        if HVPATH:
            print(f"HVPATH not found at {HVPATH}, created basic structure at: {hv_path}")
        else:
            print(f"Robin module not available, created basic HV structure at: {hv_path}")
        return hv_path

    def _get_next_file_number(self) -> int:
        """Get the next file number for incremental naming"""
        current_number = self.file_counter
        self.file_counter += 1
        return current_number

    def process_bam_file(self, bam_path: str, metadata: Dict[str, Any]) -> MGMTMetadata:
        """Process a single BAM file for MGMT analysis"""
        import logging
        logger = logging.getLogger("littlejohn.mgmt")
        
        logger.info(f"Processing BAM file: {bam_path} in MGMT analysis")
        sample_id = metadata.get('sample_id', 'unknown')
        start_time = time.time()
        
        # Get next file number for incremental naming
        file_number = self._get_next_file_number()
        
        mgmt_result = MGMTMetadata(
            sample_id=sample_id,
            bam_path=bam_path,
            analysis_timestamp=start_time
        )
        
        logger.info(f"Starting MGMT analysis for sample: {sample_id}")
        logger.debug(f"BAM file: {os.path.basename(bam_path)}")
        logger.debug(f"Work directory: {self.work_dir}")
        
        try:
            # Step 1: Check if BAM file exists and is readable
            if not os.path.exists(bam_path):
                raise FileNotFoundError(f"BAM file not found: {bam_path}")
            
            mgmt_result.processing_steps.append("file_validation")
            
            # Step 2: Check if BAM file contains MGMT spanning reads
            has_mgmt_reads = has_reads(bam_path, "chr10", 129467242, 129467244)
            
            if not has_mgmt_reads:
                mgmt_result.processing_steps.append("no_mgmt_reads")
                mgmt_result.error_message = "No MGMT spanning reads found"
                return mgmt_result
            
            mgmt_result.processing_steps.append("mgmt_reads_found")
            
            # Step 3: Extract MGMT region using bedtools
            MGMT_BED: str = self.mgmt_bed_path
            
            # Create sample-specific directory
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            
            # Create temporary BAM file for this extraction
            temp_bamfile = tempfile.NamedTemporaryFile(
                suffix=".bam",
                delete=False
            )
            temp_bamfile.close()
            
            # Extract MGMT region to temporary file
            bedtools_cmd = [
                "bedtools", "intersect",
                "-a", bam_path,
                "-b", MGMT_BED
            ]
            
            try:
                with open(temp_bamfile.name, 'wb') as f:
                    result = subprocess.run(
                        bedtools_cmd,
                        stdout=f,
                        stderr=subprocess.PIPE,
                        timeout=300  # 5 minute timeout
                    )
                
                if result.returncode != 0:
                    # Check if the output file was created and has content
                    if not os.path.exists(temp_bamfile.name) or os.path.getsize(temp_bamfile.name) == 0:
                        raise RuntimeError("Bedtools failed to create valid BAM file")
                
                # Create index for the BAM file
                try:
                    pysam.index(temp_bamfile.name, f"{temp_bamfile.name}.bai")
                except Exception as e:
                    print(f"Failed to create BAM index: {e}")
                    # Continue without index if it fails
                
                # Check if the extracted BAM has reads
                if pysam.AlignmentFile(temp_bamfile.name, "rb").count(until_eof=True) > 0:
                    # Accumulate reads in the sample's mgmt.bam file
                    mgmt_bam_output = os.path.join(sample_dir, "mgmt.bam")
                    
                    if os.path.exists(mgmt_bam_output):
                        # Concatenate with existing mgmt.bam
                        temp_holder = tempfile.NamedTemporaryFile(
                            suffix=".bam",
                            delete=False
                        )
                        temp_holder.close()
                        
                        pysam.cat("-o", temp_holder.name, mgmt_bam_output, temp_bamfile.name)
                        shutil.copy2(temp_holder.name, mgmt_bam_output)
                        try:
                            os.remove(f"{temp_holder.name}.bai")
                            os.remove(f"{temp_bamfile.name}.bai")
                        except FileNotFoundError:
                            pass
                    else:
                        # Copy to create new mgmt.bam
                        shutil.copy2(temp_bamfile.name, mgmt_bam_output)
                        try:
                            os.remove(f"{temp_bamfile.name}.bai")
                        except FileNotFoundError:
                            pass
                    
                    mgmt_result.processing_steps.append("mgmt_region_extraction")
                else:
                    mgmt_result.processing_steps.append("no_mgmt_region_reads")
                    #mgmt_result.error_message = "No reads found in extracted MGMT region"
                    return mgmt_result
                
            except subprocess.TimeoutExpired:
                raise RuntimeError("Bedtools extraction timed out")
            except FileNotFoundError:
                raise RuntimeError("bedtools not found in PATH")
            finally:
                # Clean up temporary files
                temp_files_to_clean = [temp_bamfile.name]
                if 'temp_holder' in locals():
                    temp_files_to_clean.append(temp_holder.name)
                
                for temp_file in temp_files_to_clean:
                    try:
                        if os.path.exists(temp_file):
                            os.remove(temp_file)
                    except FileNotFoundError:
                        pass
            
            # Step 4: Run matkit analysis (if available)
            try:
                print("Running matkit")
                # Use the proper run_matkit function with robin package
                if self.run_matkit(self.work_dir, mgmt_bam_output, self.hv_path, sample_dir, file_number):
                    mgmt_result.processing_steps.append("matkit_analysis")
                    print(f"Matkit analysis complete")
                else:
                    print("Matkit analysis failed or robin package not available")
                    
            except Exception as e:
                print(f"Matkit analysis error: {e}")
            
            # Step 5: Generate visualization (if methylartist is available)
            temp_files = []  # Track temporary files for cleanup
            
            try:
                # Sort and index the accumulated mgmt.bam file for methylartist
                sorted_mgmt_bam = os.path.join(sample_dir, "mgmt_sorted.bam")
                temp_files.append(sorted_mgmt_bam)
                temp_files.append(f"{sorted_mgmt_bam}.bai")
                
                try:
                    pysam.sort("-o", sorted_mgmt_bam, mgmt_bam_output)
                    pysam.index(sorted_mgmt_bam, f"{sorted_mgmt_bam}.bai")
                    print(f"Sorted and indexed accumulated MGMT BAM: {sorted_mgmt_bam}")
                    mgmt_bam_for_plot = sorted_mgmt_bam
                except Exception as e:
                    print(f"Failed to sort/index accumulated MGMT BAM: {e}")
                    mgmt_bam_for_plot = mgmt_bam_output
                
                plot_out = os.path.join(sample_dir, f"{file_number}_mgmt.png")
                methylartist_cmd = f"methylartist locus -i chr10:129466536-129467536 -b {mgmt_bam_for_plot} -o {plot_out} --motif CG --mods m"
                
                result = subprocess.run(
                    methylartist_cmd,
                    shell=True,
                    capture_output=True,
                    text=True,
                    timeout=300  # 5 minute timeout
                )
                
                if result.returncode == 0:
                    mgmt_result.processing_steps.append("visualization")
                    print(f"Visualization complete")
                else:
                    print(f"Methylartist visualization failed: {result.stderr}")
                    
            except (subprocess.TimeoutExpired, FileNotFoundError):
                print("Methylartist not available, skipping visualization")
            finally:
                # Clean up temporary files
                for temp_file in temp_files:
                    try:
                        if os.path.exists(temp_file):
                            os.remove(temp_file)
                    except Exception as e:
                        print(f"Could not remove temporary file {temp_file}: {e}")
            
            # Step 6: Collect results and metadata
            analysis_time = time.time() - start_time
            mgmt_result.results = {
                "analysis_time": analysis_time,
                "mgmt_bam_file": mgmt_bam_output,
                "processing_steps": mgmt_result.processing_steps.copy(),
                "tools_available": {
                    "bedtools": shutil.which("bedtools") is not None,
                    "robin_package": HVPATH is not None,
                    "matkit": HVPATH is not None,  # matkit comes from robin package
                    "r": shutil.which("Rscript") is not None,
                    "methylartist": shutil.which("methylartist") is not None
                },
                "hv_path": self.hv_path,
                "hv_path_source": "robin_module" if HVPATH and self.hv_path == HVPATH else "fallback"
            }
            
            mgmt_result.processing_steps.append("analysis_complete")
            print(f"MGMT analysis completed for {sample_id} in {analysis_time:.2f}s")
            
        except Exception as e:
            error_msg = f"MGMT analysis failed for {sample_id}: {str(e)}"
            print(f"{error_msg}")
            mgmt_result.error_message = error_msg
            mgmt_result.processing_steps.append("analysis_failed")
        
        return mgmt_result
    
    def convert_to_mixed_delim(self, input_file: str, output_file: str) -> bool:
        """
        Converts standard BEDMethyl format to mixed-delim format for R script compatibility.

        The R script expects column V10 to contain a space-separated string like "10 0.5"
        instead of separate columns for coverage and fraction.
        """
        try:
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in infile:
                    if line.strip():
                        fields = line.strip().split('\t')
                        if len(fields) >= 18:
                            # Standard BEDMethyl format has:
                            # V10: coverage (column 10, 0-indexed)
                            # V11: fraction modified (column 11, 0-indexed)

                            # Create mixed-delim format by combining V10 and V11 into V10
                            coverage = fields[9]  # V10 (0-indexed)
                            fraction = fields[10]  # V11 (0-indexed)

                            # Create the mixed-delim format: "coverage fraction"
                            mixed_delim = f"{coverage} {fraction}"

                            # Replace V10 with the mixed-delim string
                            fields[9] = mixed_delim

                            # Write the converted line
                            outfile.write('\t'.join(fields) + '\n')
                        else:
                            # If the line doesn't have enough columns, write it as-is
                            outfile.write(line)
            
            print(f"Converted BEDMethyl format: {input_file} → {output_file}")
            return True
            
        except Exception as e:
            print(f"Error converting BEDMethyl format: {e}")
            return False

    def run_matkit(self, temp_dir: str, mgmt_bamfile: str, hv_path: str, sample_dir: str, file_number: int) -> bool:
        """Processes BAM file with matkit and runs R script for MGMT prediction."""
        temp_files = []  # Track temporary files for cleanup
        
        try:
            print(f"Running matkit for methylation analysis")
            
            # Sort and index the BAM file
            sorted_bam = os.path.join(temp_dir, "mgmt_sorted.bam")
            temp_files.append(sorted_bam)
            temp_files.append(f"{sorted_bam}.bai")
            
            pysam.sort("-o", sorted_bam, mgmt_bamfile)
            pysam.index(sorted_bam, f"{sorted_bam}.bai")
            
            # Use matkit from robin package
            try:
                from robin.utils import run_matkit as run_matkit_util
                temp_bed_file = os.path.join(temp_dir, "temp_mgmt.bed")
                temp_files.append(temp_bed_file)
                run_matkit_util(sorted_bam, temp_bed_file)
                
                # Convert to mixed delimiter format for R script
                output_mgmt_bed = os.path.join(temp_dir, "mgmt_mixed.bed")
                temp_files.append(output_mgmt_bed)
                self.convert_to_mixed_delim(temp_bed_file, output_mgmt_bed)
                
                # Copy BED file to sample directory with incremental naming
                final_bed_file = os.path.join(sample_dir, f"{file_number}_mgmt.bed")
                shutil.copy2(output_mgmt_bed, final_bed_file)
                
                # Run R script for MGMT prediction
                r_script_path = os.path.join(hv_path, "bin", "mgmt_pred_v0.3.R")
                if os.path.exists(r_script_path):
                    cmd = f"Rscript {r_script_path} --input={output_mgmt_bed} --out_dir={sample_dir} --probes={hv_path}/bin/mgmt_probes.Rdata --model={hv_path}/bin/mgmt_137sites_mean_model.Rdata --sample=live_analysis"
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                    
                    if result.returncode == 0:
                        # Rename the output CSV file to use incremental naming
                        default_csv = os.path.join(sample_dir, "live_analysis_mgmt_status.csv")
                        final_csv = os.path.join(sample_dir, f"{file_number}_mgmt.csv")
                        if os.path.exists(default_csv):
                            shutil.move(default_csv, final_csv)
                    else:
                        print(f"R script failed: {result.stderr}")
                        # Continue anyway as we still have the BED file
                else:
                    print(f"R script not found at {r_script_path}")
                
                print(f"Matkit processing completed")
                return True
                
            except ImportError:
                print("Robin package not available, skipping matkit processing")
                return False
                
        except Exception as e:
            print(f"Error running matkit: {e}")
            return False
        finally:
            # Clean up temporary files
            for temp_file in temp_files:
                try:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
                except Exception as e:
                    print(f"Could not remove temporary file {temp_file}: {e}")


def mgmt_handler(job, work_dir=None):
    """
    Handler function for MGMT analysis jobs.
    This function processes BAM files for MGMT methylation analysis.
    
    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to BAM file directory)
    """
    try:
        bam_path = job.context.filepath
        
        # Get job-specific logger
        logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
        logger.info(f"Starting MGMT analysis for: {os.path.basename(bam_path)}")
        
        # Get metadata from preprocessing
        bam_metadata = job.context.metadata.get('bam_metadata', {})
        
        # Determine work directory
        if work_dir is None:
            # Default to BAM file directory
            work_dir = os.path.dirname(bam_path)
        else:
            # Use specified work directory, create if it doesn't exist
            os.makedirs(work_dir, exist_ok=True)
            logger.debug(f"Using specified work directory: {work_dir}")
        
        # Create MGMT analysis instance
        mgmt_analyzer = MGMTAnalysis(work_dir=work_dir)
        
        # Process the BAM file
        mgmt_result = mgmt_analyzer.process_bam_file(bam_path, bam_metadata)
        
        # Store results in job context
        job.context.add_metadata('mgmt_analysis', mgmt_result.results)
        job.context.add_metadata('mgmt_processing_steps', mgmt_result.processing_steps)
        
        if mgmt_result.error_message:
            job.context.add_error('mgmt_analysis', mgmt_result.error_message)
            logger.error(f"MGMT analysis failed: {mgmt_result.error_message}")
        else:
            job.context.add_result('mgmt_analysis', {
                'status': 'success',
                'sample_id': mgmt_result.sample_id,
                'analysis_time': mgmt_result.results.get('analysis_time', 0),
                'mgmt_bam_file': mgmt_result.results.get('mgmt_bam_file', ''),
                'processing_steps': mgmt_result.processing_steps,
                'tools_available': mgmt_result.results.get('tools_available', {})
            })
            logger.info(f"MGMT analysis complete for {os.path.basename(bam_path)}")
            logger.info(f"Sample ID: {mgmt_result.sample_id}")
            logger.info(f"Analysis time: {mgmt_result.results.get('analysis_time', 0):.2f}s")
            logger.debug(f"Processing steps: {', '.join(mgmt_result.processing_steps)}")
            logger.debug(f"Output directory: {os.path.dirname(mgmt_result.results.get('mgmt_bam_file', ''))}")
            
    except Exception as e:
        job.context.add_error('mgmt_analysis', str(e))
        logger.error(f"Error in MGMT analysis for {job.context.filepath}: {e}") 