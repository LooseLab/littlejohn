#!/usr/bin/env python3
"""
Fusion Analysis Module for LittleJohn

This module provides automated fusion analysis following the exact architecture
specified in the Fusion Analysis Documentation. It integrates with LittleJohn's 
workflow system and processes files for fusion detection.

Features:
- Automated fusion detection using BAM files
- Target panel and genome-wide fusion analysis
- Integration with LittleJohn's workflow system
- Sample-specific output directories
- Comprehensive metadata extraction and logging
- Error handling and result tracking
- State persistence and incremental processing
- Memory-optimized processing

Classes
-------
FusionMetadata
    Container for fusion analysis metadata and results.

FusionAnalysis
    Main analysis class that processes files for fusion detection.

Dependencies
-----------
- pandas: Data manipulation and analysis
- numpy: Numerical computations
- logging: Logging for debugging and monitoring
- typing: Type hints
- tempfile: Temporary file creation
- pathlib: File system paths
- os: Operating system interface
- time: Time utilities
- json: JSON serialization
- pickle: Python object serialization
- gc: Garbage collection
- pysam: BAM file processing
- asyncio: Asynchronous processing support
- subprocess: External command execution

Example Usage
-----------
.. code-block:: python

    from littlejohn.fusion_analysis import FusionAnalysis

    # Initialize analysis
    fusion_analysis = FusionAnalysis(
        work_dir="output/",
        config_path="fusion_config.json"
    )

    # Process files
    fusion_analysis.process_file("sample.bam")

Notes
-----
The module follows the LittleJohn framework patterns for:
- Integration with workflow system
- Worker process management
- State tracking and persistence
- Error handling and logging
- Output generation and file management

Authors
-------
Matt Loose
"""

import os
import sys
import tempfile
import logging
import time
import pickle
import gc
import json
import asyncio
import subprocess
import shutil
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple, Set
from dataclasses import dataclass, field
from collections import Counter, defaultdict
from io import StringIO
import numpy as np
import pandas as pd
import pysam
from littlejohn.logging_config import get_job_logger


def json_serializable(obj):
    """Convert numpy types to JSON serializable types."""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, pd.Series):
        return obj.tolist()
    elif isinstance(obj, pd.DataFrame):
        return obj.to_dict('records')
    else:
        return obj


@dataclass
class GeneRegion:
    """Represents a gene region with start, end, and name."""
    start: int
    end: int
    name: str

    def overlaps_with(self, other_start: int, other_end: int, min_overlap: int = 100) -> bool:
        """Check if this region overlaps with another region."""
        overlap_start = max(self.start, other_start)
        overlap_end = min(self.end, other_end)
        return overlap_end > overlap_start and (overlap_end - overlap_start) > min_overlap


@dataclass
class FusionMetadata:
    """Container for fusion analysis metadata and results"""
    sample_id: str
    file_path: str
    analysis_timestamp: float
    target_fusion_path: Optional[str] = None
    genome_wide_fusion_path: Optional[str] = None
    analysis_results: Optional[Dict] = None
    processing_steps: List[str] = None
    error_message: Optional[str] = None
    fusion_data: Optional[Dict] = None
    target_panel: Optional[str] = None
    
    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.analysis_results is None:
            self.analysis_results = {}
        if self.fusion_data is None:
            self.fusion_data = {}


class FusionProcessor:
    """
    Core processor for fusion detection from BAM files.
    
    Responsibilities:
    - Load and cache gene region data
    - Process BAM files for fusion candidates
    - Memory-optimized processing
    """

    def __init__(self, target_panel: str = "rCNS2"):
        self.target_panel = target_panel
        self.gene_regions_cache: Dict[str, List[GeneRegion]] = {}
        self.all_gene_regions_cache: Dict[str, List[GeneRegion]] = {}
        
        # Initialize file paths
        self._setup_file_paths()
        
        # Load gene regions into memory
        self._load_gene_regions()

    def _setup_file_paths(self) -> None:
        """Setup file paths for gene data."""
        # Try to import robin resources to get the correct path
        try:
            from robin import resources
            resources_dir = os.path.dirname(resources.__file__)
            
            # Look for BED files in robin resources first
            if self.target_panel == "rCNS2":
                rCNS2_path = os.path.join(resources_dir, "rCNS2_panel_name_uniq.bed")
                if os.path.exists(rCNS2_path):
                    self.gene_bed = rCNS2_path
                else:
                    self.gene_bed = "rCNS2_panel_name_uniq.bed"
            elif self.target_panel == "AML":
                aml_path = os.path.join(resources_dir, "AML_panel_name_uniq.bed")
                if os.path.exists(aml_path):
                    self.gene_bed = aml_path
                else:
                    self.gene_bed = "AML_panel_name_uniq.bed"
            else:
                # Default to rCNS2
                self.target_panel = "rCNS2"
                rCNS2_path = os.path.join(resources_dir, "rCNS2_panel_name_uniq.bed")
                if os.path.exists(rCNS2_path):
                    self.gene_bed = rCNS2_path
                else:
                    self.gene_bed = "rCNS2_panel_name_uniq.bed"
            
            # Look for all genes BED file
            all_genes_path = os.path.join(resources_dir, "all_genes2.bed")
            if os.path.exists(all_genes_path):
                self.all_gene_bed = all_genes_path
            else:
                self.all_gene_bed = "all_genes2.bed"
                
        except ImportError:
            # Fallback to local files if robin is not available
            if self.target_panel == "rCNS2":
                self.gene_bed = "rCNS2_panel_name_uniq.bed"
            elif self.target_panel == "AML":
                self.gene_bed = "AML_panel_name_uniq.bed"
            else:
                # Default to rCNS2
                self.target_panel = "rCNS2"
                self.gene_bed = "rCNS2_panel_name_uniq.bed"
            
            self.all_gene_bed = "all_genes2.bed"
        
        # Log the file paths for debugging
        logging.info(f"FusionProcessor initialized with target_panel: {self.target_panel}")
        logging.info(f"Gene bed path: {self.gene_bed}")
        logging.info(f"All gene bed path: {self.all_gene_bed}")

    def _load_gene_regions(self) -> None:
        """Load gene regions into memory for efficient lookup."""
        # Load target panel gene regions
        self.gene_regions_cache = self._load_bed_regions(self.gene_bed)
        
        # Load genome-wide gene regions
        self.all_gene_regions_cache = self._load_bed_regions(self.all_gene_bed)
        
        logging.info(f"Gene regions loaded - target cache: {len(self.gene_regions_cache)} regions")
        logging.info(f"Gene regions loaded - genome-wide cache: {len(self.all_gene_regions_cache)} regions")

    def _load_bed_regions(self, bed_file: str) -> Dict[str, List[GeneRegion]]:
        """
        Load BED file regions into memory for efficient lookup.
        
        Args:
            bed_file: Path to BED file
            
        Returns:
            Dictionary mapping chromosome to list of GeneRegion objects
        """
        regions = defaultdict(list)
        
        try:
            with open(bed_file, "r") as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split("\t")
                    if len(parts) < 4:
                        continue
                    
                    try:
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        gene_name = parts[3]
                        
                        regions[chrom].append(GeneRegion(start, end, gene_name))
                    except ValueError:
                        continue
                        
        except FileNotFoundError:
            # If file not found, return empty dict
            pass
        except Exception:
            # If any other error, return empty dict
            pass
            
        return dict(regions)

    def has_supplementary_alignments(self, bamfile: str) -> bool:
        """
        Quickly check if a BAM file has supplementary alignments.
        
        Args:
            bamfile: Path to the BAM file
            
        Returns:
            True if supplementary alignment is found, else False
        """
        try:
            with pysam.AlignmentFile(bamfile, "rb") as bam:
                read_count = 0
                for read in bam:
                    read_count += 1
                    if read.has_tag('SA'):
                        return True
                    if read.is_supplementary:
                        return True
                    # Limit the number of reads to check for performance
                    if read_count > 10000:
                        break
            return False
        except Exception as e:
            logging.error(f"Error checking supplementary alignments in {bamfile}: {str(e)}")
            return False

    def find_reads_with_supplementary(self, bamfile: str) -> Set[str]:
        """
        Find all reads that have supplementary alignments.
        
        Args:
            bamfile: Path to the BAM file
            
        Returns:
            Set of read names with supplementary alignments
        """
        reads_with_supp = set()
        
        try:
            with pysam.AlignmentFile(bamfile, "rb") as bam:
                for read in bam:
                    if read.has_tag('SA') or read.is_supplementary:
                        reads_with_supp.add(read.query_name)
        except Exception:
            pass
            
        return reads_with_supp

    def process_bam_for_fusions(self, bamfile: str) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
        """
        Process BAM file to find fusion candidates.
        
        Args:
            bamfile: Path to BAM file
            
        Returns:
            Tuple of (target_panel_candidates, genome_wide_candidates) DataFrames
        """
        try:
            # Find reads with supplementary alignments
            reads_with_supp = self.find_reads_with_supplementary(bamfile)
            
            if not reads_with_supp:
                return None, None
            
            # Process reads for target panel fusions
            target_candidates = self._process_reads_for_fusions(
                bamfile, reads_with_supp, self.gene_regions_cache
            )
            
            # Process reads for genome-wide fusions
            genome_wide_candidates = self._process_reads_for_fusions(
                bamfile, reads_with_supp, self.all_gene_regions_cache
            )
            
            return target_candidates, genome_wide_candidates
            
        except Exception as e:
            import traceback
            logging.error(f"Error processing BAM file for fusions: {str(e)}")
            logging.error(f"Traceback: {traceback.format_exc()}")
            return None, None

    def _process_reads_for_fusions(
        self,
        bamfile: str,
        reads_with_supp: Set[str],
        gene_regions: Dict[str, List[GeneRegion]],
    ) -> Optional[pd.DataFrame]:
        """
        Process reads to find gene intersections and create fusion candidates.
        
        Args:
            bamfile: Path to BAM file
            reads_with_supp: Set of read names with supplementary alignments
            gene_regions: Dictionary of gene regions by chromosome
            
        Returns:
            DataFrame with fusion candidates or None if no candidates found
        """
        rows = []
        
        try:
            with pysam.AlignmentFile(bamfile, "rb") as bam:
                for read in bam:
                    # Skip if not in our target reads
                    if read.query_name not in reads_with_supp:
                        continue
                    
                    # Skip secondary alignments and unmapped reads
                    if read.is_secondary or read.is_unmapped:
                        continue
                    
                    # Get reference information
                    ref_name = bam.get_reference_name(read.reference_id) if read.reference_id >= 0 else None
                    if not ref_name or ref_name == "chrM":
                        continue
                    
                    ref_start = read.reference_start
                    ref_end = read.reference_end
                    
                    # Check if this read intersects with any gene regions
                    if ref_name in gene_regions:
                        read_rows = self._find_gene_intersections(
                            read, ref_name, ref_start, ref_end, gene_regions[ref_name]
                        )
                        rows.extend(read_rows)
                        
        except Exception as e:
            logging.error(f"Error processing reads for fusions: {str(e)}")
            return None
            
        if not rows:
            return None
            
        # Create DataFrame
        df = pd.DataFrame(rows)
        return self._optimize_fusion_dataframe(df)

    def _find_gene_intersections(
        self,
        read: pysam.AlignedSegment,
        ref_name: str,
        ref_start: int,
        ref_end: int,
        gene_regions: List[GeneRegion],
    ) -> List[Dict]:
        """
        Find gene intersections for a read.
        
        Args:
            read: Pysam AlignedSegment object
            ref_name: Reference chromosome name
            ref_start: Reference start position
            ref_end: Reference end position
            gene_regions: List of gene regions for this chromosome
            
        Returns:
            List of dictionaries with fusion candidate data
        """
        intersections = []
        
        for gene_region in gene_regions:
            if gene_region.overlaps_with(ref_start, ref_end):
                intersection = {
                    'read_id': read.query_name,
                    'gene_name': gene_region.name,
                    'chromosome': ref_name,
                    'start': ref_start,
                    'end': ref_end,
                    'strand': '-' if read.is_reverse else '+',
                    'mapping_quality': read.mapping_quality,
                    'read_start': read.query_alignment_start,
                    'read_end': read.query_alignment_end,
                    'mapping_span': ref_end - ref_start,
                }
                intersections.append(intersection)
                
        return intersections

    def _optimize_fusion_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Optimize DataFrame memory usage.
        
        Args:
            df: Input DataFrame
            
        Returns:
            Optimized DataFrame
        """
        if df.empty:
            return df
            
        # Optimize dtypes
        for col in df.columns:
            if df[col].dtype == 'object':
                if col in ['read_id', 'gene_name', 'chromosome', 'strand']:
                    df[col] = df[col].astype('category')
                    
        return df


class FusionAnalysis:
    """
    Main fusion analysis class that processes files for fusion detection.
    
    This class follows the LittleJohn framework patterns for:
    - Integration with workflow system
    - Worker process management
    - State tracking and persistence
    - Error handling and logging
    - Output generation and file management
    """

    def __init__(self, work_dir=None, config_path=None, threads=4):
        """
        Initialize the fusion analysis.
        
        Args:
            work_dir: Working directory for output files
            config_path: Path to configuration file
            threads: Number of threads to use
        """
        self.work_dir = work_dir or "fusion_output"
        self.config_path = config_path
        self.threads = threads
        self.logger = logging.getLogger("littlejohn.fusion")
        
        # Ensure work directory exists
        os.makedirs(self.work_dir, exist_ok=True)
        
        # Load configuration
        self.config = self._load_config()
        
        # Initialize processor
        self.target_panel = self.config.get('target_panel', 'rCNS2')
        self.processor = FusionProcessor(self.target_panel)

    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from file or use defaults."""
        if self.config_path and os.path.exists(self.config_path):
            try:
                with open(self.config_path, 'r') as f:
                    return json.load(f)
            except Exception as e:
                self.logger.warning(f"Could not load config file: {e}")
                
        # Default configuration
        return {
            'target_panel': 'rCNS2',
            'min_overlap': 100,
            'min_mapping_quality': 10,
        }

    def _get_next_file_number(self) -> int:
        """Get the next file number for output files."""
        return int(time.time() * 1000)

    def _check_and_create_folder(self, base_dir: str, sample_id: str) -> str:
        """Check and create folder for sample output."""
        sample_dir = os.path.join(base_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        return sample_dir

    def process_file(self, file_path: str, metadata: Dict[str, Any], timestamp: Optional[float] = None) -> FusionMetadata:
        """
        Process a single file for fusion analysis.
        
        Args:
            file_path: Path to the file to process
            metadata: Metadata about the file
            timestamp: Optional timestamp for the analysis
            
        Returns:
            FusionMetadata object with results
        """
        if timestamp is None:
            timestamp = time.time()
            
        # Extract sample ID from metadata
        sample_id = metadata.get('sample_id', 'unknown')
        
        # Create metadata object
        fusion_metadata = FusionMetadata(
            sample_id=sample_id,
            file_path=file_path,
            analysis_timestamp=timestamp,
            target_panel=self.target_panel
        )
        
        try:
            self.logger.info(f"Starting fusion analysis for {file_path}")
            fusion_metadata.processing_steps.append("started")
            
            # Check if file exists
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File not found: {file_path}")
                
            # Check if file is readable
            if not os.access(file_path, os.R_OK):
                raise PermissionError(f"File not readable: {file_path}")
                
            # Check if file is a BAM file
            if not file_path.lower().endswith('.bam'):
                self.logger.warning(f"File does not have .bam extension: {file_path}")
                
            # Try to open BAM file to check if it's valid
            try:
                with pysam.AlignmentFile(file_path, "rb") as bam:
                    # Just check if we can read the header
                    header = bam.header
                    self.logger.info(f"BAM file header read successfully, {len(header.get('SQ', []))} sequences")
            except Exception as bam_error:
                raise ValueError(f"Invalid or corrupted BAM file {file_path}: {str(bam_error)}")
                
            # Create temporary directory
            with tempfile.TemporaryDirectory() as temp_dir:
                fusion_metadata.processing_steps.append("temp_dir_created")
                
                # Perform fusion analysis
                analysis_results = self._perform_fusion_analysis(
                    file_path, temp_dir, metadata, fusion_metadata
                )
                
                # Generate output files
                output_files = self._generate_output_files(
                    sample_id, analysis_results, fusion_metadata
                )
                
                # Update metadata with results
                fusion_metadata.analysis_results = analysis_results
                fusion_metadata.target_fusion_path = output_files.get('target_fusion')
                fusion_metadata.genome_wide_fusion_path = output_files.get('genome_wide_fusion')
                fusion_metadata.processing_steps.append("completed")
                
                self.logger.info(f"Fusion analysis completed for {file_path}")
                
        except Exception as e:
            import traceback
            error_msg = f"Error in fusion analysis: {str(e)}"
            self.logger.error(error_msg)
            self.logger.error(f"Traceback: {traceback.format_exc()}")
            fusion_metadata.error_message = error_msg
            fusion_metadata.processing_steps.append("failed")
            
        return fusion_metadata

    def _perform_fusion_analysis(
        self, file_path: str, temp_dir: str, metadata: Dict[str, Any], fusion_metadata: FusionMetadata
    ) -> Dict[str, Any]:
        """
        Perform the actual fusion analysis.
        
        Args:
            file_path: Path to the file to process
            temp_dir: Temporary directory for intermediate files
            metadata: File metadata
            fusion_metadata: Fusion metadata object
            
        Returns:
            Dictionary with analysis results
        """
        fusion_metadata.processing_steps.append("analysis_started")
        
        # Check for supplementary alignments
        has_supp = self.processor.has_supplementary_alignments(file_path)
        
        if not has_supp:
            self.logger.info("No supplementary alignments found, skipping fusion analysis")
            return {
                'has_supplementary': False,
                'target_candidates_count': 0,
                'genome_wide_candidates_count': 0,
                'target_candidates': None,
                'genome_wide_candidates': None,
            }
        
        fusion_metadata.processing_steps.append("supplementary_found")
        
        # Process for fusions
        target_candidates, genome_wide_candidates = self.processor.process_bam_for_fusions(file_path)
        
        fusion_metadata.processing_steps.append("fusion_processing_complete")
        
        # Prepare results
        results = {
            'has_supplementary': True,
            'target_candidates_count': len(target_candidates) if target_candidates is not None else 0,
            'genome_wide_candidates_count': len(genome_wide_candidates) if genome_wide_candidates is not None else 0,
            'target_candidates': target_candidates,
            'genome_wide_candidates': genome_wide_candidates,
        }
        
        # Store fusion data in metadata
        fusion_metadata.fusion_data = {
            'target_candidates': target_candidates.to_dict('records') if target_candidates is not None else [],
            'genome_wide_candidates': genome_wide_candidates.to_dict('records') if genome_wide_candidates is not None else [],
        }
        
        return results

    def _generate_output_files(
        self, sample_id: str, analysis_results: Dict[str, Any], fusion_metadata: FusionMetadata
    ) -> Dict[str, str]:
        """
        Generate output files for the analysis.
        
        Args:
            sample_id: Sample ID
            analysis_results: Analysis results
            fusion_metadata: Fusion metadata object
            
        Returns:
            Dictionary mapping file type to file path
        """
        # Ensure work directory exists
        os.makedirs(self.work_dir, exist_ok=True)
        
        # Create sample-specific directory
        sample_dir = self._check_and_create_folder(self.work_dir, sample_id)
        
        output_files = {}
        
        # Save target panel fusion candidates (master) to sample directory
        if analysis_results.get('target_candidates') is not None:
            target_path = os.path.join(sample_dir, "fusion_candidates_master.csv")
            analysis_results['target_candidates'].to_csv(target_path, index=False)
            output_files['target_fusion'] = target_path
            fusion_metadata.target_fusion_path = target_path
            
        # Save genome-wide fusion candidates (all) to sample directory
        if analysis_results.get('genome_wide_candidates') is not None:
            genome_path = os.path.join(sample_dir, "fusion_candidates_all.csv")
            analysis_results['genome_wide_candidates'].to_csv(genome_path, index=False)
            output_files['genome_wide_fusion'] = genome_path
            fusion_metadata.genome_wide_fusion_path = genome_path
            
        # Save metadata to sample-specific directory
        metadata_path = os.path.join(sample_dir, f"{sample_id}_fusion_metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(fusion_metadata.__dict__, f, default=json_serializable, indent=2)
            
        return output_files


def process_single_file(file_path: str, metadata: Dict[str, Any], work_dir: str, logger=None) -> Dict[str, Any]:
    """
    Process a single file for fusion analysis (standalone function).
    
    Args:
        file_path: Path to the file to process
        metadata: Metadata about the file
        work_dir: Working directory for output
        logger: Logger instance
        
    Returns:
        Dictionary with processing results
    """
    try:
        if logger:
            logger.info(f"Initializing FusionAnalysis with work_dir: {work_dir}")
        
        # Initialize analysis
        fusion_analysis = FusionAnalysis(work_dir=work_dir)
        
        if logger:
            logger.info(f"Processing file: {file_path}")
            logger.info(f"File metadata: {metadata}")
        
        # Process file
        fusion_metadata = fusion_analysis.process_file(file_path, metadata)
        
        # Convert to dictionary
        result = {
            'success': fusion_metadata.error_message is None,
            'sample_id': fusion_metadata.sample_id,
            'file_path': fusion_metadata.file_path,
            'analysis_timestamp': fusion_metadata.analysis_timestamp,
            'target_fusion_path': fusion_metadata.target_fusion_path,
            'genome_wide_fusion_path': fusion_metadata.genome_wide_fusion_path,
            'processing_steps': fusion_metadata.processing_steps,
            'error_message': fusion_metadata.error_message,
            'analysis_results': fusion_metadata.analysis_results,
        }
        
        if fusion_metadata.error_message:
            if logger:
                logger.error(f"Fusion analysis failed: {fusion_metadata.error_message}")
                logger.error(f"Processing steps completed: {fusion_metadata.processing_steps}")
        else:
            if logger:
                logger.info(f"Fusion analysis completed successfully for {file_path}")
                logger.info(f"Processing steps: {fusion_metadata.processing_steps}")
            
        return result
        
    except Exception as e:
        import traceback
        error_msg = f"Error in fusion analysis: {str(e)}"
        if logger:
            logger.error(error_msg)
            logger.error(f"Traceback: {traceback.format_exc()}")
        return {
            'success': False,
            'error_message': error_msg,
            'file_path': file_path,
        }


def fusion_handler(job, work_dir=None):
    """
    Handler function for fusion analysis jobs in the workflow system.
    
    Args:
        job: Job object from the workflow system
        work_dir: Working directory for output
    """
    try:
        # Get logger with proper parameters
        logger = get_job_logger(str(job.job_id), "fusion", job.context.filepath)
        
        # Extract file path and metadata from job
        file_path = job.context.filepath
        metadata = job.context.metadata.get('bam_metadata', {})
        
        logger.info(f"Starting fusion analysis for {file_path}")
        logger.info(f"Metadata: {metadata}")
        
        # Set default work directory if not provided
        if work_dir is None:
            work_dir = "fusion_output"
            
        logger.info(f"Using work directory: {work_dir}")
        
        # Process the file
        result = process_single_file(file_path, metadata, work_dir, logger)
        
        # Add result to job context
        job.context.add_result("fusion_analysis", result)
        
        if result['success']:
            logger.info(f"Fusion analysis completed successfully for {file_path}")
            logger.info(f"Results: {result}")
        else:
            error_msg = result.get('error_message', 'Unknown error')
            logger.error(f"Fusion analysis failed for {file_path}: {error_msg}")
            logger.error(f"Full result: {result}")
            job.context.add_error("fusion_analysis", error_msg)
            
    except Exception as e:
        import traceback
        error_msg = f"Error in fusion handler: {str(e)}"
        
        # Try to get a logger, but don't fail if we can't
        try:
            logger = get_job_logger(str(job.job_id), "fusion", job.context.filepath)
            logger.error(error_msg)
            logger.error(f"Traceback: {traceback.format_exc()}")
            logger.error(f"Job details - ID: {job.job_id}, Type: {job.job_type}, File: {job.context.filepath}")
        except Exception as logger_error:
            # If we can't get a logger, print to stderr
            import sys
            print(f"FUSION HANDLER ERROR: {error_msg}", file=sys.stderr)
            print(f"LOGGER ERROR: {logger_error}", file=sys.stderr)
            print(f"TRACEBACK: {traceback.format_exc()}", file=sys.stderr)
        
        job.context.add_error("fusion_analysis", error_msg) 