"""BAM file preprocessing and metadata extraction for LittleJohn."""

import os
import time
import re
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path

import pysam
from dateutil import parser
from littlejohn.master_csv_manager import MasterCSVManager
from littlejohn.logging_config import get_job_logger


@dataclass
class BamMetadata:
    """Container for BAM file metadata and extracted data"""
    file_path: str
    file_size: int
    creation_time: float
    extracted_data: Dict[str, Any] = field(default_factory=dict)
    processing_steps: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        if self.extracted_data is None:
            self.extracted_data = {}
        if self.processing_steps is None:
            self.processing_steps = []


class ReadBam:
    """Class for reading and analyzing BAM files to extract comprehensive statistics"""
    
    def __init__(self, bam_file: str):
        self.bam_file = bam_file
        self.sam_file = None
        self.state = "pass" if bam_file and "pass" in bam_file else "fail"
        
        # Initialize counters
        self.mapped_reads = 0
        self.unmapped_reads = 0
        self.yield_tracking = 0
        self.mapped_reads_num = 0
        self.unmapped_reads_num = 0
        self.pass_mapped_reads_num = 0
        self.fail_mapped_reads_num = 0
        self.pass_unmapped_reads_num = 0
        self.fail_unmapped_reads_num = 0
        self.mapped_bases = 0
        self.unmapped_bases = 0
        self.pass_mapped_bases = 0
        self.fail_mapped_bases = 0
        self.pass_unmapped_bases = 0
        self.fail_unmapped_bases = 0
        
    def get_rg_tags(self) -> Optional[Tuple[Optional[str], ...]]:
        """
        Extracts RG (Read Group) tags from the BAM file header.

        Returns:
            Optional[Tuple[Optional[str], ...]]: A tuple containing RG tag information
                                                 or None if no RG tags are found.
        """
        if not self.sam_file:
            return None

        rg_tags = self.sam_file.header.get("RG", [])
        if not rg_tags:
            return None

        for rg_tag in rg_tags:
            id_tag = rg_tag.get("ID")
            dt_tag = rg_tag.get("DT")
            ds_tag = rg_tag.get("DS", "")
            ds_tags = ds_tag.split(" ")
            basecall_model_tag = (
                ds_tags[1].replace("basecall_model=", "") if len(ds_tags) > 1 else None
            )
            runid_tag = ds_tags[0].replace("runid=", "") if ds_tags else None
            lb_tag = rg_tag.get("LB")
            pl_tag = rg_tag.get("PL")
            pm_tag = rg_tag.get("PM")
            pu_tag = rg_tag.get("PU")
            al_tag = rg_tag.get("al")

            return (
                id_tag,
                dt_tag,
                basecall_model_tag,
                runid_tag,
                lb_tag,
                pl_tag,
                pm_tag,
                pu_tag,
                al_tag,
            )

        return None
        
    def process_reads(self) -> Optional[Dict[str, Any]]:
        """
        Processes the reads in the BAM file and aggregates information.
        Also extracts and prints RG tags for each read.

        Returns:
            Optional[Dict[str, Any]]: A dictionary containing aggregated read information,
                                      or None if processing fails.
        """
        if not self.bam_file:
            return None

        index_file = f"{self.bam_file}.bai"
        if not os.path.isfile(index_file):
            try:
                pysam.index(self.bam_file)
            except Exception:
                pass  # Continue without index file

        try:
            self.sam_file = pysam.AlignmentFile(self.bam_file, "rb", check_sq=False)
        except (IOError, ValueError):
            return None

        rg_tags = self.get_rg_tags()
        
        if not rg_tags:
            return None

        # Ensure sample_id is not None - use a fallback if needed
        sample_id = (
            rg_tags[4]
            if rg_tags[4]
            else f"unknown_sample_{os.path.basename(self.bam_file)}"
        )

        bam_read = {
            'ID': rg_tags[0],
            'time_of_run': rg_tags[1],
            'sample_id': sample_id,
            'basecall_model': rg_tags[2],
            'runid': rg_tags[3],
            'platform': rg_tags[5],
            'flow_cell_id': rg_tags[7],
            'device_position': rg_tags[6],
            'al': rg_tags[8],
            'state': self.state,
            'last_start': None,
            'elapsed_time': None,
        }

        readset = set()
        mapped_readset = set()
        unmapped_readset = set()
        barcode_found = False

        for read in self.sam_file.fetch(until_eof=True):
            # Extract RG tag and check for barcode
            if read.has_tag("RG") and not barcode_found:
                rg_tag = read.get_tag("RG")
                # Check if RG tag ends with _barcodeNN where NN is 01-96
                barcode_match = re.search(r"_barcode(\d{1,2})$", rg_tag)
                if barcode_match and bam_read['sample_id']:
                    barcode_num = int(barcode_match.group(1))
                    if 1 <= barcode_num <= 96:
                        barcode_str = f"_barcode{barcode_num:02d}"
                        if not bam_read['sample_id'].endswith(barcode_str):
                            bam_read['sample_id'] = (
                                f"{bam_read['sample_id']}{barcode_str}"
                            )
                            barcode_found = True

            read_length = read.query_length if read.query_length else 0

            if not read.is_secondary:  # Only process primary alignments
                if read.query_name not in readset:
                    readset.add(read.query_name)

                    if not read.is_unmapped:
                        mapped_readset.add(read.query_name)
                        self.mapped_bases += read_length
                        if self.state == "pass":
                            self.pass_mapped_bases += read_length
                            self.pass_mapped_reads_num += 1
                        else:
                            self.fail_mapped_bases += read_length
                            self.fail_mapped_reads_num += 1
                    else:
                        unmapped_readset.add(read.query_name)
                        self.unmapped_bases += read_length
                        if self.state == "pass":
                            self.pass_unmapped_bases += read_length
                            self.pass_unmapped_reads_num += 1
                        else:
                            self.fail_unmapped_bases += read_length
                            self.fail_unmapped_reads_num += 1

                    if read_length > 0:
                        self.yield_tracking += read_length
            try:
                if not bam_read['last_start']:
                    bam_read['last_start'] = read.get_tag("st")
                elif read.get_tag("st") > bam_read['last_start']:
                    bam_read['last_start'] = read.get_tag("st")
            except KeyError:
                bam_read['last_start'] = bam_read['time_of_run']

        self.mapped_reads = len(mapped_readset)
        self.unmapped_reads = len(unmapped_readset)
        self.mapped_reads_num = len(mapped_readset)
        self.unmapped_reads_num = len(unmapped_readset)
        
        if bam_read['last_start'] and bam_read['time_of_run']:
            bam_read['elapsed_time'] = parser.parse(bam_read['last_start']) - parser.parse(
                bam_read['time_of_run']
            )

        return bam_read
    
    def summary(self) -> Dict[str, Any]:
        """
        Returns a summary of the BAM file.

        Returns:
            dict: A dictionary containing read statistics including counts and mean lengths.
        """
        # Calculate mean lengths
        mean_mapped_length = (
            self.mapped_bases / self.mapped_reads_num
            if self.mapped_reads_num > 0
            else 0
        )
        mean_unmapped_length = (
            self.unmapped_bases / self.unmapped_reads_num
            if self.unmapped_reads_num > 0
            else 0
        )
        mean_pass_mapped_length = (
            self.pass_mapped_bases / self.pass_mapped_reads_num
            if self.pass_mapped_reads_num > 0
            else 0
        )
        mean_fail_mapped_length = (
            self.fail_mapped_bases / self.fail_mapped_reads_num
            if self.fail_mapped_reads_num > 0
            else 0
        )
        mean_pass_unmapped_length = (
            self.pass_unmapped_bases / self.pass_unmapped_reads_num
            if self.pass_unmapped_reads_num > 0
            else 0
        )
        mean_fail_unmapped_length = (
            self.fail_unmapped_bases / self.fail_unmapped_reads_num
            if self.fail_unmapped_reads_num > 0
            else 0
        )

        return {
            "mapped_reads": self.mapped_reads,
            "unmapped_reads": self.unmapped_reads,
            "yield_tracking": self.yield_tracking,
            "state": self.state,
            # Add read numbers
            "mapped_reads_num": self.mapped_reads_num,
            "unmapped_reads_num": self.unmapped_reads_num,
            "pass_mapped_reads_num": self.pass_mapped_reads_num,
            "fail_mapped_reads_num": self.fail_mapped_reads_num,
            "pass_unmapped_reads_num": self.pass_unmapped_reads_num,
            "fail_unmapped_reads_num": self.fail_unmapped_reads_num,
            # Add bases
            "mapped_bases": self.mapped_bases,
            "unmapped_bases": self.unmapped_bases,
            "pass_mapped_bases": self.pass_mapped_bases,
            "fail_mapped_bases": self.fail_mapped_bases,
            "pass_unmapped_bases": self.pass_unmapped_bases,
            "fail_unmapped_bases": self.fail_unmapped_bases,
            # Add mean lengths
            "mean_mapped_length": mean_mapped_length,
            "mean_unmapped_length": mean_unmapped_length,
            "mean_pass_mapped_length": mean_pass_mapped_length,
            "mean_fail_mapped_length": mean_fail_mapped_length,
            "mean_pass_unmapped_length": mean_pass_unmapped_length,
            "mean_fail_unmapped_length": mean_fail_unmapped_length,
        }


def extract_bam_metadata(bam_path: str) -> BamMetadata:
    """
    Extract comprehensive metadata from a BAM file.
    
    Args:
        bam_path: Path to the BAM file
        
    Returns:
        BamMetadata: Object containing all extracted metadata
    """
    # Get basic file info
    stat = os.stat(bam_path)
    file_size = stat.st_size
    creation_time = stat.st_ctime
    
    # Initialize metadata
    metadata = BamMetadata(
        file_path=bam_path,
        file_size=file_size,
        creation_time=creation_time
    )
    
    # Use ReadBam class to extract comprehensive data
    read_bam = ReadBam(bam_path)
    
    # Extract basic metadata
    bam_info = read_bam.process_reads()
    if bam_info is None:
        # Fallback to basic extraction if ReadBam fails
        sample_id = _extract_sample_id_from_bam(bam_path)
        state = 'pass' if 'pass' in bam_path else 'fail'
        filename = os.path.basename(bam_path)
        directory = os.path.dirname(bam_path)
        
        metadata.extracted_data = {
            'sample_id': sample_id,
            'file_size': file_size,
            'state': state,
            'filename': filename,
            'directory': directory,
        }
        # Note: sample_id extraction result for debugging
        # (logger not available in this context)
    else:
        # Use comprehensive data from ReadBam
        filename = os.path.basename(bam_path)
        directory = os.path.dirname(bam_path)
        
        metadata.extracted_data = {
            'sample_id': bam_info.get('sample_id', 'unknown'),
            'file_size': file_size,
            'state': bam_info.get('state', 'unknown'),
            'filename': filename,
            'directory': directory,
            'runid': bam_info.get('runid'),
            'platform': bam_info.get('platform'),
            'device_position': bam_info.get('device_position'),
            'basecall_model': bam_info.get('basecall_model'),
            'flow_cell_id': bam_info.get('flow_cell_id'),
            'time_of_run': bam_info.get('time_of_run'),
            'file_path': bam_path
        }
        # Note: ReadBam extraction result for debugging
        # (logger not available in this context)
    
    # Extract comprehensive statistics
    bam_stats = read_bam.summary()
    if bam_stats:
        metadata.extracted_data.update(bam_stats)
    
    metadata.processing_steps.append('bam_preprocessing_complete')
    
    return metadata


def _extract_sample_id_from_bam(bam_path: str) -> str:
    """Extract sample ID from BAM file read group tags"""
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            # Get read group tags
            rg_tags = bam.header.get("RG", [])
            if rg_tags:
                for rg_tag in rg_tags:
                    # Look for LB (Library) tag which contains the sample ID
                    lb_tag = rg_tag.get("LB")
                    if lb_tag:
                        return lb_tag
            
            # If no RG tags or LB tag found, return unknown
            return "unknown"
            
    except Exception:
        # If pysam fails, return unknown
        return "unknown"


def bam_preprocessing_handler(job):
    """
    Handler function for BAM preprocessing jobs.
    This function extracts metadata from BAM files and stores it in the job context.
    """
    try:
        bam_path = job.context.filepath
        
        # Get job-specific logger
        logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
        logger.info(f"Starting preprocessing for: {os.path.basename(bam_path)}")
        
        # Extract metadata from BAM file
        logger.debug(f"Extracting metadata from: {bam_path}")
        metadata = extract_bam_metadata(bam_path)
        logger.debug(f"Extracted metadata: {metadata.extracted_data}")
        
        # Store metadata in job context
        job.context.add_metadata('bam_metadata', metadata.extracted_data)
        job.context.add_metadata('file_size', metadata.file_size)
        job.context.add_metadata('creation_time', metadata.creation_time)
        job.context.add_metadata('processing_steps', metadata.processing_steps)
        
        # Update master.csv if we have comprehensive data
        if metadata.extracted_data and 'sample_id' in metadata.extracted_data:
            sample_id = metadata.extracted_data['sample_id']
            
            # Determine work directory - check job context first, then use BAM file directory as default
            work_dir = job.context.metadata.get('work_dir', os.path.dirname(bam_path))
            
            # Create MasterCSVManager and update master.csv
            try:
                csv_manager = MasterCSVManager(work_dir)
                
                # Extract BAM statistics for CSV update
                bam_stats = {}
                for key in ['mapped_reads', 'unmapped_reads', 'yield_tracking', 
                           'mapped_reads_num', 'unmapped_reads_num',
                           'pass_mapped_reads_num', 'fail_mapped_reads_num',
                           'pass_unmapped_reads_num', 'fail_unmapped_reads_num',
                           'mapped_bases', 'unmapped_bases',
                           'pass_mapped_bases', 'fail_mapped_bases',
                           'pass_unmapped_bases', 'fail_unmapped_bases']:
                    bam_stats[key] = metadata.extracted_data.get(key, 0)
                
                # Update master.csv
                csv_manager.update_master_csv(sample_id, bam_stats, metadata.extracted_data)
                
            except Exception as e:
                logger.warning(f"Could not update master.csv for {sample_id}: {e}")
        
        # Add result to context
        job.context.add_result('preprocessing', {
            'status': 'success',
            'sample_id': metadata.extracted_data.get('sample_id', 'unknown'),
            'state': metadata.extracted_data.get('state', 'unknown'),
            'mapped_reads': metadata.extracted_data.get('mapped_reads', 0),
            'unmapped_reads': metadata.extracted_data.get('unmapped_reads', 0),
            'total_yield': metadata.extracted_data.get('yield_tracking', 0),
            'file_size': metadata.file_size
        })
        
        logger.info(f"Preprocessing complete for {os.path.basename(bam_path)}")
        logger.info(f"Sample ID: {metadata.extracted_data.get('sample_id', 'unknown')}")
        logger.info(f"State: {metadata.extracted_data.get('state', 'unknown')}")
        logger.info(f"Mapped reads: {metadata.extracted_data.get('mapped_reads', 0):,}")
        logger.info(f"Unmapped reads: {metadata.extracted_data.get('unmapped_reads', 0):,}")
        logger.info(f"Total yield: {metadata.extracted_data.get('yield_tracking', 0):,}")
        
        # Log sample ID extraction details at debug level
        logger.debug(f"Sample ID extraction details:")
        logger.debug(f"  - Extracted sample_id: {metadata.extracted_data.get('sample_id', 'unknown')}")
        logger.debug(f"  - Filename: {os.path.basename(bam_path)}")
        logger.debug(f"  - Full path: {bam_path}")
        
        # Log comprehensive metadata at debug level
        logger.debug(f"BAM File Metadata:")
        logger.debug(f"File path: {metadata.file_path}")
        logger.debug(f"File size: {metadata.file_size:,} bytes")
        logger.debug(f"Creation time: {time.ctime(metadata.creation_time)}")
        
        # Log run information if available
        if metadata.extracted_data.get('runid'):
            logger.debug(f"Run ID: {metadata.extracted_data.get('runid')}")
        if metadata.extracted_data.get('platform'):
            logger.debug(f"Platform: {metadata.extracted_data.get('platform')}")
        if metadata.extracted_data.get('device_position'):
            logger.debug(f"Device: {metadata.extracted_data.get('device_position')}")
        if metadata.extracted_data.get('flow_cell_id'):
            logger.debug(f"Flow cell: {metadata.extracted_data.get('flow_cell_id')}")
        if metadata.extracted_data.get('basecall_model'):
            logger.debug(f"Basecall model: {metadata.extracted_data.get('basecall_model')}")
        if metadata.extracted_data.get('time_of_run'):
            logger.debug(f"Run time: {metadata.extracted_data.get('time_of_run')}")
        
        # Log detailed statistics
        logger.debug(f"Read Statistics:")
        logger.debug(f"Mapped reads: {metadata.extracted_data.get('mapped_reads', 0):,}")
        logger.debug(f"Unmapped reads: {metadata.extracted_data.get('unmapped_reads', 0):,}")
        logger.debug(f"Total reads: {metadata.extracted_data.get('mapped_reads', 0) + metadata.extracted_data.get('unmapped_reads', 0):,}")
        logger.debug(f"Mapped bases: {metadata.extracted_data.get('mapped_bases', 0):,}")
        logger.debug(f"Unmapped bases: {metadata.extracted_data.get('unmapped_bases', 0):,}")
        logger.debug(f"Total yield: {metadata.extracted_data.get('yield_tracking', 0):,}")
        
        # Log mean lengths if available
        if metadata.extracted_data.get('mean_mapped_length', 0) > 0:
            logger.debug(f"Mean mapped read length: {metadata.extracted_data.get('mean_mapped_length', 0):.1f}")
        if metadata.extracted_data.get('mean_unmapped_length', 0) > 0:
            logger.debug(f"Mean unmapped read length: {metadata.extracted_data.get('mean_unmapped_length', 0):.1f}")
        
        # Log pass/fail breakdown if available
        if metadata.extracted_data.get('pass_mapped_reads_num', 0) > 0 or metadata.extracted_data.get('fail_mapped_reads_num', 0) > 0:
            logger.debug(f"Pass mapped reads: {metadata.extracted_data.get('pass_mapped_reads_num', 0):,}")
            logger.debug(f"Fail mapped reads: {metadata.extracted_data.get('fail_mapped_reads_num', 0):,}")
            logger.debug(f"Pass unmapped reads: {metadata.extracted_data.get('pass_unmapped_reads_num', 0):,}")
            logger.debug(f"Fail unmapped reads: {metadata.extracted_data.get('fail_unmapped_reads_num', 0):,}")
        
        logger.debug(f"Processing steps: {', '.join(metadata.processing_steps)}")
        
    except Exception as e:
        job.context.add_error('preprocessing', str(e))
        logger.error(f"Error preprocessing {job.context.filepath}: {e}") 