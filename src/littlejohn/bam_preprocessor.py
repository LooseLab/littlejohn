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


# ============================================================================
# CONSTANTS AND CONFIGURATION
# ============================================================================

# Compile regex pattern once for barcode detection (module-level constant)
_BARCODE_PATTERN = re.compile(r"_barcode(\d{1,2})$")

# Pre-compile common string operations
_PASS_STR = "pass"
_FAIL_STR = "fail"
_RG_TAG = "RG"
_ST_TAG = "st"
_BASECALL_MODEL_PREFIX = "basecall_model="
_RUNID_PREFIX = "runid="


# ============================================================================
# DATA STRUCTURES
# ============================================================================

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


# ============================================================================
# CORE EXTRACTION FUNCTIONS
# ============================================================================

def get_rg_tags_from_bam(sam_file) -> Optional[Tuple[Optional[str], ...]]:
    """
    Extracts RG (Read Group) tags from the BAM file header.
    This is the first step in metadata extraction.

    Args:
        sam_file: pysam AlignmentFile object

    Returns:
        Optional[Tuple[Optional[str], ...]]: A tuple containing RG tag information
                                             or None if no RG tags are found.
    """
    if not sam_file:
        return None

    rg_tags = sam_file.header.get(_RG_TAG, [])
    if not rg_tags:
        return None

    # Process only the first RG tag for efficiency
    rg_tag = rg_tags[0]
    id_tag = rg_tag.get("ID")
    dt_tag = rg_tag.get("DT")
    ds_tag = rg_tag.get("DS", "")
    
    # Optimize string splitting and processing
    if ds_tag:
        ds_tags = ds_tag.split(" ")
        ds_tags_len = len(ds_tags)
        basecall_model_tag = (
            ds_tags[1].replace(_BASECALL_MODEL_PREFIX, "") if ds_tags_len > 1 else None
        )
        runid_tag = ds_tags[0].replace(_RUNID_PREFIX, "") if ds_tags else None
    else:
        basecall_model_tag = None
        runid_tag = None
    
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


def _extract_sample_id_from_bam(bam_path: str) -> str:
    """
    Extract sample ID from BAM file read group tags.
    This is a fallback function when comprehensive processing fails.
    """
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            # Get read group tags
            rg_tags = bam.header.get(_RG_TAG, [])
            if rg_tags:
                # Process only the first RG tag for efficiency
                rg_tag = rg_tags[0]
                
                # Look for LB (Library) tag which contains the sample ID
                lb_tag = rg_tag.get("LB")
                if lb_tag:
                    return lb_tag
                
                # Fallback to runid from DS tag per ONT specification
                ds_tag = rg_tag.get("DS", "")
                if ds_tag:
                    ds_tags = ds_tag.split(" ")
                    if ds_tags and ds_tags[0].startswith(_RUNID_PREFIX):
                        runid = ds_tags[0].replace(_RUNID_PREFIX, "")
                        if runid:
                            return runid
            
            # If no RG tags or no valid sample ID found, return unknown
            return "unknown"
            
    except Exception:
        # If pysam fails, return unknown
        return "unknown"


# ============================================================================
# READ PROCESSING FUNCTIONS
# ============================================================================

def process_bam_reads(bam_file: str) -> Optional[Dict[str, Any]]:
    """
    Processes the reads in the BAM file and aggregates information.
    This is the main processing step that analyzes all reads in the file.
    Uses streaming processing to minimize memory usage.

    Args:
        bam_file: Path to the BAM file

    Returns:
        Optional[Dict[str, Any]]: A dictionary containing aggregated read information,
                                  or None if processing fails.
    """
    if not bam_file:
        return None

    # Determine state from filename - use in operator for better performance
    state = _PASS_STR if _PASS_STR in bam_file else _FAIL_STR

    # Use context manager to ensure proper file cleanup
    try:
        with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as sam_file:
            # Step 1: Extract RG tags from header
            rg_tags = get_rg_tags_from_bam(sam_file)
            
            if not rg_tags:
                return None

            # Step 2: Determine sample ID
            # Ensure sample_id is not None - use runid as fallback per ONT specification
            # LB tag is optional, but runid from DS tag is always present
            sample_id = (
                rg_tags[4]  # LB (Library) tag
                if rg_tags[4]
                else rg_tags[3]  # runid from DS tag as fallback
                if rg_tags[3]
                else f"unknown_sample_{os.path.basename(bam_file)}"  # filename as final fallback
            )

            # Step 3: Initialize data structures
            # Pre-allocate dictionary with known size for better performance
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
                'state': state,
                'last_start': None,
                'elapsed_time': None,
            }

            # Initialize counters - use single dict for better performance
            counters = {
                'mapped_reads': 0,
                'unmapped_reads': 0,
                'yield_tracking': 0,
                'mapped_reads_num': 0,
                'unmapped_reads_num': 0,
                'pass_mapped_reads_num': 0,
                'fail_mapped_reads_num': 0,
                'pass_unmapped_reads_num': 0,
                'fail_unmapped_reads_num': 0,
                'mapped_bases': 0,
                'unmapped_bases': 0,
                'pass_mapped_bases': 0,
                'fail_mapped_bases': 0,
                'pass_unmapped_bases': 0,
                'fail_unmapped_bases': 0,
                'supplementary_reads': 0,
                'reads_with_supplementary': 0,
            }

            # Use sets only for deduplication, not for counting
            # This reduces memory usage significantly for large files
            seen_reads = set()
            reads_with_supplementary = set()  # Track read IDs with supplementary mappings
            barcode_found = False
            last_start = None
            
            # Step 3.5: Initialize MGMT read detection
            has_mgmt_reads = False
            mgmt_read_count = 0

            # Step 4: Process reads in streaming fashion
            for read in sam_file.fetch(until_eof=True):
                # Extract RG tag and check for barcode - early termination optimization
                if not barcode_found and read.has_tag(_RG_TAG):
                    rg_tag = read.get_tag(_RG_TAG)
                    # Use compiled regex pattern for barcode detection
                    barcode_match = _BARCODE_PATTERN.search(rg_tag)
                    if barcode_match and sample_id:
                        barcode_num = int(barcode_match.group(1))
                        if 1 <= barcode_num <= 96:
                            barcode_str = f"_barcode{barcode_num:02d}"
                            if not sample_id.endswith(barcode_str):
                                sample_id = f"{sample_id}{barcode_str}"
                                barcode_found = True

                # Cache read properties to avoid repeated attribute access
                read_length = read.query_length or 0
                is_secondary = read.is_secondary
                is_unmapped = read.is_unmapped
                is_supplementary = read.is_supplementary
                query_name = read.query_name

                # Track supplementary reads regardless of whether they're primary or secondary
                if is_supplementary:
                    counters['supplementary_reads'] += 1
                    reads_with_supplementary.add(query_name)

                if not is_secondary:  # Only process primary alignments
                    # Use read name for deduplication only
                    if query_name not in seen_reads:
                        seen_reads.add(query_name)

                        if not is_unmapped:
                            counters['mapped_reads'] += 1
                            counters['mapped_bases'] += read_length
                            if state == _PASS_STR:
                                counters['pass_mapped_bases'] += read_length
                                counters['pass_mapped_reads_num'] += 1
                            else:
                                counters['fail_mapped_bases'] += read_length
                                counters['fail_mapped_reads_num'] += 1
                                
                            # Check for MGMT reads (chr10:129466536-129467536)
                            if not has_mgmt_reads and read.reference_name == "chr10":
                                read_start = read.reference_start
                                read_end = read.reference_end
                                # Check if read overlaps with MGMT region (129466536-129467536)
                                if (read_start < 129467536 and read_end > 129466536):
                                    has_mgmt_reads = True
                                    mgmt_read_count += 1
                        else:
                            counters['unmapped_reads'] += 1
                            counters['unmapped_bases'] += read_length
                            if state == _PASS_STR:
                                counters['pass_unmapped_bases'] += read_length
                                counters['pass_unmapped_reads_num'] += 1
                            else:
                                counters['fail_unmapped_bases'] += read_length
                                counters['fail_unmapped_reads_num'] += 1

                        if read_length > 0:
                            counters['yield_tracking'] += read_length
                    
                    # Cache tag value to avoid multiple get_tag() calls
                    try:
                        st_tag = read.get_tag(_ST_TAG)
                        if not last_start or st_tag > last_start:
                            last_start = st_tag
                    except KeyError:
                        if not last_start:
                            last_start = bam_read['time_of_run']

            # Step 5: Finalize data
            # Update counters using the streaming counts
            counters['mapped_reads_num'] = counters['mapped_reads']
            counters['unmapped_reads_num'] = counters['unmapped_reads']
            counters['reads_with_supplementary'] = len(reads_with_supplementary)
            
            # Update sample_id in bam_read
            bam_read['sample_id'] = sample_id
            bam_read['last_start'] = last_start
            
            # Add supplementary read information to metadata
            bam_read['has_supplementary_reads'] = len(reads_with_supplementary) > 0
            bam_read['supplementary_read_ids'] = list(reads_with_supplementary) if reads_with_supplementary else []
            
            # Add MGMT read information to metadata
            bam_read['has_mgmt_reads'] = has_mgmt_reads
            bam_read['mgmt_read_count'] = mgmt_read_count
            
            # Parse dates once at the end instead of during the loop
            if last_start and bam_read['time_of_run']:
                try:
                    bam_read['elapsed_time'] = parser.parse(last_start) - parser.parse(
                        bam_read['time_of_run']
                    )
                except (ValueError, TypeError):
                    # Handle invalid date formats gracefully
                    bam_read['elapsed_time'] = None

            # Add statistics to the return dictionary
            bam_read.update(counters)

            return bam_read
            
    except (IOError, ValueError):
        return None


# ============================================================================
# STATISTICS AND SUMMARY FUNCTIONS
# ============================================================================

def calculate_bam_summary(bam_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculates summary statistics from BAM processing data.
    This step processes the raw data to generate meaningful statistics.

    Args:
        bam_data: Dictionary containing BAM processing results

    Returns:
        dict: A dictionary containing read statistics including counts and mean lengths.
    """
    # Extract values from bam_data with defaults
    mapped_reads_num = bam_data.get('mapped_reads_num', 0)
    unmapped_reads_num = bam_data.get('unmapped_reads_num', 0)
    pass_mapped_reads_num = bam_data.get('pass_mapped_reads_num', 0)
    fail_mapped_reads_num = bam_data.get('fail_mapped_reads_num', 0)
    pass_unmapped_reads_num = bam_data.get('pass_unmapped_reads_num', 0)
    fail_unmapped_reads_num = bam_data.get('fail_unmapped_reads_num', 0)
    
    mapped_bases = bam_data.get('mapped_bases', 0)
    unmapped_bases = bam_data.get('unmapped_bases', 0)
    pass_mapped_bases = bam_data.get('pass_mapped_bases', 0)
    fail_mapped_bases = bam_data.get('fail_mapped_bases', 0)
    pass_unmapped_bases = bam_data.get('pass_unmapped_bases', 0)
    fail_unmapped_bases = bam_data.get('fail_unmapped_bases', 0)
    
    # Calculate mean lengths with cached values - use conditional expressions for better performance
    mean_mapped_length = mapped_bases / mapped_reads_num if mapped_reads_num > 0 else 0
    mean_unmapped_length = unmapped_bases / unmapped_reads_num if unmapped_reads_num > 0 else 0
    mean_pass_mapped_length = pass_mapped_bases / pass_mapped_reads_num if pass_mapped_reads_num > 0 else 0
    mean_fail_mapped_length = fail_mapped_bases / fail_mapped_reads_num if fail_mapped_reads_num > 0 else 0
    mean_pass_unmapped_length = pass_unmapped_bases / pass_unmapped_reads_num if pass_unmapped_reads_num > 0 else 0
    mean_fail_unmapped_length = fail_unmapped_bases / fail_unmapped_reads_num if fail_unmapped_reads_num > 0 else 0

    # Pre-allocate result dictionary with known size
    result = {
        "mapped_reads": bam_data.get('mapped_reads', 0),
        "unmapped_reads": bam_data.get('unmapped_reads', 0),
        "yield_tracking": bam_data.get('yield_tracking', 0),
        "state": bam_data.get('state', 'unknown'),
        # Add read numbers
        "mapped_reads_num": mapped_reads_num,
        "unmapped_reads_num": unmapped_reads_num,
        "pass_mapped_reads_num": pass_mapped_reads_num,
        "fail_mapped_reads_num": fail_mapped_reads_num,
        "pass_unmapped_reads_num": pass_unmapped_reads_num,
        "fail_unmapped_reads_num": fail_unmapped_reads_num,
        # Add bases
        "mapped_bases": mapped_bases,
        "unmapped_bases": unmapped_bases,
        "pass_mapped_bases": pass_mapped_bases,
        "fail_mapped_bases": fail_mapped_bases,
        "pass_unmapped_bases": pass_unmapped_bases,
        "fail_unmapped_bases": fail_unmapped_bases,
        # Add mean lengths
        "mean_mapped_length": mean_mapped_length,
        "mean_unmapped_length": mean_unmapped_length,
        "mean_pass_mapped_length": mean_pass_mapped_length,
        "mean_fail_mapped_length": mean_fail_mapped_length,
        "mean_pass_unmapped_length": mean_pass_unmapped_length,
        "mean_fail_unmapped_length": mean_fail_unmapped_length,
        # Add supplementary read statistics
        "supplementary_reads": bam_data.get('supplementary_reads', 0),
        "reads_with_supplementary": bam_data.get('reads_with_supplementary', 0),
        "has_supplementary_reads": bam_data.get('has_supplementary_reads', False),
        # Add supplementary read IDs
        "supplementary_read_ids": bam_data.get('supplementary_read_ids', []),
        # Add MGMT read statistics
        "has_mgmt_reads": bam_data.get('has_mgmt_reads', False),
        "mgmt_read_count": bam_data.get('mgmt_read_count', 0),
    }

    return result


# ============================================================================
# MAIN EXTRACTION FUNCTION
# ============================================================================

def extract_bam_metadata(bam_path: str) -> BamMetadata:
    """
    Extract comprehensive metadata from a BAM file.
    This is the main orchestration function that coordinates all processing steps.
    
    Processing order:
    1. Extract basic file information (size, creation time)
    2. Process BAM reads and extract metadata
    3. Calculate summary statistics
    4. Compile final metadata object
    
    Args:
        bam_path: Path to the BAM file
        
    Returns:
        BamMetadata: Object containing all extracted metadata
    """
    # Step 1: Extract basic file information
    # Cache path components to avoid repeated operations
    bam_path_obj = Path(bam_path)
    filename = bam_path_obj.name
    directory = str(bam_path_obj.parent)
    
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
    
    # Step 2: Process BAM reads and extract comprehensive data
    bam_info = process_bam_reads(bam_path)
    if bam_info is None:
        # Fallback to basic extraction if processing fails
        sample_id = _extract_sample_id_from_bam(bam_path)
        state = _PASS_STR if _PASS_STR in bam_path else _FAIL_STR
        
        # Pre-allocate dictionary for better performance
        metadata.extracted_data = {
            'sample_id': sample_id,
            'file_size': file_size,
            'state': state,
            'filename': filename,
            'directory': directory,
            'has_mgmt_reads': False,
            'mgmt_read_count': 0,
        }
    else:
        # Use comprehensive data from processing - pre-allocate with known size
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
            'file_path': bam_path,
            'has_mgmt_reads': bam_info.get('has_mgmt_reads', False),
            'mgmt_read_count': bam_info.get('mgmt_read_count', 0),
        }
    
    # Step 3: Calculate summary statistics
    bam_stats = calculate_bam_summary(bam_info) if bam_info else {}
    if bam_stats:
        metadata.extracted_data.update(bam_stats)
    
    # Step 4: Mark processing as complete
    metadata.processing_steps.append('bam_preprocessing_complete')
    
    return metadata


# ============================================================================
# JOB HANDLER FUNCTION
# ============================================================================

def bam_preprocessing_handler(job):
    """
    Handler function for BAM preprocessing jobs.
    This is the main entry point that orchestrates the entire preprocessing workflow.
    
    Processing order:
    1. Initialize logging
    2. Extract metadata from BAM file
    3. Store metadata in job context
    4. Update master.csv with extracted data
    5. Add results to job context
    6. Log completion and debug information
    
    This function extracts metadata from BAM files and stores it in the job context.
    """
    try:
        bam_path = job.context.filepath
        
        # Step 1: Initialize logging
        logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
        logger.info(f"Starting preprocessing for: {os.path.basename(bam_path)}")
        
        # Step 2: Extract metadata from BAM file
        logger.debug(f"Extracting metadata from: {bam_path}")
        metadata = extract_bam_metadata(bam_path)
        logger.debug(f"Extracted metadata: {metadata.extracted_data}")
        
        # Step 3: Store metadata in job context
        job.context.add_metadata('bam_metadata', metadata.extracted_data)
        job.context.add_metadata('file_size', metadata.file_size)
        job.context.add_metadata('creation_time', metadata.creation_time)
        job.context.add_metadata('processing_steps', metadata.processing_steps)
        
        # Step 4: Update master.csv if we have comprehensive data
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
                           'pass_unmapped_bases', 'fail_unmapped_bases',
                           'supplementary_reads', 'reads_with_supplementary',
                           'has_mgmt_reads', 'mgmt_read_count']:
                    bam_stats[key] = metadata.extracted_data.get(key, 0)
                
                # Update master.csv
                csv_manager.update_master_csv(sample_id, bam_stats, metadata.extracted_data)
                
            except Exception as e:
                logger.warning(f"Could not update master.csv for {sample_id}: {e}")
        
        # Step 5: Add result to context
        job.context.add_result('preprocessing', {
            'status': 'success',
            'sample_id': metadata.extracted_data.get('sample_id', 'unknown'),
            'state': metadata.extracted_data.get('state', 'unknown'),
            'mapped_reads': metadata.extracted_data.get('mapped_reads', 0),
            'unmapped_reads': metadata.extracted_data.get('unmapped_reads', 0),
            'total_yield': metadata.extracted_data.get('yield_tracking', 0),
            'file_size': metadata.file_size,
            'has_supplementary_reads': metadata.extracted_data.get('has_supplementary_reads', False),
            'supplementary_reads': metadata.extracted_data.get('supplementary_reads', 0),
            'reads_with_supplementary': metadata.extracted_data.get('reads_with_supplementary', 0),
            'has_mgmt_reads': metadata.extracted_data.get('has_mgmt_reads', False),
            'mgmt_read_count': metadata.extracted_data.get('mgmt_read_count', 0)
        })
        
        # Step 6: Log completion and summary
        logger.info(f"Preprocessing complete for {os.path.basename(bam_path)}")
        logger.info(f"Sample ID: {metadata.extracted_data.get('sample_id', 'unknown')}")
        logger.info(f"State: {metadata.extracted_data.get('state', 'unknown')}")
        logger.info(f"Mapped reads: {metadata.extracted_data.get('mapped_reads', 0):,}")
        logger.info(f"Unmapped reads: {metadata.extracted_data.get('unmapped_reads', 0):,}")
        logger.info(f"Total yield: {metadata.extracted_data.get('yield_tracking', 0):,}")
        
        # Log supplementary read information
        has_supplementary = metadata.extracted_data.get('has_supplementary_reads', False)
        supplementary_count = metadata.extracted_data.get('supplementary_reads', 0)
        reads_with_supplementary = metadata.extracted_data.get('reads_with_supplementary', 0)
        
        if has_supplementary:
            logger.info(f"Supplementary reads found: {supplementary_count:,} supplementary alignments from {reads_with_supplementary:,} unique reads")
        else:
            logger.info("No supplementary reads found")
        
        # Log MGMT read information
        has_mgmt_reads = metadata.extracted_data.get('has_mgmt_reads', False)
        mgmt_read_count = metadata.extracted_data.get('mgmt_read_count', 0)
        
        if has_mgmt_reads:
            logger.info(f"MGMT reads found: {mgmt_read_count:,} reads spanning MGMT promoter region (chr10:129466536-129467536)")
        else:
            logger.info("No MGMT reads found")
        
        # Step 7: Detailed debug logging (only if debug is enabled)
        if logger.logger.isEnabledFor(10):  # DEBUG level
            # Cache frequently accessed values to avoid repeated dict lookups
            extracted_data = metadata.extracted_data
            sample_id = extracted_data.get('sample_id', 'unknown')
            mapped_reads = extracted_data.get('mapped_reads', 0)
            unmapped_reads = extracted_data.get('unmapped_reads', 0)
            
            # Log sample ID extraction details at debug level
            logger.debug(f"Sample ID extraction details:")
            logger.debug(f"  - Extracted sample_id: {sample_id}")
            logger.debug(f"  - Filename: {os.path.basename(bam_path)}")
            logger.debug(f"  - Full path: {bam_path}")
            
            # Log comprehensive metadata at debug level
            logger.debug(f"BAM File Metadata:")
            logger.debug(f"File path: {metadata.file_path}")
            logger.debug(f"File size: {metadata.file_size:,} bytes")
            logger.debug(f"Creation time: {time.ctime(metadata.creation_time)}")
            
            # Log run information if available - use cached values
            runid = extracted_data.get('runid')
            if runid:
                logger.debug(f"Run ID: {runid}")
            platform = extracted_data.get('platform')
            if platform:
                logger.debug(f"Platform: {platform}")
            device_position = extracted_data.get('device_position')
            if device_position:
                logger.debug(f"Device: {device_position}")
            flow_cell_id = extracted_data.get('flow_cell_id')
            if flow_cell_id:
                logger.debug(f"Flow cell: {flow_cell_id}")
            basecall_model = extracted_data.get('basecall_model')
            if basecall_model:
                logger.debug(f"Basecall model: {basecall_model}")
            time_of_run = extracted_data.get('time_of_run')
            if time_of_run:
                logger.debug(f"Run time: {time_of_run}")
            
            # Log detailed statistics - use cached values
            logger.debug(f"Read Statistics:")
            logger.debug(f"Mapped reads: {mapped_reads:,}")
            logger.debug(f"Unmapped reads: {unmapped_reads:,}")
            logger.debug(f"Total reads: {mapped_reads + unmapped_reads:,}")
            
            mapped_bases = extracted_data.get('mapped_bases', 0)
            unmapped_bases = extracted_data.get('unmapped_bases', 0)
            yield_tracking = extracted_data.get('yield_tracking', 0)
            logger.debug(f"Mapped bases: {mapped_bases:,}")
            logger.debug(f"Unmapped bases: {unmapped_bases:,}")
            logger.debug(f"Total yield: {yield_tracking:,}")
            
            # Log mean lengths if available
            mean_mapped_length = extracted_data.get('mean_mapped_length', 0)
            if mean_mapped_length > 0:
                logger.debug(f"Mean mapped read length: {mean_mapped_length:.1f}")
            mean_unmapped_length = extracted_data.get('mean_unmapped_length', 0)
            if mean_unmapped_length > 0:
                logger.debug(f"Mean unmapped read length: {mean_unmapped_length:.1f}")
            
            # Log pass/fail breakdown if available
            pass_mapped_reads_num = extracted_data.get('pass_mapped_reads_num', 0)
            fail_mapped_reads_num = extracted_data.get('fail_mapped_reads_num', 0)
            if pass_mapped_reads_num > 0 or fail_mapped_reads_num > 0:
                logger.debug(f"Pass mapped reads: {pass_mapped_reads_num:,}")
                logger.debug(f"Fail mapped reads: {fail_mapped_reads_num:,}")
                logger.debug(f"Pass unmapped reads: {extracted_data.get('pass_unmapped_reads_num', 0):,}")
                logger.debug(f"Fail unmapped reads: {extracted_data.get('fail_unmapped_reads_num', 0):,}")
            
            # Log supplementary read information if available
            has_supplementary = extracted_data.get('has_supplementary_reads', False)
            if has_supplementary:
                supplementary_count = extracted_data.get('supplementary_reads', 0)
                reads_with_supplementary = extracted_data.get('reads_with_supplementary', 0)
                supplementary_read_ids = extracted_data.get('supplementary_read_ids', [])
                logger.debug(f"Supplementary reads: {supplementary_count:,} supplementary alignments")
                logger.debug(f"Reads with supplementary mappings: {reads_with_supplementary:,}")
                if supplementary_read_ids:
                    logger.debug(f"First 10 supplementary read IDs: {supplementary_read_ids[:10]}")
                    if len(supplementary_read_ids) > 10:
                        logger.debug(f"... and {len(supplementary_read_ids) - 10} more")
            else:
                logger.debug("No supplementary reads found")
            
            # Log MGMT read information if available
            has_mgmt_reads = extracted_data.get('has_mgmt_reads', False)
            if has_mgmt_reads:
                mgmt_read_count = extracted_data.get('mgmt_read_count', 0)
                logger.debug(f"MGMT reads: {mgmt_read_count:,} reads spanning MGMT promoter region")
                logger.debug(f"MGMT region: chr10:129466536-129467536")
            else:
                logger.debug("No MGMT reads found")
            
            logger.debug(f"Processing steps: {', '.join(metadata.processing_steps)}")
        
    except Exception as e:
        job.context.add_error('preprocessing', str(e))
        logger.error(f"Error preprocessing {job.context.filepath}: {e}") 