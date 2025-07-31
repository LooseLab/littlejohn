#!/usr/bin/env python3
"""
Subprocess script for running cnv_from_bam analysis.
This isolates the cnv_from_bam module from the main process to prevent signal handling issues.
"""

import sys
import os
import json
import pickle
import logging
import argparse
from pathlib import Path

def run_cnv_analysis(bam_path, copy_numbers_path, ref_cnv_dict_path, output_dir, threads=1, mapq_filter=60):
    """
    Run CNV analysis using cnv_from_bam in an isolated subprocess.
    
    Args:
        bam_path: Path to BAM file
        copy_numbers_path: Path to pickle file with copy numbers
        ref_cnv_dict_path: Path to pickle file with reference CNV data
        output_dir: Directory to save results
        threads: Number of threads to use
        mapq_filter: Mapping quality filter
    
    Returns:
        Dictionary with analysis results
    """
    try:
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Verify input files exist
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
        if not os.path.exists(copy_numbers_path):
            raise FileNotFoundError(f"Copy numbers file not found: {copy_numbers_path}")
        if not os.path.exists(ref_cnv_dict_path):
            raise FileNotFoundError(f"Reference CNV dict file not found: {ref_cnv_dict_path}")
        
        # Import cnv_from_bam only in this subprocess
        import cnv_from_bam
        
        # Load copy numbers
        with open(copy_numbers_path, 'rb') as f:
            copy_numbers = pickle.load(f)
        
        # Load reference CNV data
        with open(ref_cnv_dict_path, 'rb') as f:
            ref_cnv_dict = pickle.load(f)
        
        # First pass: process sample with accumulated copy numbers
        result = cnv_from_bam.iterate_bam_file(
            bam_path,
            _threads=threads,
            mapq_filter=mapq_filter,
            copy_numbers=copy_numbers,
            log_level=int(logging.ERROR),
        )
        
        # Second pass: process against reference using the same bin width
        result2 = cnv_from_bam.iterate_bam_file(
            bam_path,
            _threads=threads,
            mapq_filter=mapq_filter,
            copy_numbers=ref_cnv_dict,
            log_level=int(logging.ERROR),
            bin_width=result.bin_width,  # Use the same bin width as the sample
        )
        
        # Prepare results
        analysis_results = {
            'success': True,
            'r_cnv': result.cnv,
            'r_bin': result.bin_width,
            'r_var': result.variance,
            'genome_length': result.genome_length,
            'r2_cnv': result2.cnv,
            'updated_copy_numbers': copy_numbers  # The mutated copy_numbers
        }
        
        # Save results to output directory
        results_path = os.path.join(output_dir, 'cnv_analysis_results.pkl')
        with open(results_path, 'wb') as f:
            pickle.dump(analysis_results, f)
        
        return analysis_results
        
    except Exception as e:
        error_result = {
            'success': False,
            'error': str(e),
            'error_type': type(e).__name__
        }
        
        # Ensure output directory exists before trying to write error file
        try:
            os.makedirs(output_dir, exist_ok=True)
            # Save error to output directory
            error_path = os.path.join(output_dir, 'cnv_analysis_error.json')
            with open(error_path, 'w') as f:
                json.dump(error_result, f)
        except Exception as write_error:
            # If we can't write the error file, print to stderr
            print(f"Error writing error file: {write_error}", file=sys.stderr)
            print(f"Original error: {error_result}", file=sys.stderr)
        
        return error_result

def main():
    parser = argparse.ArgumentParser(description='Run CNV analysis in subprocess')
    parser.add_argument('--bam-path', required=True, help='Path to BAM file')
    parser.add_argument('--copy-numbers-path', required=True, help='Path to copy numbers pickle file')
    parser.add_argument('--ref-cnv-dict-path', required=True, help='Path to reference CNV dict pickle file')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--mapq-filter', type=int, default=60, help='Mapping quality filter')
    
    args = parser.parse_args()
    
    # Print debug information
    print(f"Subprocess working directory: {os.getcwd()}", file=sys.stderr)
    print(f"BAM path: {args.bam_path}", file=sys.stderr)
    print(f"Copy numbers path: {args.copy_numbers_path}", file=sys.stderr)
    print(f"Ref CNV dict path: {args.ref_cnv_dict_path}", file=sys.stderr)
    print(f"Output dir: {args.output_dir}", file=sys.stderr)
    
    # Check if files exist
    for path_name, path in [
        ('BAM file', args.bam_path),
        ('Copy numbers file', args.copy_numbers_path),
        ('Reference CNV dict file', args.ref_cnv_dict_path)
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {path_name} does not exist: {path}", file=sys.stderr)
            sys.exit(1)
    
    # Run the analysis
    result = run_cnv_analysis(
        args.bam_path,
        args.copy_numbers_path,
        args.ref_cnv_dict_path,
        args.output_dir,
        args.threads,
        args.mapq_filter
    )
    
    # Exit with appropriate code
    sys.exit(0 if result['success'] else 1)

if __name__ == '__main__':
    main() 