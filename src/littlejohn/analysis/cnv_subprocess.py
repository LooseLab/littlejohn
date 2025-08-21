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

# Import cnv_from_bam only in this subprocess
import cnv_from_bam

def run_cnv_analysis(
    bam_path,
    copy_numbers_path=None,
    ref_cnv_dict_path=None,
    output_dir=None,
    threads=1,
    mapq_filter=60,
    update_cnv_dict_path=None,
    sample_id=None,
):
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
    print ("CNV Subprocess started")
    print(f"BAM path: {bam_path}")
    print(f"Copy numbers path: {copy_numbers_path}")
    print(f"Update CNV dict path: {update_cnv_dict_path}")
    print(f"Sample ID: {sample_id}")
    print(f"Ref CNV dict path: {ref_cnv_dict_path}")
    print(f"Output dir: {output_dir}")
    print(f"Threads: {threads}")
    print(f"Mapq filter: {mapq_filter}")
    
    try:
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Verify input files exist
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
        if ref_cnv_dict_path is None or not os.path.exists(ref_cnv_dict_path):
            raise FileNotFoundError(f"Reference CNV dict file not found: {ref_cnv_dict_path}")
        
        # Load copy numbers either from per-sample file or from multi-sample dict
        if copy_numbers_path is not None:
            with open(copy_numbers_path, 'rb') as f:
                copy_numbers = pickle.load(f)
        else:
            if update_cnv_dict_path is None or sample_id is None:
                raise ValueError("Either copy_numbers_path must be provided or both update_cnv_dict_path and sample_id must be provided")
            if not os.path.exists(update_cnv_dict_path):
                # If the multi-sample file doesn't exist yet, start with empty dict
                copy_numbers = {}
            else:
                with open(update_cnv_dict_path, 'rb') as f:
                    multi_sample_dict = pickle.load(f)
                copy_numbers = multi_sample_dict.get(sample_id, {})
        
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
    parser.add_argument('--copy-numbers-path', required=False, help='Path to per-sample copy numbers pickle file')
    parser.add_argument('--update-cnv-dict-path', required=False, help='Path to multi-sample update_cnv_dict.pkl file')
    parser.add_argument('--sample-id', required=False, help='Sample ID to extract from multi-sample copy numbers file')
    parser.add_argument('--ref-cnv-dict-path', required=True, help='Path to reference CNV dict pickle file')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--mapq-filter', type=int, default=60, help='Mapping quality filter')
    
    args = parser.parse_args()
    
    # Print debug information
    print(f"Subprocess working directory: {os.getcwd()}", file=sys.stderr)
    print(f"BAM path: {args.bam_path}", file=sys.stderr)
    print(f"Copy numbers path: {args.copy_numbers_path}", file=sys.stderr)
    print(f"Update CNV dict path: {args.update_cnv_dict_path}", file=sys.stderr)
    print(f"Sample ID: {args.sample_id}", file=sys.stderr)
    print(f"Ref CNV dict path: {args.ref_cnv_dict_path}", file=sys.stderr)
    print(f"Output dir: {args.output_dir}", file=sys.stderr)
    
    # Check if required files exist
    if not os.path.exists(args.bam_path):
        print(f"ERROR: BAM file does not exist: {args.bam_path}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.ref_cnv_dict_path):
        print(f"ERROR: Reference CNV dict file does not exist: {args.ref_cnv_dict_path}", file=sys.stderr)
        sys.exit(1)
    if args.copy_numbers_path is None:
        # Using multi-sample mode: update_cnv_dict_path and sample_id must be provided
        if args.update_cnv_dict_path is not None and args.sample_id is not None:
            # ok even if file doesn't exist; we'll start empty
            pass
        else:
            print("ERROR: Provide either --copy-numbers-path or both --update-cnv-dict-path and --sample-id", file=sys.stderr)
            sys.exit(1)
    
    # Run the analysis
    result = run_cnv_analysis(
        args.bam_path,
        copy_numbers_path=args.copy_numbers_path,
        ref_cnv_dict_path=args.ref_cnv_dict_path,
        output_dir=args.output_dir,
        threads=args.threads,
        mapq_filter=args.mapq_filter,
        update_cnv_dict_path=args.update_cnv_dict_path,
        sample_id=args.sample_id,
    )
    
    # Exit with appropriate code
    sys.exit(0 if result['success'] else 1)

if __name__ == '__main__':
    main() 