"""Command-line interface for LittleJohn."""

import os
import sys
from pathlib import Path
from typing import Optional, List

import click

from littlejohn.watcher import FileWatcher
from littlejohn.workflow_simple import WorkflowRunner, default_file_classifier, Job
from littlejohn.bam_preprocessor import bam_preprocessing_handler
from littlejohn.mgmt_analysis import mgmt_handler
from littlejohn.cnv_analysis import cnv_handler
from littlejohn.bed_conversion import bed_conversion_handler
from littlejohn.sturgeon_analysis import sturgeon_handler
from littlejohn.nanodx_analysis import nanodx_handler, pannanodx_handler
from littlejohn.random_forest_analysis import random_forest_handler
from littlejohn.target_analysis import target_handler
from littlejohn.fusion_analysis import fusion_handler
from littlejohn.logging_config import configure_logging, set_debug_level, set_info_level, set_warning_level, set_error_level


@click.group()
@click.version_option()
def main() -> None:
    """Little John - Robin Hoods second in command who kept the merry men in line."""
    pass


@main.command()
def list_job_types() -> None:
    """List all available job types organized by queue category."""
    job_types = {
        "preprocessing": [
            "preprocessing - Extract metadata from BAM files"
        ],
        "bed_conversion": [
            "bed_conversion - Convert BAM files to BED format"
        ],
        "mgmt": [
            "mgmt - MGMT methylation analysis"
        ],
        "cnv": [
            "cnv - Copy number variation analysis"
        ],
        "target": [
            "target - Target analysis"
        ],
        "fusion": [
            "fusion - Fusion detection analysis"
        ],
        "classification": [
            "sturgeon - Sturgeon classification analysis",
            "nanodx - NanoDX analysis",
            "pannanodx - PanNanoDX analysis"
        ],
        "slow": [
            "random_forest - Random Forest analysis"
        ]
    }
    
    click.echo("Available job types in LittleJohn:\n")
    
    for queue, jobs in job_types.items():
        click.echo(f"{queue.upper()} QUEUE:")
        for job in jobs:
            click.echo(f"  • {job}")
        click.echo()
    
    click.echo("Usage examples:")
    click.echo("  • Simplified format (recommended): 'mgmt,sturgeon' (bed_conversion auto-added)")
    click.echo("  • Full pipeline (simplified): 'mgmt,cnv,target,fusion,sturgeon,nanodx,pannanodx,random_forest' (bed_conversion auto-added)")
    click.echo("  • Legacy format with queue prefixes: 'preprocessing:bed_conversion,mgmt:mgmt,classification:sturgeon'")
    click.echo("\nNote: 'preprocessing' is automatically added as the first step if not specified.")
    click.echo("Note: 'bed_conversion' is automatically added when needed for sturgeon, nanodx, pannanodx, or random_forest jobs.")
    click.echo("Note: Each analysis type (mgmt, cnv, target, fusion) has its own queue for parallel processing.")
    click.echo("Note: The system automatically determines the appropriate queue for each job type in simplified format.")





@main.command()
@click.argument("path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--workflow", "-w", required=True, help="Workflow plan. Can be specified in two formats:\n1. With queue prefixes: 'preprocessing:bed_conversion,mgmt:mgmt,classification:sturgeon'\n2. Simplified (auto-queue): 'mgmt,sturgeon' - system automatically determines appropriate queue for each job type and adds bed_conversion when needed"
)
@click.option(
    "--commands", "-c", multiple=True, help="Command mappings (e.g., 'index:samtools index {file}')"
)
@click.option(
    "--verbose", "-v", is_flag=True, help="Enable verbose output"
)
@click.option(
    "--no-process-existing", is_flag=True, help="Skip processing existing files, only watch for new changes"
)
@click.option(
    "--work-dir", "-d", type=click.Path(path_type=Path), help="Base output directory for analysis results"
)
@click.option(
    "--log-level", default="INFO", type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]), 
    help="Global log level (default: INFO)"
)
@click.option(
    "--job-log-level", multiple=True, 
    help="Set log level for specific job (e.g., 'preprocessing:DEBUG', 'mgmt:WARNING')"
)
@click.option(
    "--deduplicate-jobs", multiple=True,
    help="Job types to deduplicate by sample ID (e.g., 'sturgeon', 'mgmt'). Jobs of these types will only run once per sample, even if multiple upstream jobs complete simultaneously."
)
@click.option(
    "--no-progress", is_flag=True, help="Disable progress bars for file processing"
)
@click.option(
    "--analysis-workers", type=int, default=1, help="Number of analysis worker threads (default: 1, only queue that supports multiple workers)"
)
@click.option(
    "--legacy-analysis-queue", is_flag=True, help="Use legacy single analysis queue instead of separate queues per analysis type"
)
def workflow(path: Path, workflow: str, commands: tuple[str, ...], verbose: bool, no_process_existing: bool, work_dir: Optional[Path], log_level: str, job_log_level: tuple[str, ...], deduplicate_jobs: tuple[str, ...], no_progress: bool, analysis_workers: int, legacy_analysis_queue: bool) -> None:
    """Run various operations on BAM files in a directory. Preprocessing is automatically included as the first step."""
    try:
        # Configure logging
        job_levels = {}
        for job_level_spec in job_log_level:
            if ':' in job_level_spec:
                job_type, level = job_level_spec.split(':', 1)
                job_levels[job_type.strip()] = level.strip().upper()
        
        configure_logging(global_level=log_level, job_levels=job_levels)
        
        # Parse workflow plan and convert to standard format if needed
        workflow_steps = [step.strip() for step in workflow.split(",")]
        
        # Check if workflow uses simplified format (no queue prefixes)
        uses_simplified_format = all(':' not in step for step in workflow_steps)
        
        if uses_simplified_format:
            # Convert simplified format to standard format with automatic queue assignment
            workflow_steps = _convert_simplified_workflow(workflow_steps)
        
        # Ensure preprocessing is the first step to extract metadata
        if not workflow_steps or not workflow_steps[0].endswith(":preprocessing"):
            workflow_steps.insert(0, "preprocessing:preprocessing")
        
        # Parse command mappings
        command_map = {}
        for cmd_mapping in commands:
            if ":" in cmd_mapping:
                job_type, command = cmd_mapping.split(":", 1)
                command_map[job_type.strip()] = command.strip()
        
        # Create workflow runner with analysis worker configuration
        runner = WorkflowRunner(
            verbose=verbose,
            analysis_workers=analysis_workers,
            use_separate_analysis_queues=not legacy_analysis_queue
        )
        
        # Configure job deduplication
        if deduplicate_jobs:
            for job_type in deduplicate_jobs:
                runner.manager.add_deduplication_job_type(job_type.strip())
            click.echo(f"Job deduplication enabled for: {list(deduplicate_jobs)}")
        
        # Define handler configurations
        handler_configs = [
            # (queue_type, job_type, handler_func, legacy_queue_type, needs_work_dir)
            ("preprocessing", "preprocessing", bam_preprocessing_handler, None, False),
            ("bed_conversion", "bed_conversion", bed_conversion_handler, None, True),
            ("mgmt", "mgmt", mgmt_handler, "analysis", True),
            ("cnv", "cnv", cnv_handler, "analysis", True),
            ("target", "target", target_handler, "analysis", True),
            ("fusion", "fusion", fusion_handler, "analysis", True),
            ("classification", "sturgeon", sturgeon_handler, None, True),
            ("classification", "nanodx", nanodx_handler, None, True),
            ("classification", "pannanodx", pannanodx_handler, None, True),
            ("slow", "random_forest", random_forest_handler, None, True),
        ]
        
        # Register handlers
        for queue_type, job_type, handler_func, legacy_queue_type, needs_work_dir in handler_configs:
            # Determine the actual queue type based on legacy mode
            actual_queue_type = legacy_queue_type if legacy_analysis_queue and legacy_queue_type else queue_type
            
            # Create handler with work directory if specified and needed
            if work_dir and needs_work_dir:
                def create_handler_with_work_dir(handler, work_dir_path):
                    return lambda job: handler(job, work_dir=str(work_dir_path))
                final_handler = create_handler_with_work_dir(handler_func, work_dir)
            else:
                final_handler = handler_func
            
            runner.register_handler(actual_queue_type, job_type, final_handler)
        
        # Register command handlers
        for job_type, command in command_map.items():
            if ':' in job_type:
                queue_type, actual_job_type = job_type.split(':', 1)
                
                # Handle legacy analysis queue mapping
                if legacy_analysis_queue and queue_type in ['mgmt', 'cnv', 'target', 'fusion']:
                    actual_queue_type = 'analysis'
                else:
                    actual_queue_type = queue_type
                
                runner.register_command_handler(actual_queue_type, actual_job_type, command)
            else:
                # Default to preprocessing queue
                runner.register_command_handler("preprocessing", job_type, command)
        
        # Store work directory in job context for preprocessing
        if work_dir:
            def classifier_with_work_dir(filepath: str) -> List[Job]:
                jobs = default_file_classifier(filepath, workflow_steps)
                for job in jobs:
                    job.context.add_metadata('work_dir', str(work_dir))
                return jobs
            classifier_func = classifier_with_work_dir
        else:
            classifier_func = None
        
        # Display configuration
        _display_workflow_config(path, work_dir, workflow_steps, command_map, log_level, 
                               job_levels, deduplicate_jobs, legacy_analysis_queue, analysis_workers, 
                               no_process_existing, uses_simplified_format, workflow)
        
        # Run the workflow
        runner.run_workflow(
            watch_dir=str(path),
            workflow_plan=workflow_steps,
            recursive=True,
            patterns=["*.bam"],
            ignore_patterns=None,
            classifier_func=classifier_func,
            process_existing=not no_process_existing,
            show_progress=not no_progress
        )
        
    except KeyboardInterrupt:
        print("Stopping workflow...")
        runner.manager.stop(timeout=5.0)
        click.echo("\nWorkflow stopped by user")
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


def _convert_simplified_workflow(job_types: List[str]) -> List[str]:
    """Convert simplified job types to standard workflow format with queue prefixes."""
    # Queue mapping for job types
    queue_mapping = {
        "preprocessing": "preprocessing",
        "bed_conversion": "bed_conversion", 
        "mgmt": "mgmt",
        "cnv": "cnv",
        "target": "target",
        "fusion": "fusion",
        "sturgeon": "classification",
        "nanodx": "classification",
        "pannanodx": "classification",
        "random_forest": "slow",
    }
    
    # Jobs that require bed_conversion as a dependency
    jobs_requiring_bed_conversion = {"sturgeon", "nanodx", "pannanodx", "random_forest"}
    
    # Check if any jobs require bed_conversion
    needs_bed_conversion = any(job in jobs_requiring_bed_conversion for job in job_types)
    
    # Add bed_conversion if needed and not already present
    if needs_bed_conversion and "bed_conversion" not in job_types:
        job_types = ["bed_conversion"] + job_types
    
    converted_steps = []
    for job_type in job_types:
        queue_type = queue_mapping.get(job_type)
        if queue_type:
            converted_steps.append(f"{queue_type}:{job_type}")
        else:
            # Unknown job type, default to slow queue
            converted_steps.append(f"slow:{job_type}")
    
    return converted_steps


def _display_workflow_config(path: Path, work_dir: Optional[Path], workflow_steps: List[str], 
                           command_map: dict, log_level: str, job_levels: dict, 
                           deduplicate_jobs: tuple, legacy_analysis_queue: bool, 
                           analysis_workers: int, no_process_existing: bool, 
                           uses_simplified_format: bool = False, original_workflow: str = "") -> None:
    """Display workflow configuration information."""
    if no_process_existing:
        click.echo(f"Starting workflow on {path} for BAM files (skipping existing files)...")
    else:
        click.echo(f"Starting workflow on {path} for BAM files (will process existing files first)...")
    
    if work_dir:
        click.echo(f"Output directory: {work_dir}")
    
    if uses_simplified_format:
        # Extract job types from the final workflow steps
        final_job_types = [step.split(':')[1] for step in workflow_steps if ':' in step]
        click.echo(f"Workflow plan (simplified): {final_job_types}")
        click.echo(f"Auto-assigned queues: {workflow_steps}")
        
        # Check if bed_conversion was auto-added
        original_jobs = [step.strip() for step in original_workflow.split(",")]
        if "bed_conversion" in final_job_types and "bed_conversion" not in original_jobs:
            click.echo("Note: bed_conversion was automatically added as it's required for other jobs")
    else:
        click.echo(f"Workflow plan: {workflow_steps}")
    
    click.echo(f"Commands: {command_map}")
    click.echo(f"Log_level: {log_level}")
    
    if job_levels:
        click.echo(f"Job log levels: {job_levels}")
    if deduplicate_jobs:
        click.echo(f"Job deduplication: {list(deduplicate_jobs)}")
    
    click.echo(f"Worker configuration:")
    if legacy_analysis_queue:
        click.echo(f"  - Analysis queue mode: Legacy (single queue for all analysis types)")
        click.echo(f"  - Analysis workers: {analysis_workers}")
    else:
        click.echo(f"  - Analysis queue mode: Separate queues per analysis type")
        click.echo(f"  - Analysis workers per type: {analysis_workers} (MGMT, CNV, Target, Fusion each get {analysis_workers} workers)")
    
    click.echo(f"  - Other queues: 1 worker each (fixed)")
    click.echo("Press Ctrl+C to stop")
    click.echo("Running The Workflow!")


if __name__ == "__main__":
    main() 