"""
Command-line interface for LittleJohn.

This module provides a comprehensive CLI for running bioinformatics workflows on BAM files.
It supports both simplified and legacy workflow formats, with options for distributed
computing using Ray and traditional threading.

Key Features:
- Simplified workflow syntax: 'mgmt,sturgeon' (auto-queue assignment)
- Legacy workflow syntax: 'preprocessing:bed_conversion,mgmt:mgmt,classification:sturgeon'
- Automatic dependency management (e.g., bed_conversion for sturgeon analysis)
- Ray-based distributed computing support
- Configurable logging levels per job type
- Job deduplication by sample ID
- Progress tracking and verbose output options

Code Quality Improvements:
- Modular function design with single responsibilities
- Comprehensive input validation and error handling
- Constants for configuration values
- Better type hints and documentation
- Graceful error handling with informative messages
- Consistent error reporting to stderr
- Proper cleanup on interruption
"""

import os
import sys
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Any

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


# Constants
VALID_LOG_LEVELS = {"DEBUG", "INFO", "WARNING", "ERROR"}
VALID_JOB_TYPES = {
    "preprocessing", "bed_conversion", "mgmt", "cnv", "target", 
    "fusion", "sturgeon", "nanodx", "pannanodx", "random_forest"
}
DEFAULT_LOG_LEVEL = "INFO"
DEFAULT_ANALYSIS_WORKERS = 1
DEFAULT_TIMEOUT = 5.0

# Queue mapping for job types
QUEUE_MAPPING = {
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
JOBS_REQUIRING_BED_CONVERSION = {"sturgeon", "nanodx", "pannanodx", "random_forest"}

# Handler configurations
HANDLER_CONFIGS = [
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


def _parse_job_log_levels(job_log_level: Tuple[str, ...]) -> Dict[str, str]:
    """Parse job-specific log level specifications."""
    job_levels = {}
    
    for job_level_spec in job_log_level:
        if ':' in job_level_spec:
            job_type, level = job_level_spec.split(':', 1)
            job_type = job_type.strip()
            level = level.strip().upper()
            
            if level not in VALID_LOG_LEVELS:
                click.echo(f"Warning: Invalid log level '{level}' for job '{job_type}'. Valid levels: {', '.join(VALID_LOG_LEVELS)}", err=True)
                continue
                
            job_levels[job_type] = level
        else:
            click.echo(f"Warning: Invalid job log level format '{job_level_spec}'. Use format 'job_type:level'", err=True)
    return job_levels


def _parse_command_mappings(commands: Tuple[str, ...]) -> Dict[str, str]:
    """Parse command mapping specifications."""
    command_map = {}
    for cmd_mapping in commands:
        if ":" in cmd_mapping:
            job_type, command = cmd_mapping.split(":", 1)
            job_type = job_type.strip()
            command = command.strip()
            
            if not job_type or not command:
                click.echo(f"Warning: Empty job type or command in '{cmd_mapping}'", err=True)
                continue
                
            command_map[job_type] = command
        else:
            click.echo(f"Warning: Invalid command mapping format '{cmd_mapping}'. Use format 'job_type:command'", err=True)
    return command_map


def _validate_workflow_steps(workflow_steps: List[str]) -> List[str]:
    """Validate and clean workflow steps."""
    cleaned_steps = []
    for step in workflow_steps:
        step = step.strip()
        if not step:
            continue
            
        if ':' in step:
            queue_type, job_type = step.split(':', 1)
            queue_type = queue_type.strip()
            job_type = job_type.strip()
            
            if job_type not in VALID_JOB_TYPES:
                click.echo(f"Warning: Unknown job type '{job_type}' in step '{step}'", err=True)
                continue
                
            cleaned_steps.append(f"{queue_type}:{job_type}")
        else:
            # Simplified format - job type only
            if step not in VALID_JOB_TYPES:
                click.echo(f"Warning: Unknown job type '{step}'", err=True)
                continue
            cleaned_steps.append(step)
    
    return cleaned_steps


def _ensure_preprocessing_step(workflow_steps: List[str]) -> List[str]:
    """Ensure preprocessing is the first step in the workflow."""
    if not workflow_steps or not workflow_steps[0].endswith(":preprocessing"):
        workflow_steps.insert(0, "preprocessing:preprocessing")
    return workflow_steps


def _convert_simplified_workflow(job_types: List[str]) -> List[str]:
    """Convert simplified job types to standard workflow format with queue prefixes."""
    # Check if any jobs require bed_conversion as a dependency
    needs_bed_conversion = any(job in JOBS_REQUIRING_BED_CONVERSION for job in job_types)
    
    # Add bed_conversion if needed and not already present
    if needs_bed_conversion and "bed_conversion" not in job_types:
        job_types = ["bed_conversion"] + job_types
    
    converted_steps = []
    for job_type in job_types:
        queue_type = QUEUE_MAPPING.get(job_type)
        if queue_type:
            converted_steps.append(f"{queue_type}:{job_type}")
        else:
            # Unknown job type, default to slow queue
            converted_steps.append(f"slow:{job_type}")
    
    return converted_steps


def _create_ray_workflow_runner(verbose: bool, analysis_workers: int, 
                               legacy_analysis_queue: bool, log_level: str) -> Any:
    """Create and configure Ray-based workflow runner."""
    try:
        from littlejohn.workflow_ray import RayWorkflowRunner
        runner = RayWorkflowRunner(
            verbose=verbose,
            analysis_workers=analysis_workers,
            use_separate_analysis_queues=not legacy_analysis_queue,
            log_level=log_level
        )
        click.echo("Using Ray-based distributed computing workflow")
        return runner
    except ImportError as e:
        click.echo(f"Warning: Ray not available ({e}). Falling back to threading-based workflow.", err=True)
        return None


def _initialize_ray(num_cpus: Optional[int]) -> None:
    """Initialize Ray with specified CPU count."""
    if num_cpus is not None:
        if num_cpus <= 0:
            click.echo(f"Warning: Invalid CPU count {num_cpus}. Must be positive.", err=True)
            return
            
        import ray
        if not ray.is_initialized():
            # Suppress Ray logging to prevent interference with progress bars
            import logging
            ray_logger = logging.getLogger("ray")
            ray_logger.setLevel(logging.ERROR)
            logging.getLogger("ray.worker").setLevel(logging.ERROR)
            logging.getLogger("ray.remote").setLevel(logging.ERROR)
            logging.getLogger("ray.actor").setLevel(logging.ERROR)
            logging.getLogger("ray.util").setLevel(logging.ERROR)
            
            try:
                ray.init(num_cpus=num_cpus, ignore_reinit_error=True)
                click.echo(f"Ray initialized with {num_cpus} CPUs")
            except Exception as e:
                click.echo(f"Warning: Failed to initialize Ray: {e}", err=True)
        else:
            click.echo(f"Ray already initialized with {ray.available_resources().get('CPU', 'unknown')} CPUs")


def _configure_queue_priorities(runner: Any, queue_priority: Tuple[str, ...]) -> None:
    """Configure queue priorities for Ray workflow."""
    if not queue_priority:
        return
        
    click.echo("Configuring queue priorities:")
    for priority_spec in queue_priority:
        if ':' in priority_spec:
            queue_type, priority_str = priority_spec.split(':', 1)
            queue_type = queue_type.strip()
            priority_str = priority_str.strip()
            
            try:
                priority = int(priority_str)
                if priority < 0:
                    click.echo(f"Warning: Priority {priority} for queue '{queue_type}' is negative. Using 0 instead.", err=True)
                    priority = 0
                    
                runner.manager.set_queue_priority(queue_type, priority)
                click.echo(f"  {queue_type}: {priority}")
            except ValueError:
                click.echo(f"Warning: Invalid priority value '{priority_str}' for queue '{queue_type}'. Must be an integer.", err=True)
        else:
            click.echo(f"Warning: Invalid priority specification '{priority_spec}'. Use format 'queue:priority'", err=True)


def _create_workflow_runner(use_ray: bool, verbose: bool, analysis_workers: int, 
                           legacy_analysis_queue: bool, log_level: str) -> Any:
    """Create the appropriate workflow runner based on configuration."""
    if analysis_workers < 1:
        click.echo(f"Warning: Invalid analysis_workers value {analysis_workers}. Using {DEFAULT_ANALYSIS_WORKERS} instead.", err=True)
        analysis_workers = DEFAULT_ANALYSIS_WORKERS
        
    if use_ray:
        runner = _create_ray_workflow_runner(verbose, analysis_workers, legacy_analysis_queue, log_level)
        if runner is None:
            # Fallback to threading-based workflow
            from littlejohn.workflow_simple import WorkflowRunner
            runner = WorkflowRunner(
                verbose=verbose,
                analysis_workers=analysis_workers,
                use_separate_analysis_queues=not legacy_analysis_queue
            )
        return runner
    else:
        from littlejohn.workflow_simple import WorkflowRunner
        return WorkflowRunner(
            verbose=verbose,
            analysis_workers=analysis_workers,
            use_separate_analysis_queues=not legacy_analysis_queue
        )


def _register_handlers(runner: Any, legacy_analysis_queue: bool, work_dir: Optional[Path]) -> None:
    """Register all workflow handlers with the runner."""
    for queue_type, job_type, handler_func, legacy_queue_type, needs_work_dir in HANDLER_CONFIGS:
        # Determine the actual queue type based on legacy mode
        actual_queue_type = legacy_queue_type if legacy_analysis_queue and legacy_queue_type else queue_type
        
        # Create handler with work directory if specified and needed
        if work_dir and needs_work_dir:
            def create_handler_with_work_dir(handler, work_dir_path):
                return lambda job: handler(job, work_dir=str(work_dir_path))
            final_handler = create_handler_with_work_dir(handler_func, work_dir)
        else:
            final_handler = handler_func
        
        try:
            runner.register_handler(actual_queue_type, job_type, final_handler)
        except Exception as e:
            click.echo(f"Warning: Failed to register handler for {queue_type}:{job_type}: {e}", err=True)


def _register_command_handlers(runner: Any, command_map: Dict[str, str], 
                             legacy_analysis_queue: bool) -> None:
    """Register command handlers with the runner."""
    for job_type, command in command_map.items():
        if ':' in job_type:
            queue_type, actual_job_type = job_type.split(':', 1)
            
            # Handle legacy analysis queue mapping
            if legacy_analysis_queue and queue_type in ['mgmt', 'cnv', 'target', 'fusion']:
                actual_queue_type = 'analysis'
            else:
                actual_queue_type = queue_type
            
            try:
                runner.register_command_handler(actual_queue_type, actual_job_type, command)
            except Exception as e:
                click.echo(f"Warning: Failed to register command handler for {queue_type}:{actual_job_type}: {e}", err=True)
        else:
            # Default to preprocessing queue
            try:
                runner.register_command_handler("preprocessing", job_type, command)
            except Exception as e:
                click.echo(f"Warning: Failed to register command handler for preprocessing:{job_type}: {e}", err=True)


def _create_classifier_with_work_dir(work_dir: Path, workflow_steps: List[str]):
    """Create a classifier function that includes work directory in job context."""
    def classifier_with_work_dir(filepath: str) -> List[Job]:
        jobs = default_file_classifier(filepath, workflow_steps)
        for job in jobs:
            job.context.add_metadata('work_dir', str(work_dir))
        return jobs
    return classifier_with_work_dir


def _validate_inputs(path: Path, workflow: str, analysis_workers: int, ray_num_cpus: Optional[int]) -> None:
    """Validate input parameters and provide helpful error messages."""
    if not path.exists():
        raise click.BadParameter(f"Path '{path}' does not exist")
    
    if not path.is_dir():
        raise click.BadParameter(f"Path '{path}' is not a directory")
    
    if not workflow.strip():
        raise click.BadParameter("Workflow cannot be empty")
    
    if analysis_workers < 1:
        raise click.BadParameter("analysis_workers must be at least 1")
    
    if ray_num_cpus is not None and ray_num_cpus < 1:
        raise click.BadParameter("ray_num_cpus must be at least 1")


def _display_workflow_config(path: Path, work_dir: Optional[Path], workflow_steps: List[str], 
                           command_map: dict, log_level: str, job_levels: dict, 
                           deduplicate_jobs: tuple, legacy_analysis_queue: bool, 
                           analysis_workers: int, no_process_existing: bool, 
                           uses_simplified_format: bool = False, original_workflow: str = "", 
                           use_ray: bool = False, ray_num_cpus: Optional[int] = None,
                           queue_priority: tuple = ()) -> None:
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
    
    # Display Ray configuration if enabled
    if use_ray:
        click.echo(f"Distributed computing: Ray (experimental)")
        if ray_num_cpus:
            click.echo(f"  - Ray CPUs: {ray_num_cpus}")
        else:
            click.echo(f"  - Ray CPUs: auto-detect")
        click.echo(f"  - Analysis workers per type: {analysis_workers} (distributed across Ray cluster)")
        click.echo(f"  - Log level: {log_level} (applied to all Ray actors)")
        
        # Display priority configuration
        if queue_priority:
            click.echo(f"  - Queue priorities: {list(queue_priority)}")
        else:
            click.echo(f"  - Queue priorities: using defaults (preprocessing=10, bed_conversion=9, analysis=5, classification=3, slow=1)")
    else:
        click.echo(f"Distributed computing: Disabled (using threading)")
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


@main.command()
@click.argument("path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--workflow", "-w", required=True, 
    help="Workflow plan. Can be specified in two formats:\n1. With queue prefixes: 'preprocessing:bed_conversion,mgmt:mgmt,classification:sturgeon'\n2. Simplified (auto-queue): 'mgmt,sturgeon' - system automatically determines appropriate queue for each job type and adds bed_conversion when needed"
)
@click.option(
    "--commands", "-c", multiple=True, 
    help="Command mappings (e.g., 'index:samtools index {file}')"
)
@click.option(
    "--verbose", "-v", is_flag=True, 
    help="Enable verbose output and detailed error traces"
)
@click.option(
    "--no-process-existing", is_flag=True, 
    help="Skip processing existing files, only watch for new changes"
)
@click.option(
    "--work-dir", "-d", type=click.Path(path_type=Path), 
    help="Base output directory for analysis results"
)
@click.option(
    "--log-level", default=DEFAULT_LOG_LEVEL, type=click.Choice(list(VALID_LOG_LEVELS)), 
    help=f"Global log level (default: {DEFAULT_LOG_LEVEL})"
)
@click.option(
    "--job-log-level", multiple=True, 
    help=f"Set log level for specific job (e.g., 'preprocessing:DEBUG', 'mgmt:WARNING'). Valid levels: {', '.join(VALID_LOG_LEVELS)}"
)
@click.option(
    "--deduplicate-jobs", multiple=True,
    help="Job types to deduplicate by sample ID (e.g., 'sturgeon', 'mgmt'). Jobs of these types will only run once per sample, even if multiple upstream jobs complete simultaneously."
)
@click.option(
    "--no-progress", is_flag=True, 
    help="Disable progress bars for file processing"
)
@click.option(
    "--analysis-workers", type=int, default=DEFAULT_ANALYSIS_WORKERS, 
    help=f"Number of analysis worker threads (default: {DEFAULT_ANALYSIS_WORKERS}, only queue that supports multiple workers)"
)
@click.option(
    "--legacy-analysis-queue", is_flag=True, 
    help="Use legacy single analysis queue instead of separate queues per analysis type"
)
@click.option(
    "--use-ray", is_flag=True, 
    help="Use Ray for distributed computing (experimental). This will distribute jobs across multiple CPU cores and potentially multiple machines for improved performance."
)
@click.option(
    "--ray-num-cpus", type=int, default=None, 
    help="Number of CPUs to use for Ray (default: auto-detect). Only used when --use-ray is specified."
)
@click.option(
    "--queue-priority", multiple=True, 
    help="Set priority for a queue (e.g., 'preprocessing:10', 'bed_conversion:9', 'mgmt:5'). Only used when --use-ray is specified."
)
@click.option(
    "--show-priorities", is_flag=True, 
    help="Show current queue priorities and exit. Only used when --use-ray is specified."
)
def workflow(path: Path, workflow: str, commands: tuple[str, ...], verbose: bool, 
            no_process_existing: bool, work_dir: Optional[Path], log_level: str, 
            job_log_level: tuple[str, ...], deduplicate_jobs: tuple[str, ...], 
            no_progress: bool, analysis_workers: int, legacy_analysis_queue: bool, 
            use_ray: bool, ray_num_cpus: Optional[int], queue_priority: tuple[str, ...], 
            show_priorities: bool) -> None:
    """Run various operations on BAM files in a directory. Preprocessing is automatically included as the first step."""
    try:
        # Validate input parameters
        _validate_inputs(path, workflow, analysis_workers, ray_num_cpus)
        
        # Parse and validate inputs
        job_levels = _parse_job_log_levels(job_log_level)
        command_map = _parse_command_mappings(commands)
        
        # Configure logging
        configure_logging(global_level=log_level, job_levels=job_levels)
        
        # Parse workflow plan and convert to standard format if needed
        workflow_steps = [step.strip() for step in workflow.split(",")]
        
        # Validate and clean workflow steps
        workflow_steps = _validate_workflow_steps(workflow_steps)
        
        # Check if workflow uses simplified format (no queue prefixes)
        uses_simplified_format = all(':' not in step for step in workflow_steps)
        
        if uses_simplified_format:
            # Convert simplified format to standard format with automatic queue assignment
            workflow_steps = _convert_simplified_workflow(workflow_steps)
        
        # Ensure preprocessing is the first step to extract metadata
        workflow_steps = _ensure_preprocessing_step(workflow_steps)
        
        # Validate that we have at least some workflow steps
        if not workflow_steps:
            raise click.BadParameter("No valid workflow steps found after validation")
        
        # Create workflow runner
        runner = _create_workflow_runner(use_ray, verbose, analysis_workers, 
                                       legacy_analysis_queue, log_level)
        
        # Handle Ray-specific configuration
        if use_ray and hasattr(runner, 'manager'):
            # Handle show-priorities option
            if show_priorities:
                click.echo("Current Queue Priorities:")
                priority_info = runner.manager.get_priority_info()
                for queue_type, priority in sorted(priority_info['queue_priorities'].items(), key=lambda x: x[1], reverse=True):
                    click.echo(f"  {queue_type}: {priority}")
                return
            
            # Configure queue priorities
            _configure_queue_priorities(runner, queue_priority)
            
            # Initialize Ray with specified CPU count if provided
            if ray_num_cpus is not None:
                _initialize_ray(ray_num_cpus)
        
        # Configure job deduplication
        if deduplicate_jobs:
            for job_type in deduplicate_jobs:
                job_type = job_type.strip()
                if job_type:  # Skip empty job types
                    try:
                        runner.manager.add_deduplication_job_type(job_type)
                    except Exception as e:
                        click.echo(f"Warning: Failed to enable deduplication for '{job_type}': {e}", err=True)
            valid_dedup_jobs = [j.strip() for j in deduplicate_jobs if j.strip()]
            if valid_dedup_jobs:
                click.echo(f"Job deduplication enabled for: {valid_dedup_jobs}")
        
        # Register handlers and command handlers
        _register_handlers(runner, legacy_analysis_queue, work_dir)
        _register_command_handlers(runner, command_map, legacy_analysis_queue)
        
        # For Ray workflow, reinitialize processors after handlers are registered
        if use_ray and hasattr(runner, 'manager'):
            try:
                runner.manager._reinitialize_processors()
            except Exception as e:
                click.echo(f"Warning: Failed to reinitialize Ray processors: {e}", err=True)
        
        # Create classifier function
        classifier_func = None
        if work_dir:
            classifier_func = _create_classifier_with_work_dir(work_dir, workflow_steps)
        
        # Display configuration
        _display_workflow_config(path, work_dir, workflow_steps, command_map, log_level, 
                               job_levels, deduplicate_jobs, legacy_analysis_queue, analysis_workers, 
                               no_process_existing, uses_simplified_format, workflow, use_ray, ray_num_cpus,
                               queue_priority)
        
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
        if 'runner' in locals() and hasattr(runner, 'manager'):
            try:
                runner.manager.stop(timeout=DEFAULT_TIMEOUT)
            except Exception as e:
                click.echo(f"Warning: Error during shutdown: {e}", err=True)
        click.echo("\nWorkflow stopped by user")
    except click.BadParameter as e:
        click.echo(f"Parameter error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Unexpected error: {e}", err=True)
        if verbose:
            import traceback
            click.echo(f"Traceback: {traceback.format_exc()}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main() 