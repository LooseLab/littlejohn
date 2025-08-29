"""
Command-line interface for robin.

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

from robin.workflow_simple import default_file_classifier, Job
from robin.analysis.bam_preprocessor import bam_preprocessing_handler
from robin.analysis.mgmt_analysis import mgmt_handler
from robin.analysis.cnv_analysis import cnv_handler
from robin.analysis.bed_conversion import bed_conversion_handler
from robin.analysis.sturgeon_analysis import sturgeon_handler
from robin.analysis.nanodx_analysis import nanodx_handler, pannanodx_handler
from robin.analysis.random_forest_analysis import random_forest_handler
from robin.analysis.target_analysis import target_handler
from robin.analysis.fusion_analysis import fusion_handler
from robin.logging_config import (
    configure_logging,
)


# Constants
VALID_LOG_LEVELS = {"DEBUG", "INFO", "WARNING", "ERROR"}
VALID_JOB_TYPES = {
    "preprocessing",
    "bed_conversion",
    "mgmt",
    "cnv",
    "target",
    "fusion",
    "sturgeon",
    "nanodx",
    "pannanodx",
    "random_forest",
}
DEFAULT_LOG_LEVEL = "ERROR"
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

from robin.reporting.sections.disclaimer_text import EXTENDED_DISCLAIMER_TEXT
# Disclaimer text for user acknowledgment
DISCLAIMER_TEXT = EXTENDED_DISCLAIMER_TEXT

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


def _get_user_acknowledgment() -> bool:
    """Display a disclaimer and require explicit user acknowledgment.

    Returns:
        bool: True if the user types 'I agree' exactly, otherwise False.
    """
    click.echo("\nDISCLAIMER:")
    click.echo(DISCLAIMER_TEXT)
    click.echo("\nTo proceed, please type 'I agree' (exactly as shown):")
    try:
        response = input().strip()
    except (KeyboardInterrupt, EOFError):
        click.echo("\nAcknowledgment interrupted. Exiting.", err=True)
        return False

    if response == "I agree":
        return True

    click.echo(
        "Incorrect acknowledgment. Please run the command again and type 'I agree' to acknowledge.",
        err=True,
    )
    return False


@click.group()
@click.version_option()
def main() -> None:
    """Robin now uses Little John - his second in command who kept the merry men in line."""
    pass


@main.command()
def list_job_types() -> None:
    """List all available job types organized by queue category."""
    if not _get_user_acknowledgment():
        sys.exit(1)

    job_types = {
        "preprocessing": ["preprocessing - Extract metadata from BAM files"],
        "bed_conversion": ["bed_conversion - Convert BAM files to BED format"],
        "mgmt": ["mgmt - MGMT methylation analysis"],
        "cnv": ["cnv - Copy number variation analysis"],
        "target": ["target - Target analysis"],
        "fusion": ["fusion - Fusion detection analysis"],
        "classification": [
            "sturgeon - Sturgeon classification analysis",
            "nanodx - NanoDX analysis",
            "pannanodx - PanNanoDX analysis",
        ],
        "slow": ["random_forest - Random Forest analysis"],
    }

    click.echo("Available job types in robin:\n")

    for queue, jobs in job_types.items():
        click.echo(f"{queue.upper()} QUEUE:")
        for job in jobs:
            click.echo(f"  • {job}")
        click.echo()

    click.echo("Usage examples:")
    click.echo(
        "  • Simplified format (recommended): 'mgmt,sturgeon' (bed_conversion auto-added)"
    )
    click.echo(
        "  • Full pipeline (simplified): 'mgmt,cnv,target,fusion,sturgeon,nanodx,pannanodx,random_forest' (bed_conversion auto-added)"
    )
    click.echo(
        "  • Legacy format with queue prefixes: 'preprocessing:bed_conversion,mgmt:mgmt,classification:sturgeon'"
    )
    click.echo(
        "\nNote: 'preprocessing' is automatically added as the first step if not specified."
    )
    click.echo(
        "Note: 'bed_conversion' is automatically added when needed for sturgeon, nanodx, pannanodx, or random_forest jobs."
    )
    click.echo(
        "Note: Each analysis type (mgmt, cnv, target, fusion) has its own queue for parallel processing."
    )
    click.echo(
        "Note: The system automatically determines the appropriate queue for each job type in simplified format."
    )


def _parse_job_log_levels(job_log_level: Tuple[str, ...]) -> Dict[str, str]:
    """Parse job-specific log level specifications."""
    job_levels = {}

    for job_level_spec in job_log_level:
        if ":" in job_level_spec:
            job_type, level = job_level_spec.split(":", 1)
            job_type = job_type.strip()
            level = level.strip().upper()

            if level not in VALID_LOG_LEVELS:
                click.echo(
                    f"Warning: Invalid log level '{level}' for job '{job_type}'. Valid levels: {', '.join(VALID_LOG_LEVELS)}",
                    err=True,
                )
                continue

            job_levels[job_type] = level
        else:
            click.echo(
                f"Warning: Invalid job log level format '{job_level_spec}'. Use format 'job_type:level'",
                err=True,
            )
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
                click.echo(
                    f"Warning: Empty job type or command in '{cmd_mapping}'", err=True
                )
                continue

            command_map[job_type] = command
        else:
            click.echo(
                f"Warning: Invalid command mapping format '{cmd_mapping}'. Use format 'job_type:command'",
                err=True,
            )
    return command_map


def _validate_workflow_steps(workflow_steps: List[str]) -> List[str]:
    """Validate and clean workflow steps."""
    cleaned_steps = []
    for step in workflow_steps:
        step = step.strip()
        if not step:
            continue

        if ":" in step:
            queue_type, job_type = step.split(":", 1)
            queue_type = queue_type.strip()
            job_type = job_type.strip()

            if job_type not in VALID_JOB_TYPES:
                click.echo(
                    f"Warning: Unknown job type '{job_type}' in step '{step}'", err=True
                )
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
    needs_bed_conversion = any(
        job in JOBS_REQUIRING_BED_CONVERSION for job in job_types
    )

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


def _create_ray_workflow_runner(
    verbose: bool,
    analysis_workers: int,
    legacy_analysis_queue: bool,
    log_level: str,
    preprocessing_workers: int = 1,
    bed_workers: int = 1,
) -> Any:
    """Create and configure Ray-based workflow runner (Ray Core)."""
    try:
        # Prefer Ray Core implementation
        import asyncio
        from robin import workflow_ray as wrn

        class _RayCoreWrapper:
            def __init__(self):
                self.manager = type("_DummyManager", (), {"get_priority_info": lambda _self: {}})()
                self.coordinator = None  # Will be set when workflow starts

            def register_handler(self, *args, **kwargs):
                # Handlers are wired inside workflow_ray; keep API parity
                return None

            def register_command_handler(self, *args, **kwargs):
                return None

            def submit_sample_job(self, sample_dir: str, job_type: str, sample_id: str = None, force_regenerate: bool = False) -> bool:
                """Submit a job for an existing sample directory using Ray workflow."""
                try:
                    # Try to get coordinator with retries
                    import time
                    from robin import workflow_ray as wrn
                    
                    max_retries = 5
                    retry_delay = 1.0
                    
                    for attempt in range(max_retries):
                        # Try to get coordinator
                        if self.coordinator is None:
                            try:
                                self.coordinator = wrn.get_coordinator_sync()
                            except Exception:
                                pass
                        
                        if self.coordinator is not None:
                            break
                        
                        if attempt < max_retries - 1:
                            print(f"[Ray] Coordinator not ready, retrying in {retry_delay}s... (attempt {attempt + 1}/{max_retries})")
                            time.sleep(retry_delay)
                            retry_delay *= 1.5  # Exponential backoff
                    
                    if self.coordinator is None:
                        print(f"[Ray] No coordinator available for {job_type} job after {max_retries} attempts")
                        return False
                    
                    # Submit the job using the coordinator
                    # Ray remote calls are synchronous from the caller's perspective
                    result = self.coordinator.submit_sample_job.remote(sample_dir, job_type, sample_id, force_regenerate)
                    
                    # Wait for the result (this is the actual async part)
                    import ray
                    try:
                        final_result = ray.get(result, timeout=30.0)
                        print(f"[Ray] Successfully submitted {job_type} job for sample {sample_id or 'unknown'}")
                        return final_result
                    except Exception as e:
                        print(f"[Ray] Job submission failed: {e}")
                        return False
                        
                except Exception as e:
                    print(f"[Ray] Failed to submit {job_type} job: {e}")
                    return False

            def run_workflow(
                self,
                path: Path,
                workflow_steps: List[str],
                no_process_existing: bool,
                no_progress: bool,
                no_watch: bool,
                log_level_local: str,
                analysis_workers_local: int,
                preset: Optional[str] = None,
            ) -> None:
                # Store reference to self for coordinator access
                self_ref = self
                
                async def _run_with_coordinator():
                    import asyncio
                    from robin import workflow_ray as wrn
                    
                    # Run the workflow first to create the coordinator
                    await wrn.run(
                        plan=workflow_steps,
                        paths=[str(path)],
                        analysis_workers=analysis_workers_local,
                        process_existing=not no_process_existing,
                        monitor=not no_progress,
                        watch=(not no_watch),
                        patterns=["*.bam"],
                        ignore_patterns=None,
                        recursive=True,
                        work_dir=None,
                        log_level=log_level_local,
                        preset=preset,
                        workflow_runner=self_ref,
                    )
                    
                    # After the workflow starts, try to get the coordinator reference
                    try:
                        coord = wrn.get_coordinator_sync()
                        if coord:
                            self_ref.coordinator = coord
                            print(f"[Ray] Coordinator reference set for job submission")
                    except Exception as e:
                        print(f"[Ray] Warning: Could not set coordinator reference: {e}")
                
                asyncio.run(_run_with_coordinator())

        click.echo("Using Ray Core distributed workflow (workflow_ray)")
        return _RayCoreWrapper()
    except ImportError as e:
        click.echo(
            f"Warning: Ray Core not available ({e}). Falling back to threading-based workflow.",
            err=True,
        )
        return None


def _initialize_ray(num_cpus: Optional[int]) -> None:
    """Initialize Ray with specified CPU count."""
    if num_cpus is not None:
        if num_cpus <= 0:
            click.echo(
                f"Warning: Invalid CPU count {num_cpus}. Must be positive.", err=True
            )
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
                # Bind dashboard to 0.0.0.0 when supported so it's reachable off-host
                try:
                    ray.init(
                        num_cpus=num_cpus,
                        ignore_reinit_error=True,
                        include_dashboard=True,
                        dashboard_host=os.environ.get("RAY_DASHBOARD_HOST", "0.0.0.0"),
                    )
                except TypeError:
                    # Older Ray versions may not support dashboard args
                    ray.init(num_cpus=num_cpus, ignore_reinit_error=True)
                click.echo(f"Ray initialized with {num_cpus} CPUs")
            except Exception as e:
                click.echo(f"Warning: Failed to initialize Ray: {e}", err=True)
        else:
            click.echo(
                f"Ray already initialized with {ray.available_resources().get('CPU', 'unknown')} CPUs"
            )


def _configure_queue_priorities(runner: Any, queue_priority: Tuple[str, ...]) -> None:
    """Configure queue priorities for Ray workflow."""
    if not queue_priority:
        return

    click.echo("Configuring queue priorities:")
    for priority_spec in queue_priority:
        if ":" in priority_spec:
            queue_type, priority_str = priority_spec.split(":", 1)
            queue_type = queue_type.strip()
            priority_str = priority_str.strip()

            try:
                priority = int(priority_str)
                if priority < 0:
                    click.echo(
                        f"Warning: Priority {priority} for queue '{queue_type}' is negative. Using 0 instead.",
                        err=True,
                    )
                    priority = 0

                runner.manager.set_queue_priority(queue_type, priority)
                click.echo(f"  {queue_type}: {priority}")
            except ValueError:
                click.echo(
                    f"Warning: Invalid priority value '{priority_str}' for queue '{queue_type}'. Must be an integer.",
                    err=True,
                )
        else:
            click.echo(
                f"Warning: Invalid priority specification '{priority_spec}'. Use format 'queue:priority'",
                err=True,
            )


def _create_workflow_runner(
    use_ray: bool,
    verbose: bool,
    analysis_workers: int,
    legacy_analysis_queue: bool,
    log_level: str,
    preprocessing_workers: int = 1,
    bed_workers: int = 1,
) -> Any:
    """Create the appropriate workflow runner based on configuration."""
    if analysis_workers < 1:
        click.echo(
            f"Warning: Invalid analysis_workers value {analysis_workers}. Using {DEFAULT_ANALYSIS_WORKERS} instead.",
            err=True,
        )
        analysis_workers = DEFAULT_ANALYSIS_WORKERS

    if use_ray:
        runner = _create_ray_workflow_runner(
            verbose,
            analysis_workers,
            legacy_analysis_queue,
            log_level,
            preprocessing_workers=preprocessing_workers,
            bed_workers=bed_workers,
        )
        if runner is None:
            # Fallback to threading-based workflow
            from robin.workflow_simple import WorkflowRunner

            runner = WorkflowRunner(
                verbose=verbose,
                analysis_workers=analysis_workers,
                use_separate_analysis_queues=not legacy_analysis_queue,
                preprocessing_workers=preprocessing_workers,
                bed_workers=bed_workers,
            )
        return runner
    else:
        from robin.workflow_simple import WorkflowRunner

        return WorkflowRunner(
            verbose=verbose,
            analysis_workers=analysis_workers,
            use_separate_analysis_queues=not legacy_analysis_queue,
            preprocessing_workers=preprocessing_workers,
            bed_workers=bed_workers,
        )


def _register_handlers(
    runner: Any, legacy_analysis_queue: bool, work_dir: Optional[Path]
) -> None:
    """Register all workflow handlers with the runner."""
    for (
        queue_type,
        job_type,
        handler_func,
        legacy_queue_type,
        needs_work_dir,
    ) in HANDLER_CONFIGS:
        # Determine the actual queue type based on legacy mode
        actual_queue_type = (
            legacy_queue_type
            if legacy_analysis_queue and legacy_queue_type
            else queue_type
        )

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
            click.echo(
                f"Warning: Failed to register handler for {queue_type}:{job_type}: {e}",
                err=True,
            )


def _register_command_handlers(
    runner: Any, command_map: Dict[str, str], legacy_analysis_queue: bool
) -> None:
    """Register command handlers with the runner."""
    for job_type, command in command_map.items():
        if ":" in job_type:
            queue_type, actual_job_type = job_type.split(":", 1)

            # Handle legacy analysis queue mapping
            if legacy_analysis_queue and queue_type in [
                "mgmt",
                "cnv",
                "target",
                "fusion",
            ]:
                actual_queue_type = "analysis"
            else:
                actual_queue_type = queue_type

            try:
                runner.register_command_handler(
                    actual_queue_type, actual_job_type, command
                )
            except Exception as e:
                click.echo(
                    f"Warning: Failed to register command handler for {queue_type}:{actual_job_type}: {e}",
                    err=True,
                )
        else:
            # Default to preprocessing queue
            try:
                runner.register_command_handler("preprocessing", job_type, command)
            except Exception as e:
                click.echo(
                    f"Warning: Failed to register command handler for preprocessing:{job_type}: {e}",
                    err=True,
                )


def _create_classifier_with_work_dir(work_dir: Path, workflow_steps: List[str]):
    """Create a classifier function that includes work directory in job context."""

    def classifier_with_work_dir(filepath: str) -> List[Job]:
        jobs = default_file_classifier(filepath, workflow_steps)
        for job in jobs:
            job.context.add_metadata("work_dir", str(work_dir))
        return jobs

    return classifier_with_work_dir


def _validate_inputs(
    path: Path, workflow: str, analysis_workers: int, ray_num_cpus: Optional[int]
) -> None:
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


def _display_workflow_config(
    path: Path,
    work_dir: Optional[Path],
    workflow_steps: List[str],
    command_map: dict,
    log_level: str,
    job_levels: dict,
    deduplicate_jobs: tuple,
    legacy_analysis_queue: bool,
    analysis_workers: int,
    no_process_existing: bool,
    uses_simplified_format: bool = False,
    original_workflow: str = "",
    use_ray: bool = False,
    ray_num_cpus: Optional[int] = None,
    queue_priority: tuple = (),
) -> None:
    """Display workflow configuration information."""
    if no_process_existing:
        click.echo(
            f"Starting workflow on {path} for BAM files (skipping existing files)..."
        )
    else:
        click.echo(
            f"Starting workflow on {path} for BAM files (will process existing files first)..."
        )

    if work_dir:
        click.echo(f"Output directory: {work_dir}")

    if uses_simplified_format:
        # Extract job types from the final workflow steps
        final_job_types = [step.split(":")[1] for step in workflow_steps if ":" in step]
        click.echo(f"Workflow plan (simplified): {final_job_types}")
        click.echo(f"Auto-assigned queues: {workflow_steps}")

        # Check if bed_conversion was auto-added
        original_jobs = [step.strip() for step in original_workflow.split(",")]
        if (
            "bed_conversion" in final_job_types
            and "bed_conversion" not in original_jobs
        ):
            click.echo(
                "Note: bed_conversion was automatically added as it's required for other jobs"
            )
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
        click.echo("Distributed computing: Ray (experimental)")
        if ray_num_cpus:
            click.echo(f"  - Ray CPUs: {ray_num_cpus}")
        else:
            click.echo("  - Ray CPUs: auto-detect")
        click.echo(
            f"  - Analysis workers per type: {analysis_workers} (distributed across Ray cluster)"
        )
        click.echo(f"  - Log level: {log_level} (applied to all Ray actors)")

        # Display priority configuration
        if queue_priority:
            click.echo(f"  - Queue priorities: {list(queue_priority)}")
    else:
        click.echo("Distributed computing: Disabled (using threading)")
        click.echo("Worker configuration:")
        if legacy_analysis_queue:
            click.echo(
                "  - Analysis queue mode: Legacy (single queue for all analysis types)"
            )
            click.echo(f"  - Analysis workers: {analysis_workers}")
        else:
            click.echo("  - Analysis queue mode: Separate queues per analysis type")
            click.echo(
                f"  - Analysis workers per type: {analysis_workers} (MGMT, CNV, Target, Fusion each get {analysis_workers} workers)"
            )

        click.echo("  - Other queues: 1 worker each (fixed)")

    click.echo("Press Ctrl+C to stop")
    click.echo("Running The Workflow!")


@main.command()
@click.argument("path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--workflow",
    "-w",
    required=True,
    help="Workflow plan. Can be specified in two formats:\n1. With queue prefixes: 'preprocessing:bed_conversion,mgmt:mgmt,classification:sturgeon'\n2. Simplified (auto-queue): 'mgmt,sturgeon' - system automatically determines appropriate queue for each job type and adds bed_conversion when needed",
)
@click.option(
    "--commands",
    "-c",
    multiple=True,
    help="Command mappings (e.g., 'index:samtools index {file}')",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Enable verbose output and detailed error traces",
)
@click.option(
    "--no-process-existing",
    is_flag=True,
    help="Skip processing existing files, only watch for new changes",
)
@click.option(
    "--work-dir",
    "-d",
    type=click.Path(path_type=Path),
    help="Base output directory for analysis results",
)
@click.option(
    "--log-level",
    default=DEFAULT_LOG_LEVEL,
    type=click.Choice(list(VALID_LOG_LEVELS)),
    help=f"Global log level (default: {DEFAULT_LOG_LEVEL})",
)
@click.option(
    "--job-log-level",
    multiple=True,
    help=f"Set log level for specific job (e.g., 'preprocessing:DEBUG', 'mgmt:WARNING'). Valid levels: {', '.join(VALID_LOG_LEVELS)}",
)
@click.option(
    "--deduplicate-jobs",
    multiple=True,
    help="Job types to deduplicate by sample ID (e.g., 'sturgeon', 'mgmt'). Jobs of these types will only run once per sample, even if multiple upstream jobs complete simultaneously.",
)
@click.option(
    "--no-progress", is_flag=True, help="Disable progress bars for file processing"
)
@click.option(
    "--analysis-workers",
    type=int,
    default=DEFAULT_ANALYSIS_WORKERS,
    help=f"Number of analysis workers per analysis queue (default: {DEFAULT_ANALYSIS_WORKERS})",
)
@click.option(
    "--preprocessing-workers",
    type=int,
    default=1,
    help="Number of preprocessing workers (default: 1)",
)
@click.option(
    "--bed-workers",
    type=int,
    default=1,
    help="Number of bed_conversion workers (default: 1)",
)
@click.option(
    "--legacy-analysis-queue",
    is_flag=True,
    help="Use legacy single analysis queue instead of separate queues per analysis type",
)
@click.option(
    "--use-ray/--no-use-ray",
    default=True,
    help="Enable Ray distributed computing (default: on). Disable with --no-use-ray.",
)
@click.option(
    "--use-ray-core",
    is_flag=True,
    help="Deprecated: Ray Core is now the default when --use-ray is provided.",
)
@click.option(
    "--ray-num-cpus",
    type=int,
    default=None,
    help="Number of CPUs to use for Ray (default: auto-detect). Only used when --use-ray is specified.",
)
@click.option(
    "--queue-priority",
    multiple=True,
    help="Set priority for a queue (e.g., 'preprocessing:10', 'bed_conversion:9', 'mgmt:5'). Only used when --use-ray is specified.",
)
@click.option(
    "--show-priorities",
    is_flag=True,
    help="Show current queue priorities and exit. Only used when --use-ray is specified.",
)
@click.option(
    "--no-watch",
    is_flag=True,
    help="Do not watch directories for new files (default: watch enabled).",
)
@click.option(
    "--with-gui/--no-gui",
    default=True,
    help="Launch NiceGUI workflow monitor (default: on). Disable with --no-gui.",
)
@click.option(
    "--gui-host",
    default="0.0.0.0",
    show_default=True,
    help="Host interface for the GUI server (e.g., 0.0.0.0 to listen on all interfaces).",
)
@click.option(
    "--gui-port",
    type=int,
    default=8081,
    show_default=True,
    help="Port for the GUI server.",
)
@click.option(
    "--preset",
    type=click.Choice(["p2i", "standard", "high"]),
    default="standard",
    help="Execution preset for Ray Core: 'p2i' (2 CPUs cap, 4 grouped actors, concurrency 1), 'standard' (default; 4 CPUs cap, grouped actors, analysis uses analysis_workers), 'high' (per-job-type actors).",
)
def workflow(
    path: Path,
    workflow: str,
    commands: tuple[str, ...],
    verbose: bool,
    no_process_existing: bool,
    work_dir: Optional[Path],
    log_level: str,
    job_log_level: tuple[str, ...],
    deduplicate_jobs: tuple[str, ...],
    no_progress: bool,
    analysis_workers: int,
    preprocessing_workers: int,
    bed_workers: int,
    legacy_analysis_queue: bool,
    use_ray: bool,
    use_ray_core: bool,
    ray_num_cpus: Optional[int],
    queue_priority: tuple[str, ...],
    show_priorities: bool,
    with_gui: bool,
    gui_host: str,
    gui_port: int,
    no_watch: bool,
    preset: Optional[str],
) -> None:
    """Run various operations on BAM files in a directory. Preprocessing is automatically included as the first step."""
    try:
        # Require user acknowledgment before proceeding
        if not _get_user_acknowledgment():
            sys.exit(1)

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
        uses_simplified_format = all(":" not in step for step in workflow_steps)

        if uses_simplified_format:
            # Convert simplified format to standard format with automatic queue assignment
            workflow_steps = _convert_simplified_workflow(workflow_steps)

        # Ensure preprocessing is the first step to extract metadata
        workflow_steps = _ensure_preprocessing_step(workflow_steps)

        # Validate that we have at least some workflow steps
        if not workflow_steps:
            raise click.BadParameter("No valid workflow steps found after validation")

        # If using Ray and CPU count specified, initialize Ray BEFORE creating runner
        # so the runner respects the requested CPU resources
        if use_ray and ray_num_cpus is not None:
            _initialize_ray(ray_num_cpus)

        # Ray Core engine path (default when --use-ray is used)
        if use_ray:
            # Ensure Ray is initialized if CPUs not specified
            try:
                import ray

                if not ray.is_initialized():
                    # Apply preset CPU caps for Ray Core if provided
                    init_kwargs = {
                        "ignore_reinit_error": True,
                        "include_dashboard": True,
                        "dashboard_host": os.environ.get(
                            "RAY_DASHBOARD_HOST", "0.0.0.0"
                        ),
                    }
                    if preset in {"p2i", "standard"}:
                        init_kwargs["num_cpus"] = 2 if preset == "p2i" else 4
                    try:
                        ray.init(**init_kwargs)
                    except TypeError:
                        # Older Ray versions may not support dashboard args
                        init_kwargs.pop("include_dashboard", None)
                        init_kwargs.pop("dashboard_host", None)
                        ray.init(**init_kwargs)
            except Exception:
                pass

            # Display configuration
            _display_workflow_config(
                path,
                work_dir,
                workflow_steps,
                command_map,
                log_level,
                job_levels,
                deduplicate_jobs,
                legacy_analysis_queue,
                analysis_workers,
                no_process_existing,
                uses_simplified_format,
                workflow,
                True,
                ray_num_cpus,
                queue_priority,
            )

            # Create workflow runner for Ray workflow (needed for GUI integration)
            runner = _create_workflow_runner(
                True,  # use_ray=True
                verbose,
                analysis_workers,
                legacy_analysis_queue,
                log_level,
                preprocessing_workers=preprocessing_workers,
                bed_workers=bed_workers,
            )

            # Run Ray Core implementation
            try:
                import asyncio
                from robin import workflow_ray as wrn

                asyncio.run(
                    wrn.run(
                        plan=workflow_steps,
                        paths=[str(path)],
                        analysis_workers=analysis_workers,
                        process_existing=not no_process_existing,
                        monitor=not no_progress,
                        watch=(not no_watch),
                        patterns=["*.bam"],
                        ignore_patterns=None,
                        recursive=True,
                        work_dir=str(work_dir) if work_dir else None,
                        log_level=log_level,
                        preset=preset,
                        workflow_runner=runner,
                    )
                )
            except KeyboardInterrupt:
                print("Stopping workflow...")
                # Attempt to shutdown GUI server if running
                try:
                    from robin.gui.app import get_gui_launcher as _get  # type: ignore

                    gl = _get()
                    if gl is not None:
                        try:
                            from nicegui import app as ng_app  # type: ignore

                            ng_app.shutdown()
                        except Exception:
                            pass
                except Exception:
                    pass
            return

        # Create workflow runner
        runner = _create_workflow_runner(
            use_ray,
            verbose,
            analysis_workers,
            legacy_analysis_queue,
            log_level,
            preprocessing_workers=preprocessing_workers,
            bed_workers=bed_workers,
        )

        # Handle Ray-specific configuration
        if use_ray and hasattr(runner, "manager"):
            # Handle show-priorities option
            if show_priorities:
                click.echo("Current Queue Priorities:")
                priority_info = runner.manager.get_priority_info()
                for queue_type, priority in sorted(
                    priority_info["queue_priorities"].items(),
                    key=lambda x: x[1],
                    reverse=True,
                ):
                    click.echo(f"  {queue_type}: {priority}")
                return

            # Configure queue priorities
            _configure_queue_priorities(runner, queue_priority)

            # Ray already initialized above when ray_num_cpus is provided

        # Configure job deduplication
        if deduplicate_jobs:
            for job_type in deduplicate_jobs:
                job_type = job_type.strip()
                if job_type:  # Skip empty job types
                    try:
                        runner.manager.add_deduplication_job_type(job_type)
                    except Exception as e:
                        click.echo(
                            f"Warning: Failed to enable deduplication for '{job_type}': {e}",
                            err=True,
                        )
            valid_dedup_jobs = [j.strip() for j in deduplicate_jobs if j.strip()]
            if valid_dedup_jobs:
                click.echo(f"Job deduplication enabled for: {valid_dedup_jobs}")

        # Register handlers and command handlers
        _register_handlers(runner, legacy_analysis_queue, work_dir)
        _register_command_handlers(runner, command_map, legacy_analysis_queue)

        # For Ray workflow, reinitialize processors after handlers are registered
        if use_ray and hasattr(runner, "manager"):
            try:
                runner.manager._reinitialize_processors()
            except Exception as e:
                click.echo(
                    f"Warning: Failed to reinitialize Ray processors: {e}", err=True
                )

        # Create classifier function
        classifier_func = None
        if work_dir:
            classifier_func = _create_classifier_with_work_dir(work_dir, workflow_steps)

        # Display configuration
        _display_workflow_config(
            path,
            work_dir,
            workflow_steps,
            command_map,
            log_level,
            job_levels,
            deduplicate_jobs,
            legacy_analysis_queue,
            analysis_workers,
            no_process_existing,
            uses_simplified_format,
            workflow,
            use_ray,
            ray_num_cpus,
            queue_priority,
        )

        # Launch GUI if requested
        gui_launcher = None
        if with_gui:
            try:
                print("🚀 Launching NiceGUI workflow monitor with vertical tabs...")

                # Import the new GUI launcher
                try:
                    from robin.gui.app import launch_gui

                    # Launch GUI using the new launcher FIRST so the global sender is ready
                    gui_launcher = launch_gui(
                        host=gui_host,
                        port=gui_port,
                        show=False,
                        workflow_runner=runner,
                        workflow_steps=workflow_steps,
                        monitored_directory=str(work_dir) if work_dir else str(path),
                    )

                    # Now install workflow hooks for real-time monitoring
                    try:
                        from robin.workflow_hooks import install_workflow_hooks

                        install_workflow_hooks(runner, workflow_steps, str(path))
                        click.echo(
                            "✅ Workflow state hooks installed for real-time monitoring"
                        )
                    except Exception as e:
                        click.echo(
                            f"⚠️  Failed to install workflow hooks: {e}. GUI will show static information only."
                        )

                    base_url = f"http://{gui_host}:{gui_port}"
                    print(f"✅ GUI launched successfully on {base_url}")
                    print(f"   🏠 Welcome page: {base_url}/")
                    print(f"   📊 Workflow monitor: {base_url}/robin")
                    print(f"   📋 Sample tracking: {base_url}/live_data")
                    print(f"   🔬 Individual samples: {base_url}/live_data/sampleID/")
                    print(
                        "   Open your browser to monitor the workflow with the new navigation structure"
                    )
                    print(
                        "   The GUI will run in the background while the workflow executes"
                    )

                except ImportError as e:
                    click.echo(f"❌ Failed to import GUI launcher: {e}")
                    click.echo("   Please ensure the gui_launcher module is available")

            except Exception as e:
                click.echo(
                    f"⚠️  Failed to launch GUI: {e}. Continuing with workflow only.",
                    err=True,
                )

        # Run the workflow
        runner.run_workflow(
            watch_dir=str(path),
            workflow_plan=workflow_steps,
            recursive=True,
            patterns=["*.bam"],
            ignore_patterns=None,
            classifier_func=classifier_func,
            process_existing=not no_process_existing,
            show_progress=not no_progress,
        )

    except KeyboardInterrupt:
        print("Stopping workflow...")
        if "runner" in locals() and hasattr(runner, "manager"):
            try:
                runner.manager.stop(timeout=DEFAULT_TIMEOUT)
            except Exception as e:
                click.echo(f"Warning: Error during shutdown: {e}", err=True)
        # Attempt to shutdown GUI server if running
        try:
            from robin.gui.app import get_gui_launcher as _get  # type: ignore

            gl = _get()
            if gl is not None:
                try:
                    from nicegui import app as ng_app  # type: ignore

                    ng_app.shutdown()
                except Exception:
                    pass
        except Exception:
            pass
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
