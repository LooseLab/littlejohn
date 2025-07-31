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
    """LittleJohn - A Python CLI tool with file watching capabilities."""
    pass


@main.command()
def list_job_types() -> None:
    """List all available job types organized by queue category."""
    job_types = {
        "preprocessing": [
            "preprocessing - Extract metadata from BAM files",
            "bed_conversion - Convert BAM files to BED format"
        ],
        "analysis": [
            "mgmt - MGMT methylation analysis",
            "cnv - Copy number variation analysis", 
            "target - Target analysis",
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
    click.echo("  • Workflow with preprocessing and analysis: 'preprocessing:preprocessing,analysis:mgmt'")
    click.echo("  • Full pipeline: 'preprocessing:preprocessing,analysis:mgmt,classification:sturgeon'")
    click.echo("  • Multiple analyses: 'preprocessing:preprocessing,analysis:mgmt,analysis:cnv'")
    click.echo("\nNote: 'preprocessing:preprocessing' is automatically added as the first step if not specified.")


@main.command()
@click.argument("path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--command",
    "-c",
    help="Command to run when files change",
)
@click.option(
    "--verbose", "-v", is_flag=True, help="Enable verbose output"
)
@click.option(
    "--no-process-existing", is_flag=True, help="Skip processing existing files, only watch for new changes"
)
@click.option(
    "--no-progress", is_flag=True, help="Disable progress bars for file processing"
)
def watch(
    path: Path,
    command: Optional[str],
    verbose: bool,
    no_process_existing: bool,
    no_progress: bool,
) -> None:
    """Watch a directory for BAM file changes."""
    try:
        watcher = FileWatcher(
            path=path,
            recursive=True,
            patterns=["*.bam"],
            ignore_patterns=[],
            command=command,
            verbose=verbose,
            show_progress=not no_progress,
        )
        
        if no_process_existing:
            click.echo(f"Watching {path} for BAM file changes (skipping existing files)...")
        else:
            click.echo(f"Watching {path} for BAM file changes (will process existing files first)...")
        
        click.echo("Press Ctrl+C to stop")
        watcher.start(process_existing=not no_process_existing)
    except KeyboardInterrupt:
        click.echo("\nStopping file watcher...")
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@main.command()
@click.argument("path", type=click.Path(exists=True, path_type=Path))
def list_files(path: Path) -> None:
    """List BAM files in a directory recursively."""
    try:
        if path.is_file():
            if path.suffix == '.bam':
                click.echo(str(path))
            return

        # Find all BAM files recursively
        bam_files = list(path.rglob("*.bam"))
        
        if not bam_files:
            click.echo(f"No BAM files found in {path}")
            return

        # Sort files for consistent output
        for file_path in sorted(bam_files):
            click.echo(str(file_path))

        # Show summary
        click.echo(f"\nFound {len(bam_files)} BAM file(s)")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@main.command()
@click.argument("path", type=click.Path(exists=True, path_type=Path))
def info(path: Path) -> None:
    """Display information about a file or directory."""
    try:
        stat = path.stat()
        
        click.echo(f"Path: {path}")
        click.echo(f"Type: {'Directory' if path.is_dir() else 'File'}")
        click.echo(f"Size: {stat.st_size} bytes")
        click.echo(f"Modified: {stat.st_mtime}")
        click.echo(f"Permissions: {oct(stat.st_mode)[-3:]}")
        
        if path.is_dir():
            files = list(path.iterdir())
            click.echo(f"Contents: {len(files)} items")
            
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@main.command()
@click.argument("path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--workflow", "-w", required=True, help="Workflow plan (e.g., 'preprocessing:bed_conversion,analysis:mgmt,classification:sturgeon'). BAM preprocessing is automatically added as the first step to extract metadata."
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
def workflow(path: Path, workflow: str, commands: tuple[str, ...], verbose: bool, no_process_existing: bool, work_dir: Optional[Path], log_level: str, job_log_level: tuple[str, ...], deduplicate_jobs: tuple[str, ...], no_progress: bool) -> None:
    """Run an async workflow on BAM files in a directory. Preprocessing is automatically included as the first step."""
    try:
        
        # Configure logging
        job_levels = {}
        for job_level_spec in job_log_level:
            if ':' in job_level_spec:
                job_type, level = job_level_spec.split(':', 1)
                job_levels[job_type.strip()] = level.strip().upper()
        
        configure_logging(global_level=log_level, job_levels=job_levels)
        
        # Parse workflow plan
        workflow_steps = [step.strip() for step in workflow.split(",")]
        
        # Ensure preprocessing is the first step to extract metadata
        # If the first step is not preprocessing, add it
        if not workflow_steps or not workflow_steps[0].endswith(":preprocessing"):
            workflow_steps.insert(0, "preprocessing:preprocessing")
        
        # Parse command mappings
        command_map = {}
        for cmd_mapping in commands:
            if ":" in cmd_mapping:
                job_type, command = cmd_mapping.split(":", 1)
                command_map[job_type.strip()] = command.strip()
        
        # Create workflow runner
        runner = WorkflowRunner(verbose=verbose)
        
        # Configure job deduplication
        if deduplicate_jobs:
            for job_type in deduplicate_jobs:
                runner.manager.add_deduplication_job_type(job_type.strip())
            click.echo(f"Job deduplication enabled for: {list(deduplicate_jobs)}")
        
        # Register the mandatory preprocessing handler
        runner.register_handler("preprocessing", "preprocessing", bam_preprocessing_handler)
        
        # Register the BED conversion handler with work directory if specified
        if work_dir:
            def bed_conversion_handler_with_work_dir(job):
                return bed_conversion_handler(job, work_dir=str(work_dir))
            runner.register_handler("preprocessing", "bed_conversion", bed_conversion_handler_with_work_dir)
        else:
            runner.register_handler("preprocessing", "bed_conversion", bed_conversion_handler)
        
        # Register the MGMT analysis handler with work directory if specified
        if work_dir:
            def mgmt_handler_with_work_dir(job):
                return mgmt_handler(job, work_dir=str(work_dir))
            runner.register_handler("analysis", "mgmt", mgmt_handler_with_work_dir)
        else:
            runner.register_handler("analysis", "mgmt", mgmt_handler)
        
        # Register the CNV analysis handler with work directory if specified
        if work_dir:
            def cnv_handler_with_work_dir(job):
                return cnv_handler(job, work_dir=str(work_dir))
            runner.register_handler("analysis", "cnv", cnv_handler_with_work_dir)
        else:
            runner.register_handler("analysis", "cnv", cnv_handler)
        
        # Register the Sturgeon analysis handler with work directory if specified
        if work_dir:
            def sturgeon_handler_with_work_dir(job):
                return sturgeon_handler(job, work_dir=str(work_dir))
            runner.register_handler("classification", "sturgeon", sturgeon_handler_with_work_dir)
        else:
            runner.register_handler("classification", "sturgeon", sturgeon_handler)
        
        # Register the NanoDX analysis handler with work directory if specified
        if work_dir:
            def nanodx_handler_with_work_dir(job):
                return nanodx_handler(job, work_dir=str(work_dir))
            runner.register_handler("classification", "nanodx", nanodx_handler_with_work_dir)
        else:
            runner.register_handler("classification", "nanodx", nanodx_handler)
        
        # Register the PanNanoDX analysis handler with work directory if specified
        if work_dir:
            def pannanodx_handler_with_work_dir(job):
                return pannanodx_handler(job, work_dir=str(work_dir))
            runner.register_handler("classification", "pannanodx", pannanodx_handler_with_work_dir)
        else:
            runner.register_handler("classification", "pannanodx", pannanodx_handler)
        
        # Register the Random Forest analysis handler with work directory if specified
        if work_dir:
            def random_forest_handler_with_work_dir(job):
                return random_forest_handler(job, work_dir=str(work_dir))
            runner.register_handler("slow", "random_forest", random_forest_handler_with_work_dir)
        else:
            runner.register_handler("slow", "random_forest", random_forest_handler)
        
        # Register the Target analysis handler with work directory if specified
        if work_dir:
            def target_handler_with_work_dir(job):
                return target_handler(job, work_dir=str(work_dir))
            runner.register_handler("analysis", "target", target_handler_with_work_dir)
        else:
            runner.register_handler("analysis", "target", target_handler)
        
        # Register the Fusion analysis handler with work directory if specified
        if work_dir:
            def fusion_handler_with_work_dir(job):
                return fusion_handler(job, work_dir=str(work_dir))
            runner.register_handler("analysis", "fusion", fusion_handler_with_work_dir)
        else:
            runner.register_handler("analysis", "fusion", fusion_handler)
        
        # Store work directory in job context for preprocessing
        if work_dir:
            # Create a custom classifier that includes work directory
            def classifier_with_work_dir(filepath: str) -> List[Job]:
                jobs = default_file_classifier(filepath, workflow_steps)
                for job in jobs:
                    job.context.add_metadata('work_dir', str(work_dir))
                return jobs
            
            # Use the custom classifier
            classifier_func = classifier_with_work_dir
        else:
            classifier_func = None
        
        # Register command handlers
        for job_type, command in command_map.items():
            if job_type.startswith("preprocessing:"):
                queue_type, actual_job_type = job_type.split(":", 1)
                runner.register_command_handler("preprocessing", actual_job_type, command)
            elif job_type.startswith("analysis:"):
                queue_type, actual_job_type = job_type.split(":", 1)
                runner.register_command_handler("analysis", actual_job_type, command)
            elif job_type.startswith("classification:"):
                queue_type, actual_job_type = job_type.split(":", 1)
                runner.register_command_handler("classification", actual_job_type, command)
            elif job_type.startswith("slow:"):
                queue_type, actual_job_type = job_type.split(":", 1)
                runner.register_command_handler("slow", actual_job_type, command)
            else:
                # Default to preprocessing queue
                runner.register_command_handler("preprocessing", job_type, command)
        
        if no_process_existing:
            click.echo(f"Starting workflow on {path} for BAM files (skipping existing files)...")
        else:
            click.echo(f"Starting workflow on {path} for BAM files (will process existing files first)...")
        if work_dir:
            click.echo(f"Output directory: {work_dir}")
        click.echo(f"Workflow plan: {workflow_steps}")
        click.echo(f"Commands: {command_map}")
        click.echo(f"Log level: {log_level}")
        if job_levels:
            click.echo(f"Job log levels: {job_levels}")
        if deduplicate_jobs:
            click.echo(f"Job deduplication: {list(deduplicate_jobs)}")
        click.echo("Press Ctrl+C to stop")
        
        
        print("!!!!!!!!!! Running The Workflow !!!!!!!!!!")
        
        
        # Run the workflow
        runner.run_workflow(
            watch_dir=str(path),
            workflow_plan=workflow_steps,
            recursive=True,
            patterns=["*.bam"],
            ignore_patterns=None,
            classifier_func=classifier_func if 'classifier_func' in locals() else None,
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


if __name__ == "__main__":
    main() 