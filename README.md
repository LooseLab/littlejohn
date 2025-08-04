# LittleJohn

A Python CLI tool for automated BAM file processing and bioinformatics analysis workflows.

## Overview

LittleJohn is a specialized bioinformatics workflow engine designed for processing BAM files through automated pipelines. It provides a robust, multi-threaded system for running complex analysis workflows with built-in support for methylation analysis, copy number variation detection, fusion detection, and various classification algorithms.

## Key Features

- **BAM-Focused Processing**: Specialized for BAM file workflows with automatic metadata extraction
- **Multi-Queue Workflow System**: Separate queues for preprocessing, analysis, classification, and slow jobs
- **Built-in Analysis Modules**: 
  - MGMT methylation analysis
  - Copy number variation (CNV) analysis
  - Fusion detection analysis
  - Target analysis
  - Sturgeon classification
  - NanoDX and PanNanoDX analysis
  - Random Forest classification
- **Automatic Preprocessing**: Metadata extraction and BED conversion
- **Job Deduplication**: Prevents redundant processing per sample
- **Configurable Logging**: Global and per-job log level control
- **Progress Tracking**: Real-time progress bars and worker monitoring
- **File Watching**: Monitor directories for new BAM files
- **Master CSV Tracking**: Automatic sample-level statistics aggregation

## Installation

### From Source

1. Clone the repository:
```bash
git clone https://github.com/yourusername/littlejohn.git
cd littlejohn
```

2. Install in development mode:
```bash
pip install -e .
```

### Development Setup

For development, install with dev dependencies:

```bash
pip install -e ".[dev]"
```

## Available Commands

### `list-job-types`
List all available job types organized by queue category:

```bash
littlejohn list-job-types
```

**Available Job Types:**
- **Preprocessing Queue**: `preprocessing`, `bed_conversion`
- **Analysis Queue**: `mgmt`, `cnv`, `target`, `fusion`
- **Classification Queue**: `sturgeon`, `nanodx`, `pannanodx`
- **Slow Queue**: `random_forest`

### `watch`
Watch a directory for BAM file changes:

```bash
littlejohn watch /path/to/directory [OPTIONS]
```

**Options:**
- `--command, -c`: Command to run when files change
- `--verbose, -v`: Enable verbose output
- `--no-process-existing`: Skip processing existing files, only watch for new changes
- `--no-progress`: Disable progress bars for file processing

### `workflow`
Run an async workflow on BAM files in a directory:

```bash
littlejohn workflow /path/to/directory --workflow "workflow_plan" [OPTIONS]
```

**Required Options:**
- `--workflow, -w`: Workflow plan (e.g., 'preprocessing:bed_conversion,analysis:mgmt,classification:sturgeon')

**Optional Options:**
- `--commands, -c`: Command mappings (e.g., 'index:samtools index {file}')
- `--verbose, -v`: Enable verbose output
- `--no-process-existing`: Skip processing existing files, only watch for new changes
- `--work-dir, -d`: Base output directory for analysis results
- `--log-level`: Global log level (DEBUG|INFO|WARNING|ERROR, default: INFO)
- `--job-log-level`: Set log level for specific job (e.g., 'preprocessing:DEBUG', 'mgmt:WARNING')
- `--deduplicate-jobs`: Job types to deduplicate by sample ID (e.g., 'sturgeon', 'mgmt')
- `--no-progress`: Disable progress bars for file processing

### `list-files`
List BAM files in a directory recursively:

```bash
littlejohn list-files /path/to/directory
```

### `info`
Display information about a file or directory:

```bash
littlejohn info /path/to/file_or_directory
```

## Workflow System

LittleJohn uses a sophisticated multi-queue workflow system with four specialized queues:

### Queue Types

1. **Preprocessing Queue**: Fast jobs for metadata extraction and file conversion
   - `preprocessing`: Extract metadata from BAM files
   - `bed_conversion`: Convert BAM files to BED format

2. **Analysis Queue**: Bioinformatics analysis jobs
   - `mgmt`: MGMT methylation analysis
   - `cnv`: Copy number variation analysis
   - `target`: Target analysis
   - `fusion`: Fusion detection analysis

3. **Classification Queue**: Machine learning classification jobs
   - `sturgeon`: Sturgeon classification analysis
   - `nanodx`: NanoDX analysis
   - `pannanodx`: PanNanoDX analysis

4. **Slow Queue**: Resource-intensive jobs
   - `random_forest`: Random Forest analysis

### Workflow Syntax

Workflows are specified using the format `queue:job_type`:

```bash
# Simple workflow: preprocessing + MGMT analysis
littlejohn workflow /path/to/bam/files --workflow "preprocessing:bed_conversion,analysis:mgmt"

# Full pipeline: preprocessing + multiple analyses + classification
littlejohn workflow /path/to/bam/files --workflow "preprocessing:bed_conversion,analysis:mgmt,analysis:cnv,classification:sturgeon"

# Multiple analyses for the same sample
littlejohn workflow /path/to/bam/files --workflow "preprocessing:bed_conversion,analysis:mgmt,analysis:cnv,analysis:fusion"
```

**Note**: The `preprocessing:preprocessing` step is automatically added as the first step if not specified, ensuring metadata extraction for all workflows.

## Analysis Modules

### MGMT Analysis
MGMT (O6-methylguanine-DNA methyltransferase) promoter methylation analysis:

```bash
littlejohn workflow /path/to/bam/files --workflow "analysis:mgmt" --verbose
```

**Features:**
- Extracts reads from MGMT promoter region (chr10:129466536-129467536)
- Runs methylation calling using matkit from the robin package
- Executes R scripts for methylation status prediction
- Generates methylation plots and CSV results
- Creates sample-specific output directories

### CNV Analysis
Copy number variation analysis with breakpoint detection:

```bash
littlejohn workflow /path/to/bam/files --workflow "analysis:cnv" --verbose
```

**Features:**
- Two-pass analysis: sample vs reference CNV comparison
- Dynamic bin width calculation
- Breakpoint detection using Kernel Change Point Detection
- Genetic sex estimation from CNV patterns
- State persistence for incremental processing
- Comprehensive output in multiple formats

### Fusion Analysis
Fusion detection analysis:

```bash
littlejohn workflow /path/to/bam/files --workflow "analysis:fusion" --verbose
```

**Features:**
- Structural variant detection
- Fusion candidate identification
- Comprehensive fusion analysis pipeline

### Target Analysis
Target-specific analysis:

```bash
littlejohn workflow /path/to/bam/files --workflow "analysis:target" --verbose
```

### Classification Modules

#### Sturgeon Classification
```bash
littlejohn workflow /path/to/bam/files --workflow "classification:sturgeon" --verbose
```

#### NanoDX Analysis
```bash
littlejohn workflow /path/to/bam/files --workflow "classification:nanodx" --verbose
```

#### PanNanoDX Analysis
```bash
littlejohn workflow /path/to/bam/files --workflow "classification:pannanodx" --verbose
```

#### Random Forest Analysis
```bash
littlejohn workflow /path/to/bam/files --workflow "slow:random_forest" --verbose
```

## Advanced Features

### Job Deduplication
Prevent redundant processing when multiple upstream jobs complete simultaneously:

```bash
littlejohn workflow /path/to/bam/files \
  --workflow "preprocessing:bed_conversion,classification:sturgeon" \
  --deduplicate-jobs sturgeon \
  --verbose
```

### Configurable Logging
Fine-grained control over output verbosity:

```bash
littlejohn workflow /path/to/bam/files \
  --workflow "preprocessing:bed_conversion,analysis:mgmt" \
  --log-level WARNING \
  --job-log-level preprocessing:DEBUG \
  --job-log-level mgmt:INFO \
  --verbose
```

### Custom Output Directory
Specify a custom output directory for all analysis results:

```bash
littlejohn workflow /path/to/bam/files \
  --workflow "analysis:mgmt,analysis:cnv" \
  --work-dir /path/to/output \
  --verbose
```

### Progress Tracking
Real-time progress monitoring with multiple progress bars:

```bash
littlejohn workflow /path/to/bam/files \
  --workflow "preprocessing:bed_conversion,analysis:mgmt,classification:sturgeon" \
  --verbose
```

Progress bars show:
- Queue status (Q: queue size, A: active jobs)
- Active job details with duration
- Overall progress across all jobs
- Real-time updates as jobs complete

### Master CSV Tracking
Automatic creation and updating of `master.csv` files for each sample:
- Sample directories with cumulative statistics
- Total read counts, base counts, and yield tracking
- Pass/Fail BAM file breakdowns
- Run information aggregation
- File tracking for all processed BAM files

## Usage Examples

### Basic MGMT Analysis
```bash
# Run MGMT analysis on all BAM files in a directory
littlejohn workflow /path/to/bam/files --workflow "analysis:mgmt" --verbose
```

### Full Bioinformatics Pipeline
```bash
# Run comprehensive analysis pipeline
littlejohn workflow /path/to/bam/files \
  --workflow "preprocessing:bed_conversion,analysis:mgmt,analysis:cnv,classification:sturgeon" \
  --work-dir /path/to/results \
  --deduplicate-jobs sturgeon \
  --log-level INFO \
  --verbose

# Complete workflow with all analysis modules
littlejohn workflow ~/datasets/Demo_Run \
  --work-dir test_out2 \
  --workflow "preprocessing:bed_conversion,analysis:fusion,analysis:mgmt,analysis:cnv,analysis:target,classification:sturgeon,classification:nanodx,classification:pannanodx,slow:random_forest" \
  --verbose \
  --log-level WARNING
```

### Watch for New Files
```bash
# Process existing files first, then watch for new ones
littlejohn workflow /path/to/bam/files \
  --workflow "analysis:mgmt" \
  --verbose

# Only watch for new files (skip existing)
littlejohn workflow /path/to/bam/files \
  --workflow "analysis:mgmt" \
  --no-process-existing \
  --verbose
```

### Custom Commands
```bash
# Add custom commands to the workflow
littlejohn workflow /path/to/bam/files \
  --workflow "preprocessing:bed_conversion,analysis:mgmt" \
  --commands "notify:echo 'Processing complete for {file}'" \
  --verbose
```

## Dependencies

### Core Dependencies
- `click>=8.0.0`: CLI framework
- `watchdog>=3.0.0`: File system monitoring
- `pysam>=0.21.0`: BAM file processing
- `pandas>=1.3.0`: Data manipulation
- `numpy>=1.21.0`: Numerical computations
- `scipy>=1.7.0`: Scientific computing
- `ruptures>=1.1.0`: Change point detection
- `tqdm>=4.64.0`: Progress bars

### External Dependencies
- `bedtools`: For region extraction and BED file operations
- `robin` package: For methylation calling and R scripts
- `R` and `Rscript`: For statistical analysis and classification
- `samtools`: For BAM file manipulation (optional)

## Development

### Running Tests
```bash
pytest
```

### Code Formatting
```bash
black src/ tests/
isort src/ tests/
```

### Type Checking
```bash
mypy src/
```

### Linting
```bash
flake8 src/ tests/
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Run the test suite
6. Submit a pull request

## License

The license type to develop is not yet applied.

## Acknowledgments

- [Click](https://click.palletsprojects.com/) - Python package for creating command line interfaces
- [Watchdog](https://python-watchdog.readthedocs.io/) - Python library for monitoring file system events
- [pysam](https://pysam.readthedocs.io/) - Python interface for SAM/BAM files
- [robin](https://github.com/yourusername/robin) - Bioinformatics analysis package 