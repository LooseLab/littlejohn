# LittleJohn

A Python CLI tool with file watching capabilities built using Click and Watchdog.

## Features

- **File Watching**: Monitor directories and files for changes using Watchdog
- **BAM File Processing**: Specialized support for BAM files with automatic metadata extraction
- **Configurable Logging**: Global and per-job log level configuration for optimal verbosity control
- **MGMT Analysis**: Automated methylation analysis for the MGMT promoter region
- **CNV Analysis**: Comprehensive copy number variation analysis with breakpoint detection
- **Master CSV Tracking**: Automatic creation and updating of master.csv files for each sample
- **Pattern Matching**: Watch specific file patterns and ignore others
- **Command Execution**: Run commands when files change
- **Recursive Watching**: Watch directories recursively
- **CLI Interface**: Easy-to-use command-line interface built with Click
- **Workflow Management**: Run multi-step workflows with automatic preprocessing
- **Progress Bars**: Visual progress indication for file processing operations
- **Worker Progress Tracking**: Real-time progress monitoring for multi-worker workflows

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

## Usage

### Basic File Watching

Watch a directory for changes:

```bash
littlejohn watch /path/to/directory
```

Watch a specific file:

```bash
littlejohn watch /path/to/file.txt
```

### Advanced File Watching

Watch with specific patterns:

```bash
littlejohn watch /path/to/directory --patterns "*.py" "*.txt"
```

Watch recursively:

```bash
littlejohn watch /path/to/directory --recursive
```

Ignore certain patterns:

```bash
littlejohn watch /path/to/directory --ignore-patterns "*.tmp" "*.log"
```

Run a command when files change:

```bash
littlejohn watch /path/to/directory --command "echo {file} changed"
```

Process existing files first, then watch for changes (default behavior):

```bash
littlejohn watch /path/to/directory --patterns "*.bam" --command "samtools index {file}"
```

Skip existing files, only watch for new changes:

```bash
littlejohn watch /path/to/directory --patterns "*.bam" --command "samtools index {file}" --no-process-existing
```

Enable verbose output:

```bash
littlejohn watch /path/to/directory --verbose
```

Disable progress bars:

```bash
littlejohn watch /path/to/directory --verbose --no-progress
```

### File Listing

List files in a directory:

```bash
littlejohn list-files /path/to/directory
```

List files recursively:

```bash
littlejohn list-files /path/to/directory --recursive
```

List files matching a pattern (supports wildcards):

```bash
# Find all .bam files
littlejohn list-files /path/to/directory --pattern "*.bam"

# Find files recursively
littlejohn list-files /path/to/directory --recursive --pattern "*.bam"

# Find multiple file types
littlejohn list-files /path/to/directory --patterns "*.bam" --patterns "*.fastq"

# Exclude certain patterns
littlejohn list-files /path/to/directory --pattern "*.bam" --ignore-patterns "*.tmp"
```

### File Information

Get information about a file or directory:

```bash
littlejohn info /path/to/file.txt
```

## Command Reference

### `watch`

Watch a directory or file for changes.

**Arguments:**
- `path`: Path to the directory or file to watch

**Options:**
- `--recursive, -r`: Watch directories recursively
- `--patterns, -p`: File patterns to watch (default: all files)
- `--ignore-patterns, -i`: File patterns to ignore
- `--command, -c`: Command to run when files change
- `--verbose, -v`: Enable verbose output
- `--no-process-existing`: Skip processing existing files, only watch for new changes
- `--no-progress`: Disable progress bars for file processing

### `list-files`

List files in a directory matching a pattern.

**Arguments:**
- `path`: Path to the directory or file

**Options:**
- `--recursive, -r`: Search directories recursively
- `--pattern, -p`: File pattern to search for (default: all files, supports wildcards like `*.bam`)
- `--patterns`: Multiple file patterns to search for (use multiple times)
- `--ignore-patterns, -i`: Patterns to ignore (use multiple times)

### `workflow`

Run an async workflow on files in a directory.

**Arguments:**
- `path`: Path to the directory to watch

**Options:**
- `--workflow, -w`: Workflow plan (e.g., 'fast:index,slow:qc,fast:notify'). Preprocessing is automatically added as the first step.
- `--commands, -c`: Command mappings (e.g., 'index:samtools index {file}')
- `--verbose, -v`: Enable verbose output
- `--no-process-existing`: Skip processing existing files, only watch for new changes
- `--no-progress`: Disable progress bars for file processing
- `--work-dir, -d`: Base output directory for analysis results (defaults to BAM file directory)

### `info`

Display information about a file or directory.

**Arguments:**
- `path`: Path to the file or directory

## Examples

### Development Workflow

Watch Python files and run tests when they change:

```bash
littlejohn watch src/ --patterns "*.py" --command "python -m pytest tests/"
```

### Log File Monitoring

Watch log files and notify when they change:

```bash
littlejohn watch /var/log/ --patterns "*.log" --command "echo 'Log file changed: {file}'"
```

### Build Process

Watch source files and rebuild when they change:

```bash
littlejohn watch src/ --patterns "*.py" --ignore-patterns "*.pyc" --command "python setup.py build"
```

### BAM File Processing

LittleJohn includes specialized support for BAM files with automatic metadata extraction:

```bash
# Watch for BAM files and automatically extract metadata
littlejohn watch ~/datasets/ --command "echo 'BAM file detected: {file}'"

# List all BAM files in a directory (recursive search)
littlejohn list-files ~/datasets/

# Run workflow with automatic preprocessing (processes existing files first)
littlejohn workflow ~/datasets/ \
  --workflow "fast:index,slow:qc,fast:notify" \
  --commands "index:samtools index {file}" \
  --commands "qc:fastqc {file}" \
  --commands "notify:echo 'Processing complete for {file}'"

# Run workflow only on new files (skip existing files)
littlejohn workflow ~/datasets/ \
  --workflow "fast:index,slow:qc,fast:notify" \
  --commands "index:samtools index {file}" \
  --commands "qc:fastqc {file}" \
  --commands "notify:echo 'Processing complete for {file}'" \
  --no-process-existing
```

**Automatic Preprocessing**: The workflow command automatically includes a preprocessing step that extracts:
- Sample ID from BAM headers
- Read counts (mapped/unmapped)
- Total yield
- Run information (platform, flow cell, device position)
- Basecall model
- Pass/fail state
- File statistics

This metadata is available to all downstream workflow steps through the job context.

### Progress Bars

LittleJohn includes progress bars to show the status of file processing operations:

```bash
# Enable progress bars (default when verbose is enabled)
littlejohn watch ~/datasets/ --verbose --patterns "*.bam"

# Disable progress bars
littlejohn watch ~/datasets/ --verbose --no-progress --patterns "*.bam"

# Progress bars also work with workflows
littlejohn workflow ~/datasets/ \
  --workflow "preprocessing:bed_conversion,analysis:mgmt" \
  --verbose
```

Progress bars show:
- Total number of files to process
- Current file being processed
- Processing speed and estimated time remaining
- Overall completion percentage

### Worker Progress Tracking

For workflows, LittleJohn provides detailed worker progress tracking:

```bash
# Enable worker progress tracking (default when verbose is enabled)
littlejohn workflow ~/datasets/ \
  --workflow "preprocessing:bed_conversion,analysis:mgmt,classification:sturgeon" \
  --verbose
```

Worker progress tracking shows:
- **Multiple progress bars** - One for each worker type (Preprocessing, Analysis, Classification, Slow)
- **Queue information** - Shows queue size (Q) and active jobs (A) for each worker
- **Active job details** - Shows which specific jobs are running with their duration
- **Overall progress** - Shows total progress across all jobs
- **Real-time updates** - Progress bars update in real-time as jobs complete

Example output:
```
Preprocessing (Q:4 A:1): preprocessing:test_001.bam(3s)
Analysis (Q:0 A:0): 
Classification (Q:0 A:0): 
Slow (Q:0 A:0): 
Overall Progress: 25%|██▌ | 4/16 [00:12<00:36, Active: 1 | Completed: 3 | Failed: 0]
```

### Configurable Logging

LittleJohn provides a sophisticated logging system that allows fine-grained control over output verbosity:

```bash
# Set global log level to WARNING (suppress most output)
littlejohn workflow ~/datasets/ \
  --workflow "fast:preprocessing,fast:mgmt" \
  --log-level WARNING

# Set preprocessing to DEBUG level (show detailed preprocessing info)
littlejohn workflow ~/datasets/ \
  --workflow "fast:preprocessing,fast:mgmt" \
  --log-level WARNING \
  --job-log-level preprocessing:DEBUG

# Different log levels for different jobs
littlejohn workflow ~/datasets/ \
  --workflow "fast:preprocessing,fast:mgmt,fast:cnv" \
  --log-level INFO \
  --job-log-level preprocessing:DEBUG \
  --job-log-level mgmt:WARNING \
  --job-log-level cnv:ERROR
```

**Logging Features**:
- **Global Log Level**: Set the default log level for all components (DEBUG, INFO, WARNING, ERROR)
- **Per-Job Log Levels**: Configure different verbosity levels for specific job types
- **Structured Logging**: Each log message includes job context (job type, job ID, filename)
- **Context-Aware**: Log messages automatically include relevant file and job information
- **Flexible Configuration**: Easy configuration via CLI options or programmatic API

**Log Levels**:
- **DEBUG**: Detailed diagnostic information (file sizes, processing steps, etc.)
- **INFO**: General information about job progress and completion
- **WARNING**: Warning messages for non-critical issues
- **ERROR**: Error messages for failed operations

**Master CSV Tracking**: During preprocessing, LittleJohn automatically creates and updates `master.csv` files for each sample:
- **Sample directories**: Each sample gets its own directory in the output folder
- **Cumulative statistics**: Tracks total read counts, base counts, and yield across all BAM files
- **Pass/Fail breakdowns**: Separate counters for pass and fail BAM files
- **Run information**: Aggregates device positions, basecall models, flow cell IDs, and run times
- **File tracking**: Maintains a list of all BAM files processed for each sample

The master.csv files provide a comprehensive summary of all BAM files processed for each sample and are automatically updated during preprocessing.

### MGMT Analysis

LittleJohn includes specialized MGMT (O6-methylguanine-DNA methyltransferase) promoter methylation analysis:

```bash
# Run MGMT analysis workflow (preprocessing + MGMT analysis)
littlejohn workflow ~/datasets/ \
  --workflow "fast:mgmt" \
  --verbose

# Run MGMT analysis with custom output directory
littlejohn workflow ~/datasets/ \
  --workflow "fast:mgmt" \
  --work-dir ~/analysis_results \
  --verbose

# Run MGMT analysis with additional steps
littlejohn workflow ~/datasets/ \
  --workflow "fast:mgmt,fast:notify" \
  --commands "notify:echo 'MGMT analysis complete for {file}'" \
  --verbose
```

**MGMT Analysis Features**:
- **Region Extraction**: Uses bedtools to extract reads from the MGMT promoter region (chr10:129466536-129467536)
- **Methylation Analysis**: Runs matkit from the robin package for methylation calling
- **R-based Prediction**: Executes R scripts for MGMT methylation status prediction
- **Visualization**: Generates methylation plots using methylartist
- **Sample Organization**: Creates sample-specific output directories
- **Incremental Processing**: Accumulates reads across multiple BAM files per sample
- **Comprehensive Logging**: Detailed progress tracking and error reporting

**Output Structure**:
```
work_directory/
├── sample_id_1/
│   ├── mgmt.bam              # Accumulated MGMT region reads
│   ├── 1_mgmt.bed           # Methylation data (file 1)
│   ├── 1_mgmt.csv           # MGMT prediction results (file 1)
│   ├── 1_mgmt.png           # Methylation visualization (file 1)
│   ├── 2_mgmt.bed           # Methylation data (file 2)
│   └── ...
├── sample_id_2/
│   └── ...
└── mgmt_hg38.bed            # MGMT region definition
```

**Dependencies**: The MGMT analysis requires:
- `bedtools` for region extraction
- `robin` package for matkit methylation calling
- `R` and `Rscript` for prediction models
- `methylartist` for visualization (optional)
- `pysam` for BAM file manipulation

**Existing File Processing**: By default, workflows process all existing BAM files in the directory before starting to watch for new files. Use `--no-process-existing` to skip existing files and only process new files as they appear.

### Job Deduplication

LittleJohn includes intelligent job deduplication to prevent redundant processing when multiple upstream jobs complete simultaneously:

```bash
# Enable deduplication for sturgeon jobs (prevents duplicate sturgeon analysis per sample)
littlejohn workflow ~/datasets/ \
  --workflow "fast:bed_conversion,fast:sturgeon" \
  --deduplicate-jobs sturgeon \
  --verbose

# Enable deduplication for multiple job types
littlejohn workflow ~/datasets/ \
  --workflow "fast:bed_conversion,fast:sturgeon,fast:mgmt" \
  --deduplicate-jobs sturgeon mgmt \
  --verbose
```

**Job Deduplication Features**:
- **Sample-Based Deduplication**: Jobs are deduplicated by sample ID, not by individual files
- **Configurable Job Types**: Specify which job types should be deduplicated
- **Thread-Safe**: Uses proper locking to handle concurrent job creation
- **Automatic Cleanup**: Pending job tracking is automatically cleaned up when jobs complete
- **Comprehensive Logging**: Detailed logging of deduplication decisions

**Use Cases**:
- **Sturgeon Analysis**: When multiple BAM files for the same sample complete bed_conversion simultaneously, only one sturgeon job is queued per sample
- **Downstream Analysis**: Any analysis that operates on the output of upstream jobs can benefit from deduplication
- **Resource Optimization**: Prevents unnecessary duplicate processing and resource consumption

**How It Works**:
1. When a job is enqueued, the system checks if a job of the same type is already pending for the sample
2. If a pending job exists, the new job is skipped with a log message
3. When a job completes (success or failure), it's automatically unmarked as pending
4. This allows subsequent jobs for the same sample to be processed normally

**Default Behavior**: The `sturgeon`, `nanodx`, `pannanodx`, and `random_forest` job types are automatically included in deduplication by default, as they are common use cases for downstream analysis. Note that `random_forest` runs in the slow queue but still benefits from deduplication to prevent duplicate processing per sample.

### CNV Analysis

LittleJohn includes specialized CNV (Copy Number Variation) analysis with comprehensive copy number variation detection:

```bash
# Run CNV analysis workflow (preprocessing + CNV analysis)
littlejohn workflow ~/datasets/ \
  --workflow "fast:cnv" \
  --verbose

# Run CNV analysis with custom output directory
littlejohn workflow ~/datasets/ \
  --workflow "fast:cnv" \
  --work-dir ~/cnv_results \
  --verbose

# Run CNV analysis with additional steps
littlejohn workflow ~/datasets/ \
  --workflow "fast:cnv,fast:notify" \
  --commands "notify:echo 'CNV analysis complete for {file}'" \
  --verbose

# Run CNV analysis with MGMT analysis
littlejohn workflow ~/datasets/ \
  --workflow "fast:cnv,fast:mgmt" \
  --work-dir ~/analysis_results \
  --verbose
```

**CNV Analysis Features**:
- **Two-pass Analysis**: Sample vs reference CNV comparison
- **Dynamic Bin Width**: Automatic calculation of optimal bin width
- **Breakpoint Detection**: Kernel Change Point Detection for structural variations
- **Sex Estimation**: Genetic sex determination from CNV patterns
- **State Persistence**: Incremental processing with state tracking
- **Comprehensive Output**: Multiple file formats for downstream analysis
- **Sample Organization**: Sample-specific output directories
- **Incremental Processing**: Accumulates data across multiple BAM files per sample

**Output Structure**:
```
work_directory/
├── sample_id_1/
│   ├── CNV.npy                    # Sample CNV data
│   ├── CNV2.npy                   # Reference CNV data
│   ├── CNV3.npy                   # Normalized CNV data
│   ├── CNV_dict.npy               # Processing metadata
│   ├── cnv_data_array.npy         # Breakpoint data
│   ├── XYestimate.pkl             # Sex estimation
│   ├── update_cnv_dict.pkl        # Accumulated copy numbers
│   ├── cnv_analysis_counter.txt   # Analysis counter
│   ├── cnv_detector_state/        # State directory
│   │   └── tracker_metadata.pkl   # State metadata
│   └── bed_files/                 # BED files directory
│       ├── new_file_001.bed       # CNV regions
│       └── breakpoints_001.bed    # Breakpoints
├── sample_id_2/
│   └── ...
└── master.csv                     # Master CSV tracking
```

**Dependencies**: The CNV analysis requires:
- `cnv_from_bam` for CNV extraction from BAM files
- `robin` package for reference CNV data and resources
- `ruptures` for change point detection
- `scipy` for scientific computing
- `numpy` for numerical computations
- `pandas` for data manipulation
- `pysam` for BAM file processing

### Random Forest Analysis

LittleJohn includes specialized Random Forest methylation classification analysis:

```bash
# Run Random Forest analysis workflow (preprocessing + bed_conversion + Random Forest)
littlejohn workflow ~/datasets/ \
  --workflow "preprocessing:bed_conversion,slow:random_forest" \
  --verbose

# Run Random Forest analysis with custom output directory
littlejohn workflow ~/datasets/ \
  --workflow "preprocessing:bed_conversion,slow:random_forest" \
  --work-dir ~/random_forest_results \
  --verbose

# Run Random Forest analysis with additional steps
littlejohn workflow ~/datasets/ \
  --workflow "preprocessing:bed_conversion,analysis:mgmt,slow:random_forest,classification:sturgeon" \
  --work-dir ~/analysis_results \
  --verbose
```

**Random Forest Analysis Features**:
- **Methylation Data Extraction**: Uses modkit from the robin package for methylation calling
- **Random Forest Classification**: R-based classification using methylation patterns
- **Time Series Analysis**: Classification confidence tracking over time
- **Batch Processing**: Incremental processing with batch numbering
- **Comprehensive Output**: Multiple file formats for downstream analysis
- **Sample Organization**: Sample-specific output directories
- **R Script Integration**: Executes R scripts for classification
- **Comprehensive Logging**: Detailed progress tracking and error reporting

**Output Structure**:
```
work_directory/
├── sample_id_1/
│   ├── RandomForestBed.bed        # BED file for R script input
│   ├── random_forest_scores.csv   # Classification scores
│   └── temp_*/                    # Temporary R script output
│       └── live_1_votes.tsv       # Raw R script results
├── sample_id_2/
│   └── ...
└── master.csv                     # Master CSV tracking
```

**Dependencies**: The Random Forest analysis requires:
- `robin` package for modkit methylation calling and R scripts
- `R` and `Rscript` for Random Forest classification
- `pandas` for data manipulation
- `numpy` for numerical computations
- `subprocess` for R script execution

**R Script Requirements**: The analysis requires specific R scripts and model files from the robin package:
- `methylation_classification_nanodx_v0.2.R` - Main classification script
- `top_probes_hm450.Rdata` - Probe data
- `capper_top_100k_betas_binarised.Rdata` - Training data
- `HM450.hg38.manifest.gencode.v22.Rdata` - Array manifest

### Advanced Workflow Processing

Run complex multi-step workflows on files:

```bash
# Simple workflow: index BAM files (preprocessing automatically included)
littlejohn workflow ~/datasets/ --workflow "fast:index" --commands "index:samtools index {file}"

# Multi-step workflow: index, quality check, and notify (preprocessing automatically included)
littlejohn workflow ~/datasets/ \
  --workflow "fast:index,slow:qc,fast:notify" \
  --commands "index:samtools index {file}" \
  --commands "qc:samtools flagstat {file}" \
  --commands "notify:echo 'Processing complete for {file}'"

# Complex bioinformatics pipeline (preprocessing automatically included)
littlejohn workflow ~/sequencing_runs/ \
  --workflow "fast:validate,slow:index,slow:stats,fast:notify" \
  --commands "validate:samtools quickcheck {file}" \
  --commands "index:samtools index {file}" \
  --commands "stats:samtools flagstat {file} > {file}.stats" \
  --commands "notify:echo 'Pipeline complete: {file}'"
```

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

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- [Click](https://click.palletsprojects.com/) - Python package for creating command line interfaces
- [Watchdog](https://python-watchdog.readthedocs.io/) - Python library for monitoring file system events 