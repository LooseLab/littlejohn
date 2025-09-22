# robin

A Python CLI tool for automated BAM file processing and bioinformatics analysis workflows.

## Overview

robin is a specialized bioinformatics workflow engine designed for processing BAM files through automated pipelines. It provides a robust, multi-threaded system for running complex analysis workflows with built-in support for methylation analysis, copy number variation detection, fusion detection, and various classification algorithms.

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
- **Ray Distributed Computing**: Optional distributed processing support
- **GUI Workflow Monitor**: Built-in NiceGUI-based monitoring interface

## Installation

**We strongly recommend installing robin in a conda environment to ensure all dependencies are properly managed.**

### Prerequisites

1. **Conda**: Install Miniconda or Anaconda if you haven't already.
2. **jq** (optional): For advanced asset management. Install with `brew install jq` (macOS) or `apt-get install jq` (Ubuntu).

### Installation Steps

1. **Clone the repository with submodules**:
   ```bash
   git clone --recursive https://github.com/LooseLab/littlejohn.git
   cd littlejohn
   
   # Switch to the branch with Git LFS removal
   git checkout remove_LFS
   ```

2. **Download required model assets**:
   ```bash
   # Download all models automatically
   python setup_models.py
   
   # Or download individual models
   python scripts/fetch_asset.py general_model src/robin/models/general.zip
   python scripts/fetch_asset.py capper_model src/robin/models/Capper_et_al_NN.pkl
   python scripts/fetch_asset.py pancan_model src/robin/models/pancan_devel_v5i_NN.pkl
   ```

3. **Update and initialize submodules**:
   ```bash
   git submodule update --init --recursive
   ```
   *This ensures all submodules (nanoDX, hv_rapidCNS2) are properly initialized*

4. **Create and activate conda environment**:
   ```bash
   # For Linux/Windows
   conda env create -f robin.yml
   conda activate robin_0_5
   
   # For macOS
   conda env create -f robin_osx.yml
   conda activate robin_0_5
   ```

5. **Install robin in development mode**:
   ```bash
   pip install -e .
   ```

### Alternative Installation (without conda)

If you prefer not to use conda, you can install from source, but you'll need to manually install the system dependencies:

1. **Install system dependencies**:
   ```bash
   # Ubuntu/Debian
   sudo apt-get install samtools bedtools r-base jq
   
   # macOS
   brew install samtools bedtools r jq
   ```

2. **Clone with submodules**:
   ```bash
   git clone --recursive https://github.com/LooseLab/littlejohn.git
   cd littlejohn
   
   # Switch to the branch with Git LFS removal
   git checkout remove_LFS
   ```

3. **Download required model assets**:
   ```bash
   python setup_models.py
   ```

4. **Update and initialize submodules**:
   ```bash
   git submodule update --init --recursive
   ```
   *This ensures all submodules (nanoDX, hv_rapidCNS2) are properly initialized*

5. **Install Python dependencies**:
   ```bash
   pip install -e .
   ```

## Asset Management

Robin uses a release asset system instead of Git LFS for managing large model files. This provides better performance and reliability.

### Available Assets

The following model assets are available:

- **general_model**: `general.zip` (1.7GB) - General machine learning model archive
- **capper_model**: `Capper_et_al_NN.pkl` (132MB) - Capper et al neural network model  
- **pancan_model**: `pancan_devel_v5i_NN.pkl` (194MB) - Pan-cancer development v5i neural network model

### Downloading Assets

**Automatic setup** (recommended):
```bash
python setup_models.py
```

**Manual download**:
```bash
# Download individual models
python scripts/fetch_asset.py general_model src/robin/models/general.zip
python scripts/fetch_asset.py capper_model src/robin/models/Capper_et_al_NN.pkl
python scripts/fetch_asset.py pancan_model src/robin/models/pancan_devel_v5i_NN.pkl
```

**Using shell script**:
```bash
./scripts/fetch_asset.sh general_model src/robin/models/general.zip
./scripts/fetch_asset.sh capper_model src/robin/models/Capper_et_al_NN.pkl
./scripts/fetch_asset.sh pancan_model src/robin/models/pancan_devel_v5i_NN.pkl
```

### Authentication

For private repositories, set your GitHub token:
```bash
export GITHUB_TOKEN=your_personal_access_token
python setup_models.py
```

### Asset Verification

All assets are automatically verified using SHA256 checksums. If verification fails, the download will be retried or the corrupted file will be removed.

## Main Usage

The primary command for running robin workflows is:

```bash
robin workflow <data_folder> --work-dir <output_folder> -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest --reference ~/references/hg38_simple.fa --center <center_id>
```

### Command Breakdown

- `robin workflow`: The main workflow command
- `<data_folder>`: Directory containing your BAM files
- `--work-dir <output_folder>`: Directory where results will be saved
- `-w`: Workflow specification (comma-separated list of analysis types)
- `--reference`: Path to reference genome (required for some analyses)
- `--center <center_id>`: Center ID running the analysis (e.g., 'Sherwood', 'Auckland', 'New York')

### Example Usage

```bash
# Basic workflow with all analysis types
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/hg38_simple.fa \
  --center Sherwood

# Simplified workflow with just a few analyses
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center Auckland

# With verbose output and custom logging
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,cnv,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center New_York \
  --verbose \
  --log-level INFO
```

## Available Commands

### `list-job-types`
List all available job types organized by queue category:

```bash
robin list-job-types
```

**Available Job Types:**
- **Preprocessing Queue**: `preprocessing`
- **Bed Conversion Queue**: `bed_conversion`
- **Analysis Queue**: `mgmt`, `cnv`, `target`, `fusion`
- **Classification Queue**: `sturgeon`, `nanodx`, `pannanodx`
- **Slow Queue**: `random_forest`

### `workflow`
Run an async workflow on BAM files in a directory:

```bash
robin workflow /path/to/directory --workflow "workflow_plan" [OPTIONS]
```

**Required Options:**
- `--workflow, -w`: Workflow plan (e.g., 'mgmt,sturgeon' or 'preprocessing:bed_conversion,analysis:mgmt,classification:sturgeon')
- `--center`: Center ID running the analysis (e.g., 'Sherwood', 'Auckland', 'New York')

**Optional Options:**
- `--work-dir, -d`: Base output directory for analysis results
- `--reference, -r`: Path to reference genome (FASTA format)
- `--verbose, -v`: Enable verbose output and detailed error traces
- `--no-process-existing`: Skip processing existing files, only watch for new changes
- `--log-level`: Global log level (DEBUG|INFO|WARNING|ERROR, default: ERROR)
- `--job-log-level`: Set log level for specific job (e.g., 'preprocessing:DEBUG', 'mgmt:WARNING')
- `--deduplicate-jobs`: Job types to deduplicate by sample ID (e.g., 'sturgeon', 'mgmt')
- `--no-progress`: Disable progress bars for file processing
- `--use-ray/--no-use-ray`: Enable Ray distributed computing (default: on)
- `--with-gui/--no-gui`: Launch NiceGUI workflow monitor (default: on)



<!-- 
### `list-files`
List BAM files in a directory recursively:

```bash
robin list-files /path/to/directory
```

### `info`
Display information about a file or directory:

```bash
robin info /path/to/file_or_directory
```

## Workflow System

robin uses a sophisticated multi-queue workflow system with four specialized queues:

### Queue Types

1. **Preprocessing Queue**: Fast jobs for metadata extraction
   - `preprocessing`: Extract metadata from BAM files

2. **Bed Conversion Queue**: File format conversion jobs
   - `bed_conversion`: Convert BAM files to BED format

3. **Analysis Queue**: Bioinformatics analysis jobs
   - `mgmt`: MGMT methylation analysis
   - `cnv`: Copy number variation analysis
   - `target`: Target analysis
   - `fusion`: Fusion detection analysis

4. **Classification Queue**: Machine learning classification jobs
   - `sturgeon`: Sturgeon classification analysis
   - `nanodx`: NanoDX analysis
   - `pannanodx`: PanNanoDX analysis

5. **Slow Queue**: Resource-intensive jobs
   - `random_forest`: Random Forest analysis

### Workflow Syntax

Workflows can be specified in two formats:

#### Simplified Format (Recommended)
```bash
# Simple workflow: MGMT analysis + Sturgeon classification
robin workflow /path/to/bam/files -w "mgmt,sturgeon" --center Birmingham

# Full pipeline: all analysis types
robin workflow /path/to/bam/files -w "target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest" --center Coventry
```

#### Legacy Format (with queue prefixes)
```bash
# Simple workflow: preprocessing + MGMT analysis
robin workflow /path/to/bam/files -w "preprocessing:bed_conversion,analysis:mgmt" --center Derby

# Full pipeline: preprocessing + multiple analyses + classification
robin workflow /path/to/bam/files -w "preprocessing:bed_conversion,analysis:mgmt,analysis:cnv,classification:sturgeon" --center Plymouth
```

**Note**: The `preprocessing` step is automatically added as the first step if not specified, ensuring metadata extraction for all workflows. The `bed_conversion` step is automatically added when needed for classification jobs.

## Analysis Modules

### MGMT Analysis
MGMT (O6-methylguanine-DNA methyltransferase) promoter methylation analysis:

```bash
robin workflow /path/to/bam/files -w "mgmt" --center Sherwood --verbose
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
robin workflow /path/to/bam/files -w "cnv" --center Auckland --verbose
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
robin workflow /path/to/bam/files -w "fusion" --center New_York --verbose
```

**Features:**
- Structural variant detection
- Fusion candidate identification
- Comprehensive fusion analysis pipeline

### Target Analysis
Target-specific analysis:

```bash
robin workflow /path/to/bam/files -w "target" --center Oxford --verbose
```

### Classification Modules

#### Sturgeon Classification
```bash
robin workflow /path/to/bam/files -w "sturgeon" --center Cambridge --verbose
```

#### NanoDX Analysis
```bash
robin workflow /path/to/bam/files -w "nanodx" --center London --verbose
```

#### PanNanoDX Analysis
```bash
robin workflow /path/to/bam/files -w "pannanodx" --center Edinburgh --verbose
```

#### Random Forest Analysis
```bash
robin workflow /path/to/bam/files -w "random_forest" --center Glasgow --verbose
```

## Advanced Features

### Job Deduplication
Prevent redundant processing when multiple upstream jobs complete simultaneously:

```bash
robin workflow /path/to/bam/files \
  -w "mgmt,sturgeon" \
  --center Manchester \
  --deduplicate-jobs sturgeon \
  --verbose
```

### Configurable Logging
Fine-grained control over output verbosity:

```bash
robin workflow /path/to/bam/files \
  -w "mgmt,cnv" \
  --center Bristol \
  --log-level WARNING \
  --job-log-level preprocessing:DEBUG \
  --job-log-level mgmt:INFO \
  --verbose
```

### Ray Distributed Computing
Enable distributed processing for improved performance:

```bash
robin workflow /path/to/bam/files \
  -w "mgmt,cnv,sturgeon" \
  --center Cardiff \
  --use-ray \
  --ray-num-cpus 8 \
  --verbose
```

### GUI Workflow Monitor
Launch the built-in NiceGUI workflow monitor:

```bash
robin workflow /path/to/bam/files \
  -w "mgmt,sturgeon" \
  --center Belfast \
  --with-gui \
  --gui-port 8081 \
  --verbose
```

### Progress Tracking
Real-time progress monitoring with multiple progress bars:

```bash
robin workflow /path/to/bam/files \
  -w "mgmt,cnv,sturgeon" \
  --center Dublin \
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
robin workflow /path/to/bam/files -w "mgmt" --center Nottingham --verbose
```

### Full Bioinformatics Pipeline
```bash
# Run comprehensive analysis pipeline
robin workflow /path/to/bam/files \
  -w "mgmt,cnv,sturgeon" \
  --work-dir /path/to/results \
  --center Liverpool \
  --deduplicate-jobs sturgeon \
  --log-level INFO \
  --verbose

# Complete workflow with all analysis modules
robin workflow ~/datasets/Demo_Run \
  --work-dir test_out2 \
  -w "target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest" \
  --center Sheffield \
  --verbose \
  --log-level WARNING
```

### Watch for New Files
```bash
# Process existing files first, then watch for new ones
robin workflow /path/to/bam/files \
  -w "mgmt" \
  --center York \
  --verbose

# Only watch for new files (skip existing)
robin workflow /path/to/bam/files \
  -w "mgmt" \
  --center Leeds \
  --no-process-existing \
  --verbose
```
-->

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
- `ray[default]>=2.0.0`: Distributed computing
- `nicegui>=1.4.0`: GUI framework

### External Dependencies
- `bedtools`: For region extraction and BED file operations
- `samtools`: For BAM file manipulation
- `R` and `Rscript`: For statistical analysis and classification
- `git-lfs`: For large file storage (model files)

### Git Submodules
- `nanoDX`: NanoDX analysis tools
- `hv_rapidCNS2`: Rapid CNS analysis tools


<!--
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
-->


## License

MIT License

## Acknowledgments

- [Click](https://click.palletsprojects.com/) - Python package for creating command line interfaces
- [Watchdog](https://python-watchdog.readthedocs.io/) - Python library for monitoring file system events
- [pysam](https://pysam.readthedocs.io/) - Python interface for SAM/BAM files
- [Ray](https://ray.io/) - Distributed computing framework
- [NiceGUI](https://nicegui.io/) - Web-based GUI framework 