# ![ROBIN_logo_small.png](src/robin/gui/images/ROBIN_logo_small.png) R.O.B.I.N

<!-- [![PyPI - Version](https://img.shields.io/pypi/v/methnicegui.svg)](https://pypi.org/project/methnicegui)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/methnicegui.svg)](https://pypi.org/project/methnicegui)
-->

# ***This software is provided as is for research use only.***

**Table of Contents**

- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
- [Known Issues](#known-issues)
- [Available Commands](#available-commands)
- [Performance Features](#performance-features)
- [License](#license)


-----

## About

ROBIN (Rapid nanopOre Brain intraoperatIve classificatioN) is a comprehensive bioinformatics workflow system designed for processing and analyzing BAM files in the context of human oncology analysis. 

ROBIN was published in NeuroOncology - see [ROBIN Paper](https://academic.oup.com/neuro-oncology/article/27/8/2035/8139084)

It provides automated preprocessing, multiple analysis pipelines, and real-time monitoring capabilities. 

It now incorporates LITTLE JOHN (Lightweight Infrastructure for Task Tracking and Logging with Extensible Job Orchestration for High-throughput aNalysis), which handles the heavy lifting behind the scenes.

## Overview

ROBIN processes BAM files through automated pipelines. It provides a robust, multi-threaded system for running complex analysis workflows with built-in support for methylation analysis, copy number variation detection, fusion detection, and various classification algorithms. The platform features an advanced web-based GUI for real-time workflow monitoring, progress tracking, and interactive data visualization.

## Memory usage

We have faced challenges with the amount of memory used by ROBIN on certain systems. The incorporation of LITTLE JOHN addresses this and we can now run ROBIN on two promethION flowcells simultaneously on a Nanopore P2i sequencer. These positions can be running ligation (LSK114) or our modified ultra long protocol.

## System requirements

RAM >= 64 Gb

GPU - as per ONT guidelines for running adaptive sampling

CPU - as per ONT guidelines (more is always better)

## Known Issues

1. This first release is for testing purposes only and feedback to us.
2. Variant calling in real time is currently unavailable. We will be reintroducing this in the near future.
3. All anaylses must be interpreted by an expert.

## The Future

This version will replace the code available at https://github.com/LooseLab/ROBIN/ in the near future.

---

## Installation

**We strongly recommend installing robin in a conda environment to ensure all dependencies are properly managed.**

### Prerequisites

1. **Conda**: Install Miniconda or Anaconda if you haven't already.

### Installation Steps

1. **Clone the repository with submodules**:
   ```bash
   git clone --recursive https://github.com/LooseLab/littlejohn.git
   cd littlejohn
   ```

2. **Update and initialize submodules**:
   ```bash
   git submodule update --init --recursive
   ```
   *This ensures all submodules (nanoDX, hv_rapidCNS2) are properly initialized*

3. **Create and activate conda environment**:
   ```bash
   # For Linux/Windows
   conda env create -f robin.yml
   conda activate robin_0_5
   
   # For macOS
   conda env create -f robin_osx.yml
   conda activate robin_0_5
   ```

4. **Install robin in development mode**:
   ```bash
   pip install -e .
   ```

5. **Download required model assets**:
   ```bash
   # Download all models automatically
   python setup_models.py
   
   # Or download individual models
   python scripts/fetch_asset.py general_model src/robin/models/general.zip
   python scripts/fetch_asset.py capper_model src/robin/models/Capper_et_al_NN.pkl
   python scripts/fetch_asset.py pancan_model src/robin/models/pancan_devel_v5i_NN.pkl
   ```



### Authentication

For private repositories, set your GitHub token:
```bash
export GITHUB_TOKEN=your_personal_access_token
python setup_models.py
```

### Asset Verification

All assets are automatically verified using SHA256 checksums. If verification fails, the download will be retried or the corrupted file will be removed.

## Usage

ROBIN expects to analyse BAM files generated during sequencing by an Oxford Nanopore Technologies sequencer. ROBIN presumes real time HAC basecalling (SUP is not required). ROBIN presumes data have been called with 5hmC 5mC methylation calling in MinKNOW. ROBIN presumes real time alignment is running in MinKNOW - ROBIN does not realign your reads. 

ROBIN assumes that BAM files are being output in small batches. We recommend setting file output to one bam for every 10,000 to 50,000 reads. We do not support real time processing of BAM files in 1 hour chunks (the default output). 

ROBIN does not pod5 data or fastq data from the sequencer - you can deselect these options if you wish.

***Important***

On platforms with 64Gb of RAM or less we recommend restarting your device prior to a run. As an example, if you are running a p2i and start a run on position A and later on position B we recommend you restart at the end of the run on position B (once base calling is complete). You can simply restart dorado if you know how to do this. We find dorado holds on to memory for an indefinite period of time and this can cause problems.


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
- `--target-panel`: The specific panel that is being applied. 

### Example Usage

```bash
# Basic workflow with all analysis types
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/hg38_simple.fa \
  --center Sherwood \
  --target-panel rCNS2

# Simplified workflow with just a few analyses
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center Auckland \
  --target-panel PanCan

# With verbose output and custom logging
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,cnv,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center New_York \
  --target-panel rCNS2 \
  --verbose \
  --log-level INFO
```

## Known Issues

1. Currently SNP calling is not enabled in this version of ROBIN. It will be re-enabled in the near future.
1. CNV change inference is based on extensive heuristics - every call should be checked by visual inspection.
1. If you ctrl-c to quit ROBIN it will do its best to clean up and stop gracefully but may fail.
1. CSV data export is in development but is not currently available - it will be enabled in the near future.
1. If you wish to reanalyse a data set you must remove the exsiting result in the robin output folder.
1. Many other unknown issues - please open an issue and we will resolve where possible.

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



### Panel management

Manage built-in and custom target panels used by analyses like `target`, `cnv`, and `fusion`.

- Built-in panels include: `rCNS2`, `AML`, `PanCan`.
- Custom panels are stored internally after you add them from a BED file.

#### List available panels
```bash
robin list-panels
```

#### Add a custom panel from a BED file
```bash
# Add and register a panel (BED must have ≥4 columns: chr, start, end, gene_name[s])
robin add-panel /path/to/your_panel.bed MyCustomPanel

# Optional: validate format only, without adding
robin add-panel /path/to/your_panel.bed MyCustomPanel --validate-only
```

Notes:
- Panel names cannot be empty and cannot reuse reserved names: `rCNS2`, `AML`, `PanCan`.
- BED may be 4- or 6-column; if multiple genes are in one region, use comma-separated names.

#### Remove a custom panel
```bash
# Will prompt for confirmation
robin remove-panel MyCustomPanel

# Skip confirmation
robin remove-panel MyCustomPanel --force
```

Built-in panels cannot be removed.

#### Use a panel in a workflow
```bash
robin workflow /path/to/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion \
  --target-panel MyCustomPanel \
  --center Sherwood
```

## Performance Features

### Enhanced Processing Capabilities
- **Batched Processing**: All analysis workflows now support batched processing for improved efficiency
- **Memory Optimization**: Intelligent memory management for large datasets and long-running processes
- **Multi-threading**: Configurable multi-threaded BAM processing via `LJ_BAM_THREADS` environment variable
- **Async Updates**: Non-blocking GUI updates during analysis execution
- **Smart Progress Tracking**: Streamlined progress indicators with real-time status updates

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
- `nicegui>=3.0.4`: Advanced GUI framework with enhanced performance

### External Dependencies
- `bedtools`: For region extraction and BED file operations
- `samtools`: For BAM file manipulation
- `R` and `Rscript`: For statistical analysis and classification

### Git Submodules
- `nanoDX`: NanoDX analysis tools
- `hv_rapidCNS2`: Rapid CNS analysis tools


## License

***This software is provided "as is", and is for research use only.***


robin is distributed under a CC BY-NC 4.0 license. See LICENSE for more information. This license does not override any licenses that may be present in the third party tools used by robin.

## Acknowledgments

This tool uses a range of third party tools and applications including:

[Sturgeon] https://github.com/marcpaga/sturgeon
[Radid-CNS2] https://link.springer.com/article/10.1007/s00401-022-02415-6
[Readfish] https://github.com/LooseLab/readfish
[cnv_from_bam] https://github.com/adoni5/cnv_from_bam
[methylartist] https://github.com/adamewing/methylartist

We are grateful to the authors of these tools for their work.

We also thank a lot of people who have contributed to these tools including: Graeme Fox, Simon Deacon, Rory Munro, Satrio Wibowo, Thomas Murray, Inswasti Cahyani, Nadine Holmes, Simon Paine, Stuart Smith and many others from outside Nottingham.

We are particularly grateful to Areeba Patel, Felix Sahm and colleagues for their work on Rapid-CNS2.

This list is non-exhaustive and the software is under active development.

In addition we use:

- [Click](https://click.palletsprojects.com/) - Python package for creating command line interfaces
- [Watchdog](https://python-watchdog.readthedocs.io/) - Python library for monitoring file system events
- [pysam](https://pysam.readthedocs.io/) - Python interface for SAM/BAM files
- [Ray](https://ray.io/) - Distributed computing framework
- [NiceGUI](https://nicegui.io/) - Web-based GUI framework 