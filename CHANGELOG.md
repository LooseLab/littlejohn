# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- SNP analysis pipeline via Clair3 with snpEff and SnpSift annotation.
- SNP display data generation and GUI table with IGV navigation. Current extensive logging to the command line to track progress.
- IGV-ready BAM creation for consistent genome browser loading.
- Workflow support for `igv_bam` and `snp_analysis` job types.
- SNP analysis queueing from the GUI with concurrency control.
- MNP-FLEX upload support for bedMethyl outputs for beta users with Epignostix credentials (https://epignostix.com/). Supply credentials via `MNPFLEX_USER` and `MNPFLEX_PASS` environment variables.
- `robin utils mgmt` command to summarize MGMT CpG site methylation counts from `mgmt_sorted.bam` outputs.

### Changed
- IGV viewer initialization and BAM loading prioritization.
- Target analysis now always emits `targets_exceeding_threshold.bed` for SNP analysis.
- MGMT methylation classification uses per-read probabilities (fixes aggregated max bleed).
- Version metadata aligned to `0.4.0` across packaging and app entry points.
- `robin utils mgmt` now reports methylation percent from counts and normalizes bedmethyl fraction inputs.
- MNP-Flex results panel now hides empty fields until results are available.
- `robin utils mgmt` supports recursive search and file output for TSV summaries.
- Command-line progress output supports Rich styling; updating the installation is required to enable Rich functionality.

### Dependencies
- Added snpEff and SnpSift for SNP annotation.
- Added Rich for enhanced CLI progress and styled output (requires updating installation).

## [0.0.2] - 2024-12-19

### Added
- Model setup scripts and asset management system
- Performance optimizations for reporting system
- Enhanced error handling and loading fixes

### Changed
- Migrated from Git LFS to release asset system for model files
- Updated README with correct repository URL and installation instructions
- Improved performance across multiple analysis modules

### Fixed
- Fusion detection analysis fixes
- Tag error handling improvements
- BED file parsing to allow '.' in column 5
- Minimum read count validation (ensures minimum of 3 reads)
- Various error handling improvements

### Performance
- Speed improvements in reporting system
- Multiple performance tweaks across analysis modules
- Optimized loading and processing workflows

## [0.0.1] - 2024-12-14

### Added
- Initial release of robin bioinformatics workflow engine
- BAM file processing and analysis capabilities
- Multi-queue workflow system
- Built-in analysis modules (MGMT, CNV, fusion detection, etc.)
- GUI workflow monitor using NiceGUI
- Ray distributed computing support
- Comprehensive documentation and installation guides
