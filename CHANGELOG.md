# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Dedicated "Manage watched folders" page at `/watched_folders` (replaces modal); accessible via link from the welcome page.
- Loading modal when adding a watched folder: shows "Adding folder..." with spinner and closes automatically on success or error.
- Path overlap validation in `add_watch_path`: rejects paths that equal, contain, or are inside the work directory to prevent watching already-analysed output.
- `parquet_filter.txt` as preferred CpG filter for parquet creation (falls back to `sturg_nanodx_cpgs_0125.bed.gz` when absent).
- `scripts/compare_mgmt_parquet_methylation.py` to compare MGMT methylation counts between bed and parquet for validation.
- `METHYLATION_EXTRACTION_COMPARISON.md` documenting differences between MNP-Flex (parquet) and MGMT (bedmethyl) methylation extraction, including coordinate handling and BAM set divergence.
- SNP analysis pipeline via Clair3 with snpEff and SnpSift annotation.
- SNP display data generation and GUI table with IGV navigation. Current extensive logging to the command line to track progress.
- IGV-ready BAM creation for consistent genome browser loading.
- Workflow support for `igv_bam` and `snp_analysis` job types.
- SNP analysis queueing from the GUI with concurrency control.
- MNP-FLEX upload support for bedMethyl outputs for beta users with Epignostix credentials (https://epignostix.com/). Supply credentials via `MNPFLEX_USER` and `MNPFLEX_PASS` environment variables.
- `robin utils mgmt` command to summarize MGMT CpG site methylation counts from `mgmt_sorted.bam` outputs.

### Changed
- **Report layout and typography:**
  - **Classification:** Methylation classification plots now use seaborn for improved styling; plots show top 3 classes with larger fonts and thicker lines. Detailed classification tables use a 2-column layout (Sturgeon|NanoDX, PanNanoDX|Random Forest) with compact styling and top 10 predictions.
  - **CNV:** Summary and event tables use compact styling with unified 9pt font. Chromosome Arm Events and CNV Events Containing Genes tables are displayed side-by-side in two columns. Arm/Region columns widened to prevent text wrapping (e.g. "p-arm"). Increased cell padding and line spacing for readability.
  - **Coverage:** Removed overlapping outlier labels from the "Coverage Distribution by Chromosome" box plot. Added a "Coverage outliers (potential gains/losses)" table below the plot with Chromosome, Gene, Coverage, and Type (Gain/Loss) columns.
  - **Base:** `create_table()` now supports `compact`, `font_size`, and improved leading for consistent table typography across sections.
- Parquet filter loading: `.txt` files (e.g. `parquet_filter.txt`) are treated as 1-based Illumina coordinates and converted to 0-based for BED/PyRanges; fixes off-by-one in CpG site filtering.
- Parquet filter cache: distinct cache suffix (`_1based`) for `.txt` filters to avoid stale data after coordinate conversion change.
- Parquet filter format: support for header row and whitespace-separated columns in `parquet_filter.txt`.
- **Note:** After the parquet filter coordinate fix, existing parquet files should be rebuilt (delete parquet and re-run bed_conversion) to obtain correct methylation counts.
- IGV viewer initialization and BAM loading prioritization.
- Target analysis now always emits `targets_exceeding_threshold.bed` for SNP analysis.
- MGMT methylation classification uses per-read probabilities (fixes aggregated max bleed).
- Version metadata aligned to `0.4.0` across packaging and app entry points.
- `robin utils mgmt` now reports methylation percent from counts and normalizes bedmethyl fraction inputs.
- MNP-Flex results panel now hides empty fields until results are available.
- `robin utils mgmt` supports recursive search and file output for TSV summaries.
- Command-line progress output supports Rich styling; updating the installation is required to enable Rich functionality.

### Dependencies
- Added seaborn for methylation classification plot styling.
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
