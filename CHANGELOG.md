# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
