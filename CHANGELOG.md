# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-08-05

### Added
- Initial release of FLT3 ITD Detection Pipeline
- CIGAR-based ITD detection from Oxford Nanopore sequencing data
- Multifasta validation approach for enhanced accuracy
- Enhanced CIGAR filtering to reduce false positives
- Weighted consensus generation by supporting read counts
- Comprehensive HTML reporting with size analysis
- VCF output for standard variant format
- Size distribution analysis with ITD signature detection
- Dual validation criteria (support count + allele frequency)
- Support for hg19 and hg38 genome builds
- Configurable parameters for clinical validation
- IGV-compatible output files for visualization

### Features
- Direct detection from CIGAR insertion operations
- Advanced validation using multiple reference sequences
- Filters reads with large indels (>10bp) and poor alignment quality
- Generates weighted consensus sequences based on read support
- Creates comprehensive HTML reports with visualizations
- Exports results in standard VCF format
- Includes size analysis for ITD signature detection
- Supports parallel processing for performance
- Robust error handling and logging

### Dependencies
- Python 3.8+
- pysam >= 0.21.0
- matplotlib >= 3.5.0
- seaborn >= 0.11.0
- numpy >= 1.21.0
- pandas >= 1.3.0
- scipy >= 1.7.0
- samtools >= 1.15
- minimap2 >= 2.20
