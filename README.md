# FLT3 ITD Detection Pipeline

A comprehensive pipeline for detecting FLT3 Internal Tandem Duplications (ITDs) from Oxford Nanopore sequencing data using CIGAR-based analysis.

## Features

- **CIGAR-based ITD Detection**: Direct detection from CIGAR insertion operations
- **Multifasta Validation**: Advanced validation using multiple reference sequences
- **Enhanced CIGAR Filtering**: Filters reads with large indels and poor alignment quality
- **Consensus Generation**: Weighted consensus by supporting read counts
- **Comprehensive Reporting**: HTML reports with size analysis and visualizations
- **VCF Output**: Standard variant format for downstream analysis

## Requirements

### System Dependencies
- Python 3.8+
- samtools (>= 1.15)
- minimap2 (>= 2.20)

### Python Dependencies
```bash
pip install -r requirements.txt
```

## Installation

```bash
# Clone the repository
git clone https://github.com/teoloup/flt3-itd-detection.git
cd flt3-itd-detection

# Install dependencies
pip install -r requirements.txt

# Verify external tools are installed
samtools --version
minimap2 --version
```

## Quick Start

```bash
# Basic usage
python main_module.py -b input.bam -o output_dir/ -n SampleName

# Advanced usage with custom parameters
python main_module.py \
    -b input.bam \
    -o output_dir/ \
    -n SampleName \
    --min-supporting-reads 10 \
    --min-allele-frequency 0.05 \
    --max-itd-length 150 \
    --threads 16
```

## Usage

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-b, --bam` | Input BAM file (required) | - |
| `-o, --output` | Output directory (required) | - |
| `-n, --name` | Sample name | Sample |
| `--min-supporting-reads` | Minimum supporting reads | 3 |
| `--min-allele-frequency` | Minimum allele frequency | 0.01 |
| `--max-itd-length` | Maximum ITD length | 500 |
| `--threads` | Number of threads | 4 |
| `--debug` | Enable debug logging | False |

### Example Commands

```bash
# Standard analysis
python main_module.py -b sample.bam -o results/ -n Sample001

# High-confidence detection
python main_module.py -b sample.bam -o results/ -n Sample001 \
    --min-supporting-reads 10 --min-allele-frequency 0.05

# Large-scale analysis
python main_module.py -b sample.bam -o results/ -n Sample001 \
    --threads 32 --max-itd-length 200
```

## Output Files

- `{sample}_ITDs.vcf` - Detected ITDs in VCF format
- `{sample}_ITD_report.html` - Comprehensive HTML report
- `{sample}_summary.txt` - Text summary of results
- `{sample}_itd_detection_summary.json` - Detailed JSON results

## Algorithm Overview

1. **Read Extraction**: Extract FLT3-mapped reads with primer trimming
2. **CIGAR Analysis**: Detect ITDs from CIGAR insertion operations
3. **Consensus Generation**: Generate weighted consensus sequences
4. **Multifasta Validation**: Validate candidates using multiple references
5. **CIGAR Filtering**: Filter reads with poor alignment quality
6. **Report Generation**: Create comprehensive reports and visualizations

## Citation

If you use this pipeline in your research, please cite:

```
[Your paper citation here]
```

## License

MIT License - see LICENSE file for details.

## Contributing

Please see CONTRIBUTING.md for guidelines on contributing to this project.

## Support

For issues and questions:
- GitHub Issues: [Repository Issues](https://github.com/teoloup/flt3-itd-detection/issues)
- Email: teoloup@gmail.com
