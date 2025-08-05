# FLT3 ITD Detection Pipeline Documentation

## Algorithm Overview

The FLT3 ITD Detection Pipeline uses a CIGAR-based approach to identify Internal Tandem Duplications (ITDs) in FLT3 gene sequencing data from Oxford Nanopore technologies.

### Pipeline Steps

1. **Read Extraction and Trimming**
   - Extract reads mapping to FLT3 region
   - Trim primer sequences
   - Filter by mapping quality and read length

2. **Size Distribution Analysis**
   - Analyze read length distribution
   - Identify potential ITD signatures
   - Estimate ITD frequency

3. **CIGAR-based ITD Detection**
   - Parse CIGAR strings for insertion operations
   - Group insertions by position and sequence
   - Generate weighted consensus sequences

4. **Multifasta Validation**
   - Create reference with WT + all ITD candidates
   - Re-align reads to multifasta reference
   - Apply CIGAR filtering (>10bp indels)
   - Validate with dual criteria (support + frequency)

5. **Output Generation**
   - Generate VCF files
   - Create HTML reports with visualizations
   - Export summary statistics

### Key Features

#### Weighted Consensus Generation
- Consensus sequences weighted by supporting read counts
- More accurate than simple sequence counting
- Reduces impact of sequencing errors

#### Enhanced CIGAR Filtering
- Filters reads with large indels (>10bp)
- Removes excessive soft clipping (>50bp)
- Preserves reads with small sequencing errors

#### Dual Validation Criteria
- Minimum supporting reads threshold
- Minimum allele frequency threshold
- Both criteria must be met for validation

#### Multifasta Approach
- Single alignment to comprehensive reference
- Better read assignment for multiple ITDs
- Reduces false positives from misalignment

## Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_supporting_reads` | 3 | Minimum reads supporting an ITD |
| `min_allele_frequency` | 0.01 | Minimum allele frequency (1%) |
| `max_itd_length` | 500 | Maximum ITD length to consider |
| `min_itd_length` | 6 | Minimum ITD length to consider |
| `min_mapping_quality` | 20 | Minimum mapping quality for reads |
| `min_read_length` | 100 | Minimum read length to include |

## Output Files

### VCF Format
Standard variant call format with ITD annotations:
- Position and length information
- Allele frequency calculations
- Supporting read counts
- Confidence scores

### HTML Reports
Comprehensive visualization including:
- Summary statistics
- Size distribution plots
- ITD sequence context
- Validation metrics

### Summary Files
Text-based summary with:
- Detection parameters
- Runtime statistics
- ITD details and sequences

## Performance Considerations

- **Memory Usage**: Scales with number of reads in FLT3 region
- **Runtime**: Primarily limited by alignment step
- **Parallelization**: Use `--threads` parameter for large datasets
- **Disk Space**: Temporary files created during processing

## Troubleshooting

### Common Issues

1. **No ITDs Detected**
   - Check input BAM file quality
   - Verify FLT3 region coverage
   - Adjust sensitivity parameters

2. **High False Positive Rate**
   - Increase `min_supporting_reads`
   - Increase `min_allele_frequency`
   - Check read quality metrics

3. **Performance Issues**
   - Increase `--threads` parameter
   - Ensure sufficient memory
   - Check disk space for temporary files

### Debug Mode
Use `--debug` flag for detailed logging and temporary file retention.

## References

1. Oxford Nanopore Technologies sequencing
2. FLT3 gene structure and ITD characteristics
3. CIGAR string format specification
4. VCF format specification
