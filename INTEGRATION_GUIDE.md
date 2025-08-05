# FLT3 ITD Pipeline Integration Guide

## New Architecture Implementation

This document describes how to integrate the new CIGAR-based ITD detection architecture with your existing FLT3 pipeline.

### Architecture Overview

**OLD APPROACH (Split-read based):**
```
BAM → Read pairs → Split overlaps → ITD candidates → Validation → Results
```

**NEW APPROACH (Direct CIGAR analysis):**
```
BAM → CIGAR insertions → ITD candidates → Soft-clip fallback → Dual-reference validation → Results
```

### Key Components Implemented

1. **cigar_itd_detector.py** - Primary detection method using direct CIGAR insertion analysis
2. **softclip_itd_detector.py** - Secondary detection for reads with large soft-clipped regions  
3. **dual_reference_validator.py** - Validation using dual-reference alignment approach
4. **flt3_itd_pipeline.py** - Main orchestrator integrating all components

### Integration Steps

#### Step 1: Modify Main Module

Update your existing `main_module.py` to use the new pipeline:

```python
# OLD CODE (replace this)
from itd_generator import ITDGenerator
from reference_validator import ReferenceValidator

# NEW CODE (use this instead)
from flt3_itd_pipeline import run_flt3_itd_pipeline

def detect_flt3_itds(bam_file, reference_file, reference_sequence, args):
    """Updated ITD detection using new architecture"""
    
    # Run new pipeline
    result = run_flt3_itd_pipeline(
        bam_file=bam_file,
        reference_sequence=reference_sequence,
        reference_file=reference_file,
        min_itd_length=args.min_itd_length,
        max_itd_length=args.max_itd_length,
        min_support=args.min_support,
        min_softclip_length=getattr(args, 'min_softclip_length', 50),
        enable_softclip_fallback=getattr(args, 'enable_softclip_fallback', True),
        output_prefix=args.output_prefix
    )
    
    # Convert results to expected format for downstream modules
    final_itds = []
    for candidate in result.validated_candidates:
        itd_data = {
            'sequence': candidate.sequence,
            'length': candidate.length,
            'position': candidate.position,
            'supporting_reads': candidate.supporting_reads,
            'confidence': candidate.confidence,
            'detection_method': candidate.support_type,
            'allele_frequency': 0.0  # Will be updated by validation results
        }
        final_itds.append(itd_data)
    
    # Update with validation metrics
    for i, validation in enumerate(result.validation_results):
        if i < len(final_itds) and validation.is_valid:
            final_itds[i]['allele_frequency'] = validation.allele_frequency
            final_itds[i]['segregation_quality'] = validation.segregation_quality
    
    return final_itds, result.detection_summary
```

#### Step 2: Preserve Existing Modules

**KEEP THESE MODULES** (they work well):
- `bam_extractor.py` - FLT3 read extraction and trimming
- `vcf_writer.py` - VCF output generation  
- `html_reporter.py` - HTML report generation

**REPLACE THESE MODULES:**
- `itd_generator.py` → `cigar_itd_detector.py` + `softclip_itd_detector.py`
- `reference_validator.py` → `dual_reference_validator.py`

#### Step 3: Update Command Line Arguments

Add new parameters to your argument parser:

```python
parser.add_argument('--min-softclip-length', type=int, default=50,
                   help='Minimum soft-clip length for secondary detection')
parser.add_argument('--enable-softclip-fallback', action='store_true', default=True,
                   help='Enable soft-clip analysis if CIGAR detection insufficient')
parser.add_argument('--no-softclip-fallback', action='store_true',
                   help='Disable soft-clip fallback analysis')
```

#### Step 4: Update Dependencies

Ensure you have all required dependencies:

```bash
# Required packages
pip install pysam numpy
# minimap2 should be in PATH
```

### Workflow Comparison

#### Old Workflow
```python
# 1. Extract FLT3 reads
reads = bam_extractor.extract_flt3_reads(bam_file)

# 2. Generate ITD candidates from overlaps
itd_generator = ITDGenerator(reference_sequence)
candidates = itd_generator.find_itds_from_overlaps(reads)

# 3. Validate using reference alignment
validator = ReferenceValidator(reference_file)
validated = validator.validate_candidates(candidates)

# 4. Generate outputs
vcf_writer.write_vcf(validated, output_vcf)
html_reporter.generate_report(validated, output_html)
```

#### New Workflow
```python
# 1. Extract FLT3 reads (UNCHANGED)
reads = bam_extractor.extract_flt3_reads(bam_file)

# 2. Run integrated detection pipeline
result = run_flt3_itd_pipeline(
    bam_file=trimmed_bam,  # From bam_extractor
    reference_sequence=reference_sequence,
    reference_file=reference_file,
    min_itd_length=args.min_itd_length,
    max_itd_length=args.max_itd_length,
    min_support=args.min_support
)

# 3. Generate outputs (UNCHANGED)
vcf_writer.write_vcf(result.validated_candidates, output_vcf)
html_reporter.generate_report(result.validated_candidates, output_html)
```

### Key Improvements

1. **Direct CIGAR Analysis**: No more artificial read splitting
2. **Size Accuracy**: ITD sizes from native read CIGAR operations
3. **Dual Detection**: CIGAR primary + soft-clip fallback
4. **Better Validation**: Dual-reference approach with natural read segregation
5. **Comprehensive Logging**: Detailed detection and validation metrics

### Configuration Options

#### Conservative Settings (High Precision)
```python
result = run_flt3_itd_pipeline(
    min_itd_length=20,          # Slightly higher minimum
    max_itd_length=400,         # Exclude very large ITDs
    min_support=5,              # Require more supporting reads
    min_softclip_length=60,     # Higher soft-clip threshold
    enable_softclip_fallback=False  # CIGAR-only detection
)
```

#### Sensitive Settings (High Recall)
```python
result = run_flt3_itd_pipeline(
    min_itd_length=12,          # Lower minimum size
    max_itd_length=600,         # Include larger ITDs
    min_support=2,              # Lower support requirement
    min_softclip_length=40,     # Lower soft-clip threshold
    enable_softclip_fallback=True   # Use both detection methods
)
```

### Testing the Integration

1. **Unit Testing**: Test each component individually
```bash
python cigar_itd_detector.py test_data.bam reference.fasta --debug
python softclip_itd_detector.py test_data.bam reference.fasta --debug
python dual_reference_validator.py test_candidates.json test_data.bam reference.fasta --debug
```

2. **Pipeline Testing**: Test the complete pipeline
```bash
python flt3_itd_pipeline.py test_data.bam reference.fasta --debug --output-prefix test_run
```

3. **Comparison Testing**: Compare old vs new results
```bash
# Run old pipeline
python main_module.py --old-method test_data.bam > old_results.txt

# Run new pipeline  
python main_module.py --new-method test_data.bam > new_results.txt

# Compare outputs
diff old_results.txt new_results.txt
```

### Migration Checklist

- [ ] New detection modules implemented
- [ ] Pipeline orchestrator created
- [ ] Main module updated to use new architecture
- [ ] Command line arguments updated
- [ ] Dependencies installed
- [ ] Unit tests pass
- [ ] Integration tests pass
- [ ] Comparison tests show improved accuracy
- [ ] Documentation updated

### Expected Outcomes

Based on your size distribution analysis showing a 69bp peak but old pipeline finding only >200bp ITDs:

**Before (Split-read approach):**
- ITDs found: 202bp, 316bp, etc.
- Systematic size overestimation
- Missing smaller ITDs (69bp peak)
- `--max-itd-length 200` yields 0 results

**After (CIGAR approach):**
- ITDs found: Should include 69bp ± 10bp range
- Accurate size estimation from native CIGAR
- Better detection of smaller ITDs
- Size distribution should match expected peaks

### Troubleshooting

If you encounter issues:

1. **Import Errors**: Ensure all new modules are in the same directory
2. **pysam Issues**: Install with `pip install pysam` or `conda install -c bioconda pysam`
3. **minimap2 Missing**: Install with `conda install -c bioconda minimap2`
4. **No ITDs Found**: Check `--debug` output for detailed analysis
5. **Size Mismatches**: Verify CIGAR parsing is working correctly

### Performance Notes

- **CIGAR detection**: Fast, processes reads directly
- **Soft-clip analysis**: Slower, involves realignment
- **Dual-reference validation**: Moderate, creates temporary references
- **Overall**: Should be comparable or faster than old split-read approach

The new architecture should resolve the size overestimation issue and detect the 69bp ITDs that your size distribution analysis indicated were present but missed by the old approach.
