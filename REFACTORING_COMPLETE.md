## FLT3 ITD Detection Pipeline - Refactoring Complete! 🎉

### Summary of Changes

We have successfully refactored the FLT3 ITD detection pipeline to be completely **free of legacy split-read code** and focused on the new **CIGAR-based architecture**.

### ✅ What Was Accomplished

#### 1. **Clean Configuration System** 
- **File**: `config.py`
- **Status**: ✅ **COMPLETELY REFACTORED**
- **Changes**:
  - Removed ALL legacy split-read parameters (`split_method`, `split_overlap`, `num_splits`, etc.)
  - Streamlined to essential CIGAR-based parameters only
  - Clean `FLT3Config` class with only required parameters
  - Updated argument parser to show only user-requested parameters

#### 2. **Clean Main Module**
- **File**: `main_module.py` 
- **Status**: ✅ **COMPLETELY REFACTORED**
- **Changes**:
  - Renamed class from `FLT3ITDPipeline` to `FLT3ITDDetector`
  - Removed all legacy function calls (`cleanup_intermediate_files`, `SplitConfig`, `split_reads`, etc.)
  - Clean CIGAR-based pipeline with 4 simple steps
  - No legacy imports or dependencies
  - Proper error handling and logging

#### 3. **Essential Parameters Implemented**
All user-requested parameters are now available:

| Parameter | Description | Status |
|-----------|-------------|---------|
| `--min-itd-length` | Minimum ITD size in bp | ✅ |
| `--max-itd-length` | Maximum ITD size in bp | ✅ |
| `--min-allele-frequency` | Minimum allele frequency | ✅ |
| `--min-supporting-reads` | Minimum supporting reads | ✅ |
| `--bamout` | Export BAM/FASTA for IGV | ✅ |
| `--debug` | Debug mode | ✅ |
| `--genome` / `-g` | Human genome reference | ✅ |
| `--max-candidates` | Max candidates to process | ✅ |
| `--min-mapping-quality` | Min mapping quality | ✅ |
| `--bam` / `-b` | Input BAM file | ✅ |
| `--output` / `-o` | Output directory | ✅ |
| `--sample` / `-n` | Sample name | ✅ |

#### 4. **BAM Export Functionality**
- **Feature**: `--bamout` with dual-reference approach
- **Status**: ✅ **IMPLEMENTED**
- **Details**:
  - Creates `IGV_files/` directory
  - Exports dual-reference FASTA (WT + ITD sequences)
  - Exports supporting reads BAM files
  - Proper indexing for IGV viewing

### 🏗️ New Pipeline Architecture

```
FLT3ITDDetector
├── Step 1: Extract FLT3 reads (bam_extractor)
├── Step 2: CIGAR-based ITD detection (flt3_itd_pipeline)
├── Step 3: Generate outputs (VCF + HTML)
└── Step 4: Summary + cleanup
```

### 📋 Command Line Usage

The pipeline is now called with clean parameters:

```bash
python main_module.py \\
  --bam sample.bam \\
  --output results/ \\
  --sample my_sample \\
  --min-itd-length 10 \\
  --max-itd-length 200 \\
  --min-allele-frequency 0.01 \\
  --min-supporting-reads 3 \\
  --bamout \\
  --debug
```

### 🧹 Legacy Code Removal

**Removed Dependencies:**
- ❌ All split-read modules
- ❌ `cleanup_intermediate_files` function
- ❌ `SplitConfig` class
- ❌ `split_reads` function
- ❌ `detect_itds_from_pairs` function
- ❌ `generate_itd_candidates` function
- ❌ `validate_itd_candidates` function
- ❌ Size analysis dependencies
- ❌ All deprecated parameters

**Clean Imports:**
- ✅ `config` (clean version)
- ✅ `bam_extractor`
- ✅ `flt3_itd_pipeline` (CIGAR-based)
- ✅ `vcf_writer`
- ✅ `html_reporter`

### 🎯 Key Benefits

1. **Simplified Architecture**: 4-step pipeline vs. 7-step legacy pipeline
2. **Focused Parameters**: Only essential CIGAR detection parameters
3. **Clean Code**: No legacy function calls or deprecated imports
4. **Maintainable**: Clear separation of concerns
5. **User-Friendly**: Exactly the parameters you requested
6. **IGV Ready**: Dual-reference BAM export for visualization

### 🔧 Next Steps

The refactored pipeline is ready for use! The main requirements are:

1. **Install Dependencies**: `pysam`, `numpy`, `samtools`, `minimap2`
2. **Test with Real Data**: Run with actual BAM files
3. **Validate Results**: Compare with previous detection results

### 📁 File Status

| File | Status | Description |
|------|--------|-------------|
| `config.py` | ✅ **CLEAN** | Only CIGAR parameters |
| `main_module.py` | ✅ **CLEAN** | CIGAR-based pipeline |
| `flt3_itd_pipeline.py` | ✅ **READY** | Core CIGAR detection |
| `cigar_itd_detector.py` | ✅ **READY** | CIGAR analysis |
| `softclip_itd_detector.py` | ✅ **READY** | Soft-clip fallback |
| `dual_reference_validator.py` | ✅ **READY** | Validation |

The refactoring is **COMPLETE** and the pipeline is free of all legacy split-read code! 🚀
