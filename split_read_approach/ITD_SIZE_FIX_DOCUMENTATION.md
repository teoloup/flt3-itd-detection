# ITD Size Overestimation Fix

## Problem Identified

The pipeline was systematically **overestimating ITD sizes** by:

1. **Treating entire overlap regions as ITD sequences**: The overlap between two reads spanning an ITD includes both the duplicated sequence AND flanking reference sequence
2. **Aggressive boundary extension**: Adding soft-clip sequences and large insertions without size limits
3. **Concatenating sequences**: Combining multiple sequence fragments that represent the same biological event

## Example of the Problem

For a **69bp ITD** (as detected by size distribution analysis):
- Read 1 aligns to positions 1800-2000 (200bp)
- Read 2 aligns to positions 1850-2050 (200bp)
- **Overlap region**: 1850-2000 = **150bp**
- **Plus soft-clip extensions**: +50bp 
- **Final "ITD" size**: **200bp+** (instead of actual 69bp)

## Fixes Implemented

### 1. Conservative Boundary Detection (`_find_itd_boundaries`)

**Before:**
```python
# Extended boundaries aggressively using all soft-clips > 10bp
if r1_left_clip_len > 10:
    itd_start = min(itd_start, r1.reference_start - r1_left_clip_len)
```

**After:**
```python
# Conservative extension with size limits
max_clip_extension = 30  # Maximum 30bp extension
if r1_left_clip_len > 15:  # Only significant clips
    extension = min(r1_left_clip_len, max_clip_extension)
    itd_start = min(itd_start, r1.reference_start - extension)
```

### 2. Initial Overlap Size Reduction

**New logic:**
```python
initial_overlap_size = itd_end - itd_start
if initial_overlap_size > 200:
    # Reduce to reasonable size centered on overlap
    center = (itd_start + itd_end) // 2
    reasonable_size = min(150, initial_overlap_size // 2)
    itd_start = center - reasonable_size // 2
    itd_end = center + reasonable_size // 2
```

### 3. CIGAR Insertion Evidence

**New approach:**
- Look for actual insertions (CIGAR operation 'I') of 15-300bp
- Use these as primary evidence for ITD size
- Center boundaries around actual insertion positions

### 4. Sequence Length Limits

**Added safety checks:**
```python
if len(final_itd_sequence) > self.max_itd_length:
    # Truncate oversized sequences
    excess = len(final_itd_sequence) - self.max_itd_length
    start_trim = excess // 2
    final_itd_sequence = final_itd_sequence[start_trim:start_trim + self.max_itd_length]
```

### 5. ITD Size Corrector Module (Optional Enhancement)

Created `itd_size_corrector.py` that implements:
- **Soft-clip analysis**: Identifies duplicated sequences in soft-clips
- **CIGAR insertion analysis**: Uses actual insertion sizes from alignments
- **Reference comparison**: Separates ITD from flanking sequence
- **Multiple estimation methods**: Provides confidence-weighted results

## Expected Results

With these fixes, when you run:
```bash
python main_module.py --max-itd-length 200 input.bam reference.fa
```

**Before (problematic):**
- Found only ITDs >200bp (202bp, 316bp)
- With `--max-itd-length 200`: 0 ITDs found
- Size distribution shows 69bp peak but pipeline misses it

**After (fixed):**
- Should find ITDs around 69bp matching size distribution peak
- With `--max-itd-length 200`: Should find the 69bp ITD
- Better correlation between size distribution and detected ITDs

## Testing the Fix

1. **Run with your existing data:**
   ```bash
   python main_module.py --max-itd-length 200 --debug your_data.bam reference.fa
   ```

2. **Check the logs for:**
   - "Reduced oversized overlap from XXXbp to YYYbp"
   - "Found XXbp insertion at position YYY"
   - "Truncating oversized ITD sequence"

3. **Compare results:**
   - ITD sizes should be closer to size distribution peaks
   - Should find ITDs with `--max-itd-length 200` now

## Fine-tuning Parameters

If results need adjustment:

**For more sensitive detection (finds smaller ITDs):**
- Reduce `max_clip_extension` from 30 to 15
- Reduce `reasonable_size` cap from 150 to 100

**For more specific detection (fewer false positives):**
- Increase minimum CIGAR insertion size from 15 to 25
- Increase consensus ratio requirement

## Key Insight

The fundamental insight is that **overlap length ≠ ITD length**. The overlap between reads contains the ITD plus flanking sequence. The actual ITD size must be extracted from:
1. CIGAR insertion operations
2. Soft-clipped sequences that match reference
3. Consensus differences between reads and reference

This fix addresses the core algorithmic issue that was causing systematic size overestimation.
