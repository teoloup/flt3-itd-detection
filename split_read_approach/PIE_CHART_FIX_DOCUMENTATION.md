# Pie Chart Error Fix - "Wedge sizes must be non negative"

## Issue Description
The size analysis was failing with the error: `Error creating ITD analysis plot: Wedge sizes 'x' must be non negative values`

This error occurred when creating pie charts in the ITD analysis plot due to negative values being passed to matplotlib's pie chart function.

## Root Cause Analysis

### Problem 1: Classification Overlap
The original logic had overlapping size ranges:
- **Wild-type reads**: WT-30 to WT+30 (306-366 bp for 336bp WT)
- **ITD reads**: >= WT+min_itd_size (>=366 bp for min_itd_size=30)

This created an overlap in the range WT+20 to WT+30 (356-366 bp), causing some reads to be counted in both categories.

### Problem 2: Negative "Other" Category
The "other" category was calculated as:
```python
other_reads = total_reads - wt_reads - itd_reads
```

When `wt_reads + itd_reads > total_reads` (due to overlap), `other_reads` became negative, causing the pie chart error.

## Solution Implemented

### 1. **Fixed Classification Overlap**
```python
# Old logic (overlapping)
wt_reads = np.sum((lengths >= wt_range[0]) & (lengths <= wt_range[1]))  # WT-30 to WT+30
potential_itd_reads = np.sum(lengths >= self.expected_wt_size + self.min_itd_size)  # >= WT+30

# New logic (non-overlapping)
wt_range = (self.expected_wt_size - 30, self.expected_wt_size + 30)
wt_reads = np.sum((lengths >= wt_range[0]) & (lengths <= wt_range[1]))  # WT-30 to WT+30
itd_threshold = max(wt_range[1] + 1, self.expected_wt_size + self.min_itd_size)  # >= WT+31
potential_itd_reads = np.sum(lengths >= itd_threshold)
```

### 2. **Added Overlap Detection and Handling**
```python
# Handle classification overlaps
if wt_reads + itd_reads > total_reads:
    logger.warning(f"Classification overlap detected: WT({wt_reads}) + ITD({itd_reads}) > Total({total_reads})")
    wt_reads = max(0, total_reads - itd_reads)  # Prioritize ITD reads
    other_reads = 0
```

### 3. **Implemented Safe Pie Chart Creation**
```python
# Only include non-zero categories
sizes = []
labels = []
colors = []

if wt_reads > 0:
    sizes.append(wt_reads)
    labels.append(f'Wild-type\n({wt_reads:,})')
    colors.append('green')

if itd_reads > 0:
    sizes.append(itd_reads)
    labels.append(f'Potential ITD\n({itd_reads:,})')
    colors.append('red')

if other_reads > 0:
    sizes.append(other_reads)
    labels.append(f'Other\n({other_reads:,})')
    colors.append('gray')

# Only create pie chart if we have data
if sizes:
    ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
```

### 4. **Updated Size Classes to Avoid Overlaps**
```python
size_classes = {
    'very_short': np.sum(lengths < self.expected_wt_size - 50),
    'short': np.sum((lengths >= self.expected_wt_size - 50) & (lengths < self.expected_wt_size - 30)),  # Fixed gap
    'wild_type': wt_reads,
    'small_insertion': np.sum((lengths > self.expected_wt_size + 30) & (lengths <= self.expected_wt_size + 60)),  # Fixed gap
    'medium_itd': np.sum((lengths > self.expected_wt_size + 60) & (lengths <= self.expected_wt_size + 150)),
    'large_itd': np.sum((lengths > self.expected_wt_size + 150) & (lengths <= self.expected_wt_size + 300)),
    'very_large_itd': np.sum(lengths > self.expected_wt_size + 300)
}
```

### 5. **Added Debug Logging**
```python
logger.debug(f"Creating ITD plot with {len(lengths)} reads")
logger.debug(f"WT reads: {itd_sigs['wild_type_reads']}, ITD reads: {itd_sigs['potential_itd_reads']}")
```

## Impact on Analysis

### Before Fix
- **Error**: "Wedge sizes 'x' must be non negative values"
- **Analysis**: Failed to create ITD analysis plots
- **Data**: Some reads double-counted due to overlap

### After Fix
- **Success**: ITD analysis plots generate successfully
- **Accuracy**: No double-counting of reads
- **Robustness**: Handles edge cases gracefully

## Test Results from Your Data

Based on your logs:
- **Total reads**: 48,876
- **Wild-type reads**: 34,385 (70.4%)
- **Potential ITD reads**: 14,493 (29.7%)
- **ITD frequency**: 29.7%
- **Mean ITD size**: 62.0 bp

The fix ensures these statistics are calculated correctly without overlaps and that the pie chart can display them properly.

## Files Modified

1. **`read_size_analyzer.py`**:
   - Fixed `_detect_itd_signatures()` method
   - Updated `create_itd_analysis_plot()` method
   - Added overlap detection and handling
   - Improved error logging

## Benefits

1. **Reliability**: No more pie chart errors
2. **Accuracy**: Correct classification without overlaps
3. **Robustness**: Handles edge cases and provides warnings
4. **Debugging**: Better logging for troubleshooting

The fix ensures that your size analysis will complete successfully and provide accurate ITD frequency estimates for integration into the HTML reports.
