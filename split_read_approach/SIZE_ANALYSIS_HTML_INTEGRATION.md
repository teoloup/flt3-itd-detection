# Size Analysis Integration in HTML Reports

## Overview
The HTML report generator now includes comprehensive size distribution analysis alongside ITD detection results. This provides valuable diagnostic information about read length patterns that can indicate the presence of ITDs.

## Features Added

### 1. **Enhanced HTML Report Structure**
- New "Read Size Distribution Analysis" section
- Integrated between summary and ITD details
- Comprehensive statistical analysis
- Visual plots embedded directly in HTML

### 2. **Size Statistics Display**
- Total reads count
- Mean and median read lengths
- Length range (min-max)
- Standard deviation
- Quartile information (Q25, Q75)

### 3. **ITD Signature Analysis**
- Estimated ITD frequency percentage
- Potential ITD reads count
- Wild-type reads count
- Size class breakdown:
  - Very short reads
  - Short reads  
  - Wild-type reads
  - Small insertions
  - Medium ITDs
  - Large ITDs
  - Very large ITDs

### 4. **Visual Integration**
- **Size Distribution Plot**: Shows read length histogram with density curve, expected wild-type marker, boxplot, cumulative distribution, and size class bar chart
- **ITD Analysis Plot**: Pie chart showing read classification and ITD size distribution histogram
- All plots embedded as base64-encoded images for self-contained HTML

### 5. **Detailed Text Report**
- Complete size analysis report embedded in collapsible section
- Includes peak detection results
- Statistical summaries
- ITD size analysis

## Technical Implementation

### Modified Files
1. **`html_reporter.py`**:
   - Added `_create_size_analysis_section()` method
   - Updated `generate_html_report()` to accept size analysis parameters
   - Added CSS styling for size analysis components
   - Added `_encode_image()` for plot embedding

2. **`main_module.py`**:
   - Store size analysis results in class for HTML report
   - Pass size analysis data to HTML generator
   - Handle cases where size analysis fails or is unavailable

### New CSS Classes
- `.analysis-section`: Main container for size analysis
- `.size-stats`: Statistics grid layout
- `.stat-item`: Individual statistic display
- `.itd-signatures`: ITD signature analysis container
- `.size-plots`: Plot container with responsive layout
- `.analysis-plot`: Plot styling with borders and shadows
- `.detailed-report`: Collapsible detailed report section

## Usage

The integration is automatic when size analysis is performed. The HTML report will include:

1. **When size analysis succeeds**: Full size analysis section with statistics, plots, and detailed report
2. **When size analysis fails**: HTML report generated without size analysis section
3. **When size analysis disabled**: Standard HTML report (backward compatible)

## Benefits

### For Diagnostics
- **Early ITD Detection**: Size distribution shows ITD signatures before detailed validation
- **Quality Assessment**: Read length statistics indicate sequencing quality
- **Pattern Recognition**: Visual plots help identify unusual size distributions

### for Validation
- **Confirmation**: Size analysis results should correlate with ITD detection results
- **Troubleshooting**: Discrepancies between size analysis and ITD detection indicate potential issues
- **Confidence Building**: Consistent results across methods increase confidence

### For Reporting
- **Comprehensive View**: Single HTML report contains all analysis results
- **Visual Appeal**: Professional plots enhance report presentation
- **Self-Contained**: Embedded images make reports portable and shareable

## Example Output Structure

```html
<div class="analysis-section">
    <h2>Read Size Distribution Analysis</h2>
    
    <!-- Size Statistics -->
    <div class="size-stats">...</div>
    
    <!-- ITD Signatures -->
    <div class="itd-signatures">...</div>
    
    <!-- Visual Plots -->
    <div class="size-plots">
        <div class="plot-container">
            <h3>Size Distribution Plot</h3>
            <img src="data:image/png;base64,..." class="analysis-plot">
        </div>
        <div class="plot-container">
            <h3>ITD Analysis Plot</h3>
            <img src="data:image/png;base64,..." class="analysis-plot">
        </div>
    </div>
    
    <!-- Detailed Report -->
    <div class="detailed-report">...</div>
</div>
```

## Future Enhancements

1. **Interactive Plots**: Consider adding JavaScript-based interactive plots
2. **Peak Annotation**: Highlight significant peaks directly on plots
3. **Comparison Mode**: Compare size distributions across multiple samples
4. **Export Options**: Individual plot export functionality
5. **Mobile Optimization**: Enhanced responsive design for mobile viewing

## Testing

The integration has been tested with:
- ✅ Mock size analysis data
- ✅ Plot file embedding
- ✅ CSS styling verification
- ✅ Content validation checks
- ✅ Error handling for missing files
- ✅ Backward compatibility

All tests passed successfully, confirming robust integration.
