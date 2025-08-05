#!/usr/bin/env python3
"""
Test script to verify HTML report generation with size analysis
"""

import tempfile
import os
from pathlib import Path
import sys

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def create_mock_size_analysis_results():
    """Create mock size analysis results for testing"""
    return {
        'basic_stats': {
            'total_reads': 15000,
            'mean_length': 345.7,
            'median_length': 342.0,
            'std_length': 25.3,
            'min_length': 280,
            'max_length': 450,
            'q25_length': 330.0,
            'q75_length': 355.0
        },
        'itd_signatures': {
            'estimated_itd_frequency': 0.12,  # 12%
            'potential_itd_reads': 1800,
            'wild_type_reads': 12500,
            'size_classes': {
                'very_short': 50,
                'short': 150,
                'wild_type': 12500,
                'small_insertion': 300,
                'medium_itd': 1200,
                'large_itd': 600,
                'very_large_itd': 200
            },
            'itd_size_stats': {
                'count': 1800,
                'mean_size': 45.5,
                'median_size': 39.0,
                'min_size': 21,
                'max_size': 120,
                'std_size': 18.2
            }
        },
        'peaks': {
            'peaks': [
                {
                    'size': 336.0,
                    'height': 8500,
                    'prominence': 7200,
                    'type': 'wild_type',
                    'frequency': 0.567
                },
                {
                    'size': 375.0,
                    'height': 950,
                    'prominence': 850,
                    'type': 'potential_ITD_39bp',
                    'frequency': 0.063
                },
                {
                    'size': 381.0,
                    'height': 720,
                    'prominence': 650,
                    'type': 'potential_ITD_45bp',
                    'frequency': 0.048
                }
            ]
        }
    }

def create_mock_size_analysis_files(temp_dir):
    """Create mock size analysis plot files"""
    # Create dummy plot files (empty PNG files with headers)
    png_header = b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x02\x00\x00\x00\x90wS\xde\x00\x00\x00\tpHYs\x00\x00\x0b\x13\x00\x00\x0b\x13\x01\x00\x9a\x9c\x18\x00\x00\x00\x12IDATx\x9cc```bPPP\x00\x02\xd2\x01\x01\x00\x00\x05\x00\x01\r\n-\xdb\x00\x00\x00\x00IEND\xaeB`\x82'
    
    size_dir = Path(temp_dir) / "size_analysis"
    size_dir.mkdir(exist_ok=True)
    
    # Create mock plot files
    with open(size_dir / "test_sample_size_distribution.png", "wb") as f:
        f.write(png_header)
    
    with open(size_dir / "test_sample_itd_analysis.png", "wb") as f:
        f.write(png_header)
    
    # Create mock report file
    report_content = """========================================================
READ SIZE ANALYSIS REPORT
========================================================

BASIC STATISTICS:
  Total reads: 15,000
  Mean length: 345.7 bp
  Median length: 342.0 bp
  Standard deviation: 25.3 bp
  Range: 280 - 450 bp
  Quartiles (Q25, Q75): 330.0, 355.0 bp

ITD ANALYSIS:
  Expected WT size: 336 bp
  Wild-type reads: 12,500 (83.3%)
  Potential ITD reads: 1,800 (12.0%)
  ITD size range: 21 - 120 bp
  Mean ITD size: 45.5 bp
  Median ITD size: 39.0 bp

SIZE PEAKS DETECTED:
  Peak 1: 336.0 bp (56.7%) - wild_type
  Peak 2: 375.0 bp (6.3%) - potential_ITD_39bp
  Peak 3: 381.0 bp (4.8%) - potential_ITD_45bp

========================================================"""
    
    with open(size_dir / "test_sample_size_analysis_report.txt", "w") as f:
        f.write(report_content)
    
    return str(size_dir)

def test_html_generation():
    """Test HTML report generation with size analysis"""
    try:
        from html_reporter import generate_html_report
        
        # Create temporary directory for output
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create mock size analysis results and files
            size_analysis_results = create_mock_size_analysis_results()
            size_analysis_dir = create_mock_size_analysis_files(temp_dir)
            
            # Test HTML generation
            html_file = Path(temp_dir) / "test_report.html"
            
            # Generate HTML with size analysis (no validation results for simplicity)
            generate_html_report(
                validation_results=[],  # Empty for this test
                sample_name="TEST_SAMPLE_SIZE_ANALYSIS",
                reference_sequence="A" * 336,  # Mock reference
                total_reads=15000,
                output_file=str(html_file),
                size_analysis_results=size_analysis_results,
                size_analysis_dir=size_analysis_dir
            )
            
            # Verify HTML file was created
            if html_file.exists():
                print("✅ HTML report with size analysis generated successfully!")
                print(f"📄 Report saved to: {html_file}")
                
                # Read and check content
                with open(html_file, 'r') as f:
                    content = f.read()
                
                # Check for key sections
                checks = [
                    "Read Size Distribution Analysis" in content,
                    "Size Statistics" in content,
                    "ITD Signature Analysis" in content,
                    "Estimated ITD Frequency" in content,
                    "12.0%" in content,  # ITD frequency
                    "15,000" in content,  # Total reads
                    "Size Distribution Plot" in content,
                    "ITD Analysis Plot" in content,
                    "Detailed Size Analysis Report" in content,
                    "data:image/png;base64" in content  # Embedded images
                ]
                
                passed = sum(checks)
                total = len(checks)
                
                print(f"\n📊 Content verification: {passed}/{total} checks passed")
                
                if passed == total:
                    print("🎉 All content checks passed!")
                    return True
                else:
                    print("⚠️  Some content checks failed")
                    return False
            else:
                print("❌ HTML file was not created")
                return False
                
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("Testing HTML report generation with size analysis...\n")
    
    success = test_html_generation()
    
    if success:
        print("\n✅ Test completed successfully!")
        print("The HTML reporter now supports size analysis integration.")
    else:
        print("\n❌ Test failed!")
        print("There may be issues with the HTML report generation.")
