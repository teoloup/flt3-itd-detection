#!/usr/bin/env python3
"""
Test script to verify ITD size correction fixes
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from config import Config
from main_module import FLT3ITDAnalyzer
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_itd_size_correction():
    """Test ITD detection with size correction fixes"""
    
    # Use your existing data
    bam_file = "test_data.bam"  # Replace with your actual BAM file
    reference_file = "reference.fa"  # Replace with your actual reference
    
    if not os.path.exists(bam_file):
        logger.error(f"BAM file not found: {bam_file}")
        logger.info("Please update the script with the correct path to your test BAM file")
        return
    
    # Test with strict size limits
    config = Config(
        max_itd_length=200,  # Your problematic setting
        min_itd_length=15,
        min_supporting_reads=2,
        min_allele_frequency=0.01,
        output_dir="test_size_correction",
        keep_temp=True,
        debug=True
    )
    
    logger.info("Testing ITD detection with max_itd_length=200bp")
    logger.info("This should now find smaller ITDs around 69bp instead of only large ones")
    
    try:
        analyzer = FLT3ITDAnalyzer(config)
        results = analyzer.analyze_file(bam_file, reference_file)
        
        logger.info("\n" + "="*60)
        logger.info("ITD DETECTION RESULTS WITH SIZE CORRECTION")
        logger.info("="*60)
        
        if results and len(results) > 0:
            logger.info(f"Found {len(results)} ITD candidates:")
            
            for i, result in enumerate(results, 1):
                itd_length = result.duplication_length
                allele_freq = result.allele_frequency
                support = result.itd_coverage
                confidence = result.validation_confidence
                
                logger.info(f"  ITD {i}: {itd_length}bp "
                           f"(AF={allele_freq:.3f}, support={support}, conf={confidence:.2f})")
                
                # Check if we're getting reasonable sizes
                if 50 <= itd_length <= 100:
                    logger.info(f"    ✓ ITD size in expected range for your 69bp peak!")
                elif itd_length > 200:
                    logger.warning(f"    ⚠ ITD still oversized - size correction may need tuning")
                else:
                    logger.info(f"    → ITD size looks reasonable")
        else:
            logger.warning("No ITDs found - this might indicate the correction is too aggressive")
            logger.info("Try adjusting parameters or check if input data is correct")
            
        # Check size distribution
        logger.info("\n" + "="*60)
        logger.info("SIZE DISTRIBUTION ANALYSIS")
        logger.info("="*60)
        
        # The size analyzer should show the 69bp peak
        from read_size_analyzer import ReadSizeAnalyzer
        
        size_analyzer = ReadSizeAnalyzer(bam_file)
        size_results = size_analyzer.analyze_read_lengths()
        
        if size_results:
            logger.info("Read length statistics:")
            logger.info(f"  Median: {size_results['median_length']:.1f}bp")
            logger.info(f"  Mean: {size_results['mean_length']:.1f}bp")
            
            if 'peaks' in size_results:
                logger.info("  Detected peaks:")
                for peak in size_results['peaks'][:3]:  # Top 3 peaks
                    logger.info(f"    {peak['size']}bp (strength: {peak['strength']:.2f})")
                    
                    if 60 <= peak['size'] <= 80:
                        logger.info(f"      ✓ This matches your expected 69bp ITD!")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_itd_size_correction()
