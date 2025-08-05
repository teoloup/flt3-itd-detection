#!/usr/bin/env python3
"""
Read Size Analyzer Module
Analyzes read length distribution to detect potential ITDs/insertions
"""

import logging
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from collections import Counter
import seaborn as sns
from scipy import stats
from scipy.signal import find_peaks
import pysam

logger = logging.getLogger(__name__)

class ReadSizeAnalyzer:
    """Analyze read size distribution to detect ITD signatures"""
    
    def __init__(self, expected_wt_size: int = 336, min_itd_size: int = 30):
        """
        Initialize analyzer
        
        Args:
            expected_wt_size: Expected wild-type amplicon size
            min_itd_size: Minimum ITD size to consider significant
        """
        self.expected_wt_size = expected_wt_size
        self.min_itd_size = min_itd_size
        self.read_lengths = []
        self.analysis_results = {}
        
    def analyze_fastq_reads(self, fastq_file: str) -> Dict:
        """
        Analyze read lengths from FASTQ file
        
        Args:
            fastq_file: Path to FASTQ file
            
        Returns:
            Dictionary with analysis results
        """
        logger.info(f"Analyzing read lengths from {fastq_file}")
        
        read_lengths = []
        
        try:
            # Handle both regular and gzipped FASTQ files
            if fastq_file.endswith('.gz'):
                import gzip
                opener = gzip.open
                mode = 'rt'
            else:
                opener = open
                mode = 'r'
                
            with opener(fastq_file, mode) as f:
                line_count = 0
                for line in f:
                    line_count += 1
                    if line_count % 4 == 2:  # Sequence line
                        read_lengths.append(len(line.strip()))
                        
        except Exception as e:
            logger.error(f"Error reading FASTQ file: {e}")
            return {}
        
        self.read_lengths = read_lengths
        logger.info(f"Analyzed {len(read_lengths)} reads")
        
        return self._perform_analysis()
    
    def analyze_bam_reads(self, bam_file: str) -> Dict:
        """
        Analyze read lengths from BAM file
        
        Args:
            bam_file: Path to BAM file
            
        Returns:
            Dictionary with analysis results
        """
        logger.info(f"Analyzing read lengths from {bam_file}")
        
        read_lengths = []
        
        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                for read in bam.fetch():
                    if not read.is_unmapped and not read.is_secondary:
                        read_lengths.append(read.query_length)
                        
        except Exception as e:
            logger.error(f"Error reading BAM file: {e}")
            return {}
        
        self.read_lengths = read_lengths
        logger.info(f"Analyzed {len(read_lengths)} reads")
        
        return self._perform_analysis()
    
    def _perform_analysis(self) -> Dict:
        """Perform statistical analysis on read lengths"""
        if not self.read_lengths:
            return {}
        
        lengths = np.array(self.read_lengths)
        
        # Basic statistics
        stats_dict = {
            'total_reads': len(lengths),
            'mean_length': np.mean(lengths),
            'median_length': np.median(lengths),
            'std_length': np.std(lengths),
            'min_length': np.min(lengths),
            'max_length': np.max(lengths),
            'q25_length': np.percentile(lengths, 25),
            'q75_length': np.percentile(lengths, 75)
        }
        
        # Peak detection for ITD analysis
        peaks_analysis = self._detect_size_peaks(lengths)
        
        # ITD signature detection
        itd_signatures = self._detect_itd_signatures(lengths)
        
        self.analysis_results = {
            'basic_stats': stats_dict,
            'peaks': peaks_analysis,
            'itd_signatures': itd_signatures,
            'raw_lengths': lengths.tolist()
        }
        
        return self.analysis_results
    
    def _detect_size_peaks(self, lengths: np.ndarray) -> Dict:
        """Detect peaks in size distribution"""
        # Create histogram
        hist, bin_edges = np.histogram(lengths, bins=50)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Find peaks in histogram
        peaks, properties = find_peaks(hist, height=len(lengths) * 0.01, distance=5)
        
        peak_info = []
        for i, peak_idx in enumerate(peaks):
            peak_size = bin_centers[peak_idx]
            peak_height = hist[peak_idx]
            peak_prominence = properties.get('prominences', [0])[i] if 'prominences' in properties else 0
            
            # Classify peak
            if abs(peak_size - self.expected_wt_size) <= 20:
                peak_type = "wild_type"
            elif peak_size > self.expected_wt_size + self.min_itd_size:
                itd_size = int(peak_size - self.expected_wt_size)
                peak_type = f"potential_ITD_{itd_size}bp"
            else:
                peak_type = "other"
            
            peak_info.append({
                'size': peak_size,
                'height': peak_height,
                'prominence': peak_prominence,
                'type': peak_type,
                'frequency': peak_height / len(lengths)
            })
        
        # Sort by height (most prominent first)
        peak_info.sort(key=lambda x: x['height'], reverse=True)
        
        return {
            'peaks': peak_info,
            'histogram': {
                'bins': bin_centers.tolist(),
                'counts': hist.tolist()
            }
        }
    
    def _detect_itd_signatures(self, lengths: np.ndarray) -> Dict:
        """Detect ITD signatures in read length distribution"""
        
        # Count reads in different size ranges (ensuring no overlap)
        wt_range = (self.expected_wt_size - 30, self.expected_wt_size + 30)
        wt_reads = np.sum((lengths >= wt_range[0]) & (lengths <= wt_range[1]))
        
        # Look for reads longer than WT + minimum ITD size (non-overlapping with WT range)
        itd_threshold = max(wt_range[1] + 1, self.expected_wt_size + self.min_itd_size)
        potential_itd_reads = np.sum(lengths >= itd_threshold)
        
        # Calculate ITD frequency estimation
        total_reads = len(lengths)
        estimated_itd_frequency = potential_itd_reads / total_reads if total_reads > 0 else 0
        
        # Size classes for detailed analysis (ensuring no overlaps)
        size_classes = {
            'very_short': np.sum(lengths < self.expected_wt_size - 50),
            'short': np.sum((lengths >= self.expected_wt_size - 50) & (lengths < self.expected_wt_size - 30)),
            'wild_type': wt_reads,
            'small_insertion': np.sum((lengths > self.expected_wt_size + 30) & (lengths <= self.expected_wt_size + 60)),
            'medium_itd': np.sum((lengths > self.expected_wt_size + 60) & (lengths <= self.expected_wt_size + 150)),
            'large_itd': np.sum((lengths > self.expected_wt_size + 150) & (lengths <= self.expected_wt_size + 300)),
            'very_large_itd': np.sum(lengths > self.expected_wt_size + 300)
        }
        
        # ITD size distribution (using the corrected threshold)
        itd_reads = lengths[lengths >= itd_threshold]
        itd_sizes = itd_reads - self.expected_wt_size
        
        itd_size_stats = {}
        if len(itd_sizes) > 0:
            itd_size_stats = {
                'count': len(itd_sizes),
                'mean_size': np.mean(itd_sizes),
                'median_size': np.median(itd_sizes),
                'min_size': np.min(itd_sizes),
                'max_size': np.max(itd_sizes),
                'std_size': np.std(itd_sizes)
            }
        
        return {
            'estimated_itd_frequency': estimated_itd_frequency,
            'potential_itd_reads': potential_itd_reads,
            'wild_type_reads': wt_reads,
            'size_classes': size_classes,
            'itd_size_stats': itd_size_stats
        }
    
    def create_size_distribution_plot(self, output_file: str, title: str = "Read Length Distribution") -> bool:
        """
        Create comprehensive size distribution plot
        
        Args:
            output_file: Path to save the plot
            title: Plot title
            
        Returns:
            True if successful, False otherwise
        """
        if not self.read_lengths:
            logger.error("No read length data available for plotting")
            return False
        
        try:
            # Set up the plot style
            plt.style.use('default')
            sns.set_palette("husl")
            
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            fig.suptitle(title, fontsize=16, fontweight='bold')
            
            lengths = np.array(self.read_lengths)
            
            # 1. Histogram with density curve
            ax1 = axes[0, 0]
            ax1.hist(lengths, bins=50, alpha=0.7, density=True, color='skyblue', edgecolor='black')
            
            # Add density curve
            from scipy.stats import gaussian_kde
            density = gaussian_kde(lengths)
            xs = np.linspace(lengths.min(), lengths.max(), 200)
            ax1.plot(xs, density(xs), 'r-', linewidth=2, label='Density')
            
            # Mark expected WT size
            ax1.axvline(self.expected_wt_size, color='green', linestyle='--', linewidth=2, 
                       label=f'Expected WT ({self.expected_wt_size}bp)')
            
            ax1.set_xlabel('Read Length (bp)')
            ax1.set_ylabel('Density')
            ax1.set_title('Length Distribution with Density')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # 2. Boxplot
            ax2 = axes[0, 1]
            box_plot = ax2.boxplot(lengths, patch_artist=True)
            box_plot['boxes'][0].set_facecolor('lightcoral')
            ax2.axhline(self.expected_wt_size, color='green', linestyle='--', linewidth=2,
                       label=f'Expected WT ({self.expected_wt_size}bp)')
            ax2.set_ylabel('Read Length (bp)')
            ax2.set_title('Length Distribution Summary')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            
            # 3. Cumulative distribution
            ax3 = axes[1, 0]
            sorted_lengths = np.sort(lengths)
            cumulative = np.arange(1, len(sorted_lengths) + 1) / len(sorted_lengths)
            ax3.plot(sorted_lengths, cumulative, 'b-', linewidth=2)
            ax3.axvline(self.expected_wt_size, color='green', linestyle='--', linewidth=2,
                       label=f'Expected WT ({self.expected_wt_size}bp)')
            ax3.set_xlabel('Read Length (bp)')
            ax3.set_ylabel('Cumulative Probability')
            ax3.set_title('Cumulative Distribution')
            ax3.legend()
            ax3.grid(True, alpha=0.3)
            
            # 4. Size classes bar plot
            ax4 = axes[1, 1]
            if 'itd_signatures' in self.analysis_results:
                size_classes = self.analysis_results['itd_signatures']['size_classes']
                classes = list(size_classes.keys())
                counts = list(size_classes.values())
                
                bars = ax4.bar(classes, counts, color=['red', 'orange', 'green', 'yellow', 'blue', 'purple', 'brown'])
                ax4.set_xlabel('Size Class')
                ax4.set_ylabel('Read Count')
                ax4.set_title('Reads by Size Class')
                ax4.tick_params(axis='x', rotation=45)
                
                # Add count labels on bars
                for bar, count in zip(bars, counts):
                    if count > 0:
                        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01,
                                str(count), ha='center', va='bottom', fontweight='bold')
            
            ax4.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Size distribution plot saved to {output_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error creating size distribution plot: {e}")
            return False
    
    def create_itd_analysis_plot(self, output_file: str) -> bool:
        """Create focused ITD analysis plot"""
        if not self.read_lengths or 'itd_signatures' not in self.analysis_results:
            logger.error("No ITD analysis data available for plotting")
            return False
        
        try:
            fig, axes = plt.subplots(1, 2, figsize=(15, 6))
            fig.suptitle('ITD Analysis', fontsize=16, fontweight='bold')
            
            lengths = np.array(self.read_lengths)
            itd_sigs = self.analysis_results['itd_signatures']
            
            logger.debug(f"Creating ITD plot with {len(lengths)} reads")
            logger.debug(f"WT reads: {itd_sigs['wild_type_reads']}, ITD reads: {itd_sigs['potential_itd_reads']}")
            
            # 1. WT vs ITD comparison
            ax1 = axes[0]
            wt_reads = itd_sigs['wild_type_reads']
            itd_reads = itd_sigs['potential_itd_reads']
            total_reads = len(lengths)
            other_reads = max(0, total_reads - wt_reads - itd_reads)  # Ensure non-negative
            
            # Handle the case where classifications might overlap
            if wt_reads + itd_reads > total_reads:
                # Adjust to prevent overlap (prioritize ITD reads)
                logger.warning(f"Classification overlap detected: WT({wt_reads}) + ITD({itd_reads}) > Total({total_reads})")
                wt_reads = max(0, total_reads - itd_reads)
                other_reads = 0
            
            # Only include non-zero categories in pie chart
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
                ax1.set_title('Read Classification')
            else:
                ax1.text(0.5, 0.5, 'No classification data', ha='center', va='center', 
                        transform=ax1.transAxes, fontsize=14)
                ax1.set_title('Read Classification')
            
            # 2. ITD size distribution
            ax2 = axes[1]
            if itd_sigs['itd_size_stats'] and itd_sigs['itd_size_stats']['count'] > 0:
                # Use the same threshold logic as in _detect_itd_signatures
                wt_range = (self.expected_wt_size - 30, self.expected_wt_size + 30)
                itd_threshold = max(wt_range[1] + 1, self.expected_wt_size + self.min_itd_size)
                itd_reads_full = lengths[lengths >= itd_threshold]
                itd_sizes = itd_reads_full - self.expected_wt_size
                
                if len(itd_sizes) > 0:
                    ax2.hist(itd_sizes, bins=20, alpha=0.7, color='red', edgecolor='black')
                    ax2.axvline(np.mean(itd_sizes), color='blue', linestyle='--', linewidth=2,
                               label=f'Mean: {np.mean(itd_sizes):.1f}bp')
                    ax2.set_xlabel('ITD Size (bp)')
                    ax2.set_ylabel('Count')
                    ax2.set_title('ITD Size Distribution')
                    ax2.legend()
                    ax2.grid(True, alpha=0.3)
                else:
                    ax2.text(0.5, 0.5, 'No ITDs detected', ha='center', va='center', 
                            transform=ax2.transAxes, fontsize=14)
                    ax2.set_title('ITD Size Distribution')
            
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"ITD analysis plot saved to {output_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error creating ITD analysis plot: {e}")
            return False
    
    def generate_report(self) -> str:
        """Generate text report of size analysis"""
        if not self.analysis_results:
            return "No analysis results available."
        
        report = []
        report.append("=" * 60)
        report.append("READ SIZE ANALYSIS REPORT")
        report.append("=" * 60)
        
        # Basic statistics
        stats = self.analysis_results['basic_stats']
        report.append(f"\nBASIC STATISTICS:")
        report.append(f"  Total reads: {stats['total_reads']:,}")
        report.append(f"  Mean length: {stats['mean_length']:.1f} bp")
        report.append(f"  Median length: {stats['median_length']:.1f} bp")
        report.append(f"  Standard deviation: {stats['std_length']:.1f} bp")
        report.append(f"  Range: {stats['min_length']:.0f} - {stats['max_length']:.0f} bp")
        report.append(f"  Quartiles (Q25, Q75): {stats['q25_length']:.1f}, {stats['q75_length']:.1f} bp")
        
        # ITD signatures
        if 'itd_signatures' in self.analysis_results:
            itd_sigs = self.analysis_results['itd_signatures']
            report.append(f"\nITD ANALYSIS:")
            report.append(f"  Expected WT size: {self.expected_wt_size} bp")
            report.append(f"  Wild-type reads: {itd_sigs['wild_type_reads']:,} ({itd_sigs['wild_type_reads']/stats['total_reads']*100:.1f}%)")
            report.append(f"  Potential ITD reads: {itd_sigs['potential_itd_reads']:,} ({itd_sigs['estimated_itd_frequency']*100:.1f}%)")
            
            if itd_sigs['itd_size_stats']:
                itd_stats = itd_sigs['itd_size_stats']
                report.append(f"  ITD size range: {itd_stats['min_size']:.0f} - {itd_stats['max_size']:.0f} bp")
                report.append(f"  Mean ITD size: {itd_stats['mean_size']:.1f} bp")
                report.append(f"  Median ITD size: {itd_stats['median_size']:.1f} bp")
        
        # Peaks
        if 'peaks' in self.analysis_results:
            peaks = self.analysis_results['peaks']['peaks']
            if peaks:
                report.append(f"\nSIZE PEAKS DETECTED:")
                for i, peak in enumerate(peaks[:5]):  # Top 5 peaks
                    report.append(f"  Peak {i+1}: {peak['size']:.1f} bp ({peak['frequency']*100:.1f}%) - {peak['type']}")
        
        report.append("\n" + "=" * 60)
        
        return "\n".join(report)


def analyze_read_sizes(input_file: str, output_dir: str, expected_wt_size: int = 336) -> Dict:
    """
    Main function to analyze read sizes and generate plots
    
    Args:
        input_file: Path to FASTQ or BAM file
        output_dir: Directory to save output files
        expected_wt_size: Expected wild-type amplicon size
        
    Returns:
        Analysis results dictionary
    """
    analyzer = ReadSizeAnalyzer(expected_wt_size=expected_wt_size)
    
    # Determine file type and analyze
    if input_file.lower().endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        results = analyzer.analyze_fastq_reads(input_file)
    elif input_file.lower().endswith(('.bam', '.sam')):
        results = analyzer.analyze_bam_reads(input_file)
    else:
        logger.error(f"Unsupported file format: {input_file}")
        return {}
    
    if not results:
        return {}
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    base_name = Path(input_file).stem
    
    # Main size distribution plot
    dist_plot = output_path / f"{base_name}_size_distribution.png"
    analyzer.create_size_distribution_plot(str(dist_plot))
    
    # ITD-focused plot
    itd_plot = output_path / f"{base_name}_itd_analysis.png"
    analyzer.create_itd_analysis_plot(str(itd_plot))
    
    # Generate text report
    report = analyzer.generate_report()
    report_file = output_path / f"{base_name}_size_analysis_report.txt"
    with open(report_file, 'w') as f:
        f.write(report)
    
    logger.info(f"Size analysis complete. Results saved to {output_dir}")
    logger.info(f"Report:\n{report}")
    
    return results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze read size distribution for ITD detection")
    parser.add_argument("input_file", help="Input FASTQ or BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--expected-wt-size", type=int, default=336, help="Expected wild-type size (default: 336)")
    parser.add_argument("--min-itd-size", type=int, default=30, help="Minimum ITD size (default: 30)")
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    results = analyze_read_sizes(args.input_file, args.output, args.expected_wt_size)
    
    if results:
        print("Analysis completed successfully!")
    else:
        print("Analysis failed!")
