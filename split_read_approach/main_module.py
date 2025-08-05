#!/usr/bin/env python3
"""
FLT3 ITD Detection Main Module
Orchestrates the complete ITD detection pipeline using split-read approach
"""

import sys
import logging
from pathlib import Path
import shutil
from datetime import datetime
from typing import List
import pysam

# Import all modules
from config import parse_arguments, FLT3Config
from bam_extractor import extract_flt3_reads
from read_splitter import split_reads, SplitConfig
from paired_aligner import detect_itds_from_pairs
from itd_generator import generate_itd_candidates, ITDCandidate
from reference_validator import validate_itd_candidates, ValidationResult, cleanup_intermediate_files
from vcf_writer import write_vcf_results
from html_reporter import generate_html_report





# Set up logging
def setup_logging(debug: bool = False):
    """Configure logging"""
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

class FLT3ITDPipeline:
    """Main pipeline for FLT3 ITD detection"""
    
    def __init__(self, config: FLT3Config):
        self.config = config
        self.logger = setup_logging(config.debug)
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        # Initialize size analysis results
        self.size_analysis_results = None
        # Initialize intermediate files tracker for centralized cleanup
        self.intermediate_files = {
            'fastq': [],
            'igv': [],
            'temp': [],
            'alignment': [],
            'igv_dirs': []  # Track IGV directories separately
        }
        
    def run(self):
        """Run the complete ITD detection pipeline"""
        self.logger.info("=" * 60)
        self.logger.info("FLT3 ITD Detection Pipeline - Split Read Approach")
        self.logger.info("=" * 60)
        self.logger.info(f"Sample: {self.config.sample_name}")
        self.logger.info(f"Input BAM: {self.config.bam_file}")
        self.logger.info(f"Output directory: {self.config.output_dir}")
        
        try:
            # Step 1: Extract and trim reads
            self.logger.info("\n[Step 1/7] Extracting FLT3 reads and trimming primers...")
            reads, fastq_file = extract_flt3_reads(
                self.config.bam_file,
                self.config.genome_build,
                self.config.min_mapping_quality,
                self.config.min_read_length,
                config=self.config
            )
            
            # Track the trimmed FASTQ file for cleanup
            if fastq_file:
                self.intermediate_files['fastq'].append(fastq_file)
                # Also pass the FASTQ file to config for validator use
                self.config.fastq_file = fastq_file
            
            if not reads:
                self.logger.warning("No reads found in FLT3 region!")
                self._write_empty_results()
                # Cleanup any intermediate files that were created
                cleanup_intermediate_files(self.config, self.intermediate_files, self.logger)
                return
            
            self.logger.info(f"Extracted {len(reads)} reads")
            
            # Step 1.5: Analyze read size distribution
            self.logger.info("\n[Step 1.5/7] Analyzing read size distribution...")
            if fastq_file and Path(fastq_file).exists():
                try:
                    from read_size_analyzer import analyze_read_sizes
                    
                    # Create size analysis output directory
                    size_analysis_dir = Path(self.config.output_dir) / "size_analysis"
                    size_analysis_dir.mkdir(parents=True, exist_ok=True)
                    
                    # Track size analysis directory for potential cleanup
                    self.intermediate_files['temp'].append(str(size_analysis_dir))
                    
                    # Perform size analysis
                    size_results = analyze_read_sizes(
                        input_file=fastq_file,
                        output_dir=str(size_analysis_dir),
                        expected_wt_size=len(self.config.reference_sequence)  # Use actual reference length
                    )
                    
                    # Store results for HTML report
                    self.size_analysis_results = size_results
                    
                    if size_results and 'itd_signatures' in size_results:
                        itd_sigs = size_results['itd_signatures']
                        estimated_itd_freq = itd_sigs['estimated_itd_frequency']
                        potential_itd_reads = itd_sigs['potential_itd_reads']
                        
                        self.logger.info(f"Size analysis results:")
                        self.logger.info(f"  Estimated ITD frequency: {estimated_itd_freq:.1%}")
                        self.logger.info(f"  Potential ITD reads: {potential_itd_reads:,}")
                        self.logger.info(f"  Size analysis plots saved to: {size_analysis_dir}")
                        
                        # Log size distribution summary
                        if 'basic_stats' in size_results:
                            stats = size_results['basic_stats']
                            self.logger.info(f"  Mean read length: {stats['mean_length']:.1f} bp")
                            self.logger.info(f"  Length range: {stats['min_length']:.0f}-{stats['max_length']:.0f} bp")
                        
                        # If ITDs are detected in size analysis, note it
                        if estimated_itd_freq > 0.01:  # More than 1%
                            self.logger.info(f"⚠️  Size distribution suggests potential ITDs present!")
                        else:
                            self.logger.info("ℹ️  Size distribution appears normal (no obvious ITD signature)")
                    
                except ImportError:
                    self.logger.warning("Size analysis skipped - missing dependencies (scipy, seaborn)")
                    self.size_analysis_results = None
                except Exception as e:
                    self.logger.warning(f"Size analysis failed: {e}")
                    self.size_analysis_results = None
                except Exception as e:
                    self.logger.warning(f"Size analysis failed: {e}")
            else:
                self.logger.warning("Size analysis skipped - no trimmed FASTQ file available")
            
            # Step 2: Split reads into artificial pairs
            self.logger.info("\n[Step 2/7] Splitting reads into artificial pairs...")
            split_config = SplitConfig(
                split_method=self.config.split_method,
                overlap=self.config.split_overlap,
                min_split_size=self.config.min_split_size,
                num_splits=self.config.num_splits
            )
            
            pairs = split_reads(reads, 
                              method=self.config.split_method,
                              overlap=self.config.split_overlap,
                              min_split_size=self.config.min_split_size)
            
            if not pairs:
                self.logger.warning("No valid read pairs generated!")
                self._write_empty_results()
                # Cleanup any intermediate files that were created
                cleanup_intermediate_files(self.config, self.intermediate_files, self.logger)
                return
            
            self.logger.info(f"Generated {len(pairs)} read pairs")
            
            # Step 3: Align pairs and detect overlaps/soft clips
            self.logger.info("\n[Step 3/7] Aligning paired reads...")
            
            # Create temporary reference file
            ref_file = self.output_dir / "temp_reference.fa"
            with open(ref_file, 'w') as f:
                f.write(">FLT3_reference\n")
                f.write(self.config.reference_sequence + "\n")
            
            detection_results = detect_itds_from_pairs(
                pairs,
                str(ref_file),
                self.config.min_overlap_length,
                self.config.min_softclip_length
            )
            
            # Track temp files from paired aligner
            if 'temp_dir' in detection_results:
                # Add the entire temp directory for cleanup
                temp_dir = detection_results['temp_dir']
                if Path(temp_dir).exists():
                    self.intermediate_files['temp'].append(temp_dir)
            
            # Step 4: Generate ITD candidates
            self.logger.info("\n[Step 4/7] Generating ITD candidates...")
            candidates = generate_itd_candidates(
                detection_results,
                self.config.reference_sequence,
                self.config.min_itd_length,
                self.config.min_supporting_reads,
                self.config.max_itd_length
            )
            
            if not candidates:
                self.logger.info("No ITD candidates found")
                self._write_empty_results()
                if ref_file.exists():
                    self.intermediate_files['temp'].append(str(ref_file))
                cleanup_intermediate_files(self.config, self.intermediate_files, self.logger)
                return
            
            self.logger.info(f"Generated {len(candidates)} ITD candidates")
            
            # Apply length filters as hard filters to save computation time
            pre_filter_count = len(candidates)
            candidates = [c for c in candidates if self.config.min_itd_length <= c.length <= self.config.max_itd_length]
            
            if len(candidates) < pre_filter_count:
                filtered_out = pre_filter_count - len(candidates)
                self.logger.info(f"Filtered out {filtered_out} candidates outside length range "
                               f"({self.config.min_itd_length}-{self.config.max_itd_length} bp)")
            
            if not candidates:
                self.logger.info("No candidates remaining after length filtering")
                self._write_empty_results()
                if ref_file.exists():
                    self.intermediate_files['temp'].append(str(ref_file))
                cleanup_intermediate_files(self.config, self.intermediate_files, self.logger)
                return
            
            self.logger.info(f"Proceeding with {len(candidates)} length-filtered candidates")
            
            # Log candidate statistics
            position_bins = {}
            for cand in candidates:
                pos_bin = cand.position // 10
                if pos_bin not in position_bins:
                    position_bins[pos_bin] = []
                position_bins[pos_bin].append(cand)
            
            self.logger.info(f"Candidates distributed across {len(position_bins)} position bins")
            for pos_bin, cands in sorted(position_bins.items())[:5]:  # Show top 5
                self.logger.debug(f"  Position ~{pos_bin*10}: {len(cands)} candidates")
            
            # Step 5: Validate candidates with original reads
            self.logger.info("\n[Step 5/7] Validating ITD candidates...")
            
            # Pre-filter candidates to reduce validation time
            filtered_candidates = self._prefilter_candidates(candidates)
            self.logger.info(f"Pre-filtered to {len(filtered_candidates)} candidates for validation")
            
            #self.fix_candidate_read_names(filtered_candidates)
            # Note: IGV/BAM files are now handled within the validator based on config flags

            # Use config values for min_supporting_reads and min_high_confidence_ratio if available
            validate_kwargs = {
                'candidates': filtered_candidates,
                'reference_sequence': self.config.reference_sequence,
                'original_reads': reads,
                'max_itd_length': self.config.max_itd_length,
                'config': self.config
            }
            # Pass min_supporting_reads if present in config
            if hasattr(self.config, 'validation_min_reads'):
                validate_kwargs['min_supporting_reads'] = self.config.validation_min_reads
            if hasattr(self.config, 'min_high_confidence_ratio'):
                validate_kwargs['min_high_confidence_ratio'] = self.config.min_high_confidence_ratio

            # Set output dir for IGV validation files if bamout is set
            if getattr(self.config, 'bamout', False):
                # Pass a directory for IGV outputs to the validator config
                self.config.validation_igv_dir = str(self.output_dir / 'itd_validations')
            validation_results, igv_files, igv_dirs = validate_itd_candidates(**validate_kwargs)
            
            # Track IGV files and directories for cleanup
            if igv_files:
                self.intermediate_files['igv'].extend(igv_files)
            if igv_dirs:
                self.intermediate_files['igv_dirs'].extend(igv_dirs)
            
            if not validation_results:
                self.logger.info("No ITDs passed validation")
                self._write_empty_results()
                if ref_file.exists():
                    self.intermediate_files['temp'].append(str(ref_file))
                cleanup_intermediate_files(self.config, self.intermediate_files, self.logger)
                return
            
            self.logger.info(f"Validated {len(validation_results)} ITDs")
            
            # Step 6: Write results
            self.logger.info("\n[Step 6/7] Writing results...")
            
            # Write VCF
            if self.config.write_vcf:
                vcf_file = self.output_dir / f"{self.config.sample_name}_FLT3_ITDs.vcf"
                write_vcf_results(
                    validation_results,
                    self.config.reference_sequence,
                    str(vcf_file),
                    self.config.sample_name
                )
                self.logger.info(f"VCF written to: {vcf_file}")
            
            # Step 7: Generate HTML report
            if self.config.write_html:
                self.logger.info("\n[Step 7/7] Generating HTML report...")
                html_file = self.output_dir / f"{self.config.sample_name}_FLT3_ITD_report.html"
                
                # Check if size analysis results are available
                size_analysis_results = None
                size_analysis_dir = None
                if hasattr(self, 'size_analysis_results') and self.size_analysis_results:
                    size_analysis_results = self.size_analysis_results
                    size_analysis_dir = str(Path(self.config.output_dir) / "size_analysis")
                
                generate_html_report(
                    validation_results,
                    self.config.sample_name,
                    self.config.reference_sequence,
                    len(reads),
                    str(html_file),
                    size_analysis_results,
                    size_analysis_dir
                )
                self.logger.info(f"HTML report written to: {html_file}")
            
            # Write summary
            self._write_summary(validation_results, len(reads))
            
            # Cleanup - remove temp reference first, then centralized cleanup
            if ref_file.exists():
                self.intermediate_files['temp'].append(str(ref_file))
            
            # Centralized cleanup of all intermediate files
            cleanup_intermediate_files(self.config, self.intermediate_files, self.logger)
            
            self.logger.info("\n" + "=" * 60)
            self.logger.info("Pipeline completed successfully!")
            self.logger.info("=" * 60)
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}", exc_info=True)
            if self.config.debug:
                import traceback
                traceback.print_exc()
            # Cleanup intermediate files even on failure
            cleanup_intermediate_files(self.config, self.intermediate_files, self.logger)
            raise

    # def fix_candidate_read_names(self, candidates):
    #     """Fix read names in candidates to match original BAM"""
    #     for candidate in candidates:
    #         original_names = []
    #         for split_name in candidate.supporting_reads:
    #             original_name = split_name.split('_R')[0]  # Remove _R1_0, _R2_1 etc
    #             if original_name not in original_names:
    #                 original_names.append(original_name)
    #         candidate.supporting_reads = original_names

    def _write_empty_results(self):
        """Write empty results when no ITDs found"""
        # Write empty VCF
        if self.config.write_vcf:
            vcf_file = self.output_dir / f"{self.config.sample_name}_FLT3_ITDs.vcf"
            with open(vcf_file, 'w') as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
                f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.config.sample_name}\n")
        
        # Write summary
        summary_file = self.output_dir / f"{self.config.sample_name}_summary.txt"
        with open(summary_file, 'w') as f:
            f.write(f"FLT3 ITD Detection Summary\n")
            f.write(f"Sample: {self.config.sample_name}\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Result: No ITDs detected\n")
    
    def _prefilter_candidates(self, candidates: List[ITDCandidate]) -> List[ITDCandidate]:
        """Pre-filter candidates to reduce validation time with better prioritization"""
        # Group candidates by position and length
        position_groups = {}
        
        for cand in candidates:
            # Create bins for position and length
            pos_bin = cand.position // 5
            len_bin = cand.length // 5
            key = (pos_bin, len_bin)
            if key not in position_groups:
                position_groups[key] = []
            position_groups[key].append(cand)
        
        # Select best candidates from each group
        filtered = []
        for group in position_groups.values():
            # IMPROVED RANKING: Prioritize supporting reads first, then confidence
            # Supporting reads are the most reliable indicator of real ITDs
            group.sort(key=lambda x: (len(x.supporting_reads), x.confidence), reverse=True)
            
            # Take top candidates from each group
            max_per_group = min(3, self.config.max_candidates_per_position)
            filtered.extend(group[:max_per_group])
        
        # IMPROVED FINAL RANKING: Support-weighted confidence score
        def calculate_priority_score(candidate):
            """Calculate priority score for validation ranking"""
            support_score = len(candidate.supporting_reads)
            confidence_score = candidate.confidence
            length_bonus = min(0.2, candidate.length / self.config.max_itd_length)  # Bonus for longer ITDs
            
            # Weighted score: supporting reads are most important
            priority_score = (support_score * 3.0) + confidence_score + length_bonus
            return priority_score
        
        # Sort by priority score (support-weighted)
        filtered.sort(key=calculate_priority_score, reverse=True)
        
        # Limit total number of candidates
        max_candidates = 50  # Reasonable limit for validation
        if len(filtered) > max_candidates:
            self.logger.warning(f"Limiting validation to top {max_candidates} candidates")
            filtered = filtered[:max_candidates]
        
        # Debug: Show top candidates being selected for validation
        self.logger.info("Top candidates selected for validation:")
        for i, cand in enumerate(filtered[:5]):
            score = calculate_priority_score(cand)
            self.logger.info(f"  {i+1}. {cand.length}bp at pos {cand.position}, support={len(cand.supporting_reads)}, conf={cand.confidence:.2f}, score={score:.1f}")
        
        return filtered
    
    def _write_summary(self, results: List[ValidationResult], total_reads: int):
        """Write summary file"""
        summary_file = self.output_dir / f"{self.config.sample_name}_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write(f"FLT3 ITD Detection Summary\n")
            f.write(f"{'=' * 50}\n")
            f.write(f"Sample: {self.config.sample_name}\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total reads analyzed: {total_reads}\n")
            f.write(f"ITDs detected: {len(results)}\n\n")
            
            if results:
                f.write("ITD Details:\n")
                f.write("-" * 50 + "\n")
                for i, result in enumerate(results, 1):
                    f.write(f"\nITD {i}:\n")
                    f.write(f"  Length: {result.duplication_length} bp\n")
                    f.write(f"  Position: {result.insertion_position}\n")
                    f.write(f"  Allele Frequency: {result.allele_frequency:.1%}\n")
                    f.write(f"  Supporting Reads: {result.itd_coverage}\n")
                    f.write(f"  Confidence: {result.validation_confidence:.1%}\n")

def main():
    """Main entry point"""
    # Parse arguments
    config = parse_arguments()
    
    # Check dependencies
    required_tools = ['minimap2', 'samtools']
    for tool in required_tools:
        if shutil.which(tool) is None:
            print(f"Error: Required tool '{tool}' not found in PATH")
            sys.exit(1)
    
    # Check input file
    if not Path(config.bam_file).exists():
        print(f"Error: Input BAM file not found: {config.bam_file}")
        sys.exit(1)
    
    # Run pipeline
    pipeline = FLT3ITDPipeline(config)
    pipeline.run()

if __name__ == "__main__":
    main()