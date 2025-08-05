#!/usr/bin/env python3
"""
FLT3 ITD Detection Main Module
Clean CIGAR-based ITD detection pipeline - no legacy code
"""

import sys
import os
import logging
import shutil
import subprocess
from pathlib import Path
from datetime import datetime

# Core modules
from config import parse_arguments, FLT3Config
from bam_extractor import extract_flt3_reads
from flt3_itd_pipeline import run_flt3_itd_pipeline
from vcf_writer import write_vcf_results
from html_reporter import generate_html_report

def setup_logging(debug: bool = False):
    """Configure logging"""
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)

class FLT3ITDDetector:
    
    def __init__(self, config: FLT3Config):
        self.config = config
        self.logger = setup_logging(config.debug)
        
        # Create output directory with mkdir -p behavior
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Output directory: {self.output_dir}")
        
        # Create size analysis subdirectory
        self.size_analysis_dir = self.output_dir / "flt3_size_analysis"
        self.size_analysis_dir.mkdir(parents=True, exist_ok=True)
        
        # Track files for cleanup (only if not in debug mode)
        self.temp_files = []
        self.keep_temp_files = config.debug
        
    def run(self):
        """Run the complete ITD detection pipeline"""
        self.logger.info("=" * 60)
        self.logger.info("FLT3 ITD Detection - CIGAR-based Pipeline")
        self.logger.info("=" * 60)
        self.logger.info(f"Sample: {self.config.sample_name}")
        self.logger.info(f"Input BAM: {self.config.bam_file}")
        self.logger.info(f"Output directory: {self.config.output_dir}")
        self.logger.info(f"ITD size range: {self.config.min_itd_length}-{self.config.max_itd_length}bp")
        self.logger.info(f"Min support: {self.config.min_supporting_reads} reads")
        self.logger.info(f"Min allele frequency: {self.config.min_allele_frequency}")
        
        try:
            # Step 1: Extract and trim FLT3 reads
            self.logger.info("\n[Step 1/4] Extracting FLT3 reads and trimming primers...")
            reads, trimmed_fastq = self._extract_reads()
            
            if not reads:
                self.logger.warning("No reads found in FLT3 region!")
                self._write_empty_results()
                return
            
            # Store total reads for reporting
            self.total_reads = len(reads)
            self.logger.info(f"Extracted {len(reads)} reads")
            
            # Step 1.5: Analyze read size distribution
            self.logger.info("\\n[Step 1.5/4] Analyzing read size distribution...")
            self.size_analysis_results = self._analyze_read_sizes(trimmed_fastq)
            
            # Step 2: Create reference file for pipeline
            reference_file = self._create_reference_file()
            
            # Step 3: Run CIGAR-based ITD detection
            self.logger.info("\n[Step 3/5] Running CIGAR-based ITD detection...")
            # Use trimmed FASTQ for re-alignment and CIGAR detection
            result = self._run_itd_detection(trimmed_fastq, reference_file)
            
            if not result.validated_candidates:
                self.logger.warning("No ITDs detected!")
                self._write_empty_results()
                return
            
            self.logger.info(f"Found {len(result.validated_candidates)} validated ITDs")
            
            # Step 4: Generate outputs
            self.logger.info("\n[Step 4/5] Generating output files...")
            self._write_outputs(result)
            
            # Step 5: Summary (removed - not needed for default output)
            # self._write_summary(result)
            
            self.logger.info("\n" + "=" * 60)
            self.logger.info("FLT3 ITD Detection Complete!")
            self.logger.info("=" * 60)
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}")
            if self.config.debug:
                import traceback
                self.logger.error(traceback.format_exc())
            raise
        
        finally:
            # Cleanup temporary files
            self._cleanup()
    
    def _extract_reads(self):
        """Extract and trim FLT3 reads"""
        try:
            reads, fastq_file = extract_flt3_reads(
                self.config.bam_file,
                self.config.genome_build,
                self.config.min_mapping_quality,
                self.config.min_read_length,
                config=self.config
            )
            
            if fastq_file:
                self.temp_files.append(fastq_file)
            
            return reads, fastq_file
            
        except Exception as e:
            self.logger.error(f"Failed to extract reads: {e}")
            raise
    
    def _create_reference_file(self) -> str:
        """Create reference FASTA file for the pipeline"""
        reference_file = self.output_dir / f"{self.config.sample_name}_reference.fasta"
        
        with open(reference_file, 'w') as f:
            f.write(f">FLT3_reference_{self.config.genome_build}\n")
            f.write(f"{self.config.reference_sequence}\n")
        
        self.temp_files.append(str(reference_file))
        self.logger.debug(f"Created reference file: {reference_file}")
        return str(reference_file)
    
    def _analyze_read_sizes(self, fastq_file: str) -> dict:
        """Analyze read size distribution to detect ITD signatures"""
        try:
            from read_size_analyzer import analyze_read_sizes
            
            # Run read size analysis
            size_results = analyze_read_sizes(
                input_file=fastq_file,
                output_dir=str(self.size_analysis_dir),  # Use size analysis subdirectory
                expected_wt_size=self.config.amplicon_length,
                wt_tolerance=self.config.wt_tolerance,
                min_allele_frequency=self.config.min_allele_frequency,
                sample_name=self.config.sample_name
            )
            
            # Extract ITD sizes from the analysis results
            detected_itd_sizes = []

            # Debug: Log peak detection results
            if 'peaks' in size_results and size_results['peaks']['peaks']:
                self.logger.info(f"Read size analysis detected {len(size_results['peaks']['peaks'])} peaks:")
                for i, peak in enumerate(size_results['peaks']['peaks'][:10]):  # Show first 10 peaks
                    self.logger.info(f"  Peak {i+1}: {peak['size']:.1f}bp ({peak['height']} reads) - {peak['type']}")
            else:
                self.logger.warning("No peaks detected in read size analysis!")
            
            # Method 1: Extract from peaks with ITD type
            if 'peaks' in size_results and 'peaks' in size_results['peaks']:
                for peak in size_results['peaks']['peaks']:
                    peak_type = peak.get('type', '')
                    if 'potential_ITD_' in peak_type:
                        # Extract size from "potential_ITD_61bp" format
                        import re
                        match = re.search(r'potential_ITD_(\d+)bp', peak_type)
                        if match:
                            detected_itd_sizes.append(int(match.group(1)))
            
            # Method 2: Extract from ITD size statistics
            if 'itd_signatures' in size_results and size_results['itd_signatures'].get('itd_size_stats'):
                itd_stats = size_results['itd_signatures']['itd_size_stats']
                if itd_stats.get('count', 0) > 0:
                    # Add mean and median ITD sizes
                    if 'mean_size' in itd_stats:
                        detected_itd_sizes.append(int(round(itd_stats['mean_size'])))
                    if 'median_size' in itd_stats:
                        detected_itd_sizes.append(int(round(itd_stats['median_size'])))
            
            # Method 3: Parse from report text if available
            report_file = None
            for file_path in self.output_dir.glob("*size_analysis_report.txt"):
                report_file = file_path
                break
            
            if report_file and report_file.exists():
                try:
                    with open(report_file, 'r') as f:
                        report_text = f.read()
                        size_results['report_text'] = report_text
                        
                        # Extract ITD sizes from the report text
                        import re
                        # Look for patterns like "potential_ITD_61bp"
                        matches = re.findall(r'potential_ITD_(\d+)bp', report_text)
                        for match in matches:
                            detected_itd_sizes.append(int(match))
                        
                        # Also look for "Mean ITD size: X bp"
                        mean_match = re.search(r'Mean ITD size:\s*(\d+(?:\.\d+)?)\s*bp', report_text)
                        if mean_match:
                            detected_itd_sizes.append(int(float(mean_match.group(1))))
                        
                        # And "Median ITD size: X bp"
                        median_match = re.search(r'Median ITD size:\s*(\d+(?:\.\d+)?)\s*bp', report_text)
                        if median_match:
                            detected_itd_sizes.append(int(float(median_match.group(1))))
                        
                except Exception as e:
                    self.logger.warning(f"Failed to parse size analysis report: {e}")
            
            # Store detected ITD sizes (remove duplicates and sort)
            if detected_itd_sizes:
                size_results['detected_itd_sizes'] = sorted(list(set(detected_itd_sizes)))
                self.logger.info(f"Size analysis detected ITDs: {size_results['detected_itd_sizes']}")
            else:
                size_results['detected_itd_sizes'] = []
                self.logger.info("No ITD sizes detected in size analysis")
            
            self.logger.info("Read size analysis completed")
            
            return size_results
            
        except Exception as e:
            self.logger.warning(f"Read size analysis failed: {e}")
            return {}
    
    def _realign_reads(self, fastq_file: str, reference_file: str) -> str:
        """Re-align trimmed FASTQ reads to reference using minimap2 and samtools"""
        try:
            # Output SAM and BAM file paths
            sam_file = f"{self.output_dir}/{self.config.sample_name}_realigned.sam"
            bam_file = f"{self.output_dir}/{self.config.sample_name}_realigned.bam"
            sorted_bam = f"{self.output_dir}/{self.config.sample_name}_realigned_sorted.bam"
            
            # Step 1: Align with minimap2
            self.logger.info(f"Aligning {fastq_file} to {reference_file} using minimap2...")
            minimap_cmd = f"minimap2 -ax sr -t {self.config.threads} {reference_file} {fastq_file} > {sam_file}"
            result = subprocess.run(minimap_cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise RuntimeError(f"minimap2 alignment failed: {result.stderr}")
            
            # Step 2: Convert SAM to BAM
            self.logger.info("Converting SAM to BAM...")
            samtools_view_cmd = f"samtools view -bS {sam_file} > {bam_file}"
            result = subprocess.run(samtools_view_cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise RuntimeError(f"samtools view failed: {result.stderr}")
            
            # Step 3: Sort BAM file
            self.logger.info("Sorting BAM file...")
            samtools_sort_cmd = f"samtools sort -@ {self.config.threads} {bam_file} -o {sorted_bam}"
            result = subprocess.run(samtools_sort_cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise RuntimeError(f"samtools sort failed: {result.stderr}")
            
            # Step 4: Index BAM file
            self.logger.info("Indexing BAM file...")
            samtools_index_cmd = f"samtools index {sorted_bam}"
            result = subprocess.run(samtools_index_cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise RuntimeError(f"samtools index failed: {result.stderr}")
            
            # Clean up intermediate files
            os.remove(sam_file)
            os.remove(bam_file)
            
            # Add the final sorted BAM and its index to temp files for cleanup
            self.temp_files.append(sorted_bam)
            self.temp_files.append(f"{sorted_bam}.bai")
            
            self.logger.info(f"Re-alignment completed: {sorted_bam}")
            return sorted_bam
            
        except Exception as e:
            self.logger.error(f"Re-alignment failed: {e}")
            raise

    def _run_itd_detection(self, trimmed_fastq: str, reference_file: str):
        """Run the CIGAR-based ITD detection pipeline"""
        try:
            # Step 1: Re-align trimmed reads to FLT3 reference to get CIGAR strings
            self.logger.info("Re-aligning trimmed reads to FLT3 reference...")
            aligned_bam = self._realign_reads(trimmed_fastq, reference_file)
            
            # Step 2: Run CIGAR-based ITD detection on aligned BAM
            result = run_flt3_itd_pipeline(
                bam_file=aligned_bam,
                reference_sequence=self.config.reference_sequence,
                reference_file=reference_file,
                trimmed_fastq=trimmed_fastq,
                min_itd_length=self.config.min_itd_length,
                max_itd_length=self.config.max_itd_length,
                min_support=self.config.min_supporting_reads,
                min_softclip_length=self.config.min_softclip_length,
                enable_softclip_fallback=self.config.enable_softclip_fallback,
                size_analysis_results=self.size_analysis_results,
                output_prefix=f"{self.output_dir}/{self.config.sample_name}",
                config=self.config
            )
            
            return result
            
        except Exception as e:
            self.logger.error(f"ITD detection failed: {e}")
            raise
    
    def _write_outputs(self, result):
        """Write default outputs: VCF, HTML, and size analysis files"""
        
        # Convert candidates to legacy format for compatibility
        legacy_results = self._convert_to_legacy_format(result)
        
        # Always write VCF (default output)
        vcf_file = self.output_dir / f"{self.config.sample_name}_ITDs.vcf"
        try:
            write_vcf_results(
                legacy_results, 
                self.config.reference_sequence, 
                str(vcf_file), 
                self.config.sample_name,
                genome_chr=self.config.vcf_chr,
                genome_start=self.config.vcf_start_pos
            )
            self.logger.info(f"VCF written to {vcf_file}")
        except Exception as e:
            self.logger.error(f"Failed to write VCF: {e}")
        
        # Always write HTML report (default output, even if negative)
        html_file = self.output_dir / f"{self.config.sample_name}_ITD_report.html"
        try:
            generate_html_report(
                validation_results=legacy_results,
                sample_name=self.config.sample_name,
                reference_sequence=self.config.reference_sequence,
                total_reads=getattr(self, 'total_reads', 0),
                output_file=str(html_file),
                size_analysis_results=self.size_analysis_results,
                size_analysis_dir=str(self.size_analysis_dir)  # Point to size analysis subdirectory
            )
            self.logger.info(f"HTML report written to {html_file}")
        except Exception as e:
            self.logger.error(f"Failed to write HTML report: {e}")
        
        # Size analysis files are already generated in flt3_size_analysis directory
        self.logger.info(f"Size analysis files written to {self.size_analysis_dir}")
        
        # Export BAM/FASTA for IGV if --bamout is specified
        if self.config.bamout:
            self._export_igv_files(result)
    
    def _export_igv_files(self, result):
        """Export BAM and FASTA files for IGV viewing when --bamout is specified"""
        try:
            # Check if bamout files were created during validation
            bamout_dir = self.output_dir / "bamout"
            if bamout_dir.exists():
                self.logger.info(f"BAM validation files created in: {bamout_dir}")
                # Files are kept in bamout directory for --bamout mode
                # Don't move them to main output directory
            else:
                self.logger.warning("BAM validation files not found. This may be due to:")
                self.logger.warning("  - Missing samtools (BAM indexing failed)")
                self.logger.warning("  - BAM file creation issues")
                self.logger.warning("  - Validation process failure")
                self.logger.info("Continuing with standard output files (VCF, HTML)")
        except Exception as e:
            self.logger.error(f"Failed to export IGV files: {e}")
    
    def _convert_to_legacy_format(self, result):
        """Convert new pipeline results to legacy format for VCF/HTML writers"""
        legacy_results = []
        
        for i, candidate in enumerate(result.validated_candidates):
            validation = result.validation_results[i] if i < len(result.validation_results) else None
            
            # Create legacy candidate object
            class LegacyCandidate:
                def __init__(self, cand):
                    self.sequence = cand.sequence
                    self.length = cand.length
                    self.position = cand.position
                    self.supporting_reads = cand.supporting_reads
                    self.confidence = cand.confidence
                    self.support_type = getattr(cand, 'support_type', 'cigar_insertion')
                    self.duplication_start = getattr(cand, 'duplication_start', None)
                    self.duplication_end = getattr(cand, 'duplication_end', None)
                    self.is_primary = getattr(cand, 'is_primary', False)
            
            # Create legacy validation result
            class LegacyValidation:
                def __init__(self, cand, val):
                    self.itd_candidate = LegacyCandidate(cand)
                    
                    # Handle supporting_reads being either int or list
                    support_count = cand.supporting_reads if isinstance(cand.supporting_reads, int) else len(cand.supporting_reads)
                    
                    # Handle both object and dict validation results
                    if isinstance(val, dict):
                        self.is_valid = val.get('is_valid', val.get('passed', True))
                        self.validation_confidence = val.get('validation_confidence', cand.confidence)
                        self.allele_frequency = val.get('allele_frequency', 0.1)
                        self.total_coverage = val.get('total_coverage', val.get('total_reads_processed', support_count * 2))
                        self.itd_coverage = val.get('itd_coverage', val.get('supporting_reads', support_count))
                    else:
                        # Original object-based validation results
                        self.is_valid = val.is_valid if val else True
                        self.validation_confidence = val.validation_confidence if val else cand.confidence
                        self.allele_frequency = val.allele_frequency if val else 0.1
                        self.total_coverage = val.total_coverage if val else support_count * 2
                        self.itd_coverage = val.itd_coverage if val else support_count
                    
                    self.insertion_position = cand.position
                    self.duplication_length = cand.length
                    self.duplication_start = getattr(cand, 'duplication_start', None)
                    self.duplication_end = getattr(cand, 'duplication_end', None)
            
            legacy_results.append(LegacyValidation(candidate, validation))
        
        return legacy_results
    
    def _create_itd_reference(self, candidate):
        """Create ITD reference sequence"""
        ref_seq = self.config.reference_sequence
        pos = candidate.position
        
        # Insert ITD sequence at position
        if pos < len(ref_seq):
            itd_ref = ref_seq[:pos] + candidate.sequence + ref_seq[pos:]
        else:
            itd_ref = ref_seq + candidate.sequence
        
        return itd_ref
    
    def _write_summary(self, result):
        """Write summary file"""
        summary_file = self.output_dir / f"{self.config.sample_name}_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write(f"FLT3 ITD Detection Summary - {self.config.sample_name}\\n")
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
            f.write(f"Pipeline: CIGAR-based ITD detection\\n")
            f.write(f"Genome Build: {self.config.genome_build}\\n")
            f.write(f"ITD Size Range: {self.config.min_itd_length}-{self.config.max_itd_length}bp\\n")
            f.write(f"Min Support: {self.config.min_supporting_reads} reads\\n")
            f.write(f"Min AF: {self.config.min_allele_frequency}\\n\\n")
            
            # Detection results
            f.write(f"CIGAR candidates: {len(result.cigar_candidates)}\\n")
            f.write(f"Soft-clip candidates: {len(result.softclip_candidates)}\\n")
            f.write(f"Validated ITDs: {len(result.validated_candidates)}\\n")
            
            # Runtime stats
            runtime = result.runtime_stats
            f.write(f"Detection time: {runtime.get('total_time', 0):.2f}s\\n")
            if 'cigar_detection_time' in runtime:
                f.write(f"  CIGAR detection: {runtime['cigar_detection_time']:.2f}s\\n")
            if 'softclip_detection_time' in runtime:
                f.write(f"  Soft-clip detection: {runtime['softclip_detection_time']:.2f}s\\n")
            if 'validation_time' in runtime:
                f.write(f"  Validation: {runtime['validation_time']:.2f}s\\n")
            
            f.write("\\n")
            
            # ITD details
            if result.validated_candidates:
                f.write("ITD Details:\\n")
                for i, candidate in enumerate(result.validated_candidates, 1):
                    validation = result.validation_results[i-1] if i-1 < len(result.validation_results) else None
                    
                    # Helper function to get support count
                    def get_support_count(supporting_reads):
                        if isinstance(supporting_reads, list):
                            return len(supporting_reads)
                        elif isinstance(supporting_reads, int):
                            return supporting_reads
                        else:
                            return 0
                    
                    support_count = get_support_count(candidate.supporting_reads)
                    
                    f.write(f"{i}. Length: {candidate.length}bp, ")
                    f.write(f"Position: {candidate.position}, ")
                    f.write(f"Method: {candidate.support_type}, ")
                    f.write(f"Support: {support_count} reads, ")
                    f.write(f"Confidence: {candidate.confidence:.3f}")
                    
                    if validation:
                        f.write(f", AF: {validation.allele_frequency:.3f}")
                    
                    # Add primary marking if available
                    if hasattr(candidate, 'is_primary') and candidate.is_primary:
                        f.write(f", PRIMARY")
                    
                    f.write("\\n")
                    f.write(f"   Sequence: {candidate.sequence}\\n")
                    
                    if hasattr(candidate, 'duplication_start') and candidate.duplication_start:
                        f.write(f"   Duplication region: {candidate.duplication_start}-{candidate.duplication_end}\\n")
                    
                    f.write("\\n")
            else:
                f.write("No ITDs detected.\\n")
        
        self.logger.info(f"Summary written to {summary_file}")
    
    def _write_empty_results(self):
        """Write default outputs when no ITDs found"""
        # Always write VCF (default output, even if empty)
        vcf_file = self.output_dir / f"{self.config.sample_name}_ITDs.vcf"
        try:
            write_vcf_results([], self.config.reference_sequence, str(vcf_file),
                             self.config.sample_name, genome_chr=self.config.vcf_chr,
                             genome_start=self.config.vcf_start_pos)
            self.logger.info(f"Empty VCF written to {vcf_file}")
        except Exception as e:
            self.logger.error(f"Failed to write empty VCF: {e}")

        # Always write HTML report (default output, even if negative)
        html_file = self.output_dir / f"{self.config.sample_name}_ITD_report.html"
        try:
            generate_html_report(
                validation_results=[],
                sample_name=self.config.sample_name,
                reference_sequence=self.config.reference_sequence,
                total_reads=getattr(self, 'total_reads', 0),
                output_file=str(html_file),
                size_analysis_results=self.size_analysis_results,
                size_analysis_dir=str(self.size_analysis_dir)
            )
            self.logger.info(f"HTML report written to {html_file}")
        except Exception as e:
            self.logger.error(f"Failed to write HTML report: {e}")

        # Summary file removed - not needed for default output
    
    def _cleanup(self):
        """Clean up temporary files"""
        if not self.config.debug:  # Keep temp files in debug mode
            for temp_file in self.temp_files:
                try:
                    if Path(temp_file).exists():
                        # If it's a file in a temp directory, also try to clean up the directory
                        temp_path = Path(temp_file)
                        temp_path.unlink()
                        self.logger.debug(f"Cleaned up: {temp_file}")
                        
                        # Try to clean up the parent directory if it's a temp directory
                        parent_dir = temp_path.parent
                        if parent_dir.name.startswith(('flt3_cutadapt_', 'tmp')) and parent_dir.exists():
                            try:
                                import shutil
                                shutil.rmtree(parent_dir)
                                self.logger.debug(f"Cleaned up temp directory: {parent_dir}")
                            except Exception as e:
                                self.logger.warning(f"Failed to cleanup temp directory {parent_dir}: {e}")
                except Exception as e:
                    self.logger.warning(f"Failed to cleanup {temp_file}: {e}")

def check_dependencies():
    """Check required dependencies"""
    missing = []
    
    # Check external tools
    for tool in ['samtools', 'minimap2']:
        if not shutil.which(tool):
            missing.append(f"{tool} (not in PATH)")
    
    # Check Python packages
    try:
        import pysam
    except ImportError:
        missing.append("pysam (pip install pysam)")
    
    try:
        import numpy
    except ImportError:
        missing.append("numpy (pip install numpy)")
    
    return missing

def main():
    """Main entry point"""
    # Parse arguments
    config = parse_arguments()
    
    # Check dependencies
    missing_deps = check_dependencies()
    if missing_deps:
        print("Error: Missing dependencies:")
        for dep in missing_deps:
            print(f"  - {dep}")
        print("\\nPlease install missing dependencies and try again.")
        sys.exit(1)
    
    # Check input file
    if not Path(config.bam_file).exists():
        print(f"Error: Input BAM file not found: {config.bam_file}")
        sys.exit(1)
    
    # Run pipeline
    detector = FLT3ITDDetector(config)
    try:
        detector.run()
        print(f"\\n✅ Analysis complete! Results saved to: {config.output_dir}")
    except Exception as e:
        print(f"\\n❌ Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
