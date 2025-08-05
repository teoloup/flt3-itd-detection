#!/usr/bin/env python3
"""
Dual-Reference Validator Module
Validates ITDs using dual-reference alignment approach
"""

import logging
import tempfile
import subprocess
import pysam
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import shutil

logger = logging.getLogger(__name__)

@dataclass
class ValidationResult:
    """Results from dual-reference ITD validation"""
    itd_candidate: 'ITDCandidate'
    supporting_reads: List[str]
    allele_frequency: float
    total_coverage: int
    itd_coverage: int
    validation_confidence: float
    insertion_position: int
    duplication_length: int
    wt_coverage: int  # Coverage on wild-type reference
    segregation_quality: float  # How well reads segregate between references
    is_valid: bool = True  # Whether this validation result is considered valid

class DualReferenceValidator:
    """Validate ITDs using dual-reference alignment"""
    
    def __init__(self, reference_sequence: str, original_reads_fastq: str = None, config=None):
        """
        Initialize dual-reference validator
        
        Args:
            reference_sequence: Original FLT3 reference sequence
            original_reads_fastq: Path to FASTQ file with original reads
            config: Configuration object
        """
        self.reference_sequence = reference_sequence
        self.original_reads_fastq = original_reads_fastq
        self.config = config
        self.temp_dir = tempfile.mkdtemp(prefix="flt3_dual_ref_")
        self.validation_files = []  # Track files for cleanup
        
        logger.info(f"Initialized dual-reference validator with temp dir: {self.temp_dir}")
    
    def validate_candidate(self, itd_candidate: 'ITDCandidate', 
                          min_supporting_reads: int = 3,
                          min_allele_frequency: float = 0.01) -> Optional[ValidationResult]:
        """
        Validate single ITD candidate using dual-reference approach
        
        Args:
            itd_candidate: ITD candidate to validate
            min_supporting_reads: Minimum supporting reads required
            min_allele_frequency: Minimum allele frequency required
            
        Returns:
            ValidationResult if validated, None otherwise
        """
        actual_length = len(itd_candidate.sequence)
        if actual_length != itd_candidate.length:
            logger.info(f"Validating ITD candidate: {itd_candidate.length}bp sequence (actual_length={actual_length}) at position {itd_candidate.position}")
        else:
            logger.info(f"Validating ITD candidate: {itd_candidate.length}bp sequence at position {itd_candidate.position}")
        
        try:
            # Step 1: Create dual-reference FASTA
            dual_ref_file = self._create_dual_reference(itd_candidate)
            
            # Step 2: Align reads to dual-reference
            alignment_bam = self._align_to_dual_reference(dual_ref_file)
            
            if not alignment_bam:
                logger.warning("Alignment to dual-reference failed")
                return None
            
            # Step 3: Analyze read segregation
            segregation_results = self._analyze_read_segregation(alignment_bam, itd_candidate)
            
            # Step 4: Calculate validation metrics
            validation_result = self._calculate_validation_metrics(
                itd_candidate, segregation_results, min_supporting_reads, min_allele_frequency
            )
            
            # Step 5: Export for IGV if requested
            if validation_result and self._should_export_igv():
                self._export_validation_files(itd_candidate, dual_ref_file, alignment_bam)
            
            return validation_result
            
        except Exception as e:
            logger.error(f"Validation failed for {itd_candidate.length}bp ITD: {e}")
            return None
    
    def _create_dual_reference(self, itd_candidate: 'ITDCandidate') -> str:
        """
        Create dual-reference FASTA with WT and ITD sequences
        
        Args:
            itd_candidate: ITD candidate
            
        Returns:
            Path to dual-reference FASTA file
        """
        # Determine insertion position
        insertion_pos = itd_candidate.insertion_site or itd_candidate.position
        
        # Create ITD-containing reference
        itd_reference = (
            self.reference_sequence[:insertion_pos] +
            itd_candidate.sequence +
            self.reference_sequence[insertion_pos:]
        )
        
        # Write dual-reference FASTA
        dual_ref_file = Path(self.temp_dir) / f"dual_ref_{itd_candidate.length}bp.fa"
        
        with open(dual_ref_file, 'w') as f:
            # Wild-type reference
            f.write(">FLT3_WT\n")
            f.write(self.reference_sequence + "\n")
            
            # ITD-containing reference
            f.write(f">FLT3_ITD_{itd_candidate.length}bp_pos{insertion_pos}\n")
            f.write(itd_reference + "\n")
        
        self.validation_files.append(str(dual_ref_file))
        
        logger.debug(f"Created dual-reference: WT({len(self.reference_sequence)}bp) + "
                    f"ITD({len(itd_reference)}bp) at position {insertion_pos}")
        
        return str(dual_ref_file)
    
    def _align_to_dual_reference(self, dual_ref_file: str) -> Optional[str]:
        """
        Align original reads to dual-reference
        
        Args:
            dual_ref_file: Path to dual-reference FASTA
            
        Returns:
            Path to alignment BAM file or None if failed
        """
        if not self.original_reads_fastq:
            logger.error("No FASTQ file provided for alignment")
            return None
        
        output_sam = Path(self.temp_dir) / "dual_ref_alignment.sam"
        
        # Run minimap2 with no secondary alignments
        threads = self.config.threads if self.config else 4
        cmd = [
            'minimap2',
            '-ax', 'map-ont',  # Nanopore preset
            '-t', str(threads),
            '--secondary=no',  # CRITICAL: No secondary alignments for clean segregation
            '--max-intron-len', '0',  # No introns expected
            dual_ref_file,
            self.original_reads_fastq
        ]
        
        logger.debug(f"Running minimap2: {' '.join(cmd)}")
        
        try:
            with open(output_sam, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
                
            if result.returncode != 0:
                logger.error(f"Minimap2 failed: {result.stderr}")
                return None
            
            # Convert SAM to BAM using samtools subprocess calls only
            output_bam = Path(self.temp_dir) / "dual_ref_alignment.bam"
            sorted_bam = Path(self.temp_dir) / "dual_ref_alignment_sorted.bam"
            
            try:
                # Convert SAM to BAM
                sam_to_bam_cmd = ['samtools', 'view', '-bS', str(output_sam)]
                with open(output_bam, 'wb') as bam_file:
                    result = subprocess.run(sam_to_bam_cmd, stdout=bam_file, stderr=subprocess.PIPE, text=False)
                    if result.returncode != 0:
                        logger.error(f"SAM to BAM conversion failed: {result.stderr.decode()}")
                        return None
                
                # Sort BAM
                sort_cmd = ['samtools', 'sort', '-o', str(sorted_bam), str(output_bam)]
                result = subprocess.run(sort_cmd, stderr=subprocess.PIPE, text=True)
                if result.returncode != 0:
                    logger.error(f"BAM sorting failed: {result.stderr}")
                    return None
                
                # Index sorted BAM
                index_cmd = ['samtools', 'index', str(sorted_bam)]
                result = subprocess.run(index_cmd, stderr=subprocess.PIPE, text=True)
                if result.returncode != 0:
                    logger.error(f"BAM indexing failed: {result.stderr}")
                    return None
                
                # Verify the BAM file was created correctly
                if not sorted_bam.exists() or sorted_bam.stat().st_size == 0:
                    logger.error("BAM file was not created or is empty")
                    return None
                
                index_file = Path(str(sorted_bam) + ".bai")
                if not index_file.exists():
                    logger.error("BAM index was not created")
                    return None
                
                # Clean up intermediate files
                if output_sam.exists():
                    output_sam.unlink()
                if output_bam.exists():
                    output_bam.unlink()
                
                self.validation_files.extend([str(sorted_bam), str(index_file)])
                
                logger.debug(f"Alignment completed successfully: {sorted_bam}")
                logger.debug(f"BAM size: {sorted_bam.stat().st_size} bytes")
                logger.debug(f"Index size: {index_file.stat().st_size} bytes")
                
                return str(sorted_bam)
                
            except FileNotFoundError:
                logger.error("samtools not found in PATH. Please install samtools.")
                return None
            except Exception as e:
                logger.error(f"BAM creation failed: {e}")
                return None
            
        except Exception as e:
            logger.error(f"Alignment failed: {e}")
            return None
    
    def _analyze_read_segregation(self, alignment_bam: str, itd_candidate: 'ITDCandidate') -> Dict:
        """
        Analyze how reads segregate between WT and ITD references
        
        Args:
            alignment_bam: Path to dual-reference alignment BAM
            itd_candidate: ITD candidate
            
        Returns:
            Dictionary with segregation analysis results
        """
        wt_reads = []
        itd_reads = []
        total_reads = 0
        unmapped_reads = 0
        
        try:
            # First verify BAM file exists and is readable
            bam_path = Path(alignment_bam)
            if not bam_path.exists():
                logger.error(f"BAM file not found: {alignment_bam}")
                return self._get_empty_segregation_result()
            
            if bam_path.stat().st_size < 100:
                logger.error(f"BAM file too small: {bam_path.stat().st_size} bytes")
                return self._get_empty_segregation_result()
            
            # Try to open and read the BAM file
            with pysam.AlignmentFile(alignment_bam, "rb") as bam:
                for read in bam.fetch():
                    total_reads += 1
                    
                    if read.is_unmapped:
                        unmapped_reads += 1
                        continue
                    
                    # Determine which reference the read aligned to
                    ref_name = bam.get_reference_name(read.reference_id)
                    
                    if ref_name == "FLT3_WT":
                        wt_reads.append({
                            'name': read.query_name,
                            'mapping_quality': read.mapping_quality,
                            'start': read.reference_start,
                            'end': read.reference_end,
                            'cigar': read.cigarstring
                        })
                    elif ref_name.startswith("FLT3_ITD"):
                        itd_reads.append({
                            'name': read.query_name,
                            'mapping_quality': read.mapping_quality,
                            'start': read.reference_start,
                            'end': read.reference_end,
                            'cigar': read.cigarstring
                        })
        
        except Exception as e:
            logger.error(f"Error analyzing read segregation: {e}")
            logger.error(f"This may be due to missing pysam or corrupted BAM file: {alignment_bam}")
            return self._get_empty_segregation_result()
        
        # Calculate segregation metrics
        mapped_reads = len(wt_reads) + len(itd_reads)
        wt_coverage = len(wt_reads)
        itd_coverage = len(itd_reads)
        
        # Calculate segregation quality (how cleanly reads separate)
        if mapped_reads > 0:
            segregation_quality = max(wt_coverage, itd_coverage) / mapped_reads
        else:
            segregation_quality = 0.0
        
        # Calculate allele frequency
        allele_frequency = itd_coverage / mapped_reads if mapped_reads > 0 else 0.0
        
        results = {
            'total_reads': total_reads,
            'mapped_reads': mapped_reads,
            'unmapped_reads': unmapped_reads,
            'wt_reads': wt_reads,
            'itd_reads': itd_reads,
            'wt_coverage': wt_coverage,
            'itd_coverage': itd_coverage,
            'allele_frequency': allele_frequency,
            'segregation_quality': segregation_quality
        }
        
        logger.debug(f"Read segregation: {wt_coverage} WT, {itd_coverage} ITD, "
                    f"AF={allele_frequency:.3f}, segregation_quality={segregation_quality:.3f}")
        
        return results
    
    def _get_empty_segregation_result(self) -> Dict:
        """Return empty segregation result for error cases"""
        return {
            'total_reads': 0,
            'mapped_reads': 0,
            'unmapped_reads': 0,
            'wt_reads': [],
            'itd_reads': [],
            'wt_coverage': 0,
            'itd_coverage': 0,
            'allele_frequency': 0.0,
            'segregation_quality': 0.0
        }
    
    def _calculate_validation_metrics(self, itd_candidate: 'ITDCandidate', 
                                    segregation_results: Dict,
                                    min_supporting_reads: int,
                                    min_allele_frequency: float) -> Optional[ValidationResult]:
        """
        Calculate final validation metrics and determine if ITD is valid
        
        Args:
            itd_candidate: ITD candidate
            segregation_results: Results from read segregation analysis
            min_supporting_reads: Minimum supporting reads required
            min_allele_frequency: Minimum allele frequency required
            
        Returns:
            ValidationResult if ITD passes validation, None otherwise
        """
        if not segregation_results:
            return None
        
        # Extract metrics
        itd_coverage = segregation_results['itd_coverage']
        total_coverage = segregation_results['mapped_reads']
        allele_frequency = segregation_results['allele_frequency']
        segregation_quality = segregation_results['segregation_quality']
        wt_coverage = segregation_results['wt_coverage']
        
        # Check validation criteria
        passes_support = itd_coverage >= min_supporting_reads
        passes_af = allele_frequency >= min_allele_frequency
        passes_segregation = segregation_quality >= 0.6  # At least 60% of reads in one group
        
        if not (passes_support and passes_af and passes_segregation):
            logger.debug(f"ITD validation failed: support={passes_support} ({itd_coverage}>={min_supporting_reads}), "
                        f"AF={passes_af} ({allele_frequency:.3f}>={min_allele_frequency}), "
                        f"segregation={passes_segregation} ({segregation_quality:.3f}>=0.6)")
            return None
        
        # Calculate validation confidence
        validation_confidence = self._calculate_validation_confidence(
            itd_candidate, segregation_results
        )
        
        # Extract supporting read names
        supporting_reads = [read['name'] for read in segregation_results['itd_reads']]
        
        result = ValidationResult(
            itd_candidate=itd_candidate,
            supporting_reads=supporting_reads,
            allele_frequency=allele_frequency,
            total_coverage=total_coverage,
            itd_coverage=itd_coverage,
            validation_confidence=validation_confidence,
            insertion_position=itd_candidate.insertion_site or itd_candidate.position,
            duplication_length=itd_candidate.length,
            wt_coverage=wt_coverage,
            segregation_quality=segregation_quality,
            is_valid=True
        )
        
        actual_length = len(itd_candidate.sequence)
        if actual_length != itd_candidate.length:
            logger.info(f"ITD validation successful: {itd_candidate.length}bp sequence (actual_length={actual_length}), "
                        f"AF={allele_frequency:.3f}, support={itd_coverage}, confidence={validation_confidence:.3f}")
        else:
            logger.info(f"ITD validation successful: {itd_candidate.length}bp sequence, "
                        f"AF={allele_frequency:.3f}, support={itd_coverage}, confidence={validation_confidence:.3f}")
        
        return result
    
    def _calculate_validation_confidence(self, itd_candidate: 'ITDCandidate', 
                                       segregation_results: Dict) -> float:
        """
        Calculate validation confidence based on multiple factors
        
        Args:
            itd_candidate: ITD candidate
            segregation_results: Results from segregation analysis
            
        Returns:
            Validation confidence score (0.0 to 1.0)
        """
        # Base confidence from original detection
        base_confidence = itd_candidate.confidence * 0.4
        
        # Confidence from allele frequency
        af = segregation_results['allele_frequency']
        af_confidence = min(0.25, af * 0.25)  # Up to 0.25 for high AF
        
        # Confidence from read segregation quality
        segregation_quality = segregation_results['segregation_quality']
        segregation_confidence = segregation_quality * 0.2  # Up to 0.2 for perfect segregation
        
        # Confidence from supporting read count
        support_count = segregation_results['itd_coverage']
        support_confidence = min(0.15, support_count * 0.01)  # Up to 0.15 for high support
        
        total_confidence = base_confidence + af_confidence + segregation_confidence + support_confidence
        
        return min(0.95, total_confidence)
    
    def _should_export_igv(self) -> bool:
        """Check if IGV files should be exported"""
        if not self.config:
            return False
        
        return (getattr(self.config, 'bamout', False) or 
                getattr(self.config, 'keep_temp', False) or 
                getattr(self.config, 'debug', False))
    
    def _export_validation_files(self, itd_candidate: 'ITDCandidate', 
                                dual_ref_file: str, alignment_bam: str):
        """
        Export validation files for IGV visualization
        
        Args:
            itd_candidate: ITD candidate
            dual_ref_file: Path to dual-reference FASTA
            alignment_bam: Path to alignment BAM
        """
        try:
            # Determine output directory
            output_dir = getattr(self.config, 'output_dir', '.')
            igv_dir = Path(output_dir) / "validation_igv"
            igv_dir.mkdir(parents=True, exist_ok=True)
            
            # Create unique prefix
            prefix = f"itd_{itd_candidate.length}bp_pos{itd_candidate.position}_dual_ref"
            
            # Copy reference file
            ref_dest = igv_dir / f"{prefix}.fa"
            shutil.copy2(dual_ref_file, ref_dest)
            
            # Copy and rename BAM file (with integrity check)
            bam_dest = igv_dir / f"{prefix}.bam"
            bai_dest = igv_dir / f"{prefix}.bam.bai"
            
            # Verify BAM file integrity before copying using samtools
            try:
                # Use samtools quickcheck to verify file integrity
                check_cmd = ['samtools', 'quickcheck', alignment_bam]
                result = subprocess.run(check_cmd, stderr=subprocess.PIPE, text=True, timeout=30)
                
                if result.returncode != 0:
                    logger.error(f"BAM file integrity check failed: {result.stderr}")
                    logger.error(f"Skipping export of corrupted BAM: {alignment_bam}")
                    return
                
                # Additional check: verify file size is reasonable
                bam_size = Path(alignment_bam).stat().st_size
                if bam_size < 100:  # Too small to be a valid BAM
                    logger.error(f"BAM file too small ({bam_size} bytes), likely corrupted")
                    return
                
                logger.debug(f"BAM file integrity verified: {bam_size} bytes")
                    
            except subprocess.TimeoutExpired:
                logger.error("BAM integrity check timed out - file may be corrupted")
                return
            except FileNotFoundError:
                logger.warning("samtools not available for integrity check, proceeding with copy...")
            except Exception as e:
                logger.warning(f"BAM integrity check failed: {e}, proceeding with copy...")
            
            # BAM file appears valid, proceed with copy
            try:
                shutil.copy2(alignment_bam, bam_dest)
                
                # Copy index if it exists
                index_src = Path(alignment_bam + ".bai")
                if index_src.exists():
                    shutil.copy2(index_src, bai_dest)
                else:
                    # Try to create index if it doesn't exist
                    try:
                        index_cmd = ['samtools', 'index', str(bam_dest)]
                        result = subprocess.run(index_cmd, stderr=subprocess.PIPE, text=True, timeout=60)
                        if result.returncode != 0:
                            logger.warning(f"Failed to create BAM index: {result.stderr}")
                    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
                        logger.warning(f"Failed to create BAM index: {e}")
                
                logger.info(f"Exported validation files to {igv_dir}: {ref_dest.name}, {bam_dest.name}")
                
            except Exception as e:
                logger.error(f"Failed to copy BAM file: {e}")
            
        except Exception as e:
            logger.warning(f"Failed to export validation files: {e}")
    
    def validate_all_candidates(self, candidates: List['ITDCandidate'],
                               min_supporting_reads: int = 3,
                               min_allele_frequency: float = 0.01) -> List[ValidationResult]:
        """
        Validate all ITD candidates
        
        Args:
            candidates: List of ITD candidates
            min_supporting_reads: Minimum supporting reads
            min_allele_frequency: Minimum allele frequency
            
        Returns:
            List of validated ITD results
        """
        validated_results = []
        
        logger.info(f"Validating {len(candidates)} ITD candidates using dual-reference approach")
        
        for i, candidate in enumerate(candidates, 1):
            actual_length = len(candidate.sequence)
            logger.info(f"Validating candidate {i}/{len(candidates)}: {actual_length}bp sequence "
                       f"(candidate.length={candidate.length})")
            
            result = self.validate_candidate(
                candidate, min_supporting_reads, min_allele_frequency
            )
            
            if result:
                validated_results.append(result)
                logger.info(f"Candidate {i} validated successfully")
            else:
                logger.info(f"Candidate {i} failed validation")
        
        # Sort by allele frequency
        validated_results.sort(key=lambda x: x.allele_frequency, reverse=True)
        
        logger.info(f"Dual-reference validation complete: {len(validated_results)} validated ITDs")
        
        return validated_results
    
    def cleanup(self):
        """Clean up temporary files"""
        try:
            if Path(self.temp_dir).exists():
                shutil.rmtree(self.temp_dir)
                logger.debug(f"Cleaned up temporary directory: {self.temp_dir}")
        except Exception as e:
            logger.warning(f"Failed to clean up temporary directory: {e}")


def validate_itd_candidates_dual_reference(candidates: List['ITDCandidate'], 
                                          reference_sequence: str,
                                          reads_fastq: str,
                                          min_supporting_reads: int = 3,
                                          min_allele_frequency: float = 0.01,
                                          config=None) -> List[ValidationResult]:
    """
    Main function to validate ITD candidates using dual-reference approach
    
    Args:
        candidates: List of ITD candidates
        reference_sequence: Original FLT3 reference sequence
        reads_fastq: Path to FASTQ file with original reads
        min_supporting_reads: Minimum supporting reads
        min_allele_frequency: Minimum allele frequency
        config: Configuration object
        
    Returns:
        List of validated ITD results
    """
    validator = DualReferenceValidator(
        reference_sequence=reference_sequence,
        original_reads_fastq=reads_fastq,
        config=config
    )
    
    try:
        results = validator.validate_all_candidates(
            candidates, min_supporting_reads, min_allele_frequency
        )
        return results
    finally:
        validator.cleanup()


if __name__ == "__main__":
    # Test functionality
    import argparse
    
    parser = argparse.ArgumentParser(description="Test dual-reference ITD validation")
    parser.add_argument("fastq_file", help="Input FASTQ file with reads")
    parser.add_argument("reference_file", help="Reference FASTA file")
    parser.add_argument("--itd-sequence", required=True, help="ITD sequence to test")
    parser.add_argument("--itd-position", type=int, required=True, help="ITD insertion position")
    parser.add_argument("--min-support", type=int, default=3, help="Minimum supporting reads")
    parser.add_argument("--min-af", type=float, default=0.01, help="Minimum allele frequency")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Read reference sequence
    with open(args.reference_file, 'r') as f:
        lines = f.readlines()
        reference_seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    
    # Create test ITD candidate
    from cigar_itd_detector import ITDCandidate
    test_candidate = ITDCandidate(
        sequence=args.itd_sequence,
        length=len(args.itd_sequence),
        position=args.itd_position,
        support_type='test',
        supporting_reads=['test_read'],
        confidence=0.8,
        insertion_site=args.itd_position
    )
    
    # Validate ITD
    results = validate_itd_candidates_dual_reference(
        candidates=[test_candidate],
        reference_sequence=reference_seq,
        reads_fastq=args.fastq_file,
        min_supporting_reads=args.min_support,
        min_allele_frequency=args.min_af
    )
    
    # Print results
    if results:
        print(f"Validation successful:")
        for result in results:
            print(f"  ITD: {result.duplication_length}bp")
            print(f"  Allele frequency: {result.allele_frequency:.3f}")
            print(f"  Supporting reads: {result.itd_coverage}")
            print(f"  Confidence: {result.validation_confidence:.3f}")
    else:
        print("Validation failed - no valid ITDs found")
