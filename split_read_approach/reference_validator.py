#!/usr/bin/env python3
"""
Reference Validator Module
Inserts ITD candidates into reference and validates with original long reads
"""

import logging
import subprocess
import tempfile
import pysam

from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from itd_generator import ITDCandidate
import shutil

logger = logging.getLogger(__name__)

def should_keep_file(filetype: str, config) -> bool:
    """Determine if a file should be kept based on config flags and file type.
    filetype: 'all', 'igv', 'fastq', 'temp', etc.
    """
    if getattr(config, 'debug', False):
        return True  # Keep everything
    if getattr(config, 'keep_temp', False):
        return True  # Keep all intermediates
    if getattr(config, 'bamout', False):
        return filetype == 'igv'  # Only keep IGV files
    return False  # Default: delete

def cleanup_intermediate_files(config, files_dict, logger=None):
    """Centralized cleanup of intermediate files based on config flags and filetype.
    files_dict: {'fastq': [...], 'igv': [...], 'temp': [...], 'igv_dirs': [...]}
    """
    for filetype, files in files_dict.items():
        # Special handling for IGV directories
        if filetype == 'igv_dirs':
            # Clean up IGV directories unless we want to keep IGV files
            if not should_keep_file('igv', config):
                for f in files:
                    try:
                        f_path = Path(f)
                        if f_path.is_dir():
                            shutil.rmtree(f_path)
                            if logger:
                                logger.info(f"Deleted IGV directory: {f}")
                    except Exception as e:
                        if logger:
                            logger.warning(f"Could not remove IGV directory {f}: {e}")
        elif not should_keep_file(filetype, config):
            for f in files:
                try:
                    f_path = Path(f)
                    if f_path.is_dir():
                        shutil.rmtree(f_path)
                        if logger:
                            logger.info(f"Deleted intermediate directory: {f}")
                    elif f_path.is_file():
                        f_path.unlink(missing_ok=True)
                        if logger:
                            logger.info(f"Deleted intermediate file: {f}")
                except Exception as e:
                    if logger:
                        logger.warning(f"Could not remove {f}: {e}")

@dataclass
class ValidationResult:
    """Results from validating an ITD candidate"""
    itd_candidate: 'ITDCandidate'
    supporting_reads: List[str]
    allele_frequency: float
    total_coverage: int
    itd_coverage: int
    validation_confidence: float
    insertion_position: int
    duplication_length: int

class ReferenceValidator:
    def _export_validation_igv(self, itd_candidate, ref_file, bam_file):
        """Export the custom reference and BAM for IGV, ensuring BAM is coordinate-sorted before indexing."""
        # Determine output directory
        outdir = getattr(self.config, 'validation_igv_dir', None)
        if not outdir:
            outdir = getattr(self.config, 'output_dir', '.')
        outdir = Path(outdir) / "debug"
        outdir.mkdir(parents=True, exist_ok=True)

        # Compose unique prefix for this ITD
        ref_prefix = f"itd_{itd_candidate.length}bp_{itd_candidate.position}_customref"
        ref_dest = outdir / f"{ref_prefix}.fa"
        bam_dest = outdir / f"{ref_prefix}.bam"
        bai_dest = outdir / f"{ref_prefix}.bam.bai"

        # Rewrite reference FASTA with correct header
        with open(ref_file, "r") as fin:
            lines = fin.readlines()
        lines[0] = f">{ref_prefix}\n"
        with open(ref_dest, "w") as fout:
            fout.writelines(lines)

        # Reheader BAM to match reference FASTA
        bam_tmp = str(outdir / f"{ref_prefix}.tmp.bam")
        with pysam.AlignmentFile(bam_file, "rb") as bam_in:
            header = bam_in.header.to_dict()
            if "SQ" in header:
                for sq in header["SQ"]:
                    sq["SN"] = ref_prefix
                    sq["LN"] = len("".join(lines[1:]).replace("\n", ""))
            with pysam.AlignmentFile(bam_tmp, "wb", header=header) as bam_out:
                for read in bam_in:
                    if not read.is_unmapped:
                        read.reference_id = 0
                        bam_out.write(read)

        # Always sort BAM by coordinate before indexing
        bam_sorted = str(outdir / f"{ref_prefix}.sorted.bam")
        try:
            pysam.sort("-o", bam_sorted, bam_tmp)
        except Exception as e:
            logger.error(f"pysam.sort failed: {e}")
            Path(bam_tmp).unlink(missing_ok=True)
            return

        # Move sorted BAM to final destination
        Path(bam_sorted).replace(bam_dest)
        Path(bam_tmp).unlink(missing_ok=True)

        # Index BAM
        try:
            pysam.index(str(bam_dest))
        except Exception as e:
            logger.error(f"pysam.index failed: {e}")
            return

        # Clean up any old BAI
        if Path(bam_dest.with_suffix(".bam.bai")).exists():
            Path(bam_dest.with_suffix(".bam.bai")).unlink()

        # Copy BAI to expected name
        if Path(str(bam_dest) + ".bai").exists():
            Path(str(bam_dest) + ".bai").rename(bai_dest)
        elif Path(str(bam_dest)[:-4] + ".bai").exists():
            Path(str(bam_dest)[:-4] + ".bai").rename(bai_dest)

        logger.info(f"Exported IGV files: {ref_dest}, {bam_dest}, {bai_dest}")
        
        # Return list of created IGV files for cleanup tracking
        return [str(ref_dest), str(bam_dest), str(bai_dest)]
    """Validate ITD candidates by inserting into reference and realigning"""
    
    def __init__(self, reference_sequence: str, original_reads: List[Dict] = None, fastq_file: str = None, config=None):
        self.reference_sequence = reference_sequence
        self.original_reads = original_reads
        self.fastq_file = fastq_file  # Path to FASTQ file if available
        self.temp_dir = tempfile.mkdtemp(prefix="flt3_itd_validation_")
        self.config = config  # Add config attribute
        self.igv_files = []  # Track IGV files for cleanup
        self.igv_dirs = set()  # Track unique IGV directories for cleanup
        
    def create_modified_reference(self, itd_candidate: 'ITDCandidate') -> Tuple[str, int]:
        """Create reference with ITD sequence inserted at appropriate position
        
        Returns:
            Tuple of (reference_file_path, actual_insertion_position_in_modified_ref)
        """
        
        # The ITD candidate sequence already represents the complete duplication
        # We just need to insert it at the right position in the reference
        
        # Determine insertion position
        if itd_candidate.insertion_site is not None:
            insert_pos = itd_candidate.insertion_site
        elif itd_candidate.duplication_start is not None:
            # Insert at the start of the region that was duplicated
            insert_pos = itd_candidate.duplication_start
        else:
            # Fallback: use the candidate position
            insert_pos = itd_candidate.position
        
        # Ensure insertion position is within reference bounds
        insert_pos = max(0, min(insert_pos, len(self.reference_sequence)))
        
        # For ITD validation, the candidate sequence represents the complete ITD (duplicated region)
        # We need to create a reference where this sequence is inserted at the insertion point
        # representing: original[0:pos] + ITD_sequence + original[pos:]
        # But we need to be careful about what the ITD sequence actually represents
        
        # Debug: let's see what we're working with
        logger.debug(f"ITD candidate: {len(itd_candidate.sequence)}bp at pos {insert_pos}")
        logger.debug(f"Original reference: {len(self.reference_sequence)}bp")
        logger.debug(f"ITD sequence (first 50bp): {itd_candidate.sequence[:50]}...")
        
        # For ITD validation, we need to insert the COMPLETE ITD sequence
        # The biological reality is that ITD reads contain the full duplication
        # Decomposition is useful for analysis, but for validation alignment,
        # we need the reference to match what the ITD reads actually contain
        
        logger.info(f"Creating validation reference with COMPLETE ITD sequence")
        logger.info(f"ITD length: {len(itd_candidate.sequence)}bp at position {insert_pos}")
        
        # Insert the complete ITD sequence for proper validation alignment
        modified_ref = (
            self.reference_sequence[:insert_pos] +
            itd_candidate.sequence +
            self.reference_sequence[insert_pos:]
        )
        
        itd_start_in_modified_ref = insert_pos
        itd_validation_length = len(itd_candidate.sequence)
        
        logger.info(f"Created validation reference: original {len(self.reference_sequence)}bp → "
                   f"modified {len(modified_ref)}bp, ITD region: "
                   f"{itd_start_in_modified_ref}-{itd_start_in_modified_ref + itd_validation_length}")
        
        # Store the validation length for later use
        self._last_itd_validation_length = itd_validation_length
        
        # Write to file
        ref_file = Path(self.temp_dir) / f"ref_with_itd_{itd_candidate.length}bp.fa"
        with open(ref_file, 'w') as f:
            f.write(f">FLT3_with_ITD_{itd_candidate.length}bp_pos{insert_pos}\n")
            f.write(modified_ref + "\n")
        
        return str(ref_file), itd_start_in_modified_ref
    
    def get_reads_fastq(self) -> str:
        """Return path to FASTQ file with original reads (do not rewrite if already exists)"""
        if self.fastq_file and Path(self.fastq_file).exists():
            return self.fastq_file
        # Fallback: write reads to FASTQ if only dicts are available
        fastq_file = Path(self.temp_dir) / "original_reads.fastq"
        if self.original_reads:
            with open(fastq_file, 'w') as f:
                for read in self.original_reads:
                    f.write(f"@{read['name']}\n")
                    f.write(f"{read['sequence']}\n")
                    f.write("+\n")
                    if read.get('qualities'):
                        f.write(''.join(chr(q + 33) for q in read['qualities']) + "\n")
                    else:
                        f.write('I' * len(read['sequence']) + "\n")
            return str(fastq_file)
        raise RuntimeError("No FASTQ file or original reads available for validation.")
    
    def align_to_modified_reference(self, ref_file: str, reads_file: str) -> str:
        """Align original reads to modified reference"""
        output_sam = Path(self.temp_dir) / "validation_alignment.sam"
        
        # Run minimap2
        cmd = [
            'minimap2',
            '-ax', 'map-ont',  # Nanopore preset
            '-t', '4',
            '--secondary=no',
            ref_file,
            reads_file
        ]
        
        with open(output_sam, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
            if result.returncode != 0:
                logger.warning(f"Minimap2 validation failed: {result.stderr.decode()}")
                return None
        
        # Convert to sorted BAM
        output_bam = Path(self.temp_dir) / "validation_alignment.bam"
        pysam.sort("-o", str(output_bam), str(output_sam))
        pysam.index(str(output_bam))
        
        return str(output_bam)
    
    def calculate_coverage(self, bam_file: str, itd_position: int, itd_length: int) -> Dict:
        """Calculate coverage by checking ALL reads against the ITD region, using CIGAR analysis for high-confidence support"""
        total_reads = 0
        itd_supporting_reads = []
        high_confidence_reads = []
        spanning_reads = 0
        non_spanning_reads = 0
        frame_shift_reads = 0
        wt_reads = 0
        ambiguous_reads = 0

        # Define ITD region in the modified reference
        itd_start = itd_position
        itd_end = itd_position + itd_length
        
        logger.info(f"Analyzing ITD region in modified reference: {itd_start}-{itd_end} (length {itd_length})")

        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                for read in bam.fetch():
                    if read.is_unmapped or read.is_duplicate or read.is_secondary:
                        continue

                    total_reads += 1

                    # Check if read spans the ITD region
                    read_start = read.reference_start
                    read_end = read.reference_end
                    
                    # Debug: log read coverage for first few reads
                    if total_reads <= 5:
                        logger.debug(f"Read {read.query_name}: {read_start}-{read_end}, "
                                   f"spans ITD: {read_start <= itd_start and read_end >= itd_end}")

                    if read_start <= itd_start and read_end >= itd_end:
                        spanning_reads += 1
                        # Use robust CIGAR analysis
                        cigar_result = self._analyze_cigar_for_itd(read, itd_start, itd_end, itd_length)
                        
                        # Debug: log CIGAR analysis for first few spanning reads
                        if spanning_reads <= 3:
                            logger.debug(f"CIGAR analysis for {read.query_name}: {cigar_result}")
                            # Additional frame-shift debugging
                            if cigar_result.get('status') == 'frame_shift':
                                logger.warning(f"Frame-shift detected in {read.query_name}: "
                                             f"soft_clips_outside={cigar_result.get('soft_clips_outside', [])}, "
                                             f"outside_match_ratio={cigar_result.get('outside_match_ratio', 0):.3f}")
                        
                        if cigar_result['status'] == 'itd_confirmed':
                            itd_supporting_reads.append(read.query_name)
                            if cigar_result.get('high_confidence', False):
                                high_confidence_reads.append(read.query_name)
                        elif cigar_result['status'] == 'frame_shift':
                            frame_shift_reads += 1
                        elif cigar_result['status'] == 'wt':
                            wt_reads += 1
                        else:
                            ambiguous_reads += 1
                    else:
                        non_spanning_reads += 1

        except Exception as e:
            logger.error(f"Error calculating coverage: {e}")
            return {
                'total_coverage': 0,
                'itd_coverage': 0,
                'allele_frequency': 0.0,
                'high_confidence_ratio': 0.0,
                'itd_supporting_reads': []
            }

        itd_coverage = len(itd_supporting_reads)
        allele_frequency = itd_coverage / total_reads if total_reads > 0 else 0.0
        high_confidence_ratio = (len(high_confidence_reads) / itd_coverage) if itd_coverage > 0 else 0.0

        logger.info(f"Coverage analysis results:")
        logger.info(f"  Total reads: {total_reads}")
        logger.info(f"  Spanning ITD region: {spanning_reads}")
        logger.info(f"  Non-spanning: {non_spanning_reads}")
        logger.info(f"  ITD supporting: {itd_coverage}")
        logger.info(f"  High confidence: {len(high_confidence_reads)}")
        logger.info(f"  Frame-shift rejected: {frame_shift_reads}")
        logger.info(f"  Wild-type: {wt_reads}")
        logger.info(f"  Ambiguous: {ambiguous_reads}")
        logger.info(f"  Allele frequency: {allele_frequency:.3f}")

        return {
            'total_coverage': total_reads,
            'itd_coverage': itd_coverage,
            'allele_frequency': allele_frequency,
            'high_confidence_ratio': high_confidence_ratio,
            'itd_supporting_reads': itd_supporting_reads
        }
    
    def _analyze_cigar_for_itd(self, read: pysam.AlignedRead, 
                               itd_start: int, itd_end: int, 
                               itd_length: int) -> Dict:
        """Analyze CIGAR string to confirm ITD presence and check for frame-shift artifacts"""
        if not read.cigartuples:
            return {'status': 'ambiguous', 'high_confidence': False}
        
        # Track position in reference and query
        ref_pos = read.reference_start
        query_pos = 0
        
        # Track if we've properly spanned the ITD region
        entered_itd = False
        exited_itd = False
        bases_in_itd = 0
        total_matches_in_itd = 0
        total_errors_in_itd = 0
        
        # Track alignment quality outside ITD region (to detect frame-shifts)
        soft_clips_outside_itd = []
        total_matches_outside_itd = 0
        total_bases_outside_itd = 0
        
        # Debug: Track major CIGAR operations in ITD region
        large_insertions_in_itd = []
        large_deletions_in_itd = []
        
        # Tolerance for nanopore errors
        max_error_rate = 0.15  # 15% error rate tolerance
        min_match_ratio = 0.7   # At least 70% of bases should match
        
        # Frame-shift detection thresholds
        max_soft_clip_outside_itd = 50  # Max soft-clipping allowed outside ITD
        min_match_ratio_outside_itd = 0.8  # High match ratio required outside ITD
        
        for op, length in read.cigartuples:
            # Check if we're in ITD region
            in_itd_region = ref_pos < itd_end and (ref_pos + length if op in [0, 2] else ref_pos) > itd_start
            
            if ref_pos >= itd_start and not entered_itd:
                entered_itd = True
            
            if ref_pos >= itd_end and entered_itd and not exited_itd:
                exited_itd = True
            
            # Process CIGAR operations
            if op == 0:  # Match/mismatch
                if in_itd_region:
                    # We're in the ITD region
                    overlap_start = max(ref_pos, itd_start)
                    overlap_end = min(ref_pos + length, itd_end)
                    overlap_length = max(0, overlap_end - overlap_start)
                    
                    bases_in_itd += overlap_length
                    total_matches_in_itd += overlap_length
                    
                    # Count bases outside ITD region in this operation
                    bases_outside_this_op = length - overlap_length
                    if bases_outside_this_op > 0:
                        total_matches_outside_itd += bases_outside_this_op
                        total_bases_outside_itd += bases_outside_this_op
                else:
                    # Entirely outside ITD region
                    total_matches_outside_itd += length
                    total_bases_outside_itd += length
                
                ref_pos += length
                query_pos += length
                
            elif op == 1:  # Insertion
                if in_itd_region:
                    # Track large insertions in ITD region
                    if length > 3:
                        large_insertions_in_itd.append(length)
                    # Small insertions are OK (nanopore errors)
                    if length <= 3:
                        total_errors_in_itd += length
                    else:
                        # Large insertion might indicate problem
                        total_errors_in_itd += length * 2
                # Note: insertions don't consume reference, so don't affect outside ITD tracking
                query_pos += length
                
            elif op == 2:  # Deletion
                if in_itd_region:
                    # Track deletions in ITD region
                    large_deletions_in_itd.append(length)
                    
                    # ANY deletion >5bp in ITD region indicates missing ITD (too strict for nanopore errors)
                    # For complex ITDs with novel insertions, any significant deletion suggests WT read
                    if length > 5:
                        logger.debug(f"Large deletion in ITD region: {length}bp - classifying as wild-type")
                        return {'status': 'wt', 'high_confidence': True, 'large_deletions': large_deletions_in_itd, 'reason': f'{length}bp_deletion_in_itd'}
                    else:
                        # Small deletions count as errors
                        total_errors_in_itd += length
                else:
                    # Deletion outside ITD region - count as error in outside region
                    total_bases_outside_itd += length
                ref_pos += length
                
            elif op == 4:  # Soft clip
                if in_itd_region and length > 10:
                    # Large soft clips in ITD region suggest alignment issues
                    return {'status': 'ambiguous', 'high_confidence': False, 'reason': 'large_soft_clip_in_itd'}
                elif not in_itd_region:
                    # Track soft clips outside ITD region (potential frame-shift indicator)
                    soft_clips_outside_itd.append(length)
                query_pos += length
                
            elif op == 5:  # Hard clip
                pass
        
        # Determine ITD support based on coverage and frame-shift analysis
        if entered_itd and exited_itd:
            # Check for frame-shift artifacts
            total_soft_clips_outside = sum(soft_clips_outside_itd)
            outside_match_ratio = (total_matches_outside_itd / total_bases_outside_itd) if total_bases_outside_itd > 0 else 1.0
            
            # Frame-shift detection
            is_frame_shift = (
                total_soft_clips_outside > max_soft_clip_outside_itd or
                (total_bases_outside_itd > 50 and outside_match_ratio < min_match_ratio_outside_itd)
            )
            
            if is_frame_shift:
                return {
                    'status': 'frame_shift', 
                    'high_confidence': False,
                    'soft_clips_outside': soft_clips_outside_itd,
                    'outside_match_ratio': outside_match_ratio,
                    'reason': 'poor_alignment_outside_itd'
                }
            
            # CRITICAL: Check if there were significant deletions in the ITD region
            # If there are large deletions, this is likely a wild-type read
            total_deletion_length = sum(large_deletions_in_itd)
            if total_deletion_length > itd_length * 0.15:  # >15% of ITD length deleted
                logger.debug(f"Large deletions totaling {total_deletion_length}bp in ITD region - classifying as wild-type")
                return {
                    'status': 'wt', 
                    'high_confidence': True, 
                    'large_deletions': large_deletions_in_itd,
                    'reason': f'total_{total_deletion_length}bp_deletions_in_itd'
                }
            
            # Calculate match ratio in ITD region
            if bases_in_itd > 0:
                match_ratio = total_matches_in_itd / bases_in_itd
                error_rate = total_errors_in_itd / bases_in_itd
                
                # For ITD support, we need good coverage WITHOUT large deletions
                # AND the read should show evidence of the insertion (not just spanning it)
                
                # High confidence if good coverage with low errors AND good outside alignment AND no large deletions
                high_confidence = (
                    bases_in_itd >= itd_length * 0.8 and  # Covers 80% of ITD
                    match_ratio >= min_match_ratio and
                    error_rate <= max_error_rate and
                    not is_frame_shift and  # No frame-shift artifacts
                    total_deletion_length <= 3  # Minimal deletions in ITD region
                )
                
                # Determine status - be more strict about what counts as ITD support
                if bases_in_itd >= itd_length * 0.7:  # At least 70% coverage (more strict)
                    if error_rate <= max_error_rate and total_deletion_length <= itd_length * 0.1:  # Max 10% deletion
                        return {
                            'status': 'itd_confirmed', 
                            'high_confidence': high_confidence,
                            'match_ratio_itd': match_ratio,
                            'match_ratio_outside': outside_match_ratio,
                            'soft_clips_outside': total_soft_clips_outside,
                            'total_deletions': total_deletion_length
                        }
                    else:
                        return {'status': 'ambiguous', 'high_confidence': False, 'reason': 'high_error_or_deletion_rate_in_itd'}
                else:
                    # Poor coverage of ITD region
                    return {'status': 'wt', 'high_confidence': False, 'reason': 'poor_itd_coverage'}
            else:
                # No bases aligned to ITD region
                return {'status': 'wt', 'high_confidence': True, 'reason': 'no_itd_coverage'}
        else:
            # Didn't properly span ITD region
            return {'status': 'ambiguous', 'high_confidence': False, 'reason': 'incomplete_spanning'}
    
    def validate_candidate(self, itd_candidate: 'ITDCandidate', 
                          min_supporting_reads: int = 3,
                          min_high_confidence_ratio: float = 0.5,
                          max_itd_length: int = 500) -> Optional[ValidationResult]:
        """Validate a single ITD candidate with strict CIGAR checking"""
        logger.info(f"Validating ITD candidate: {itd_candidate.length}bp at position {itd_candidate.position}")
        
        # Apply length filter as hard filter to save computation time
        if itd_candidate.length > max_itd_length:
            logger.info(f"Skipping candidate: length {itd_candidate.length}bp exceeds maximum {max_itd_length}bp")
            return None
        
        # Create modified reference
        ref_file, itd_position_in_modified_ref = self.create_modified_reference(itd_candidate)
        
        # Use existing FASTQ file for reads
        reads_file = self.get_reads_fastq()
        
        # Align to modified reference
        bam_file = self.align_to_modified_reference(ref_file, reads_file)
        if not bam_file:
            return None
        
        # Calculate coverage with CIGAR validation using CORRECT ITD coordinates
        # Use the stored validation length which accounts for novel + duplicated content
        validation_length = getattr(self, '_last_itd_validation_length', itd_candidate.length)
        
        coverage_stats = self.calculate_coverage(
            bam_file, 
            itd_position_in_modified_ref,  # Use position in modified reference
            validation_length  # Use the decomposed ITD length
        )

        # Export for IGV only if bamout, keep_temp, or debug is set
        if (getattr(self.config, 'bamout', False) or getattr(self.config, 'keep_temp', False) or getattr(self.config, 'debug', False)):
            igv_files = self._export_validation_igv(itd_candidate, ref_file, bam_file)
            if igv_files:
                self.igv_files.extend(igv_files)
                # Track the debug directory
                outdir = getattr(self.config, 'validation_igv_dir', None)
                if not outdir:
                    outdir = getattr(self.config, 'output_dir', '.')
                debug_dir = Path(outdir) / "debug"
                self.igv_dirs.add(str(debug_dir))
        
        # More stringent validation criteria
        passes_validation = (
            coverage_stats['itd_coverage'] >= min_supporting_reads and
            coverage_stats['high_confidence_ratio'] >= min_high_confidence_ratio and
            coverage_stats['allele_frequency'] >= self.config.min_allele_frequency
        )
        
        if passes_validation:
            # Calculate validation confidence based on multiple factors
            validation_confidence = self._calculate_validation_confidence(
                itd_candidate, coverage_stats
            )
            
            result = ValidationResult(
                itd_candidate=itd_candidate,
                supporting_reads=coverage_stats['itd_supporting_reads'],
                allele_frequency=coverage_stats['allele_frequency'],
                total_coverage=coverage_stats['total_coverage'],
                itd_coverage=coverage_stats['itd_coverage'],
                validation_confidence=validation_confidence,
                insertion_position=itd_candidate.insertion_site or itd_candidate.position,
                duplication_length=itd_candidate.length
            )
            
            logger.info(f"Validation successful: AF={result.allele_frequency:.3f}, "
                       f"Support={result.itd_coverage}/{result.total_coverage}, "
                       f"High-confidence ratio={coverage_stats['high_confidence_ratio']:.2f}")
            
            return result
        else:
            logger.info(f"Validation failed: AF={coverage_stats['allele_frequency']:.3f}, "
                       f"Support={coverage_stats['itd_coverage']}, "
                       f"High-confidence ratio={coverage_stats['high_confidence_ratio']:.2f}")
            return None
    
    def _calculate_validation_confidence(self, candidate: 'ITDCandidate', 
                                       coverage_stats: Dict) -> float:
        """Calculate validation confidence based on multiple factors"""
        # Base confidence from original detection
        base_conf = candidate.confidence * 0.3
        
        # Confidence from allele frequency
        af_conf = min(0.3, coverage_stats['allele_frequency'] * 0.3)
        
        # Confidence from high-quality read ratio
        hq_conf = coverage_stats['high_confidence_ratio'] * 0.2
        
        # Confidence from read support
        support_conf = min(0.2, (coverage_stats['itd_coverage'] / 10) * 0.2)
        
        total_confidence = base_conf + af_conf + hq_conf + support_conf
        
        return min(0.95, total_confidence)
    
    def validate_all_candidates(self, candidates: List['ITDCandidate'],
                               min_supporting_reads: int = 3,
                               min_high_confidence_ratio: float = 0.5,
                               max_itd_length: int = 500) -> List[ValidationResult]:
        """Validate all ITD candidates with enhanced biological event deduplication"""
        validated_results = []
        
        # RESPECT the pre-filtering order from main module
        # DO NOT re-sort - the main module has already ranked candidates optimally
        # by supporting reads count which is the best predictor of real ITDs
        sorted_candidates = candidates  # Use the order provided by main module
        
        for candidate in sorted_candidates:
            # Validate this candidate
            result = self.validate_candidate(candidate, min_supporting_reads, 
                                           min_high_confidence_ratio, max_itd_length)
            
            if result:
                # Check if this represents the same biological event as existing results
                is_duplicate_event = False
                for existing in validated_results:
                    if self._are_itds_same_biological_event(existing, result):
                        logger.info(f"Skipping validation - same biological event already found: "
                                   f"{candidate.length}bp at {candidate.position} "
                                   f"(similar to {existing.duplication_length}bp, AF={existing.allele_frequency:.3f})")
                        is_duplicate_event = True
                        break
                
                if not is_duplicate_event:
                    validated_results.append(result)
                    
                    # Early stopping if we have enough high-confidence results
                    if len(validated_results) >= 10 and result.validation_confidence < 0.7:
                        logger.info("Stopping validation - sufficient high-confidence ITDs found")
                        break
        
        # Sort by allele frequency
        validated_results.sort(key=lambda x: x.allele_frequency, reverse=True)
        
        # Final deduplication of results (in case any were missed in the loop)
        final_results = self._deduplicate_results(validated_results)
        
        return final_results
    
    def _calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity between two ITD sequences"""
        if not seq1 or not seq2:
            return 0.0
        
        # Handle length differences
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))
        
        if min_len == 0:
            return 0.0
        
        # Check for substring relationships (one ITD is contained in another)
        if seq1 in seq2 or seq2 in seq1:
            return min_len / max_len  # Similarity based on overlap
        
        # Check for significant overlap at start or end
        overlap_threshold = min(min_len * 0.7, 30)  # At least 70% of shorter sequence or 30bp
        
        # Check if sequences overlap at start
        for i in range(int(overlap_threshold), min_len + 1):
            if seq1[:i] == seq2[:i]:
                return i / max_len
        
        # Check if sequences overlap at end
        for i in range(int(overlap_threshold), min_len + 1):
            if seq1[-i:] == seq2[-i:]:
                return i / max_len
        
        # Check for internal overlaps (one sequence contained within another with gaps)
        longer_seq = seq1 if len(seq1) > len(seq2) else seq2
        shorter_seq = seq2 if len(seq1) > len(seq2) else seq1
        
        if shorter_seq in longer_seq:
            return len(shorter_seq) / len(longer_seq)
        
        return 0.0
    
    def _are_itds_same_biological_event(self, result1: ValidationResult, result2: ValidationResult) -> bool:
        """Determine if two ITD results represent the same biological event"""
        
        # Check allele frequency similarity (should be very close for same biological event)
        af_diff = abs(result1.allele_frequency - result2.allele_frequency)
        if af_diff > 0.05:  # More than 5% difference suggests different events
            return False
        
        # Check supporting read count similarity
        support_ratio = min(result1.itd_coverage, result2.itd_coverage) / max(result1.itd_coverage, result2.itd_coverage)
        if support_ratio < 0.8:  # Less than 80% similarity in read support
            return False
        
        # Check sequence similarity
        seq_similarity = self._calculate_sequence_similarity(
            result1.itd_candidate.sequence, 
            result2.itd_candidate.sequence
        )
        
        if seq_similarity >= 0.7:  # 70% sequence similarity
            logger.info(f"ITDs likely represent same biological event: "
                       f"AF_diff={af_diff:.3f}, support_ratio={support_ratio:.3f}, "
                       f"seq_similarity={seq_similarity:.3f}")
            return True
        
        # Check for overlapping genomic regions with high AF similarity
        pos1_start = getattr(result1.itd_candidate, 'duplication_start', None)
        if pos1_start is None:
            pos1_start = result1.insertion_position
        
        pos1_end = getattr(result1.itd_candidate, 'duplication_end', None)
        if pos1_end is None:
            pos1_end = result1.insertion_position + result1.duplication_length
            
        pos2_start = getattr(result2.itd_candidate, 'duplication_start', None)
        if pos2_start is None:
            pos2_start = result2.insertion_position
            
        pos2_end = getattr(result2.itd_candidate, 'duplication_end', None)
        if pos2_end is None:
            pos2_end = result2.insertion_position + result2.duplication_length
        
        # Additional safety check - ensure all positions are valid integers
        if any(pos is None for pos in [pos1_start, pos1_end, pos2_start, pos2_end]):
            logger.warning("Cannot compare ITDs - missing position information")
            return False
        
        # Calculate overlap
        overlap_start = max(pos1_start, pos2_start)
        overlap_end = min(pos1_end, pos2_end)
        overlap_length = max(0, overlap_end - overlap_start)
        
        if overlap_length > 0:
            total_span = max(pos1_end, pos2_end) - min(pos1_start, pos2_start)
            overlap_ratio = overlap_length / total_span
            
            # If significant genomic overlap AND very similar AF, likely same event
            if overlap_ratio >= 0.5 and af_diff <= 0.02:
                logger.info(f"ITDs likely represent same biological event: "
                           f"genomic_overlap={overlap_ratio:.3f}, AF_diff={af_diff:.3f}")
                return True
        
        return False

    def _deduplicate_results(self, results: List[ValidationResult]) -> List[ValidationResult]:
        """Remove duplicate validation results using enhanced biological event detection"""
        if not results:
            return []
        
        # Sort by allele frequency descending to prioritize higher AF ITDs
        results_sorted = sorted(results, key=lambda x: x.allele_frequency, reverse=True)
        
        deduped = []
        
        for result in results_sorted:
            is_duplicate = False
            
            # Check against all existing results
            for existing in deduped:
                if self._are_itds_same_biological_event(existing, result):
                    # Merge information from the duplicate
                    logger.info(f"Merging duplicate ITD: {result.duplication_length}bp "
                               f"(AF={result.allele_frequency:.3f}) into "
                               f"{existing.duplication_length}bp (AF={existing.allele_frequency:.3f})")
                    
                    # Keep the one with higher allele frequency as primary
                    if result.allele_frequency > existing.allele_frequency:
                        # Replace existing with this better result
                        idx = deduped.index(existing)
                        # Merge supporting reads
                        result.supporting_reads.extend(existing.supporting_reads)
                        result.supporting_reads = list(set(result.supporting_reads))
                        deduped[idx] = result
                    else:
                        # Merge supporting reads into existing
                        existing.supporting_reads.extend(result.supporting_reads)
                        existing.supporting_reads = list(set(existing.supporting_reads))
                    
                    is_duplicate = True
                    break
            
            if not is_duplicate:
                deduped.append(result)
        
        # Final sort by allele frequency
        deduped.sort(key=lambda x: x.allele_frequency, reverse=True)
        
        logger.info(f"Deduplication: {len(results)} → {len(deduped)} ITDs after biological event merging")
        
        return deduped
    
    def cleanup(self):
        """Clean up temporary validation files only. FASTQ cleanup is handled centrally."""
        # Remove temp_dir
        if Path(self.temp_dir).exists():
            shutil.rmtree(self.temp_dir)

def validate_itd_candidates(candidates: List['ITDCandidate'], 
                           reference_sequence: str,
                           original_reads: List[Dict],
                           min_supporting_reads: int = 3,
                           min_high_confidence_ratio: float = 0.5,
                           max_itd_length: int = 500,
                           config=None) -> Tuple[List[ValidationResult], List[str], List[str]]:
    """Main function to validate ITD candidates. Returns results, IGV files list, and IGV directories list."""
    # Accept fastq_file as part of config if available
    fastq_file = getattr(config, 'fastq_file', None) if config else None
    validator = ReferenceValidator(reference_sequence, original_reads, fastq_file=fastq_file, config=config)
    try:
        results = validator.validate_all_candidates(candidates, min_supporting_reads, 
                                                   min_high_confidence_ratio, max_itd_length)
        return results, validator.igv_files, list(validator.igv_dirs)
    finally:
        validator.cleanup()