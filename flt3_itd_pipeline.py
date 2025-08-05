#!/usr/bin/env python3
"""
FLT3 ITD Detection Pipeline Orchestrator
Integrates CIGAR-based detection, soft-clip analysis, and dual-reference validation
"""

import logging
import time
import subprocess
import tempfile
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
from dataclasses import dataclass
import json
import pysam
import numpy as np

# Import our detection modules
from cigar_itd_detector import CigarITDDetector, ITDCandidate, detect_cigar_itds
from softclip_itd_detector import SoftClipITDDetector, SoftClipITDCandidate, detect_softclip_itds
from dual_reference_validator import DualReferenceValidator, ValidationResult, validate_itd_candidates_dual_reference

logger = logging.getLogger(__name__)

@dataclass
class ITDDetectionResult:
    """Comprehensive ITD detection result"""
    cigar_candidates: List[ITDCandidate]
    softclip_candidates: List[SoftClipITDCandidate]
    validated_candidates: List[Union[ITDCandidate, SoftClipITDCandidate]]
    validation_results: List[ValidationResult]
    detection_summary: Dict
    runtime_stats: Dict

class FLT3ITDPipeline:
    """
    Main pipeline orchestrator for FLT3 ITD detection
    Implements the new architecture: CIGAR detection → Soft-clip analysis → Dual-reference validation
    """
    
    def __init__(self, reference_sequence: str, reference_file: str,
                 trimmed_fastq: str = None,
                 min_itd_length: int = 15, max_itd_length: int = 500,
                 min_support: int = 3, min_softclip_length: int = 50,
                 enable_softclip_fallback: bool = True,
                 size_analysis_results: dict = None, config=None):
        """
        Initialize the FLT3 ITD detection pipeline
        
        Args:
            reference_sequence: FLT3 reference sequence string
            reference_file: Path to reference FASTA file
            trimmed_fastq: Path to trimmed FASTQ file for re-alignment
            min_itd_length: Minimum ITD size to detect
            max_itd_length: Maximum ITD size to detect
            min_support: Minimum supporting reads required
            min_softclip_length: Minimum soft-clip length for secondary analysis
            enable_softclip_fallback: Whether to use soft-clip analysis if CIGAR detection is insufficient
            size_analysis_results: Results from read size analysis for primary marking
        """
        self.reference_sequence = reference_sequence
        self.reference_file = reference_file
        self.trimmed_fastq = trimmed_fastq
        self.min_itd_length = min_itd_length
        self.max_itd_length = max_itd_length
        self.min_support = min_support
        self.min_softclip_length = min_softclip_length
        self.enable_softclip_fallback = enable_softclip_fallback
        self.size_analysis_results = size_analysis_results or {}
        self.config = config
        
        # Get threads from config or use default
        self.threads = config.threads if config else 4
        
        logger.info(f"Initialized FLT3 ITD Pipeline with new architecture:")
        logger.info(f"  ITD size range: {min_itd_length}-{max_itd_length}bp")
        logger.info(f"  Minimum support: {min_support} reads")
        logger.info(f"  Soft-clip fallback: {'enabled' if enable_softclip_fallback else 'disabled'}")
        logger.info(f"  Threads: {self.threads}")
        
    def detect_itds(self, bam_file: str, output_prefix: Optional[str] = None) -> ITDDetectionResult:
        """
        Main ITD detection pipeline
        
        Args:
            bam_file: Path to BAM file with FLT3 reads
            output_prefix: Optional prefix for output files
            
        Returns:
            Comprehensive detection results
        """
        start_time = time.time()
        runtime_stats = {}
        
        # Store output_prefix for use in validation
        self.output_prefix = output_prefix
        
        logger.info(f"Starting FLT3 ITD detection pipeline on {bam_file}")
        
        # Step 1: Primary detection using CIGAR analysis
        logger.info("=== Step 1: CIGAR-based ITD Detection ===")
        cigar_start = time.time()
        
        cigar_candidates = self._detect_cigar_itds(bam_file)
        runtime_stats['cigar_detection_time'] = time.time() - cigar_start
        runtime_stats['cigar_candidates'] = len(cigar_candidates)
        
        logger.info(f"CIGAR detection found {len(cigar_candidates)} candidates")
        
        # Step 1.5: Left-align ITD positions to normalize alignment ambiguity
        logger.info("=== Step 1.5: Left-alignment Normalization ===")
        left_aligned_candidates = self._left_align_candidates(cigar_candidates)
        logger.info(f"Left-aligned {len(cigar_candidates)} candidates")
        
        # Step 1.6: Merge similar CIGAR candidates into consensus ITDs
        logger.info("=== Step 1.6: Consensus Merging ===")
        merged_candidates = self._merge_similar_candidates(
            left_aligned_candidates,
            position_tolerance=50,  # Increased tolerance for tandem repeat regions
            length_tolerance=10     # Increased tolerance for length differences
        )
        logger.info(f"Merged {len(left_aligned_candidates)} candidates into {len(merged_candidates)} consensus ITDs")
        
        # Step 1.7: Mark primary candidates based on size analysis
        if self.size_analysis_results and merged_candidates:
            logger.info("=== Step 1.7: Primary Candidate Marking ===")
            self._mark_primary_candidates(merged_candidates)
            primary_count = sum(1 for c in merged_candidates if c.is_primary)
            logger.info(f"Marked {primary_count} candidates as primary based on size analysis")
        
        # Step 2: Secondary detection using soft-clip analysis (if needed)
        softclip_candidates = []
        if self._should_run_softclip_analysis(merged_candidates):
            logger.info("=== Step 2: Soft-clip ITD Detection ===")
            softclip_start = time.time()
            
            softclip_candidates = self._detect_softclip_itds(bam_file)
            runtime_stats['softclip_detection_time'] = time.time() - softclip_start
            runtime_stats['softclip_candidates'] = len(softclip_candidates)
            
            logger.info(f"Soft-clip detection found {len(softclip_candidates)} candidates")
        else:
            logger.info("=== Step 2: Soft-clip Detection Skipped ===")
            logger.info("Sufficient CIGAR candidates found, skipping soft-clip analysis")
            runtime_stats['softclip_detection_time'] = 0
            runtime_stats['softclip_candidates'] = 0
        
        # Step 3: Multifasta validation with CIGAR filtering
        logger.info("=== Step 3: Multifasta Validation with CIGAR Filtering ===")
        validation_start = time.time()
        
        # NEW MULTIFASTA VALIDATION PROCESS:
        # 1. Create multifasta reference with WT + all ITD candidate sequences
        # 2. Align all reads to this multifasta reference once
        # 3. Parse alignments to count reads mapped to each reference (WT, ITD_1, ITD_2, etc.)
        # 4. Filter reads based on CIGAR: exclude reads with large indels (>10bp) or excessive soft clipping
        # 5. Filter reads with insertion sizes that don't match the ITD candidate size
        # 6. Calculate allele frequency and validate candidates based on sufficient support
        # 7. Report validated ITDs with proper read support counts
        
        all_candidates = merged_candidates + softclip_candidates
        validated_candidates, validation_results = self._validate_candidates(all_candidates, bam_file)
        runtime_stats['validation_time'] = time.time() - validation_start
        runtime_stats['validated_candidates'] = len(validated_candidates)
        
        logger.info(f"Validation: {len(all_candidates)} → {len(validated_candidates)} candidates")
        
        # Step 4: Generate detection summary
        detection_summary = self._generate_detection_summary(
            merged_candidates, softclip_candidates, validated_candidates, validation_results
        )
        
        runtime_stats['total_time'] = time.time() - start_time
        
        # Step 5: Save results if output prefix provided
        if output_prefix:
            self._save_results(output_prefix, detection_summary, runtime_stats)
        
        result = ITDDetectionResult(
            cigar_candidates=merged_candidates,
            softclip_candidates=softclip_candidates,
            validated_candidates=validated_candidates,
            validation_results=validation_results,
            detection_summary=detection_summary,
            runtime_stats=runtime_stats
        )
        
        logger.info(f"Pipeline completed in {runtime_stats['total_time']:.2f}s")
        return result
    
    def _detect_cigar_itds(self, bam_file: str) -> List[ITDCandidate]:
        """
        Run CIGAR-based ITD detection
        
        Args:
            bam_file: Path to BAM file
            
        Returns:
            List of CIGAR ITD candidates
        """
        try:
            candidates = detect_cigar_itds(
                bam_file=bam_file,
                reference_sequence=self.reference_sequence,
                min_itd_length=self.min_itd_length,
                max_itd_length=self.max_itd_length,
                min_support=self.min_support
            )
            
            logger.info(f"CIGAR detection: {len(candidates)} candidates found")
            for i, candidate in enumerate(candidates[:5]):  # Log first 5
                logger.debug(f"  CIGAR candidate {i+1}: {candidate.length}bp at pos {candidate.position}, "
                           f"support={len(candidate.supporting_reads)}")
            
            return candidates
            
        except Exception as e:
            logger.error(f"CIGAR detection failed: {e}")
            return []
    
    def _detect_softclip_itds(self, bam_file: str) -> List[SoftClipITDCandidate]:
        """
        Run soft-clip based ITD detection
        
        Args:
            bam_file: Path to BAM file
            
        Returns:
            List of soft-clip ITD candidates
        """
        try:
            candidates = detect_softclip_itds(
                bam_file=bam_file,
                reference_sequence=self.reference_sequence,
                reference_file=self.reference_file,
                min_itd_length=self.min_itd_length,
                max_itd_length=self.max_itd_length,
                min_support=self.min_support,
                min_softclip_length=self.min_softclip_length,
                threads=self.threads
            )
            
            logger.info(f"Soft-clip detection: {len(candidates)} candidates found")
            for i, candidate in enumerate(candidates[:5]):  # Log first 5
                logger.debug(f"  Soft-clip candidate {i+1}: {candidate.length}bp at pos {candidate.position}, "
                           f"support={len(candidate.supporting_reads)}")
            
            return candidates
            
        except Exception as e:
            logger.error(f"Soft-clip detection failed: {e}")
            return []
    
    def _should_run_softclip_analysis(self, cigar_candidates: List[ITDCandidate]) -> bool:
        """
        Decide whether to run soft-clip analysis based on CIGAR results
        
        Args:
            cigar_candidates: CIGAR detection results
            
        Returns:
            True if soft-clip analysis should be run
        """
        if not self.enable_softclip_fallback:
            return False
        
        # Run soft-clip analysis if:
        # 1. No CIGAR candidates found
        # 2. Only low-confidence CIGAR candidates
        # 3. Very few CIGAR candidates relative to expected
        
        if not cigar_candidates:
            logger.info("No CIGAR candidates found, running soft-clip analysis")
            return True
        
        high_confidence_candidates = [c for c in cigar_candidates if c.confidence > 0.8]
        if not high_confidence_candidates:
            logger.info("No high-confidence CIGAR candidates, running soft-clip analysis")
            return True
        
        if len(cigar_candidates) < 2:  # Very few candidates
            logger.info("Few CIGAR candidates found, running soft-clip analysis for completeness")
            return True
        
        logger.info("Sufficient CIGAR candidates found, skipping soft-clip analysis")
        return False
    
    def _merge_similar_candidates(self, candidates: List[ITDCandidate], 
                                position_tolerance: int = 10, 
                                length_tolerance: int = 3) -> List[ITDCandidate]:
        """
        Merge similar ITD candidates into consensus ITDs
        
        Args:
            candidates: List of ITD candidates to merge
            position_tolerance: Maximum position difference for merging (bp)
            length_tolerance: Maximum length difference for merging (bp)
            
        Returns:
            List of merged consensus ITD candidates
        """
        if not candidates:
            return []
        
        merged_candidates = []
        used_indices = set()
        
        for i, candidate in enumerate(candidates):
            if i in used_indices:
                continue
                
            # Find all similar candidates
            similar_candidates = [candidate]
            similar_indices = {i}
            
            for j, other_candidate in enumerate(candidates):
                if j <= i or j in used_indices:
                    continue
                    
                # Check if candidates are similar
                position_diff = abs(candidate.position - other_candidate.position)
                length_diff = abs(candidate.length - other_candidate.length)
                
                # Also check sequence similarity for tandem repeat regions
                sequence_similarity = self._calculate_sequence_similarity(
                    candidate.sequence, other_candidate.sequence
                )
                
                # Candidates are similar if they meet position AND length criteria
                # OR if they have high sequence similarity (indicating same biological event)
                position_similar = (position_diff <= position_tolerance and 
                                  length_diff <= length_tolerance)
                sequence_similar = sequence_similarity >= 0.7  # 70% sequence similarity
                
                if position_similar or sequence_similar:
                    similar_candidates.append(other_candidate)
                    similar_indices.add(j)
                    if sequence_similar and not position_similar:
                        logger.debug(f"Merged by sequence similarity: {sequence_similarity:.2f} "
                                   f"(pos diff: {position_diff}bp, len diff: {length_diff}bp)")
            
            # If we found similar candidates, merge them
            if len(similar_candidates) > 1:
                consensus_candidate = self._create_consensus_candidate(similar_candidates)
                merged_candidates.append(consensus_candidate)
                logger.debug(f"Merged {len(similar_candidates)} similar candidates into consensus: "
                           f"{consensus_candidate.length}bp at pos {consensus_candidate.position}")
            else:
                merged_candidates.append(candidate)
            
            used_indices.update(similar_indices)
        
        return merged_candidates
    
    def _calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate sequence similarity between two ITD sequences
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Similarity score between 0.0 and 1.0
        """
        if not seq1 or not seq2:
            return 0.0
        
        if seq1 == seq2:
            return 1.0
        
        # Method 1: Longest common subsequence similarity
        lcs_length = self._longest_common_subsequence(seq1, seq2)
        lcs_similarity = (2.0 * lcs_length) / (len(seq1) + len(seq2))
        
        # Method 2: Check for substantial overlap (one sequence contained in the other)
        if seq1 in seq2 or seq2 in seq1:
            overlap_similarity = min(len(seq1), len(seq2)) / max(len(seq1), len(seq2))
            return max(lcs_similarity, overlap_similarity)
        
        # Method 3: Check for shared k-mers (for fragmented matches)
        kmer_similarity = self._kmer_similarity(seq1, seq2, k=8)
        
        return max(lcs_similarity, kmer_similarity)

    def _validate_candidates_dual_reference(self, candidates: List[ITDCandidate], bam_file: str) -> Tuple[List[ITDCandidate], List[Dict]]:
        """
        Validate ITD candidates using multifasta approach with CIGAR filtering
        
        Args:
            candidates: List of ITD candidates to validate
            bam_file: Path to BAM file with aligned reads
            
        Returns:
            Tuple of (validated_candidates, validation_results)
        """
        try:
            logger.info(f"Starting multifasta validation for {len(candidates)} candidates")
            
            # Write multifasta for all candidates and WT
            output_prefix = getattr(self, 'output_prefix', None)
            multifasta_path = self._write_multifasta_for_validation(candidates, output_prefix)
            
            # Use trimmed FASTQ if available, otherwise extract reads from BAM
            reads_file = getattr(self, 'trimmed_fastq', None)
            if not reads_file:
                reads_file = self._extract_reads_from_bam(bam_file, output_prefix)
            
            # Align reads to multifasta using minimap2
            sam_output = self._align_reads_to_multifasta(reads_file, multifasta_path, output_prefix)
            
            # Parse alignments and validate with CIGAR filtering
            validated_candidates, validation_results = self._parse_multifasta_alignments(
                sam_output, candidates, multifasta_path
            )
            
            logger.info(f"Validation complete: {len(validated_candidates)}/{len(candidates)} candidates passed")
            return validated_candidates, validation_results
            
        except Exception as e:
            logger.error(f"Multifasta validation failed: {e}")
            return candidates, []  # Return unvalidated candidates if validation fails
    
    def _extract_reads_from_bam(self, bam_file: str, output_prefix: Optional[str] = None) -> str:
        """
        Extract reads from BAM file to FASTQ for alignment
        
        Args:
            bam_file: Path to input BAM file
            output_prefix: Optional prefix for output files
            
        Returns:
            Path to extracted FASTQ file
        """
        try:
            # Generate output FASTQ filename
            if output_prefix:
                fastq_file = f"{output_prefix}_reads_for_validation.fastq"
            else:
                fastq_file = "reads_for_validation.fastq"
            
            # Use samtools to convert BAM to FASTQ
            cmd = [
                "samtools", "fastq", 
                "-F", "2304",  # Exclude secondary and supplementary alignments
                "-o", fastq_file,
                bam_file
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info(f"Extracted reads to {fastq_file}")
            return fastq_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to extract reads from BAM: {e.stderr}")
            raise
        except Exception as e:
            logger.error(f"Error extracting reads: {e}")
            raise
    
    def _align_reads_to_multifasta(self, reads_file: str, multifasta_path: str, 
                                 output_prefix: Optional[str] = None) -> str:
        """
        Align reads to multifasta reference using minimap2
        
        Args:
            reads_file: Path to FASTQ file with reads
            multifasta_path: Path to multifasta reference
            output_prefix: Optional prefix for output files
            
        Returns:
            Path to SAM output file
        """
        try:
            # Generate output SAM filename
            if output_prefix:
                sam_file = f"{output_prefix}_multifasta_alignment.sam"
            else:
                sam_file = "multifasta_alignment.sam"
            
            # Run minimap2 alignment
            cmd = [
                "minimap2", 
                "-ax", "map-ont",  # Oxford Nanopore preset
                "--secondary=no",  # No secondary alignments
                "--max-intron-len", "0",  # No introns expected
                multifasta_path,
                reads_file
            ]
            
            with open(sam_file, 'w') as sam_out:
                result = subprocess.run(cmd, stdout=sam_out, stderr=subprocess.PIPE, 
                                      text=True, check=True)
            
            logger.info(f"Alignment complete: {sam_file}")
            return sam_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Minimap2 alignment failed: {e.stderr}")
            raise
        except Exception as e:
            logger.error(f"Error during alignment: {e}")
            raise
    
    def _parse_multifasta_alignments(self, sam_file: str, candidates: List[ITDCandidate], 
                                   multifasta_path: str) -> Tuple[List[ITDCandidate], List[Dict]]:
        """
        Parse SAM alignments and validate candidates with CIGAR filtering
        
        Args:
            sam_file: Path to SAM alignment file
            candidates: Original ITD candidates
            multifasta_path: Path to multifasta reference
            
        Returns:
            Tuple of (validated_candidates, validation_results)
        """
        try:
            # Create mapping from reference names to candidates
            ref_to_candidate = {}
            ref_to_index = {}
            for i, candidate in enumerate(candidates):
                ref_name = f"ITD_{i+1}"
                ref_to_candidate[ref_name] = candidate
                ref_to_index[ref_name] = i
            
            # Track supporting reads for each candidate using index
            candidate_support = {}
            validation_results = []
            
            # Parse SAM file - try pysam first, fallback to manual parsing
            try:
                import pysam
                alignment_file = pysam.AlignmentFile(sam_file, "r")
                use_pysam = True
            except ImportError:
                logger.warning("pysam not available, using manual SAM parsing")
                use_pysam = False
                alignment_file = open(sam_file, 'r')
            
            total_reads_processed = 0
            filtered_reads = {'unmapped': 0, 'secondary': 0, 'large_indels': 0, 'low_quality': 0, 'soft_clipped': 0}
            
            try:
                if use_pysam:
                    # Use pysam for parsing
                    for read in alignment_file:
                        total_reads_processed += 1
                        
                        if read.is_unmapped or read.is_secondary or read.is_supplementary:
                            filtered_reads['unmapped' if read.is_unmapped else 'secondary'] += 1
                            continue
                        
                        ref_name = read.reference_name
                        
                        # Skip WT alignments for now (focus on ITD validation)
                        if ref_name == "WT":
                            continue
                        
                        if ref_name not in ref_to_candidate:
                            continue
                        
                        # Enhanced CIGAR filtering
                        if self._should_filter_read_by_cigar(read.cigartuples, ref_to_candidate[ref_name]):
                            filtered_reads['large_indels'] += 1
                            continue
                        
                        # Check mapping quality
                        if read.mapping_quality < 20:
                            filtered_reads['low_quality'] += 1
                            continue
                        
                        # Count as supporting read
                        candidate = ref_to_candidate[ref_name]
                        candidate_index = ref_to_index[ref_name]
                        if candidate_index not in candidate_support:
                            candidate_support[candidate_index] = []
                        candidate_support[candidate_index].append(read.query_name)
                else:
                    # Manual SAM parsing
                    for line in alignment_file:
                        if line.startswith('@'):
                            continue  # Skip header
                        
                        total_reads_processed += 1
                        fields = line.strip().split('\t')
                        if len(fields) < 11:
                            continue
                        
                        flag = int(fields[1])
                        ref_name = fields[2]
                        mapq = int(fields[4])
                        cigar = fields[5]
                        
                        # Check flags
                        if flag & 4:  # Unmapped
                            filtered_reads['unmapped'] += 1
                            continue
                        if flag & 256:  # Secondary
                            filtered_reads['secondary'] += 1
                            continue
                        
                        # Skip WT alignments
                        if ref_name == "WT":
                            continue
                        
                        if ref_name not in ref_to_candidate:
                            continue
                        
                        # Parse CIGAR and filter
                        cigar_tuples = self._parse_cigar_string(cigar)
                        if self._should_filter_read_by_cigar(cigar_tuples, ref_to_candidate[ref_name]):
                            filtered_reads['large_indels'] += 1
                            continue
                        
                        # Check mapping quality
                        if mapq < 20:
                            filtered_reads['low_quality'] += 1
                            continue
                        
                        # Count as supporting read
                        candidate = ref_to_candidate[ref_name]
                        candidate_index = ref_to_index[ref_name]
                        if candidate_index not in candidate_support:
                            candidate_support[candidate_index] = []
                        candidate_support[candidate_index].append(fields[0])  # Query name
            
            finally:
                alignment_file.close()
            
            logger.info(f"Processed {total_reads_processed} reads, filtered: {filtered_reads}")
            logger.info(f"Reads mapping to ITD candidates: {sum(len(reads) for reads in candidate_support.values())}")
            for i, candidate in enumerate(candidates):
                support_count = len(candidate_support.get(i, []))
                logger.info(f"  ITD_{i+1} ({candidate.length}bp): {support_count} supporting reads")
            
            # Filter candidates based on support AND allele frequency
            validated_candidates = []
            for i, candidate in enumerate(candidates):
                support_count = len(candidate_support.get(i, []))
                
                # Calculate allele frequency for this candidate
                allele_frequency = support_count / total_reads_processed if total_reads_processed > 0 else 0.0
                
                # Get minimum thresholds from config
                min_support = max(5, self.min_support)  # At least 5 reads or configured minimum
                min_af = getattr(self.config, 'min_allele_frequency', 0.01) if self.config else 0.01
                
                # Check both support and allele frequency criteria
                passes_support = support_count >= min_support
                passes_af = allele_frequency >= min_af
                
                if passes_support and passes_af:
                    validated_candidates.append(candidate)
                    validation_results.append({
                        'candidate': candidate,
                        'supporting_reads': support_count,
                        'total_reads_processed': total_reads_processed,  # Include total reads for AF calculation
                        'passed': True,
                        'reason': f'Sufficient support: {support_count} reads (>= {min_support}) and AF: {allele_frequency:.4f} (>= {min_af:.4f})'
                    })
                    logger.info(f"Validated ITD candidate: {candidate.length}bp at pos {candidate.position}, support={support_count}, AF={allele_frequency:.4f}")
                else:
                    # Determine reason for failure
                    if not passes_support and not passes_af:
                        reason = f'Insufficient support: {support_count} < {min_support} reads AND low AF: {allele_frequency:.4f} < {min_af:.4f}'
                    elif not passes_support:
                        reason = f'Insufficient support: {support_count} < {min_support} reads'
                    else:  # not passes_af
                        reason = f'Low allele frequency: {allele_frequency:.4f} < {min_af:.4f}'
                    
                    validation_results.append({
                        'candidate': candidate,
                        'supporting_reads': support_count,
                        'total_reads_processed': total_reads_processed,  # Include total reads for AF calculation
                        'passed': False,
                        'reason': reason
                    })
                    logger.info(f"Filtered ITD candidate: {candidate.length}bp at pos {candidate.position}, support={support_count}, AF={allele_frequency:.4f} - {reason}")
            
            return validated_candidates, validation_results
            
        except Exception as e:
            logger.error(f"Error parsing alignments: {e}")
            return candidates, []
    
    def _should_filter_read_by_cigar(self, cigar_tuples: List[Tuple[int, int]], 
                                   candidate: ITDCandidate) -> bool:
        """
        Enhanced CIGAR filtering to check if read should be excluded
        
        In multifasta alignment, reads should map to their best matching reference.
        We simply filter reads with large indels (>10bp) which indicate poor alignment
        or sequencing artifacts, regardless of which reference they mapped to.
        
        Args:
            cigar_tuples: List of (operation, length) tuples from pysam or parsed manually
            candidate: ITD candidate the read is aligned to (not used in filtering logic)
            
        Returns:
            True if read should be filtered out, False otherwise
        """
        if not cigar_tuples:
            return True
        
        total_soft_clipped = 0
        large_indel_threshold = 10  # Filter reads with any indel > 10bp (nanopore noise threshold)
        
        for operation, length in cigar_tuples:
            # CIGAR operations: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7=X, 8=
            if operation == 1:  # Insertion
                if length > large_indel_threshold:
                    logger.debug(f"Filtering read with large insertion: {length}bp > {large_indel_threshold}bp")
                    return True
                    
            elif operation == 2:  # Deletion
                if length > large_indel_threshold:
                    logger.debug(f"Filtering read with large deletion: {length}bp > {large_indel_threshold}bp")
                    return True
                    
            elif operation == 4:  # Soft clipping
                total_soft_clipped += length
        
        # Filter reads with excessive soft clipping (might indicate poor alignment)
        if total_soft_clipped > 50:  # More than 50bp soft clipped
            logger.debug(f"Filtering read with excessive soft clipping: {total_soft_clipped}bp")
            return True
        
        return False
    
    def _parse_cigar_string(self, cigar: str) -> List[Tuple[int, int]]:
        """
        Parse CIGAR string into list of (operation, length) tuples
        
        Args:
            cigar: CIGAR string (e.g., "100M50I86M")
            
        Returns:
            List of (operation, length) tuples
        """
        if not cigar or cigar == "*":
            return []
        
        import re
        # Parse CIGAR operations
        operations = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
        
        cigar_tuples = []
        for match in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
            length = int(match.group(1))
            op_char = match.group(2)
            operation = operations.get(op_char, 0)
            cigar_tuples.append((operation, length))
        
        return cigar_tuples
    
    def _longest_common_subsequence(self, seq1: str, seq2: str) -> int:
        """
        Calculate the length of the longest common subsequence
        """
        m, n = len(seq1), len(seq2)
        # Create DP table
        dp = [[0] * (n + 1) for _ in range(m + 1)]
        
        # Fill the DP table
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if seq1[i-1] == seq2[j-1]:
                    dp[i][j] = dp[i-1][j-1] + 1
                else:
                    dp[i][j] = max(dp[i-1][j], dp[i][j-1])
        
        return dp[m][n]

    def _kmer_similarity(self, seq1: str, seq2: str, k: int = 8) -> float:
        """
        Calculate k-mer based similarity between sequences
        """
        if len(seq1) < k or len(seq2) < k:
            return 0.0
        
        kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
        kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
        
        if not kmers1 or not kmers2:
            return 0.0
        
        intersection = len(kmers1 & kmers2)
        union = len(kmers1 | kmers2)
        
        return intersection / union if union > 0 else 0.0

    def _write_multifasta_for_validation(self, candidates: List[ITDCandidate], 
                                       output_prefix: Optional[str] = None) -> str:
        """
        Write multifasta file with WT sequence and all ITD candidates
        
        Args:
            candidates: List of ITD candidates
            output_prefix: Optional prefix for output files
            
        Returns:
            Path to created multifasta file
        """
        try:
            # Generate output filename
            if output_prefix:
                multifasta_file = f"{output_prefix}_multifasta_validation.fa"
            else:
                multifasta_file = "multifasta_validation.fa"
            
            with open(multifasta_file, 'w') as f:
                # Write WT reference sequence
                f.write(">WT\n")
                f.write(f"{self.reference_sequence}\n")
                
                # Write each ITD candidate sequence
                for i, candidate in enumerate(candidates):
                    # Create ITD reference by inserting the ITD sequence at the specified position
                    itd_ref = (self.reference_sequence[:candidate.position] + 
                              candidate.sequence + 
                              self.reference_sequence[candidate.position:])
                    
                    f.write(f">ITD_{i+1}\n")
                    f.write(f"{itd_ref}\n")
            
            logger.info(f"Created multifasta reference: {multifasta_file} with WT + {len(candidates)} ITD sequences")
            return multifasta_file
            
        except Exception as e:
            logger.error(f"Failed to create multifasta file: {e}")
            raise
    
    def _create_consensus_candidate(self, candidates: List[ITDCandidate]) -> ITDCandidate:
        """
        Create consensus ITD candidate from similar candidates
        
        Args:
            candidates: List of similar ITD candidates
            
        Returns:
            Consensus ITD candidate
        """
        # Helper function to get supporting reads count
        def get_support_count(candidate):
            if isinstance(candidate.supporting_reads, list):
                return len(candidate.supporting_reads)
            elif isinstance(candidate.supporting_reads, int):
                return candidate.supporting_reads
            else:
                return 0

        # Build weighted pools for consensus
        sequence_weights = []
        position_weights = []
        all_supporting_reads = []
        max_confidence = 0.0

        for c in candidates:
            support = get_support_count(c)
            sequence_weights.extend([c.sequence] * support)
            position_weights.extend([c.position] * support)
            if isinstance(c.supporting_reads, list):
                all_supporting_reads.extend(c.supporting_reads)
            elif isinstance(c.supporting_reads, int):
                all_supporting_reads.extend([None] * c.supporting_reads)
            max_confidence = max(max_confidence, c.confidence)

        # Consensus sequence: most supported sequence by read count
        from collections import Counter
        sequence_counts = Counter(sequence_weights)
        consensus_sequence = sequence_counts.most_common(1)[0][0]

        # Consensus position: median of all supporting reads' positions
        import numpy as np
        consensus_position = int(np.median(position_weights)) if position_weights else candidates[0].position

        # Create consensus candidate
        consensus_candidate = ITDCandidate(
            sequence=consensus_sequence,
            length=len(consensus_sequence),
            position=consensus_position,
            supporting_reads=all_supporting_reads,
            confidence=max_confidence,
            support_type=f"consensus_of_{len(candidates)}"
        )

        return consensus_candidate
    
    def _left_align_candidates(self, candidates: List[ITDCandidate]) -> List[ITDCandidate]:
        """
        Left-align ITD candidates to normalize alignment ambiguity
        
        This function shifts ITD positions as far left as possible while maintaining
        the same sequence. This helps consolidate ITDs that represent the same
        biological event but appear at different positions due to alignment ambiguity.
        
        Args:
            candidates: List of ITD candidates to left-align
            
        Returns:
            List of left-aligned ITD candidates
        """
        aligned_candidates = []
        
        for candidate in candidates:
            try:
                # Left-align the ITD position
                aligned_position = self._left_align_position(
                    position=candidate.position,
                    sequence=candidate.sequence,
                    reference=self.reference_sequence
                )
                
                # Create new candidate with aligned position
                aligned_candidate = ITDCandidate(
                    sequence=candidate.sequence,
                    length=candidate.length,
                    position=aligned_position,
                    support_type=candidate.support_type,
                    supporting_reads=candidate.supporting_reads,
                    confidence=candidate.confidence,
                    insertion_site=aligned_position,
                    duplication_start=getattr(candidate, 'duplication_start', None),
                    duplication_end=getattr(candidate, 'duplication_end', None),
                    is_primary=getattr(candidate, 'is_primary', False)
                )
                
                aligned_candidates.append(aligned_candidate)
                
                if aligned_position != candidate.position:
                    logger.debug(f"Left-aligned ITD: {candidate.position} → {aligned_position} "
                               f"({candidate.length}bp)")
                
            except Exception as e:
                logger.warning(f"Failed to left-align candidate at position {candidate.position}: {e}")
                # Keep original candidate if left-alignment fails
                aligned_candidates.append(candidate)
        
        return aligned_candidates
    
    def _left_align_position(self, position: int, sequence: str, reference: str) -> int:
        """
        Left-align a specific ITD position
        
        Args:
            position: Original insertion position
            sequence: ITD sequence
            reference: Reference sequence
            
        Returns:
            Left-aligned position
        """
        if not sequence or position < 0 or position >= len(reference):
            return position
        
        # Start from the original position and shift left
        current_pos = position
        
        # Keep shifting left while the sequence can still be "inserted" at that position
        # This means checking if inserting the sequence there would create the same
        # local context as inserting it at the original position
        while current_pos > 0:
            # Check if we can shift one position to the left
            if self._can_shift_left(current_pos, sequence, reference):
                current_pos -= 1
            else:
                break
        
        return current_pos
    
    def _can_shift_left(self, position: int, sequence: str, reference: str) -> bool:
        """
        Check if an ITD can be shifted one position to the left
        
        This uses a more sophisticated approach that checks if shifting left
        would still represent the same biological duplication event.
        
        Args:
            position: Current insertion position
            sequence: ITD sequence
            reference: Reference sequence
            
        Returns:
            True if the ITD can be shifted left
        """
        if position <= 0 or not sequence:
            return False
        
        # Method 1: Simple base matching
        left_base = reference[position - 1]
        last_itd_base = sequence[-1]
        
        if left_base == last_itd_base:
            return True
        
        # Method 2: Check for sequence overlap
        # If the ITD sequence starts with bases that match the reference
        # immediately to the left of the insertion, we can shift left
        max_check = min(len(sequence), position, 10)  # Check up to 10bp
        
        for i in range(1, max_check + 1):
            # Check if the first i bases of the ITD match the i bases to the left
            ref_left = reference[position - i:position]
            itd_start = sequence[:i]
            
            if ref_left == itd_start:
                return True
        
        return False
    
    def _mark_primary_candidates(self, candidates: List[ITDCandidate]):
        """
        Mark candidates as primary if their length matches size analysis results
        Uses dynamic tolerance that scales with ITD size for better accuracy
        
        Args:
            candidates: List of ITD candidates to mark
        """
        if not self.size_analysis_results or not candidates:
            return
        
        # Get detected ITD sizes from size analysis
        detected_sizes = self.size_analysis_results.get('detected_itd_sizes', [])
        
        # If no detected_itd_sizes key, try to extract from other fields
        if not detected_sizes:
            detected_sizes = self._extract_itd_sizes_from_analysis()
        
        if not detected_sizes:
            logger.info("No ITD sizes detected in size analysis")
            return
        
        logger.info(f"Size analysis detected ITDs: {detected_sizes}")
        
        # Mark candidates that match detected sizes with dynamic tolerance
        marked_count = 0
        for candidate in candidates:
            # Use actual sequence length for comparison
            actual_length = len(candidate.sequence)
            
            for detected_size in detected_sizes:
                # Calculate dynamic tolerance based on ITD size
                # Larger ITDs have proportionally larger tolerance
                # Base tolerance: 5bp for small ITDs
                # Additional tolerance: 5% of detected size for larger ITDs
                base_tolerance = 5
                size_based_tolerance = max(3, int(detected_size * 0.05))  # Min 3bp, 5% of size
                dynamic_tolerance = base_tolerance + size_based_tolerance
                
                if abs(actual_length - detected_size) <= dynamic_tolerance:
                    candidate.is_primary = True
                    marked_count += 1
                    if actual_length != candidate.length:
                        logger.info(f"Marked candidate as primary: {actual_length}bp sequence "
                                    f"(candidate.length={candidate.length}, matches size analysis: {detected_size}bp, "
                                    f"tolerance: ±{dynamic_tolerance}bp)")
                    else:
                        logger.info(f"Marked candidate as primary: {candidate.length}bp sequence "
                                    f"(matches size analysis: {detected_size}bp, tolerance: ±{dynamic_tolerance}bp)")
                    break  # Only mark once per candidate
        
        if marked_count == 0:
            logger.info("No candidates matched size analysis results")
    
    def _extract_itd_sizes_from_analysis(self) -> List[int]:
        """
        Extract ITD sizes from size analysis results using ReadSizeAnalyzer structure
        
        Returns:
            List of detected ITD sizes in bp
        """
        detected_sizes = []
        
        # Method 1: Extract from peaks with ITD classification
        if 'peaks' in self.size_analysis_results and isinstance(self.size_analysis_results['peaks'], dict):
            peaks = self.size_analysis_results['peaks'].get('peaks', [])
            if isinstance(peaks, list):
                for peak in peaks:
                    if isinstance(peak, dict):
                        peak_type = peak.get('type', '').lower()
                        
                        # Look for "potential_itd_XXbp" pattern
                        if 'potential_itd_' in peak_type:
                            import re
                            match = re.search(r'potential_itd_(\d+)bp', peak_type)
                            if match:
                                detected_sizes.append(int(match.group(1)))
        
        # Method 2: Extract from ITD signature statistics
        if 'itd_signatures' in self.size_analysis_results:
            itd_sigs = self.size_analysis_results['itd_signatures']
            if isinstance(itd_sigs, dict) and 'itd_size_stats' in itd_sigs:
                itd_stats = itd_sigs['itd_size_stats']
                if isinstance(itd_stats, dict) and itd_stats.get('count', 0) > 0:
                    # Add mean and median ITD sizes if available
                    if 'mean_size' in itd_stats:
                        detected_sizes.append(int(round(itd_stats['mean_size'])))
                    if 'median_size' in itd_stats:
                        detected_sizes.append(int(round(itd_stats['median_size'])))
        
        # Method 3: Extract from basic stats if available
        if 'basic_stats' in self.size_analysis_results:
            stats = self.size_analysis_results['basic_stats']
            if isinstance(stats, dict):
                # If there are any custom ITD size fields
                if 'mean_itd_size' in stats:
                    detected_sizes.append(int(round(stats['mean_itd_size'])))
                if 'median_itd_size' in stats:
                    detected_sizes.append(int(round(stats['median_itd_size'])))
        
        # Method 4: Parse from report text if available
        if 'report_text' in self.size_analysis_results:
            import re
            text = self.size_analysis_results['report_text']
            
            # Look for patterns like "potential_ITD_61bp" in the report text
            matches = re.findall(r'potential_ITD_(\d+)bp', text)
            for match in matches:
                detected_sizes.append(int(match))
            
            # Also look for "Mean ITD size: X bp" patterns
            mean_match = re.search(r'Mean ITD size:\s*(\d+(?:\.\d+)?)\s*bp', text)
            if mean_match:
                detected_sizes.append(int(float(mean_match.group(1))))
            
            # And "Median ITD size: X bp" patterns
            median_match = re.search(r'Median ITD size:\s*(\d+(?:\.\d+)?)\s*bp', text)
            if median_match:
                detected_sizes.append(int(float(median_match.group(1))))
        
        # Remove duplicates and return sorted list
        return sorted(list(set(detected_sizes)))
    
    def _validate_candidates(self, candidates: List[Union[ITDCandidate, SoftClipITDCandidate]], 
                           bam_file: str) -> Tuple[List[Union[ITDCandidate, SoftClipITDCandidate]], 
                                                  List[ValidationResult]]:
        """
        Validate candidates using multifasta approach with CIGAR filtering
        
        Args:
            candidates: List of all candidates to validate
            bam_file: Original BAM file for validation
            
        Returns:
            Tuple of (validated candidates, validation results)
        """
        if not candidates:
            logger.info("No candidates to validate")
            return [], []
        
        try:
            # Convert SoftClipITDCandidate to ITDCandidate if needed for consistency
            itd_candidates = []
            for candidate in candidates:
                if hasattr(candidate, 'sequence') and hasattr(candidate, 'position'):
                    if not isinstance(candidate, ITDCandidate):
                        # Convert SoftClipITDCandidate to ITDCandidate
                        itd_candidate = ITDCandidate(
                            sequence=candidate.sequence,
                            length=candidate.length,
                            position=candidate.position,
                            supporting_reads=candidate.supporting_reads,
                            confidence=getattr(candidate, 'confidence', 0.5),
                            support_type=getattr(candidate, 'support_type', 'softclip')
                        )
                        itd_candidates.append(itd_candidate)
                    else:
                        itd_candidates.append(candidate)
            
            # Use our new multifasta validation approach
            validated_candidates, validation_results = self._validate_candidates_dual_reference(
                itd_candidates, bam_file
            )
            
            # Convert validation results to ValidationResult format for compatibility
            validation_result_objects = []
            for result in validation_results:
                if result.get('passed', False):
                    supporting_reads = result.get('supporting_reads', 0)
                    total_reads = result.get('total_reads_processed', 30000)  # Get actual total from validation
                    
                    # Calculate proper allele frequency: supporting reads / total reads after CIGAR filtering
                    allele_frequency = supporting_reads / total_reads if total_reads > 0 else 0.0
                    
                    candidate = result.get('candidate')
                    logger.info(f"ITD {candidate.length}bp: AF = {supporting_reads}/{total_reads} = {allele_frequency:.4f}")
                    
                    # Create a proper ValidationResult object
                    mock_result = type('ValidationResult', (), {
                        'is_valid': True,
                        'allele_frequency': allele_frequency,
                        'segregation_quality': 0.8,  # Default value
                        'supporting_reads': supporting_reads,
                        'total_coverage': total_reads,  # Add missing total_coverage attribute
                        'itd_coverage': supporting_reads,  # ITD coverage = supporting reads for this ITD
                        'confidence': result.get('candidate').confidence if result.get('candidate') else 0.5,
                        'validation_confidence': result.get('candidate').confidence if result.get('candidate') else 0.5
                    })()
                    validation_result_objects.append(mock_result)
            
            logger.info(f"Multifasta validation completed: {len(validated_candidates)}/{len(candidates)} candidates passed")
            return validated_candidates, validation_result_objects
            
        except Exception as e:
            logger.error(f"Multifasta validation failed, falling back to old method: {e}")
            # Fallback to old validation method if multifasta fails
            reads_file = self.trimmed_fastq if self.trimmed_fastq else bam_file
            validation_results = validate_itd_candidates_dual_reference(candidates=candidates,
                reads_fastq=reads_file,
                reference_sequence=self.reference_sequence,
                min_supporting_reads=self.min_support,
                min_allele_frequency=0.01,
                config=self.config)
            
            validated_candidates = []
            for i, result in enumerate(validation_results):
                validated_candidates.append(candidates[i])
                logger.debug(f"Validated candidate {i+1}: {candidates[i].length}bp "
                           f"(allele_freq={result.allele_frequency:.3f})")
            
            return validated_candidates, validation_results
    
    def _generate_detection_summary(self, cigar_candidates: List[ITDCandidate],
                                  softclip_candidates: List[SoftClipITDCandidate],
                                  validated_candidates: List[Union[ITDCandidate, SoftClipITDCandidate]],
                                  validation_results: List[ValidationResult]) -> Dict:
        """
        Generate comprehensive detection summary
        
        Args:
            cigar_candidates: CIGAR detection results
            softclip_candidates: Soft-clip detection results
            validated_candidates: Final validated candidates
            validation_results: Validation results
            
        Returns:
            Detection summary dictionary
        """
        summary = {
            'total_candidates_found': len(cigar_candidates) + len(softclip_candidates),
            'cigar_candidates': len(cigar_candidates),
            'softclip_candidates': len(softclip_candidates),
            'validated_candidates': len(validated_candidates),
            'validation_rate': len(validated_candidates) / max(1, len(cigar_candidates) + len(softclip_candidates)),
            'detection_methods_used': [],
            'itd_characteristics': [],
            'validation_metrics': {}
        }
        
        # Track detection methods used
        if cigar_candidates:
            summary['detection_methods_used'].append('CIGAR_analysis')
        if softclip_candidates:
            summary['detection_methods_used'].append('soft_clip_analysis')
        
        # Analyze ITD characteristics
        if validated_candidates:
            lengths = [c.length for c in validated_candidates]
            positions = [c.position for c in validated_candidates]
            confidences = [c.confidence for c in validated_candidates]
            
            summary['itd_characteristics'] = {
                'count': len(validated_candidates),
                'length_range': f"{min(lengths)}-{max(lengths)}bp",
                'mean_length': sum(lengths) / len(lengths),
                'position_range': f"{min(positions)}-{max(positions)}",
                'mean_confidence': sum(confidences) / len(confidences)
            }
        
        # Validation metrics
        if validation_results:
            valid_results = [r for r in validation_results if r.is_valid]
            if valid_results:
                allele_freqs = [r.allele_frequency for r in valid_results]
                segregation_scores = [r.segregation_quality for r in valid_results]
                
                summary['validation_metrics'] = {
                    'mean_allele_frequency': sum(allele_freqs) / len(allele_freqs),
                    'mean_segregation_quality': sum(segregation_scores) / len(segregation_scores),
                    'validation_success_rate': len(valid_results) / len(validation_results)
                }
        
        return summary
    
    def _save_results(self, output_prefix: str, detection_summary: Dict, runtime_stats: Dict):
        """
        Save detection results to files
        
        Args:
            output_prefix: Prefix for output files
            detection_summary: Detection summary data
            runtime_stats: Runtime statistics
        """
        try:
            # Save summary as JSON
            summary_file = f"{output_prefix}_itd_detection_summary.json"
            combined_data = {
                'detection_summary': detection_summary,
                'runtime_stats': runtime_stats,
                'pipeline_config': {
                    'min_itd_length': self.min_itd_length,
                    'max_itd_length': self.max_itd_length,
                    'min_support': self.min_support,
                    'min_softclip_length': self.min_softclip_length,
                    'enable_softclip_fallback': self.enable_softclip_fallback
                }
            }
            
            with open(summary_file, 'w') as f:
                json.dump(combined_data, f, indent=2)
            
            logger.info(f"Saved detection summary to {summary_file}")
            
        except Exception as e:
            logger.warning(f"Failed to save results: {e}")


def run_flt3_itd_pipeline(bam_file: str, reference_sequence: str, reference_file: str,
                         trimmed_fastq: str = None,
                         min_itd_length: int = 15, max_itd_length: int = 500,
                         min_support: int = 3, min_softclip_length: int = 50,
                         enable_softclip_fallback: bool = True,
                         size_analysis_results: dict = None,
                         output_prefix: Optional[str] = None, config=None) -> ITDDetectionResult:
    """
    Main function to run the complete FLT3 ITD detection pipeline
    
    Args:
        bam_file: Path to BAM file with FLT3 reads
        reference_sequence: FLT3 reference sequence string
        reference_file: Path to reference FASTA file
        trimmed_fastq: Path to trimmed FASTQ file for re-alignment
        min_itd_length: Minimum ITD size to detect
        max_itd_length: Maximum ITD size to detect
        min_support: Minimum supporting reads required
        min_softclip_length: Minimum soft-clip length for secondary analysis
        enable_softclip_fallback: Whether to use soft-clip analysis if CIGAR detection is insufficient
        size_analysis_results: Results from read size analysis for primary marking
        output_prefix: Optional prefix for output files
        
    Returns:
        Comprehensive detection results
    """
    pipeline = FLT3ITDPipeline(
        reference_sequence=reference_sequence,
        reference_file=reference_file,
        trimmed_fastq=trimmed_fastq,
        min_itd_length=min_itd_length,
        max_itd_length=max_itd_length,
        min_support=min_support,
        min_softclip_length=min_softclip_length,
        enable_softclip_fallback=enable_softclip_fallback,
        size_analysis_results=size_analysis_results,
        config=config
    )
    
    return pipeline.detect_itds(bam_file, output_prefix)


if __name__ == "__main__":
    # Test the complete pipeline
    import argparse
    
    parser = argparse.ArgumentParser(description="FLT3 ITD Detection Pipeline")
    parser.add_argument("bam_file", help="Input BAM file with FLT3 reads")
    parser.add_argument("reference_file", help="Reference FASTA file")
    parser.add_argument("--min-itd-length", type=int, default=15, help="Minimum ITD length")
    parser.add_argument("--max-itd-length", type=int, default=500, help="Maximum ITD length")
    parser.add_argument("--min-support", type=int, default=3, help="Minimum supporting reads")
    parser.add_argument("--min-softclip", type=int, default=50, help="Minimum soft-clip length")
    parser.add_argument("--no-softclip-fallback", action="store_true", 
                       help="Disable soft-clip fallback analysis")
    parser.add_argument("--output-prefix", help="Output file prefix")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Read reference sequence
    try:
        with open(args.reference_file, 'r') as f:
            lines = f.readlines()
            reference_seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    except Exception as e:
        logger.error(f"Failed to read reference file: {e}")
        exit(1)
    
    # Run pipeline
    try:
        result = run_flt3_itd_pipeline(
            bam_file=args.bam_file,
            reference_sequence=reference_seq,
            reference_file=args.reference_file,
            min_itd_length=args.min_itd_length,
            max_itd_length=args.max_itd_length,
            min_support=args.min_support,
            min_softclip_length=args.min_softclip,
            enable_softclip_fallback=not args.no_softclip_fallback,
            output_prefix=args.output_prefix
        )
        
        # Print results
        print(f"\n=== FLT3 ITD Detection Results ===")
        print(f"CIGAR candidates: {len(result.cigar_candidates)}")
        print(f"Soft-clip candidates: {len(result.softclip_candidates)}")
        print(f"Validated ITDs: {len(result.validated_candidates)}")
        print(f"Runtime: {result.runtime_stats['total_time']:.2f}s")
        
        if result.validated_candidates:
            print(f"\nValidated ITDs:")
            for i, candidate in enumerate(result.validated_candidates, 1):
                print(f"  {i}. {candidate.length}bp at position {candidate.position}")
                print(f"     Support: {len(candidate.supporting_reads)} reads")
                print(f"     Confidence: {candidate.confidence:.3f}")
                print(f"     Method: {candidate.support_type}")
        else:
            print("\nNo validated ITDs found.")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        exit(1)
