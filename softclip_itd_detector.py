#!/usr/bin/env python3
"""
Soft-Clip ITD Detector Module
Secondary ITD detection from soft-clipped regions in reads
"""

import logging
import tempfile
import subprocess
import pysam
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Set
from collections import defaultdict, Counter
from dataclasses import dataclass
import shutil

logger = logging.getLogger(__name__)

@dataclass  
class SoftClipITDCandidate:
    """ITD candidate from soft-clip analysis"""
    sequence: str
    length: int
    position: int
    support_type: str = 'soft_clip'
    supporting_reads: List[str] = None
    confidence: float = 0.0
    insertion_site: Optional[int] = None
    overlap_evidence: Optional[Dict] = None  # Evidence of overlap between mapped and soft-clip regions

    def __post_init__(self):
        if self.supporting_reads is None:
            self.supporting_reads = []

class SoftClipITDDetector:
    """Detect ITDs from soft-clipped regions that realign to overlapping regions"""
    
    def __init__(self, reference_sequence: str, reference_file: str,
                 min_itd_length: int = 15, min_support: int = 3,
                 max_itd_length: int = 500, min_softclip_length: int = 50,
                 threads: int = 4):
        """
        Initialize soft-clip ITD detector
        
        Args:
            reference_sequence: FLT3 reference sequence string
            reference_file: Path to reference FASTA file
            min_itd_length: Minimum ITD size to consider
            min_support: Minimum supporting reads required
            max_itd_length: Maximum ITD size to consider
            min_softclip_length: Minimum soft-clip length to analyze
            threads: Number of threads for minimap2
        """
        self.reference_sequence = reference_sequence
        self.reference_file = reference_file
        self.min_itd_length = min_itd_length
        self.max_itd_length = max_itd_length
        self.min_support = min_support
        self.min_softclip_length = min_softclip_length
        self.threads = threads
        self.temp_dir = tempfile.mkdtemp(prefix="flt3_softclip_")
        
        logger.info(f"Initialized soft-clip ITD detector: min_softclip={min_softclip_length}bp, "
                   f"ITD_range={min_itd_length}-{max_itd_length}bp, threads={threads}")
    
    def detect_itds_from_bam(self, bam_file: str) -> List[SoftClipITDCandidate]:
        """
        Main function to detect ITDs from soft-clipped regions
        
        Args:
            bam_file: Path to BAM file with FLT3 reads
            
        Returns:
            List of soft-clip ITD candidates
        """
        logger.info(f"Starting soft-clip ITD detection from {bam_file}")
        
        # Step 1: Extract soft-clipped reads
        softclip_reads = self._extract_softclip_reads(bam_file)
        logger.info(f"Found {len(softclip_reads)} reads with significant soft-clips")
        
        if not softclip_reads:
            logger.info("No reads with significant soft-clips found")
            return []
        
        # Step 2: Extract and realign soft-clipped sequences
        realignment_results = self._realign_softclips(softclip_reads)
        logger.info(f"Realigned {len(realignment_results)} soft-clip sequences")
        
        # Step 3: Detect overlaps between original mapping and soft-clip realignment
        overlap_evidence = self._detect_overlaps(softclip_reads, realignment_results)
        logger.info(f"Found overlap evidence in {len(overlap_evidence)} reads")
        
        # Step 4: Group overlapping evidence and generate ITD candidates
        candidates = self._generate_itd_candidates(overlap_evidence)
        logger.info(f"Generated {len(candidates)} ITD candidates from soft-clip analysis")
        
        # Step 5: Validate candidates
        validated_candidates = self._validate_candidates(candidates)
        logger.info(f"Validated soft-clip candidates: {len(candidates)} → {len(validated_candidates)}")
        
        return validated_candidates
    
    def _extract_softclip_reads(self, bam_file: str) -> List[Dict]:
        """
        Extract reads with significant soft-clips
        
        Args:
            bam_file: Path to BAM file
            
        Returns:
            List of reads with soft-clip information
        """
        softclip_reads = []
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                # Skip unmapped, secondary, or supplementary reads
                if (read.is_unmapped or read.is_secondary or 
                    read.is_supplementary or read.is_duplicate):
                    continue
                
                # Skip reads without CIGAR or sequence
                if not read.cigartuples or not read.query_sequence:
                    continue
                
                # Extract soft-clip information
                softclips = self._parse_softclips(read)
                
                # Check if any soft-clip is significant
                significant_clips = [clip for clip in softclips 
                                   if clip['length'] >= self.min_softclip_length]
                
                if significant_clips:
                    softclip_reads.append({
                        'read_name': read.query_name,
                        'read_sequence': read.query_sequence,
                        'read_length': len(read.query_sequence),
                        'mapping_start': read.reference_start,
                        'mapping_end': read.reference_end,
                        'mapping_quality': read.mapping_quality,
                        'cigar': read.cigarstring,
                        'softclips': significant_clips,
                        'is_reverse': read.is_reverse
                    })
        
        logger.debug(f"Extracted {len(softclip_reads)} reads with soft-clips ≥{self.min_softclip_length}bp")
        return softclip_reads
    
    def _parse_softclips(self, read: pysam.AlignedRead) -> List[Dict]:
        """
        Parse soft-clips from a read's CIGAR string
        
        Args:
            read: Aligned read
            
        Returns:
            List of soft-clip dictionaries
        """
        softclips = []
        query_pos = 0
        
        for i, (operation, length) in enumerate(read.cigartuples):
            if operation == 4:  # Soft clip
                clip_sequence = read.query_sequence[query_pos:query_pos + length]
                
                # Determine position relative to mapping
                if i == 0:  # Left soft-clip
                    clip_type = 'left'
                    ref_position = read.reference_start  # Where it would map
                elif i == len(read.cigartuples) - 1:  # Right soft-clip
                    clip_type = 'right'
                    ref_position = read.reference_end  # Where it would map
                else:
                    clip_type = 'internal'  # Rare, but possible
                    ref_position = read.reference_start  # Approximate
                
                softclips.append({
                    'type': clip_type,
                    'sequence': clip_sequence,
                    'length': length,
                    'query_start': query_pos,
                    'query_end': query_pos + length,
                    'ref_position': ref_position
                })
            
            if operation in [0, 1, 4]:  # Operations that advance query
                query_pos += length
        
        return softclips
    
    def _realign_softclips(self, softclip_reads: List[Dict]) -> Dict[str, List[Dict]]:
        """
        Realign soft-clipped sequences to reference
        
        Args:
            softclip_reads: List of reads with soft-clips
            
        Returns:
            Dictionary mapping read names to realignment results
        """
        # Create FASTQ file with soft-clipped sequences
        softclip_fastq = self._create_softclip_fastq(softclip_reads)
        
        if not softclip_fastq:
            return {}
        
        # Realign soft-clips to reference
        realignment_sam = self._align_softclips(softclip_fastq)
        
        if not realignment_sam:
            return {}
        
        # Parse realignment results
        realignment_results = self._parse_realignment_results(realignment_sam)
        
        # Clean up temporary files
        if Path(softclip_fastq).exists():
            Path(softclip_fastq).unlink()
        if Path(realignment_sam).exists():
            Path(realignment_sam).unlink()
        
        return realignment_results
    
    def _create_softclip_fastq(self, softclip_reads: List[Dict]) -> Optional[str]:
        """
        Create FASTQ file with soft-clipped sequences for realignment
        
        Args:
            softclip_reads: List of reads with soft-clips
            
        Returns:
            Path to FASTQ file or None if failed
        """
        fastq_file = Path(self.temp_dir) / "softclips.fastq"
        
        try:
            with open(fastq_file, 'w') as f:
                for read in softclip_reads:
                    for i, clip in enumerate(read['softclips']):
                        # Create unique identifier for each soft-clip
                        clip_id = f"{read['read_name']}_clip_{clip['type']}_{i}"
                        
                        f.write(f"@{clip_id}\n")
                        f.write(f"{clip['sequence']}\n")
                        f.write("+\n")
                        f.write("I" * clip['length'] + "\n")  # Dummy quality scores
            
            return str(fastq_file)
            
        except Exception as e:
            logger.error(f"Failed to create soft-clip FASTQ: {e}")
            return None
    
    def _align_softclips(self, softclip_fastq: str) -> Optional[str]:
        """
        Align soft-clipped sequences to reference
        
        Args:
            softclip_fastq: Path to FASTQ with soft-clips
            
        Returns:
            Path to alignment SAM file or None if failed
        """
        output_sam = Path(self.temp_dir) / "softclip_realignment.sam"
        
        # Use minimap2 with sensitive settings for short sequences
        cmd = [
            'minimap2',
            '-ax', 'sr',  # Short read mode for better sensitivity
            '-t', str(self.threads),
            '--secondary=yes',  # Allow secondary alignments for multi-mapping
            '-k', '10',  # Smaller k-mer for short sequences
            '-w', '5',   # Smaller window
            self.reference_file,
            softclip_fastq
        ]
        
        try:
            with open(output_sam, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
                
            if result.returncode != 0:
                logger.warning(f"Soft-clip realignment failed: {result.stderr}")
                return None
            
            return str(output_sam)
            
        except Exception as e:
            logger.error(f"Soft-clip alignment failed: {e}")
            return None
    
    def _parse_realignment_results(self, realignment_sam: str) -> Dict[str, List[Dict]]:
        """
        Parse realignment results and group by original read
        
        Args:
            realignment_sam: Path to realignment SAM file
            
        Returns:
            Dictionary mapping read names to realignment results
        """
        realignment_results = defaultdict(list)
        
        try:
            with pysam.AlignmentFile(realignment_sam, "r") as sam:
                for alignment in sam:
                    if alignment.is_unmapped:
                        continue
                    
                    # Extract original read name from clip ID
                    clip_id = alignment.query_name
                    read_name = clip_id.split('_clip_')[0]
                    clip_info = clip_id.split('_clip_')[1] if '_clip_' in clip_id else 'unknown'
                    
                    realignment_results[read_name].append({
                        'clip_id': clip_id,
                        'clip_info': clip_info,
                        'ref_start': alignment.reference_start,
                        'ref_end': alignment.reference_end,
                        'mapping_quality': alignment.mapping_quality,
                        'cigar': alignment.cigarstring,
                        'sequence': alignment.query_sequence,
                        'is_reverse': alignment.is_reverse
                    })
        
        except Exception as e:
            logger.error(f"Failed to parse realignment results: {e}")
            
        return dict(realignment_results)
    
    def _detect_overlaps(self, softclip_reads: List[Dict], 
                        realignment_results: Dict[str, List[Dict]]) -> List[Dict]:
        """
        Detect overlaps between original mapping and soft-clip realignments
        
        Args:
            softclip_reads: Original reads with soft-clips
            realignment_results: Realignment results for soft-clips
            
        Returns:
            List of overlap evidence
        """
        overlap_evidence = []
        
        for read in softclip_reads:
            read_name = read['read_name']
            
            if read_name not in realignment_results:
                continue
            
            # Check each realigned soft-clip for overlap with original mapping
            for realignment in realignment_results[read_name]:
                overlap_info = self._calculate_overlap(read, realignment)
                
                if overlap_info and overlap_info['overlap_length'] > 0:
                    # This suggests an ITD: the soft-clip realigns to a region 
                    # that overlaps with the original mapping
                    evidence = {
                        'read_name': read_name,
                        'original_mapping': {
                            'start': read['mapping_start'],
                            'end': read['mapping_end']
                        },
                        'softclip_realignment': {
                            'start': realignment['ref_start'],
                            'end': realignment['ref_end'],
                            'sequence': realignment['sequence'],
                            'clip_info': realignment['clip_info']
                        },
                        'overlap_info': overlap_info,
                        'full_read': read
                    }
                    
                    overlap_evidence.append(evidence)
                    
                    logger.debug(f"Found overlap evidence in {read_name}: "
                                f"original={read['mapping_start']}-{read['mapping_end']}, "
                                f"softclip={realignment['ref_start']}-{realignment['ref_end']}, "
                                f"overlap={overlap_info['overlap_length']}bp")
        
        return overlap_evidence
    
    def _calculate_overlap(self, read: Dict, realignment: Dict) -> Optional[Dict]:
        """
        Calculate overlap between original mapping and soft-clip realignment
        
        Args:
            read: Original read information
            realignment: Soft-clip realignment information
            
        Returns:
            Overlap information dictionary or None
        """
        orig_start = read['mapping_start']
        orig_end = read['mapping_end']
        realign_start = realignment['ref_start']
        realign_end = realignment['ref_end']
        
        # Calculate overlap
        overlap_start = max(orig_start, realign_start)
        overlap_end = min(orig_end, realign_end)
        overlap_length = max(0, overlap_end - overlap_start)
        
        if overlap_length > 0:
            # Calculate potential ITD characteristics
            # The ITD length would be the non-overlapping part of the soft-clip
            if 'left' in realignment['clip_info']:
                # Left soft-clip: ITD = part before original mapping
                itd_start = realign_start
                itd_end = orig_start
                itd_insertion_point = orig_start
            elif 'right' in realignment['clip_info']:
                # Right soft-clip: ITD = part after original mapping  
                itd_start = orig_end
                itd_end = realign_end
                itd_insertion_point = orig_end
            else:
                # Default case
                itd_start = min(orig_start, realign_start)
                itd_end = max(orig_end, realign_end)
                itd_insertion_point = overlap_start
            
            itd_length = max(0, itd_end - itd_start - overlap_length)
            
            return {
                'overlap_start': overlap_start,
                'overlap_end': overlap_end,
                'overlap_length': overlap_length,
                'itd_start': itd_start,
                'itd_end': itd_end,
                'itd_length': itd_length,
                'itd_insertion_point': itd_insertion_point,
                'itd_sequence': realignment['sequence']  # Approximate
            }
        
        return None
    
    def _generate_itd_candidates(self, overlap_evidence: List[Dict]) -> List[SoftClipITDCandidate]:
        """
        Generate ITD candidates from overlap evidence
        
        Args:
            overlap_evidence: List of overlap evidence
            
        Returns:
            List of ITD candidates
        """
        # Group evidence by similar ITD characteristics
        itd_groups = self._group_overlap_evidence(overlap_evidence)
        
        candidates = []
        
        for group_id, evidence_list in itd_groups.items():
            if len(evidence_list) < self.min_support:
                continue
            
            # Generate consensus ITD from group
            consensus = self._generate_consensus_itd(evidence_list)
            
            if (consensus and 
                self.min_itd_length <= consensus['length'] <= self.max_itd_length):
                
                candidate = SoftClipITDCandidate(
                    sequence=consensus['sequence'],
                    length=consensus['length'],
                    position=consensus['insertion_point'],
                    support_type='soft_clip',
                    supporting_reads=[ev['read_name'] for ev in evidence_list],
                    confidence=consensus['confidence'],
                    insertion_site=consensus['insertion_point'],
                    overlap_evidence=evidence_list
                )
                
                candidates.append(candidate)
                
                logger.debug(f"Generated soft-clip ITD candidate: {candidate.length}bp "
                            f"at pos {candidate.position}, support={len(candidate.supporting_reads)}")
        
        return candidates
    
    def _group_overlap_evidence(self, overlap_evidence: List[Dict]) -> Dict[str, List[Dict]]:
        """
        Group overlap evidence by similar ITD characteristics
        
        Args:
            overlap_evidence: List of overlap evidence
            
        Returns:
            Dictionary mapping group IDs to evidence lists
        """
        groups = defaultdict(list)
        
        for evidence in overlap_evidence:
            overlap_info = evidence['overlap_info']
            
            # Create grouping key based on ITD characteristics
            insertion_point = overlap_info['itd_insertion_point']
            itd_length = overlap_info['itd_length']
            
            # Group by position (±20bp) and length (±10bp)
            position_bin = insertion_point // 20
            length_bin = itd_length // 10
            
            group_key = f"pos_{position_bin}_len_{length_bin}"
            groups[group_key].append(evidence)
        
        return dict(groups)
    
    def _generate_consensus_itd(self, evidence_list: List[Dict]) -> Optional[Dict]:
        """
        Generate consensus ITD from grouped evidence
        
        Args:
            evidence_list: List of evidence for similar ITDs
            
        Returns:
            Consensus ITD information or None
        """
        if not evidence_list:
            return None
        
        # Extract ITD characteristics
        insertion_points = []
        itd_lengths = []
        sequences = []
        
        for evidence in evidence_list:
            overlap_info = evidence['overlap_info']
            insertion_points.append(overlap_info['itd_insertion_point'])
            itd_lengths.append(overlap_info['itd_length'])
            sequences.append(overlap_info['itd_sequence'])
        
        # Calculate consensus values
        import numpy as np
        consensus_insertion_point = int(np.median(insertion_points))
        consensus_length = int(np.median(itd_lengths))
        
        # Most common sequence
        sequence_counts = Counter(sequences)
        consensus_sequence = sequence_counts.most_common(1)[0][0]
        
        # Calculate confidence
        support_count = len(evidence_list)
        position_consistency = 1.0 - (np.std(insertion_points) / 20.0) if len(insertion_points) > 1 else 1.0
        length_consistency = 1.0 - (np.std(itd_lengths) / 10.0) if len(itd_lengths) > 1 else 1.0
        
        base_confidence = 0.6  # Base confidence for soft-clip detection
        support_bonus = min(0.2, support_count * 0.03)
        consistency_bonus = (position_consistency + length_consistency) * 0.1
        
        confidence = min(0.9, base_confidence + support_bonus + consistency_bonus)
        
        return {
            'insertion_point': consensus_insertion_point,
            'length': consensus_length,
            'sequence': consensus_sequence,
            'confidence': confidence,
            'support_count': support_count
        }
    
    def _validate_candidates(self, candidates: List[SoftClipITDCandidate]) -> List[SoftClipITDCandidate]:
        """
        Validate soft-clip ITD candidates
        
        Args:
            candidates: List of candidates to validate
            
        Returns:
            List of validated candidates
        """
        validated = []
        
        for candidate in candidates:
            if self._validate_single_candidate(candidate):
                validated.append(candidate)
                logger.debug(f"Validated soft-clip candidate: {candidate.length}bp")
            else:
                logger.debug(f"Rejected soft-clip candidate: {candidate.length}bp")
        
        return validated
    
    def _validate_single_candidate(self, candidate: SoftClipITDCandidate) -> bool:
        """
        Validate a single soft-clip ITD candidate
        
        Args:
            candidate: Candidate to validate
            
        Returns:
            True if candidate is valid
        """
        # Check sequence complexity
        if len(set(candidate.sequence)) < 3:
            return False
        
        # Check for excessive homopolymers
        max_homopolymer = max(len(list(g)) for k, g in 
                             __import__('itertools').groupby(candidate.sequence))
        if max_homopolymer > min(8, len(candidate.sequence) * 0.4):
            return False
        
        # Check support
        if len(candidate.supporting_reads) < self.min_support:
            return False
        
        # Check confidence
        if candidate.confidence < 0.5:
            return False
        
        return True
    
    def cleanup(self):
        """Clean up temporary files"""
        try:
            if Path(self.temp_dir).exists():
                shutil.rmtree(self.temp_dir)
                logger.debug(f"Cleaned up soft-clip temp directory: {self.temp_dir}")
        except Exception as e:
            logger.warning(f"Failed to clean up soft-clip temp directory: {e}")


def detect_softclip_itds(bam_file: str, reference_sequence: str, reference_file: str,
                        min_itd_length: int = 15, min_support: int = 3,
                        max_itd_length: int = 500, 
                        min_softclip_length: int = 50, threads: int = 4) -> List[SoftClipITDCandidate]:
    """
    Main function to detect ITDs from soft-clipped regions
    
    Args:
        bam_file: Path to BAM file with FLT3 reads
        reference_sequence: FLT3 reference sequence string
        reference_file: Path to reference FASTA file
        min_itd_length: Minimum ITD size
        min_support: Minimum supporting reads
        max_itd_length: Maximum ITD size
        min_softclip_length: Minimum soft-clip length to analyze
        
    Returns:
        List of soft-clip ITD candidates
    """
    detector = SoftClipITDDetector(
        reference_sequence=reference_sequence,
        reference_file=reference_file,
        min_itd_length=min_itd_length,
        min_support=min_support,
        max_itd_length=max_itd_length,
        min_softclip_length=min_softclip_length,
        threads=threads
    )
    
    try:
        return detector.detect_itds_from_bam(bam_file)
    finally:
        detector.cleanup()


if __name__ == "__main__":
    print("This module is part of the FLT3 ITD detection pipeline.")
    print("Please use main_module.py for command-line execution.")
    print("This module should be imported and used programmatically.")
