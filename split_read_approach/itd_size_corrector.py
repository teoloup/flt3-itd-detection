#!/usr/bin/env python3
"""
ITD Size Corrector Module
Correctly estimates ITD sizes by separating actual duplication from flanking sequence
"""

import logging
import pysam
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)

class ITDSizeCorrector:
    """Corrects ITD size estimation by properly analyzing read alignment patterns"""
    
    def __init__(self, reference_sequence: str):
        self.reference_sequence = reference_sequence
        
    def correct_itd_size(self, r1: pysam.AlignedRead, r2: pysam.AlignedRead, 
                        overlap: Dict) -> Optional[Dict]:
        """
        Correctly estimate ITD size from overlapping read pairs
        
        Returns:
            Dict with corrected ITD information:
            - actual_itd_length: The true ITD size
            - itd_sequence: The actual duplicated sequence
            - insertion_point: Where the ITD is inserted
            - confidence: Estimation confidence
        """
        
        # Method 1: Analyze soft-clips to identify ITD boundaries
        soft_clip_result = self._analyze_soft_clips_for_itd_size(r1, r2, overlap)
        if soft_clip_result:
            return soft_clip_result
            
        # Method 2: Analyze CIGAR insertions 
        cigar_result = self._analyze_cigar_insertions_for_itd_size(r1, r2, overlap)
        if cigar_result:
            return cigar_result
            
        # Method 3: Reference comparison method
        ref_result = self._analyze_reference_comparison_for_itd_size(r1, r2, overlap)
        if ref_result:
            return ref_result
            
        return None
    
    def _analyze_soft_clips_for_itd_size(self, r1: pysam.AlignedRead, r2: pysam.AlignedRead, 
                                        overlap: Dict) -> Optional[Dict]:
        """Estimate ITD size by analyzing soft-clipped sequences"""
        
        # Get soft-clipped sequences from both reads
        r1_clips = self._extract_soft_clips(r1)
        r2_clips = self._extract_soft_clips(r2)
        
        # Look for soft-clips that match reference sequence (indicating duplication)
        for read_name, clips in [("r1", r1_clips), ("r2", r2_clips)]:
            for clip_type, clip_info in clips.items():
                if not clip_info or len(clip_info['sequence']) < 15:
                    continue
                    
                # Find where this soft-clip matches in the reference
                ref_matches = self._find_sequence_in_reference(clip_info['sequence'])
                
                for match in ref_matches:
                    # Check if this match indicates a duplication
                    if self._is_likely_duplication(clip_info, match, overlap):
                        itd_size = len(clip_info['sequence'])
                        
                        # Validate that this size makes biological sense
                        if 15 <= itd_size <= 300:  # Reasonable ITD size range
                            logger.info(f"ITD size estimated from {read_name} {clip_type} soft-clip: {itd_size}bp")
                            
                            return {
                                'actual_itd_length': itd_size,
                                'itd_sequence': clip_info['sequence'],
                                'insertion_point': match['position'],
                                'confidence': 0.8,
                                'method': f'soft_clip_{read_name}_{clip_type}'
                            }
        
        return None
    
    def _analyze_cigar_insertions_for_itd_size(self, r1: pysam.AlignedRead, r2: pysam.AlignedRead, 
                                              overlap: Dict) -> Optional[Dict]:
        """Estimate ITD size by analyzing large insertions in CIGAR strings"""
        
        for read in [r1, r2]:
            if not read.cigartuples:
                continue
                
            ref_pos = read.reference_start
            query_pos = 0
            
            for op, length in read.cigartuples:
                if op == 1 and length >= 15:  # Large insertion
                    # Extract the inserted sequence
                    insertion_seq = read.query_sequence[query_pos:query_pos + length]
                    
                    # Check if this insertion matches reference sequence (duplication)
                    ref_matches = self._find_sequence_in_reference(insertion_seq)
                    
                    for match in ref_matches:
                        # Check if insertion position makes sense relative to match
                        if abs(ref_pos - match['position']) < 200:  # Within reasonable distance
                            logger.info(f"ITD size estimated from CIGAR insertion: {length}bp")
                            
                            return {
                                'actual_itd_length': length,
                                'itd_sequence': insertion_seq,
                                'insertion_point': ref_pos,
                                'confidence': 0.9,
                                'method': 'cigar_insertion'
                            }
                
                # Update positions
                if op in [0, 1, 4]:  # Operations that consume query
                    query_pos += length
                if op in [0, 2]:  # Operations that consume reference
                    ref_pos += length
        
        return None
    
    def _analyze_reference_comparison_for_itd_size(self, r1: pysam.AlignedRead, r2: pysam.AlignedRead, 
                                                  overlap: Dict) -> Optional[Dict]:
        """
        Estimate ITD size by comparing read sequences to reference
        This method identifies the actual duplicated portion
        """
        
        # Extract sequences from the overlap region
        overlap_start = overlap['overlap_start']
        overlap_end = overlap['overlap_end']
        overlap_length = overlap_end - overlap_start
        
        # Get reference sequence for the overlap region
        ref_seq = self.reference_sequence[overlap_start:overlap_end]
        
        # Extract read sequences that align to this region
        r1_seq = self._extract_sequence_for_region(r1, overlap_start, overlap_end)
        r2_seq = self._extract_sequence_for_region(r2, overlap_start, overlap_end)
        
        if not r1_seq or not r2_seq:
            return None
        
        # Method: Find the longest common subsequence that doesn't match reference
        # This represents the duplicated (ITD) portion
        
        # Compare read sequences to find consensus ITD
        consensus_itd = self._find_consensus_itd_sequence(r1_seq, r2_seq, ref_seq)
        
        if consensus_itd and 15 <= len(consensus_itd) <= 300:
            # Find where this sequence exists in the reference (original location)
            ref_matches = self._find_sequence_in_reference(consensus_itd)
            
            if ref_matches:
                # Calculate insertion point (where the duplication occurred)
                insertion_point = self._estimate_insertion_point(consensus_itd, ref_matches, overlap)
                
                logger.info(f"ITD size estimated from reference comparison: {len(consensus_itd)}bp")
                
                return {
                    'actual_itd_length': len(consensus_itd),
                    'itd_sequence': consensus_itd,
                    'insertion_point': insertion_point,
                    'confidence': 0.7,
                    'method': 'reference_comparison'
                }
        
        return None
    
    def _extract_soft_clips(self, read: pysam.AlignedRead) -> Dict:
        """Extract soft-clipped sequences from a read"""
        clips = {'left': None, 'right': None}
        
        if not read.cigartuples:
            return clips
        
        query_seq = read.query_sequence
        if not query_seq:
            return clips
        
        # Left soft clip
        if read.cigartuples[0][0] == 4:  # Soft clip
            clip_len = read.cigartuples[0][1]
            clips['left'] = {
                'sequence': query_seq[:clip_len],
                'length': clip_len,
                'ref_position': read.reference_start
            }
        
        # Right soft clip  
        if read.cigartuples[-1][0] == 4:  # Soft clip
            clip_len = read.cigartuples[-1][1]
            clips['right'] = {
                'sequence': query_seq[-clip_len:],
                'length': clip_len,
                'ref_position': read.reference_end
            }
        
        return clips
    
    def _find_sequence_in_reference(self, sequence: str, min_match_length: int = 15) -> List[Dict]:
        """Find matches of a sequence in the reference genome"""
        matches = []
        
        if len(sequence) < min_match_length:
            return matches
        
        # Search for exact matches
        for i in range(len(self.reference_sequence) - len(sequence) + 1):
            if self.reference_sequence[i:i + len(sequence)] == sequence:
                matches.append({
                    'position': i,
                    'length': len(sequence),
                    'similarity': 1.0
                })
        
        # If no exact matches, search for partial matches
        if not matches:
            for i in range(len(self.reference_sequence) - min_match_length + 1):
                for match_len in range(min_match_length, min(len(sequence) + 1, 100)):
                    if i + match_len > len(self.reference_sequence):
                        break
                    
                    ref_subseq = self.reference_sequence[i:i + match_len]
                    if sequence.startswith(ref_subseq) or sequence.endswith(ref_subseq):
                        similarity = match_len / len(sequence)
                        if similarity >= 0.7:  # At least 70% match
                            matches.append({
                                'position': i,
                                'length': match_len,
                                'similarity': similarity
                            })
        
        # Sort by similarity and length
        matches.sort(key=lambda x: (x['similarity'], x['length']), reverse=True)
        return matches[:5]  # Return top 5 matches
    
    def _is_likely_duplication(self, clip_info: Dict, ref_match: Dict, overlap: Dict) -> bool:
        """Check if a soft-clip likely represents a duplication"""
        
        # The reference match should be within or near the overlap region
        match_pos = ref_match['position']
        overlap_start = overlap['overlap_start']
        overlap_end = overlap['overlap_end']
        
        # Check if match is within reasonable distance of overlap
        if not (overlap_start - 100 <= match_pos <= overlap_end + 100):
            return False
        
        # Check if the clip position makes sense for a duplication
        clip_ref_pos = clip_info['ref_position']
        
        # For a duplication, we expect the clip to be either:
        # 1. At the start of the duplicated region (right clip)
        # 2. At the end of the duplicated region (left clip)
        
        distance_to_match = abs(clip_ref_pos - match_pos)
        if distance_to_match > 200:  # Too far to be related
            return False
        
        return True
    
    def _extract_sequence_for_region(self, read: pysam.AlignedRead, start: int, end: int) -> Optional[str]:
        """Extract read sequence that aligns to a specific reference region"""
        
        if not read.cigartuples or not read.query_sequence:
            return None
        
        # Build mapping from reference position to query position
        ref_pos = read.reference_start
        query_pos = 0
        ref_to_query = {}
        
        for op, length in read.cigartuples:
            if op == 0:  # Match/mismatch
                for i in range(length):
                    if start <= ref_pos + i < end:
                        ref_to_query[ref_pos + i] = query_pos + i
                ref_pos += length
                query_pos += length
            elif op == 1:  # Insertion
                # For insertions, map to the reference position where they occur
                if start <= ref_pos < end:
                    for i in range(length):
                        ref_to_query[ref_pos] = query_pos + i
                query_pos += length
            elif op == 2:  # Deletion
                ref_pos += length
            elif op == 4:  # Soft clip
                query_pos += length
        
        # Extract sequence for the region
        extracted_seq = []
        for ref_pos in range(start, end):
            if ref_pos in ref_to_query:
                query_idx = ref_to_query[ref_pos]
                if 0 <= query_idx < len(read.query_sequence):
                    extracted_seq.append(read.query_sequence[query_idx])
        
        return ''.join(extracted_seq) if extracted_seq else None
    
    def _find_consensus_itd_sequence(self, r1_seq: str, r2_seq: str, ref_seq: str) -> Optional[str]:
        """Find the consensus ITD sequence from two read sequences"""
        
        # Simple approach: find the longest common subsequence between reads
        # that is NOT present in the reference at this position
        
        min_len = min(len(r1_seq), len(r2_seq))
        if min_len < 15:
            return None
        
        # Find longest common prefix
        common_prefix = ""
        for i in range(min_len):
            if r1_seq[i] == r2_seq[i]:
                common_prefix += r1_seq[i]
            else:
                break
        
        # Find longest common suffix  
        common_suffix = ""
        for i in range(1, min_len + 1):
            if r1_seq[-i] == r2_seq[-i]:
                common_suffix = r1_seq[-i] + common_suffix
            else:
                break
        
        # The ITD is likely the part that's common between reads but extra compared to reference
        if len(common_prefix) >= 15 and common_prefix != ref_seq[:len(common_prefix)]:
            return common_prefix
        elif len(common_suffix) >= 15 and common_suffix != ref_seq[-len(common_suffix):]:
            return common_suffix
        
        # If neither prefix nor suffix works, look for internal duplications
        # This is more complex and would require sequence alignment algorithms
        
        return None
    
    def _estimate_insertion_point(self, itd_sequence: str, ref_matches: List[Dict], overlap: Dict) -> int:
        """Estimate where the ITD was inserted based on reference matches"""
        
        if not ref_matches:
            return overlap['overlap_start']
        
        # Use the best reference match to estimate insertion point
        best_match = ref_matches[0]
        
        # The insertion point is typically after the original location
        # of the duplicated sequence
        insertion_point = best_match['position'] + best_match['length']
        
        # Ensure it's within reasonable bounds
        overlap_start = overlap['overlap_start']
        overlap_end = overlap['overlap_end']
        
        if not (overlap_start - 50 <= insertion_point <= overlap_end + 50):
            # Fallback to overlap midpoint
            insertion_point = (overlap_start + overlap_end) // 2
        
        return insertion_point


def correct_itd_candidates_sizes(candidates: List, reference_sequence: str, 
                                bam_file: str) -> List:
    """
    Correct ITD size estimates for all candidates
    
    Args:
        candidates: List of ITD candidates with potentially inflated sizes
        reference_sequence: Reference genome sequence
        bam_file: BAM file with aligned read pairs
        
    Returns:
        List of candidates with corrected sizes
    """
    
    corrector = ITDSizeCorrector(reference_sequence)
    corrected_candidates = []
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        reads_dict = {read.query_name: read for read in bam}
    
    for candidate in candidates:
        if candidate.support_type == 'overlap' and hasattr(candidate, 'source_overlaps'):
            # Try to correct size based on original overlap data
            for overlap in candidate.source_overlaps or []:
                r1 = reads_dict.get(overlap.get('r1_name'))
                r2 = reads_dict.get(overlap.get('r2_name'))
                
                if r1 and r2:
                    correction = corrector.correct_itd_size(r1, r2, overlap)
                    
                    if correction:
                        # Create corrected candidate
                        corrected_candidate = candidate
                        corrected_candidate.length = correction['actual_itd_length']
                        corrected_candidate.sequence = correction['itd_sequence']
                        corrected_candidate.position = correction['insertion_point']
                        
                        # Update confidence based on correction method
                        if correction['confidence'] > candidate.confidence:
                            corrected_candidate.confidence = correction['confidence']
                        
                        logger.info(f"Corrected ITD size: {candidate.length}bp → {correction['actual_itd_length']}bp "
                                   f"(method: {correction['method']})")
                        
                        corrected_candidates.append(corrected_candidate)
                        break
            else:
                # Keep original if no correction found
                corrected_candidates.append(candidate)
        else:
            # Keep soft-clip candidates as they are typically more accurate
            corrected_candidates.append(candidate)
    
    return corrected_candidates
