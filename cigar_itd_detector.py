#!/usr/bin/env python3
"""
CIGAR ITD Detector Module
Direct ITD detection from CIGAR insertion operations in whole amplicon reads
"""

import logging
import pysam
from typing import List, Dict, Optional, Tuple
from collections import defaultdict, Counter
from dataclasses import dataclass
import numpy as np

logger = logging.getLogger(__name__)

@dataclass
class ITDCandidate:
    """Represents an ITD candidate from CIGAR analysis"""
    sequence: str
    length: int
    position: int
    support_type: str = 'cigar_insertion'
    supporting_reads: List[str] = None
    confidence: float = 0.0
    insertion_site: Optional[int] = None
    duplication_start: Optional[int] = None
    duplication_end: Optional[int] = None
    is_primary: bool = False

    def __post_init__(self):
        if self.supporting_reads is None:
            self.supporting_reads = []

class CigarITDDetector:
    """Detect ITDs directly from CIGAR insertion operations"""
    
    def __init__(self, reference_sequence: str, min_itd_length: int = 15, 
                 min_support: int = 3, max_itd_length: int = 500,
                 position_tolerance: int = 10):
        """
        Initialize CIGAR ITD detector
        
        Args:
            reference_sequence: FLT3 reference sequence
            min_itd_length: Minimum ITD size to consider
            min_support: Minimum supporting reads required
            max_itd_length: Maximum ITD size to consider
            position_tolerance: Position clustering tolerance (bp)
        """
        self.reference_sequence = reference_sequence
        self.min_itd_length = min_itd_length
        self.max_itd_length = max_itd_length
        self.min_support = min_support
        self.position_tolerance = position_tolerance
        
        logger.info(f"Initialized CIGAR ITD detector: {min_itd_length}-{max_itd_length}bp, "
                   f"min_support={min_support}, position_tolerance={position_tolerance}bp")
    
    def detect_itds_from_bam(self, bam_file: str) -> List[ITDCandidate]:
        """
        Main function to detect ITDs from BAM file using CIGAR analysis
        
        Args:
            bam_file: Path to BAM file with FLT3 reads
            
        Returns:
            List of ITD candidates
        """
        logger.info(f"Starting CIGAR ITD detection from {bam_file}")
        
        # Step 1: Extract CIGAR insertions from all reads
        insertions = self._extract_cigar_insertions(bam_file)
        logger.info(f"Found {len(insertions)} insertions ≥{self.min_itd_length}bp in CIGAR strings")
        
        if not insertions:
            logger.info("No significant insertions found in CIGAR strings")
            return []
        
        # Step 2: Group insertions by position and sequence similarity
        insertion_groups = self._group_insertions(insertions)
        logger.info(f"Grouped insertions into {len(insertion_groups)} position clusters")
        
        # Step 3: Generate consensus ITD candidates from groups
        candidates = self._generate_candidates_from_groups(insertion_groups)
        logger.info(f"Generated {len(candidates)} ITD candidates from CIGAR insertions")
        
        # Step 4: Validate and filter candidates
        validated_candidates = self._validate_candidates(candidates)
        logger.info(f"Validated candidates: {len(candidates)} → {len(validated_candidates)}")
        
        return validated_candidates
    
    def _extract_cigar_insertions(self, bam_file: str) -> List[Dict]:
        """
        Extract insertion operations from CIGAR strings
        
        Returns:
            List of insertion dictionaries with sequence, position, and read info
        """
        insertions = []
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                # Skip unmapped, secondary, or supplementary reads
                if (read.is_unmapped or read.is_secondary or 
                    read.is_supplementary or read.is_duplicate):
                    continue
                
                # Skip reads without CIGAR or sequence
                if not read.cigartuples or not read.query_sequence:
                    continue
                
                # Parse CIGAR for insertions
                read_insertions = self._parse_cigar_insertions(read)
                insertions.extend(read_insertions)
        
        # Filter by size
        size_filtered = [ins for ins in insertions 
                        if self.min_itd_length <= ins['length'] <= self.max_itd_length]
        
        logger.debug(f"CIGAR insertions: {len(insertions)} total → {len(size_filtered)} size-filtered")
        
        return size_filtered
    
    def _parse_cigar_insertions(self, read: pysam.AlignedRead) -> List[Dict]:
        """
        Parse CIGAR string to extract insertion operations
        
        Args:
            read: Aligned read
            
        Returns:
            List of insertion dictionaries
        """
        insertions = []
        ref_pos = read.reference_start
        query_pos = 0
        
        for operation, length in read.cigartuples:
            if operation == 1:  # Insertion
                if length >= self.min_itd_length:
                    # Extract insertion sequence
                    insertion_seq = read.query_sequence[query_pos:query_pos + length]
                    
                    insertions.append({
                        'sequence': insertion_seq,
                        'length': length,
                        'ref_position': ref_pos,
                        'query_start': query_pos,
                        'query_end': query_pos + length,
                        'read_name': read.query_name,
                        'read_length': len(read.query_sequence),
                        'mapping_quality': read.mapping_quality,
                        'is_reverse': read.is_reverse
                    })
                
                query_pos += length
                # Note: insertions don't advance reference position
                
            elif operation in [0, 2]:  # Match/mismatch or deletion
                ref_pos += length
                if operation == 0:  # Only match/mismatch advances query
                    query_pos += length
                    
            elif operation in [4, 5]:  # Soft/hard clip
                if operation == 4:  # Soft clip advances query
                    query_pos += length
        
        return insertions
    
    def _group_insertions(self, insertions: List[Dict]) -> Dict[str, List[Dict]]:
        """
        Group insertions by genomic position and sequence similarity
        
        Args:
            insertions: List of insertion dictionaries
            
        Returns:
            Dictionary mapping group keys to insertion lists
        """
        position_groups = defaultdict(list)
        
        # First, group by position with tolerance
        for insertion in insertions:
            position_bin = insertion['ref_position'] // self.position_tolerance
            position_groups[position_bin].append(insertion)
        
        # Second, within each position group, cluster by sequence similarity
        final_groups = {}
        group_id = 0
        
        for position_bin, pos_insertions in position_groups.items():
            sequence_clusters = self._cluster_by_sequence_similarity(pos_insertions)
            
            for cluster in sequence_clusters:
                if len(cluster) >= 1:  # Keep all clusters for now, filter later
                    group_key = f"group_{group_id}"
                    final_groups[group_key] = cluster
                    group_id += 1
        
        # Log grouping statistics
        group_sizes = [len(group) for group in final_groups.values()]
        if group_sizes:
            logger.debug(f"Insertion grouping: {len(final_groups)} groups, "
                        f"sizes: min={min(group_sizes)}, max={max(group_sizes)}, "
                        f"mean={np.mean(group_sizes):.1f}")
        
        return final_groups
    
    def _cluster_by_sequence_similarity(self, insertions: List[Dict], 
                                       similarity_threshold: float = 0.85) -> List[List[Dict]]:
        """
        Cluster insertions by sequence similarity, preserving read support information
        
        Args:
            insertions: List of insertions at similar positions
            similarity_threshold: Minimum similarity to cluster together
            
        Returns:
            List of clusters (each cluster contains aggregated insertion data with support counts)
        """
        if len(insertions) <= 1:
            return [insertions] if insertions else []
        
        # First, aggregate identical sequences to get read support counts
        sequence_groups = {}
        for insertion in insertions:
            seq = insertion['sequence']
            pos = insertion['ref_position']
            
            # Create a key that includes both sequence and position for exact matching
            key = f"{seq}_{pos}"
            
            if key not in sequence_groups:
                # Create aggregated insertion with support count
                sequence_groups[key] = {
                    'sequence': seq,
                    'ref_position': pos,
                    'length': insertion['length'],
                    'original_support': 1,  # Start with 1 read
                    'supporting_reads': [insertion['read_name']],
                    'mapping_qualities': [insertion['mapping_quality']],
                    'query_positions': [insertion['query_start']] if 'query_start' in insertion else []
                }
            else:
                # Increment support count for this exact sequence+position
                sequence_groups[key]['original_support'] += 1
                sequence_groups[key]['supporting_reads'].append(insertion['read_name'])
                sequence_groups[key]['mapping_qualities'].append(insertion['mapping_quality'])
                if 'query_start' in insertion:
                    sequence_groups[key]['query_positions'].append(insertion['query_start'])
        
        # Convert to list of aggregated insertions
        aggregated_insertions = list(sequence_groups.values())
        
        # Now cluster these aggregated insertions by similarity
        clusters = []
        used_indices = set()
        
        for i, insertion in enumerate(aggregated_insertions):
            if i in used_indices:
                continue
            
            # Start new cluster
            cluster = [insertion]
            used_indices.add(i)
            
            # Find similar insertions
            for j, other_insertion in enumerate(aggregated_insertions[i+1:], i+1):
                if j in used_indices:
                    continue
                
                similarity = self._calculate_sequence_similarity(
                    insertion['sequence'], other_insertion['sequence']
                )
                
                if similarity >= similarity_threshold:
                    cluster.append(other_insertion)
                    used_indices.add(j)
            
            clusters.append(cluster)
        
        return clusters
    
    def _calculate_sequence_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate sequence similarity between two sequences
        
        Args:
            seq1, seq2: Sequences to compare
            
        Returns:
            Similarity score (0.0 to 1.0)
        """
        if not seq1 or not seq2:
            return 0.0
        
        # Quick length check
        len_diff = abs(len(seq1) - len(seq2))
        max_len = max(len(seq1), len(seq2))
        
        if len_diff > max_len * 0.2:  # >20% length difference
            return 0.0
        
        # Simple string similarity for now (can be improved with alignment)
        if seq1 == seq2:
            return 1.0
        
        # Check for substring relationships
        if seq1 in seq2 or seq2 in seq1:
            min_len = min(len(seq1), len(seq2))
            return min_len / max_len
        
        # Hamming distance for same-length sequences
        if len(seq1) == len(seq2):
            matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
            return matches / len(seq1)
        
        # For different lengths, use longest common subsequence approximation
        return self._lcs_similarity(seq1, seq2)
    
    def _lcs_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate similarity based on longest common subsequence
        
        Args:
            seq1, seq2: Sequences to compare
            
        Returns:
            Similarity score based on LCS
        """
        if not seq1 or not seq2:
            return 0.0
        
        # Simple LCS approximation - count common k-mers
        k = 3  # Use 3-mers
        if len(seq1) < k or len(seq2) < k:
            return 0.0
        
        kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
        kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
        
        common_kmers = len(kmers1 & kmers2)
        total_kmers = len(kmers1 | kmers2)
        
        return common_kmers / total_kmers if total_kmers > 0 else 0.0
    
    def _generate_candidates_from_groups(self, insertion_groups: Dict[str, List[Dict]]) -> List[ITDCandidate]:
        """
        Generate ITD candidates from insertion groups
        
        Args:
            insertion_groups: Grouped insertions
            
        Returns:
            List of ITD candidates
        """
        candidates = []
        
        for group_name, insertions in insertion_groups.items():
            # Calculate total support from aggregated insertions
            total_support = sum(ins.get('original_support', 1) for ins in insertions)
            
            if total_support < self.min_support:
                logger.debug(f"Skipping {group_name}: insufficient support ({total_support} < {self.min_support})")
                continue
            
            # Generate consensus sequence and position
            consensus_result = self._generate_consensus(insertions)
            
            if not consensus_result:
                continue
            
            # Collect all supporting reads from aggregated insertions
            all_supporting_reads = []
            for ins in insertions:
                if 'supporting_reads' in ins:
                    all_supporting_reads.extend(ins['supporting_reads'])
                else:
                    # Fallback for non-aggregated insertions
                    all_supporting_reads.append(ins.get('read_name', 'unknown'))
            
            # Create ITD candidate
            candidate = ITDCandidate(
                sequence=consensus_result['sequence'],
                length=len(consensus_result['sequence']),
                position=consensus_result['position'],
                support_type='cigar_insertion',
                supporting_reads=all_supporting_reads,
                confidence=consensus_result['confidence'],
                insertion_site=consensus_result['position']
            )
            
            candidates.append(candidate)
            
            logger.debug(f"Generated candidate from {group_name}: {candidate.length}bp at pos {candidate.position}, "
                        f"support={len(candidate.supporting_reads)}, confidence={candidate.confidence:.3f}")
        
        return candidates
    
    def _generate_consensus(self, insertions: List[Dict]) -> Optional[Dict]:
        """
        Generate consensus sequence and position from group of insertions
        
        Args:
            insertions: List of similar insertions
            
        Returns:
            Dictionary with consensus information or None
        """
        if not insertions:
            return None
        
        # Calculate weighted consensus based on original read support
        # Each insertion carries information about how many reads supported it
        
        # Weighted position and sequence consensus by supporting reads
        position_weights = []
        sequence_weights = []
        total_original_support = 0

        for ins in insertions:
            original_support = ins.get('original_support', 1)
            total_original_support += original_support
            # Weight positions and sequences by supporting reads
            position_weights.extend([ins['ref_position']] * original_support)
            sequence_weights.extend([ins['sequence']] * original_support)

        # Consensus position: median of all supporting reads' positions
        consensus_position = int(np.median(position_weights))

        # Consensus sequence: most supported sequence by read count
        sequence_counts = Counter(sequence_weights)
        consensus_sequence = sequence_counts.most_common(1)[0][0]
        most_supported_count = sequence_counts.most_common(1)[0][1]

        # Calculate confidence based on weighted support
        support_count = total_original_support
        sequence_consistency = most_supported_count / total_original_support

        # Position clustering quality (use all supporting positions)
        position_std = np.std(position_weights) if len(position_weights) > 1 else 0
        position_quality = max(0, 1 - (position_std / self.position_tolerance))

        # Overall confidence calculation
        base_confidence = 0.5
        support_bonus = min(0.3, support_count * 0.05)
        consistency_bonus = sequence_consistency * 0.15
        position_bonus = position_quality * 0.05

        confidence = min(0.95, base_confidence + support_bonus + consistency_bonus + position_bonus)

        return {
            'sequence': consensus_sequence,
            'position': consensus_position,
            'confidence': confidence,
            'support_count': support_count,
            'sequence_consistency': sequence_consistency,
            'position_std': position_std
        }
    
    def _validate_candidates(self, candidates: List[ITDCandidate]) -> List[ITDCandidate]:
        """
        Validate ITD candidates using biological and technical criteria
        
        Args:
            candidates: List of ITD candidates
            
        Returns:
            List of validated candidates
        """
        validated = []
        
        for candidate in candidates:
            validation_result = self._validate_single_candidate(candidate)
            
            if validation_result['is_valid']:
                # Update confidence based on validation
                candidate.confidence = min(0.95, candidate.confidence * validation_result['validation_factor'])
                validated.append(candidate)
                
                logger.debug(f"Validated candidate: {candidate.length}bp, "
                            f"final_confidence={candidate.confidence:.3f}")
            else:
                logger.debug(f"Rejected candidate: {candidate.length}bp, "
                            f"reason={validation_result['reason']}")
        
        return validated
    
    def _validate_single_candidate(self, candidate: ITDCandidate) -> Dict:
        """
        Validate a single ITD candidate
        
        Args:
            candidate: ITD candidate to validate
            
        Returns:
            Dictionary with validation results
        """
        # Check 1: Sequence complexity
        if not self._check_sequence_complexity(candidate.sequence):
            return {'is_valid': False, 'reason': 'low_sequence_complexity', 'validation_factor': 0.0}
        
        # Check 2: Reference similarity (ITDs should have some reference similarity)
        ref_similarity = self._check_reference_similarity(candidate.sequence)
        if ref_similarity < 0.3:  # At least 30% similarity to reference
            return {'is_valid': False, 'reason': 'low_reference_similarity', 'validation_factor': 0.0}
        
        # Check 3: Support quality
        if len(candidate.supporting_reads) < self.min_support:
            return {'is_valid': False, 'reason': 'insufficient_support', 'validation_factor': 0.0}
        
        # Check 4: Length within bounds
        if not (self.min_itd_length <= candidate.length <= self.max_itd_length):
            return {'is_valid': False, 'reason': 'length_out_of_bounds', 'validation_factor': 0.0}
        
        # Calculate validation factor based on quality metrics
        validation_factor = min(1.2, 1.0 + ref_similarity * 0.2)  # Boost for good reference similarity
        
        return {
            'is_valid': True,
            'reason': 'passed_validation',
            'validation_factor': validation_factor,
            'reference_similarity': ref_similarity
        }
    
    def _check_sequence_complexity(self, sequence: str) -> bool:
        """
        Check if sequence has sufficient complexity
        
        Args:
            sequence: Sequence to check
            
        Returns:
            True if sequence has good complexity
        """
        if len(sequence) < self.min_itd_length:
            return False
        
        # Check nucleotide diversity
        unique_bases = len(set(sequence))
        if unique_bases < 3:  # Need at least 3 different bases
            return False
        
        # Check for excessive homopolymer runs
        max_homopolymer = max(len(list(g)) for k, g in 
                             __import__('itertools').groupby(sequence))
        
        if max_homopolymer > min(10, len(sequence) * 0.3):  # Max 30% homopolymer
            return False
        
        return True
    
    def _check_reference_similarity(self, sequence: str) -> float:
        """
        Check similarity to reference sequence
        
        Args:
            sequence: ITD sequence to check
            
        Returns:
            Similarity score (0.0 to 1.0)
        """
        if not sequence:
            return 0.0
        
        # Simple approach: find best local alignment to reference
        best_similarity = 0.0
        seq_len = len(sequence)
        
        # Sliding window approach
        for i in range(len(self.reference_sequence) - seq_len + 1):
            ref_window = self.reference_sequence[i:i + seq_len]
            similarity = self._calculate_sequence_similarity(sequence, ref_window)
            best_similarity = max(best_similarity, similarity)
        
        # Also check substring matches
        for k in range(self.min_itd_length, min(seq_len + 1, 50)):  # Check up to 50bp substrings
            for i in range(seq_len - k + 1):
                subseq = sequence[i:i + k]
                if subseq in self.reference_sequence:
                    substring_similarity = k / seq_len
                    best_similarity = max(best_similarity, substring_similarity)
        
        return best_similarity


def detect_cigar_itds(bam_file: str, reference_sequence: str, 
                     min_itd_length: int = 15, min_support: int = 3,
                     max_itd_length: int = 500) -> List[ITDCandidate]:
    """
    Main function to detect ITDs from CIGAR insertions
    
    Args:
        bam_file: Path to BAM file with FLT3 reads
        reference_sequence: FLT3 reference sequence
        min_itd_length: Minimum ITD size
        min_support: Minimum supporting reads
        max_itd_length: Maximum ITD size
        
    Returns:
        List of ITD candidates
    """
    detector = CigarITDDetector(
        reference_sequence=reference_sequence,
        min_itd_length=min_itd_length,
        min_support=min_support,
        max_itd_length=max_itd_length
    )
    
    return detector.detect_itds_from_bam(bam_file)


if __name__ == "__main__":
    print("This module should be imported and used through the main FLT3 ITD detection pipeline.")
    print("Run 'python main_module.py' instead for command-line usage.")
