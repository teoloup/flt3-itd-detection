#!/usr/bin/env python3
"""
ITD Generator Module
Generates putative ITD sequences from overlaps and soft clips
"""

import logging
from typing import List, Dict, Optional, Tuple
from collections import defaultdict, Counter
import pysam
from dataclasses import dataclass
import re

# Try to import Levenshtein for fast sequence similarity
try:
    import Levenshtein
    _has_levenshtein = True
except ImportError:
    _has_levenshtein = False

logger = logging.getLogger(__name__)

@dataclass
class ITDCandidate:
    """Represents an ITD candidate"""
    sequence: str
    length: int
    position: int
    support_type: str  # 'overlap' or 'softclip'
    supporting_reads: List[str]
    confidence: float
    insertion_site: Optional[int] = None
    duplication_start: Optional[int] = None
    duplication_end: Optional[int] = None
    source_overlaps: Optional[List] = None  # Add this field for compatibility

class ITDGenerator:
    """Generate ITD candidates from alignment results"""
    
    def __init__(self, reference_sequence: str, min_itd_length: int = 15, 
                 min_support: int = 2, max_itd_length: int = 500):
        self.reference_sequence = reference_sequence
        self.min_itd_length = min_itd_length
        self.max_itd_length = max_itd_length
        self.min_support = min_support
        # Pre-compute reference k-mers for faster matching
        self._build_reference_index()
        
    def _build_reference_index(self):
        """Build k-mer index of reference for faster matching"""
        self.ref_kmers = {}
        k = 10  # k-mer size
        for i in range(len(self.reference_sequence) - k + 1):
            kmer = self.reference_sequence[i:i+k]
            if kmer not in self.ref_kmers:
                self.ref_kmers[kmer] = []
            self.ref_kmers[kmer].append(i)
    
    def extract_overlap_sequences(self, overlaps: List[Dict], bam_file: str) -> List[ITDCandidate]:
        """Extract ITD sequences from overlapping pairs (parallelized per original read)"""
        import concurrent.futures
        candidates = []
        # Pre-filter overlaps to reduce redundancy
        filtered_overlaps = self._prefilter_overlaps(overlaps)
        logger.info(f"Pre-filtered {len(overlaps)} overlaps to {len(filtered_overlaps)}")
        # Group overlaps by original read
        overlap_groups = defaultdict(list)
        for overlap in filtered_overlaps:
            overlap_groups[overlap['original_read']].append(overlap)
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Create read dictionary for quick lookup
            reads_dict = {read.query_name: read for read in bam}

            def process_group(args):
                original_read, read_overlaps = args
                best_overlap = max(read_overlaps, key=lambda x: x['length'])
                r1 = reads_dict.get(best_overlap['r1_name'])
                r2 = reads_dict.get(best_overlap['r2_name'])
                if not r1 or not r2:
                    return None
                
                # Extract complete ITD sequence with enhanced reconstruction
                overlap_seq = self._extract_complete_itd_sequence(r1, r2, best_overlap)
                
                if overlap_seq and len(overlap_seq) >= self.min_itd_length:
                    if self._validate_itd_sequence(overlap_seq):
                        # Calculate confidence based on multiple factors
                        base_confidence = 0.6
                        length_bonus = min(0.2, len(overlap_seq) / 200)  # Bonus for longer sequences
                        
                        return ITDCandidate(
                            sequence=overlap_seq,
                            length=len(overlap_seq),
                            position=best_overlap['overlap_start'],
                            support_type='overlap',
                            supporting_reads=[original_read],
                            confidence=base_confidence + length_bonus,
                            insertion_site=best_overlap['overlap_start'],
                            source_overlaps=[best_overlap]  # Initialize source_overlaps
                        )
                return None

            with concurrent.futures.ThreadPoolExecutor() as executor:
                results = list(executor.map(process_group, overlap_groups.items()))
            candidates = [c for c in results if c is not None]
        return candidates
    
 
    def _quick_validate_overlap(self, overlap: Dict) -> bool:
        """Quick validation of overlap before expensive processing"""
        # Check if overlap length is reasonable
        if overlap['length'] < self.min_itd_length or overlap['length'] > self.max_itd_length:
            return False
        
        # Check if alignment lengths suggest ITD
        r1_len = overlap.get('r1_aligned_length', 0)
        r2_len = overlap.get('r2_aligned_length', 0)
        
        # If both align to similar lengths, probably not ITD
        if abs(r1_len - r2_len) < 10:
            return False
        
        return True

    def _classify_itd_sequence(self, sequence):
        """Classify ITD sequence and assign confidence with relaxed criteria"""
        
        # Check for exact match in reference
        if sequence in self.reference_sequence:
            return "pure_duplication", 0.9
        
        # Check for partial matches (sequence contains reference parts)
        best_match_length = 0
        k = 15  # Check 15bp windows
        
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if kmer in self.reference_sequence:
                best_match_length = max(best_match_length, k)
                # Try to extend the match
                for j in range(k, min(len(sequence) - i, 100)):
                    extended = sequence[i:i+j]
                    if extended in self.reference_sequence:
                        best_match_length = max(best_match_length, j)
                    else:
                        break
        
        # Classify based on how much matches reference
        match_ratio = best_match_length / len(sequence)
        
        if match_ratio >= 0.8:
            return "duplication_with_insertion", 0.8
        elif match_ratio >= 0.5:
            return "complex_itd", 0.7
        elif match_ratio >= 0.3:
            return "partial_duplication", 0.6
        elif match_ratio >= 0.1:  # New category: minimal reference similarity
            return "novel_itd_with_reference", 0.5
        else:
            # Novel sequences with good read support should still be considered
            # ITDs can be entirely novel insertions or complex rearrangements
            return "novel_insertion", 0.4  # Reduced but still acceptable confidence
    
    def _quick_reference_check(self, sequence: str) -> bool:
        """Quick check if sequence likely exists in reference"""
        # Check first k-mer
        k = 10
        if len(sequence) >= k:
            first_kmer = sequence[:k]
            return first_kmer in self.ref_kmers
        return False
    
    def merge_candidates(self, candidates: List[ITDCandidate]) -> List[ITDCandidate]:
        """Merge similar ITD candidates - optimized version"""
        if not candidates:
            return []
        
        # Sort by position and length
        candidates.sort(key=lambda x: (x.position, x.length))
        
        # Quick merge based on position and length only
        merged = []
        position_bins = defaultdict(list)
        
        # Bin candidates
        for cand in candidates:
            key = (cand.position // 5, cand.length // 5)
            position_bins[key].append(cand)
        
        # Take best from each bin
        for bin_candidates in position_bins.values():
            if len(bin_candidates) == 1:
                merged.append(bin_candidates[0])
            else:
                # Merge by combining support
                best = max(bin_candidates, key=lambda x: len(x.supporting_reads))  # Prioritize support over confidence
                
                # Combine supporting reads
                all_reads = []
                for c in bin_candidates:
                    all_reads.extend(c.supporting_reads)
                best.supporting_reads = list(set(all_reads))
                
                # Update confidence based on supporting read count
                support_bonus = min(0.3, len(best.supporting_reads) * 0.02)  # Bonus for high support
                merge_bonus = min(0.1, (len(bin_candidates) - 1) * 0.02)  # Bonus for merging multiple candidates
                best.confidence = min(0.95, best.confidence + support_bonus + merge_bonus)
                
                merged.append(best)
        
    def _merge_similar_candidates(self, candidates: List[ITDCandidate]) -> List[ITDCandidate]:
        """Merge ITD candidates that likely represent the same biological event"""
        if len(candidates) <= 1:
            return candidates
        
        # Sort by confidence and length (prioritize better candidates)
        candidates.sort(key=lambda x: (x.confidence, x.length), reverse=True)
        
        merged = []
        used_indices = set()
        
        for i, candidate in enumerate(candidates):
            if i in used_indices:
                continue
            
            # Look for similar candidates to merge
            similar_candidates = [candidate]
            used_indices.add(i)
            
            for j, other in enumerate(candidates[i+1:], i+1):
                if j in used_indices:
                    continue
                
                # Check if candidates represent the same biological event
                if self._are_candidates_similar(candidate, other):
                    similar_candidates.append(other)
                    used_indices.add(j)
                    logger.debug(f"Merging similar ITD candidates: {candidate.length}bp vs {other.length}bp")
            
            # Create merged candidate if we found similar ones
            if len(similar_candidates) > 1:
                merged_candidate = self._create_merged_candidate(similar_candidates)
                merged.append(merged_candidate)
            else:
                merged.append(candidate)
        
        return merged
    
    def _are_candidates_similar(self, candidate1: ITDCandidate, candidate2: ITDCandidate) -> bool:
        """Check if two ITD candidates represent the same biological event"""
        
        # Check position similarity (within 50bp)
        pos_diff = abs(candidate1.position - candidate2.position)
        if pos_diff > 50:
            return False
        
        # Check sequence similarity
        seq1, seq2 = candidate1.sequence, candidate2.sequence
        
        # Check if one sequence contains the other (substring relationship)
        if seq1 in seq2 or seq2 in seq1:
            return True
        
        # Check for significant overlap at sequence level
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))
        
        # Check overlap at start or end
        for overlap_len in range(min(min_len, 30), min_len + 1):
            if (seq1[:overlap_len] == seq2[:overlap_len] or  # Same start
                seq1[-overlap_len:] == seq2[-overlap_len:]):  # Same end
                overlap_ratio = overlap_len / max_len
                if overlap_ratio >= 0.6:  # 60% overlap
                    return True
        
        return False
    
    def _create_merged_candidate(self, candidates: List[ITDCandidate]) -> ITDCandidate:
        """Create a merged ITD candidate from similar candidates"""
        
        # Choose the longest sequence as the representative
        best_candidate = max(candidates, key=lambda x: len(x.sequence))
        
        # Merge supporting reads from all candidates
        all_supporting_reads = []
        for c in candidates:
            all_supporting_reads.extend(c.supporting_reads)
        unique_supporting_reads = list(set(all_supporting_reads))
        
        # Merge confidence (weighted by length and confidence)
        total_weight = sum(c.confidence * c.length for c in candidates)
        total_length = sum(c.length for c in candidates)
        merged_confidence = min(0.95, total_weight / total_length if total_length > 0 else best_candidate.confidence)
        
        # Use the best insertion site and duplication coordinates
        insertion_sites = [c.insertion_site for c in candidates if c.insertion_site is not None]
        duplication_starts = [c.duplication_start for c in candidates if c.duplication_start is not None]
        duplication_ends = [c.duplication_end for c in candidates if c.duplication_end is not None]
        
        # Use most common or median values
        insertion_site = insertion_sites[0] if insertion_sites else None
        duplication_start = duplication_starts[0] if duplication_starts else None
        duplication_end = duplication_ends[0] if duplication_ends else None
        
        # Collect source overlaps if available
        all_source_overlaps = []
        for c in candidates:
            if hasattr(c, 'source_overlaps') and c.source_overlaps:
                all_source_overlaps.extend(c.source_overlaps)
        
        logger.info(f"Merged {len(candidates)} ITD candidates into {len(best_candidate.sequence)}bp ITD "
                   f"(confidence: {merged_confidence:.3f})")
        
        return ITDCandidate(
            sequence=best_candidate.sequence,
            length=len(best_candidate.sequence),
            position=best_candidate.position,
            support_type=best_candidate.support_type,
            supporting_reads=unique_supporting_reads,
            confidence=merged_confidence,
            insertion_site=insertion_site,
            duplication_start=duplication_start,
            duplication_end=duplication_end,
            source_overlaps=all_source_overlaps if all_source_overlaps else None
        )
    
    def _prefilter_overlaps(self, overlaps: List[Dict]) -> List[Dict]:
        """Pre-filter overlaps to remove obvious false positives"""
        filtered = []
        
        # Group by original read and overlap characteristics
        overlap_groups = defaultdict(list)
        for overlap in overlaps:
            # Create a key based on overlap characteristics
            key = (
                overlap['original_read'],
                overlap['length'] // 10,  # Group by 10bp bins
                overlap['overlap_start'] // 20  # Group by 20bp position bins
            )
            overlap_groups[key].append(overlap)
        
        # Keep only the best overlap from each group
        for group in overlap_groups.values():
            # Sort by overlap length and alignment quality indicators
            best = max(group, key=lambda x: (
                x['length'],
                -abs(x.get('r1_aligned_length', 0) - x.get('r2_aligned_length', 0))
            ))
            filtered.append(best)
        
        return filtered
    
    def _validate_itd_sequence(self, sequence: str) -> bool:
        """Validate if a sequence looks like a real ITD"""
        # Check sequence complexity (not just homopolymers)
        if len(set(sequence)) < 3:  # Less than 3 different bases (relaxed from 4)
            return False
        
        # Check for excessive homopolymer runs
        for base in 'ATGC':
            if base * 10 in sequence:  # 10 or more of the same base (relaxed from 8)
                return False
        
        # Check if sequence has some similarity to reference (relaxed validation)
        # ITDs can include novel insertions, so we don't require exact matches
        ref_matches = self._find_in_reference_fuzzy(sequence, max_mismatches=len(sequence) // 5)  # Allow more mismatches
        
        # Accept if we find partial matches or if sequence is long enough with good complexity
        if ref_matches:
            return True  # Has some reference similarity
        elif len(sequence) >= self.min_itd_length * 2 and len(set(sequence)) >= 3:
            return True  # Long enough with good complexity, could be novel ITD
        else:
            return False  # Too short and no reference similarity
    
    def _extract_overlap_sequence(self, r1: pysam.AlignedRead, r2: pysam.AlignedRead, 
                                  overlap: Dict) -> Optional[str]:
        """Extract the complete ITD sequence using enhanced reconstruction with full duplication detection"""
        try:
            # Helper to build comprehensive ref->query mapping from CIGAR
            def build_comprehensive_mapping(aligned_read):
                ref_pos = aligned_read.reference_start
                query_pos = 0
                ref_to_query = {}
                query_to_ref = {}
                
                for op, length in aligned_read.cigartuples:
                    if op == 0:  # Match/mismatch
                        for i in range(length):
                            ref_to_query[ref_pos + i] = query_pos + i
                            query_to_ref[query_pos + i] = ref_pos + i
                        ref_pos += length
                        query_pos += length
                    elif op == 1:  # Insertion - query advances but not reference
                        for i in range(length):
                            query_to_ref[query_pos + i] = ref_pos  # Insert at this ref position
                        query_pos += length
                    elif op == 2:  # Deletion - reference advances but not query
                        ref_pos += length
                    elif op == 4:  # Soft clip - query advances but not reference
                        query_pos += length
                    elif op == 5:  # Hard clip - nothing advances
                        pass
                        
                return ref_to_query, query_to_ref

            # Helper to extract soft-clipped sequences with positions
            def get_comprehensive_softclips(aligned_read):
                cigartuples = aligned_read.cigartuples
                seq = aligned_read.query_sequence
                if not cigartuples:
                    return {}, {}
                
                left_softclip = {}
                right_softclip = {}
                query_pos = 0
                
                # Left soft clip
                if cigartuples[0][0] == 4:
                    clip_len = cigartuples[0][1]
                    left_softclip = {
                        'sequence': seq[:clip_len],
                        'length': clip_len,
                        'ref_position': aligned_read.reference_start  # Where it would map
                    }
                    query_pos += clip_len
                
                # Process middle operations to find right clip position
                for op, length in cigartuples[1:-1]:
                    if op in [0, 1, 4]:  # Operations that advance query
                        query_pos += length
                
                # Right soft clip
                if cigartuples[-1][0] == 4:
                    clip_len = cigartuples[-1][1]
                    right_softclip = {
                        'sequence': seq[-clip_len:],
                        'length': clip_len,
                        'ref_position': aligned_read.reference_end  # Where it would map
                    }
                
                return left_softclip, right_softclip

            # Build mappings for both reads
            r1_ref_to_query, r1_query_to_ref = build_comprehensive_mapping(r1)
            r2_ref_to_query, r2_query_to_ref = build_comprehensive_mapping(r2)
            
            # Get soft-clipped information
            r1_left_clip, r1_right_clip = get_comprehensive_softclips(r1)
            r2_left_clip, r2_right_clip = get_comprehensive_softclips(r2)

            # Strategy: Reconstruct the complete ITD by combining all evidence
            
            # 1. Start with the overlap consensus (high confidence region)
            overlap_seq_parts = []
            consensus_count = 0
            total_bases = 0
            
            for ref_pos in range(overlap['overlap_start'], overlap['overlap_end']):
                base_r1 = None
                base_r2 = None
                
                if ref_pos in r1_ref_to_query:
                    query_pos_r1 = r1_ref_to_query[ref_pos]
                    base_r1 = r1.query_sequence[query_pos_r1]
                if ref_pos in r2_ref_to_query:
                    query_pos_r2 = r2_ref_to_query[ref_pos]
                    base_r2 = r2.query_sequence[query_pos_r2]
                
                # Build consensus
                if base_r1 and base_r2:
                    total_bases += 1
                    if base_r1 == base_r2:
                        overlap_seq_parts.append(base_r1)
                        consensus_count += 1
                    else:
                        # For disagreements, choose the base from the read with better alignment
                        overlap_seq_parts.append(base_r1)  # Default to r1
                elif base_r1:
                    overlap_seq_parts.append(base_r1)
                elif base_r2:
                    overlap_seq_parts.append(base_r2)

            overlap_consensus = ''.join(overlap_seq_parts)
            
            # 2. Extend the sequence using soft-clips and insertions to capture full ITD
            
            # Collect all sequence extensions
            extensions = []
            
            # Add overlap consensus as the core
            extensions.append({
                'sequence': overlap_consensus,
                'type': 'overlap',
                'priority': 100,  # Highest priority
                'ref_start': overlap['overlap_start'],
                'ref_end': overlap['overlap_end']
            })
            
            # Add soft-clipped sequences
            if r1_left_clip and r1_left_clip['sequence']:
                extensions.append({
                    'sequence': r1_left_clip['sequence'],
                    'type': 'left_clip',
                    'priority': 80,
                    'ref_start': r1_left_clip['ref_position'] - len(r1_left_clip['sequence']),
                    'ref_end': r1_left_clip['ref_position']
                })
            
            if r1_right_clip and r1_right_clip['sequence']:
                extensions.append({
                    'sequence': r1_right_clip['sequence'],
                    'type': 'right_clip',
                    'priority': 80,
                    'ref_start': r1_right_clip['ref_position'],
                    'ref_end': r1_right_clip['ref_position'] + len(r1_right_clip['sequence'])
                })
            
            if r2_left_clip and r2_left_clip['sequence']:
                extensions.append({
                    'sequence': r2_left_clip['sequence'],
                    'type': 'left_clip',
                    'priority': 70,
                    'ref_start': r2_left_clip['ref_position'] - len(r2_left_clip['sequence']),
                    'ref_end': r2_left_clip['ref_position']
                })
            
            if r2_right_clip and r2_right_clip['sequence']:
                extensions.append({
                    'sequence': r2_right_clip['sequence'],
                    'type': 'right_clip',
                    'priority': 70,
                    'ref_start': r2_right_clip['ref_position'],
                    'ref_end': r2_right_clip['ref_position'] + len(r2_right_clip['sequence'])
                })

            # 3. Reconstruct complete ITD by ordering and merging extensions
            
            # Sort by reference position to get proper order
            extensions.sort(key=lambda x: x['ref_start'])
            
            # Merge overlapping or adjacent extensions
            merged_extensions = []
            for ext in extensions:
                if not merged_extensions:
                    merged_extensions.append(ext)
                    continue
                
                last_ext = merged_extensions[-1]
                
                # Check if extensions can be merged (overlap or adjacent)
                if ext['ref_start'] <= last_ext['ref_end'] + 10:  # Allow small gaps
                    # Merge extensions
                    if ext['priority'] > last_ext['priority']:
                        # Replace with higher priority extension
                        merged_extensions[-1] = ext
                    elif ext['type'] == 'right_clip' and last_ext['type'] == 'overlap':
                        # Extend overlap with right clip
                        merged_extensions[-1] = {
                            'sequence': last_ext['sequence'] + ext['sequence'],
                            'type': 'extended',
                            'priority': max(last_ext['priority'], ext['priority']),
                            'ref_start': last_ext['ref_start'],
                            'ref_end': ext['ref_end']
                        }
                    elif ext['type'] == 'left_clip' and last_ext['type'] == 'overlap':
                        # Prepend left clip to overlap
                        merged_extensions[-1] = {
                            'sequence': ext['sequence'] + last_ext['sequence'],
                            'type': 'extended',
                            'priority': max(last_ext['priority'], ext['priority']),
                            'ref_start': ext['ref_start'],
                            'ref_end': last_ext['ref_end']
                        }
                else:
                    merged_extensions.append(ext)
            
            # 4. Build final ITD sequence
            final_itd_parts = []
            for ext in merged_extensions:
                final_itd_parts.append(ext['sequence'])
            
            final_itd_sequence = ''.join(final_itd_parts)
            
            # 5. Quality control and validation
            consensus_ratio = (consensus_count / total_bases) if total_bases > 0 else 0
            min_consensus_ratio = 0.3  # Relaxed threshold
            
            # Additional validation: check if sequence has reasonable complexity
            if (len(final_itd_sequence) >= self.min_itd_length and 
                consensus_ratio >= min_consensus_ratio and
                len(set(final_itd_sequence)) >= 3):  # At least 3 different bases
                
                # CRITICAL FIX: Limit ITD sequence length more aggressively
                if len(final_itd_sequence) > self.max_itd_length:
                    logger.debug(f"Truncating oversized ITD sequence: {len(final_itd_sequence)}bp → {self.max_itd_length}bp")
                    # Keep the middle portion which is most likely to be the actual ITD
                    excess = len(final_itd_sequence) - self.max_itd_length
                    start_trim = excess // 2
                    final_itd_sequence = final_itd_sequence[start_trim:start_trim + self.max_itd_length]
                
                logger.debug(f"Reconstructed ITD: {len(final_itd_sequence)}bp, "
                           f"consensus_ratio={consensus_ratio:.3f}, "
                           f"extensions={len(merged_extensions)}")
                
                return final_itd_sequence
            
            # Fallback: return the overlap consensus if full reconstruction fails
            elif len(overlap_consensus) >= self.min_itd_length and consensus_ratio >= min_consensus_ratio:
                
                # Also limit overlap consensus size
                if len(overlap_consensus) > self.max_itd_length:
                    logger.debug(f"Truncating oversized overlap consensus: {len(overlap_consensus)}bp → {self.max_itd_length}bp")
                    excess = len(overlap_consensus) - self.max_itd_length
                    start_trim = excess // 2
                    overlap_consensus = overlap_consensus[start_trim:start_trim + self.max_itd_length]
                
                logger.debug(f"Using overlap consensus: {len(overlap_consensus)}bp, "
                           f"consensus_ratio={consensus_ratio:.3f}")
                return overlap_consensus
                
        except Exception as e:
            logger.debug(f"Error extracting overlap sequence: {e}")
            return None
        
        return None

    def _find_itd_boundaries(self, r1: pysam.AlignedRead, r2: pysam.AlignedRead, overlap: Dict) -> Tuple[int, int]:
        """Determine the true ITD boundaries by analyzing soft-clips and insertions - CORRECTED VERSION"""
        
        # Start with overlap boundaries - but be more conservative about extensions
        itd_start = overlap['overlap_start']
        itd_end = overlap['overlap_end']
        
        # CRITICAL FIX: Limit the initial ITD size estimate
        initial_overlap_size = itd_end - itd_start
        if initial_overlap_size > 200:
            # The overlap is too large - likely includes flanking sequence
            # Reduce to a more reasonable ITD size centered on the overlap
            center = (itd_start + itd_end) // 2
            reasonable_size = min(150, initial_overlap_size // 2)  # Cap at 150bp
            itd_start = center - reasonable_size // 2
            itd_end = center + reasonable_size // 2
            logger.debug(f"Reduced oversized overlap from {initial_overlap_size}bp to {reasonable_size}bp")
        
        # Look for CIGAR insertions that might indicate the actual ITD size
        actual_itd_evidence = []
        
        for read in [r1, r2]:
            if not read.cigartuples:
                continue
            
            ref_pos = read.reference_start
            query_pos = 0
            
            for op, length in read.cigartuples:
                if op == 1 and 15 <= length <= 300:  # Reasonable ITD-sized insertion
                    # This insertion is within the expected ITD size range
                    if itd_start - 50 <= ref_pos <= itd_end + 50:  # Near our overlap region
                        actual_itd_evidence.append({
                            'position': ref_pos,
                            'size': length,
                            'type': 'insertion'
                        })
                        logger.debug(f"Found {length}bp insertion at position {ref_pos}")
                
                # Update reference position
                if op in [0, 2]:  # Operations that advance reference
                    ref_pos += length
                if op in [0, 1, 4]:  # Operations that advance query
                    query_pos += length
        
        # Use insertion evidence to refine ITD boundaries
        if actual_itd_evidence:
            # Use the largest insertion as the primary ITD size indicator
            best_evidence = max(actual_itd_evidence, key=lambda x: x['size'])
            
            # Center the ITD boundaries around this insertion
            itd_center = best_evidence['position']
            itd_size = best_evidence['size']
            
            # Add some padding but keep it reasonable
            padding = min(20, itd_size // 4)  # Max 20bp padding or 25% of ITD size
            
            itd_start = itd_center - padding
            itd_end = itd_center + itd_size + padding
            
            logger.debug(f"Refined ITD boundaries using {itd_size}bp insertion: {itd_start}-{itd_end}")
        
        # Conservative soft-clip extension (much smaller than before)
        r1_left_clip_len = 0
        r1_right_clip_len = 0
        r2_left_clip_len = 0
        r2_right_clip_len = 0
        
        if r1.cigartuples:
            if r1.cigartuples[0][0] == 4:  # Left soft clip
                r1_left_clip_len = r1.cigartuples[0][1]
            if r1.cigartuples[-1][0] == 4:  # Right soft clip
                r1_right_clip_len = r1.cigartuples[-1][1]
        
        if r2.cigartuples:
            if r2.cigartuples[0][0] == 4:  # Left soft clip
                r2_left_clip_len = r2.cigartuples[0][1]
            if r2.cigartuples[-1][0] == 4:  # Right soft clip
                r2_right_clip_len = r2.cigartuples[-1][1]
        
        # MUCH more conservative soft-clip extension
        max_clip_extension = 30  # Maximum 30bp extension from soft clips
        
        if r1_left_clip_len > 15:  # Only significant left clips
            extension = min(r1_left_clip_len, max_clip_extension)
            itd_start = min(itd_start, r1.reference_start - extension)
        if r2_left_clip_len > 15:  # Only significant left clips
            extension = min(r2_left_clip_len, max_clip_extension)
            itd_start = min(itd_start, r2.reference_start - extension)
        
        if r1_right_clip_len > 15:  # Only significant right clips
            extension = min(r1_right_clip_len, max_clip_extension)
            itd_end = max(itd_end, r1.reference_end + extension)
        if r2_right_clip_len > 15:  # Only significant right clips
            extension = min(r2_right_clip_len, max_clip_extension)
            itd_end = max(itd_end, r2.reference_end + extension)
        
        # Final safety check: ensure ITD size is reasonable
        final_itd_size = itd_end - itd_start
        if final_itd_size > self.max_itd_length:
            # If still too large, cap it at max_itd_length
            center = (itd_start + itd_end) // 2
            itd_start = center - self.max_itd_length // 2
            itd_end = center + self.max_itd_length // 2
            logger.debug(f"Capped ITD size at {self.max_itd_length}bp")
        
        return itd_start, itd_end

    def _extract_complete_itd_sequence(self, r1: pysam.AlignedRead, r2: pysam.AlignedRead, 
                                       overlap: Dict) -> Optional[str]:
        """Extract complete ITD sequence using boundary detection and comprehensive reconstruction"""
        
        # First, determine the true ITD boundaries
        itd_start, itd_end = self._find_itd_boundaries(r1, r2, overlap)
        
        logger.debug(f"ITD boundaries: {itd_start}-{itd_end} (extended from {overlap['overlap_start']}-{overlap['overlap_end']})")
        
        # Extract sequence using the enhanced method with extended boundaries
        extended_overlap = {
            'overlap_start': itd_start,
            'overlap_end': itd_end,
            'original_read': overlap['original_read']
        }
        
        return self._extract_overlap_sequence(r1, r2, extended_overlap)
        return None
    
    def process_soft_clips(self, soft_clips: List[Dict]) -> List[ITDCandidate]:
        """Generate ITD candidates from soft-clipped sequences (parallelized reference search)"""
        import concurrent.futures
        candidates = []
        # Group soft clips by sequence similarity (with fuzzy matching)
        clip_groups = self._cluster_soft_clips_fuzzy(soft_clips)

        def process_clip_group(args):
            representative_seq, clips = args
            if len(clips) >= self.min_support and len(representative_seq) >= self.min_itd_length:
                ref_matches = self._find_in_reference_fuzzy(representative_seq)
                if ref_matches:
                    best_match = ref_matches[0]
                    base_confidence = 0.7
                    support_bonus = min(0.2, (len(clips) - self.min_support) * 0.05)
                    match_penalty = (1 - best_match['similarity']) * 0.3
                    return ITDCandidate(
                        sequence=representative_seq,
                        length=len(representative_seq),
                        position=clips[0]['ref_position'],
                        support_type='softclip',
                        supporting_reads=[c['read_name'] for c in clips],
                        confidence=max(0.3, base_confidence + support_bonus - match_penalty),
                        duplication_start=best_match['start'],
                        duplication_end=best_match['end']
                    )
            return None

        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(process_clip_group, clip_groups.items()))
        candidates = [c for c in results if c is not None]
        return candidates
    
    def _cluster_soft_clips_fuzzy(self, soft_clips: List[Dict], similarity_threshold: float = 0.85) -> Dict[str, List[Dict]]:
        """Cluster soft-clipped sequences by similarity with fuzzy matching"""
        clusters = {}
        
        for clip in soft_clips:
            seq = clip['sequence']
            matched = False
            
            # Try to match with existing clusters
            for rep_seq in list(clusters.keys()):
                similarity = self._sequence_similarity(seq, rep_seq)
                if similarity >= similarity_threshold:
                    clusters[rep_seq].append(clip)
                    matched = True
                    break
            
            # Create new cluster if no match
            if not matched:
                clusters[seq] = [clip]
        
        return clusters
    
    def _sequence_similarity(self, seq1: str, seq2: str) -> float:
        """Calculate sequence similarity using Levenshtein distance (fast if available)"""
        if len(seq1) == 0 or len(seq2) == 0:
            return 0.0
        # Simple length-based check first
        len_diff = abs(len(seq1) - len(seq2))
        if len_diff > max(len(seq1), len(seq2)) * 0.2:  # More than 20% length difference
            return 0.0
        if _has_levenshtein:
            distance = Levenshtein.distance(seq1, seq2)
        else:
            # Fallback to slow pure Python implementation
            matrix = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
            for i in range(len(seq1) + 1):
                matrix[i][0] = i
            for j in range(len(seq2) + 1):
                matrix[0][j] = j
            for i in range(1, len(seq1) + 1):
                for j in range(1, len(seq2) + 1):
                    if seq1[i-1] == seq2[j-1]:
                        cost = 0
                    else:
                        cost = 1
                    matrix[i][j] = min(
                        matrix[i-1][j] + 1,      # deletion
                        matrix[i][j-1] + 1,      # insertion
                        matrix[i-1][j-1] + cost  # substitution
                    )
            distance = matrix[len(seq1)][len(seq2)]
        max_len = max(len(seq1), len(seq2))
        similarity = 1 - (distance / max_len)
        return similarity
    
    def _find_in_reference_fuzzy(self, sequence: str, max_mismatches: int = 2) -> List[Dict]:
        """Find sequence in reference with fuzzy matching, using k-mer prefilter for speed"""
        matches = []
        seq_len = len(sequence)
        k = 10 if seq_len >= 10 else seq_len  # Use k-mer size up to 10
        if seq_len < k:
            # Fallback to original method for very short queries
            for i in range(len(self.reference_sequence) - seq_len + 1):
                ref_substr = self.reference_sequence[i:i + seq_len]
                mismatches = sum(1 for a, b in zip(sequence, ref_substr) if a != b)
                if mismatches <= max_mismatches:
                    similarity = 1 - (mismatches / seq_len)
                    matches.append({
                        'start': i,
                        'end': i + seq_len,
                        'sequence': ref_substr,
                        'mismatches': mismatches,
                        'similarity': similarity
                    })
            matches.sort(key=lambda x: x['similarity'], reverse=True)
            return matches

        # Use k-mer prefilter
        query_kmer = sequence[:k]
        candidate_starts = self.ref_kmers.get(query_kmer, [])
        for i in candidate_starts:
            # Only consider windows where the k-mer matches
            if i + seq_len > len(self.reference_sequence):
                continue
            ref_substr = self.reference_sequence[i:i + seq_len]
            mismatches = sum(1 for a, b in zip(sequence, ref_substr) if a != b)
            if mismatches <= max_mismatches:
                similarity = 1 - (mismatches / seq_len)
                matches.append({
                    'start': i,
                    'end': i + seq_len,
                    'sequence': ref_substr,
                    'mismatches': mismatches,
                    'similarity': similarity
                })
        matches.sort(key=lambda x: x['similarity'], reverse=True)
        return matches
    
    def merge_candidates(self, candidates: List[ITDCandidate]) -> List[ITDCandidate]:
        """Merge similar ITD candidates more aggressively"""
        if not candidates:
            return []
        
        # Sort by position and length
        candidates.sort(key=lambda x: (x.position, x.length))
        
        # Group candidates by similar characteristics
        merged_groups = []
        used = set()
        
        for i, cand1 in enumerate(candidates):
            if i in used:
                continue
            
            # Start a new group
            group = [cand1]
            used.add(i)
            
            # Find all similar candidates
            for j, cand2 in enumerate(candidates[i+1:], i+1):
                if j in used:
                    continue
                
                # Check if candidates are similar enough to merge
                if self._are_similar_aggressive(cand1, cand2):
                    group.append(cand2)
                    used.add(j)
            
            merged_groups.append(group)
        
        # Merge each group into a single candidate
        merged = []
        for group in merged_groups:
            if len(group) == 1:
                merged.append(group[0])
            else:
                merged_candidate = self._merge_cluster(group)
                # Only keep high-confidence merged candidates
                if merged_candidate.confidence >= 0.5 or len(merged_candidate.supporting_reads) >= self.min_support * 2:
                    merged.append(merged_candidate)
        
        # Final deduplication based on position and length
        final_candidates = []
        seen_positions = {}
        
        for candidate in sorted(merged, key=lambda x: x.confidence, reverse=True):
            key = (candidate.position // 5, candidate.length // 5)  # Group by 5bp bins
            if key not in seen_positions:
                seen_positions[key] = candidate
                final_candidates.append(candidate)
            else:
                # Merge supporting reads if very similar
                existing = seen_positions[key]
                if abs(existing.position - candidate.position) <= 2 and abs(existing.length - candidate.length) <= 2:
                    existing.supporting_reads.extend(candidate.supporting_reads)
                    existing.supporting_reads = list(set(existing.supporting_reads))
                    existing.confidence = max(existing.confidence, candidate.confidence)
        
        return final_candidates
    
    def _are_similar_aggressive(self, cand1: ITDCandidate, cand2: ITDCandidate) -> bool:
        """Check if two candidates are similar enough to merge (more aggressive)"""
        # Position proximity (within 10bp)
        if abs(cand1.position - cand2.position) > 10:
            return False
        
        # Length similarity (within 5bp or 10%)
        length_diff = abs(cand1.length - cand2.length)
        if length_diff > 5 and length_diff > max(cand1.length, cand2.length) * 0.1:
            return False
        
        # Sequence similarity if both have sequences
        if cand1.sequence and cand2.sequence:
            similarity = self._sequence_similarity(cand1.sequence, cand2.sequence)
            if similarity < 0.8:  # 80% similarity threshold
                return False
        
        return True
    
    def _merge_cluster(self, cluster: List[ITDCandidate]) -> ITDCandidate:
        """Merge a cluster of similar candidates"""
        # Use the most common sequence
        sequences = [c.sequence for c in cluster]
        seq_counts = Counter(sequences)
        best_seq = seq_counts.most_common(1)[0][0]
        
        # Combine supporting reads
        all_reads = []
        for c in cluster:
            all_reads.extend(c.supporting_reads)
        unique_reads = list(set(all_reads))
        
        # Average position
        avg_position = sum(c.position for c in cluster) // len(cluster)
        
        # Combined confidence based on support
        max_confidence = max(c.confidence for c in cluster)
        support_bonus = min(0.3, len(unique_reads) * 0.015)  # Strong bonus for high support
        cluster_bonus = min(0.1, len(cluster) * 0.02)  # Bonus for multiple similar candidates
        
        return ITDCandidate(
            sequence=best_seq,
            length=len(best_seq),
            position=avg_position,
            support_type='merged',
            supporting_reads=unique_reads,
            confidence=min(0.95, max_confidence + support_bonus + cluster_bonus),
            insertion_site=cluster[0].insertion_site,
            duplication_start=cluster[0].duplication_start,
            duplication_end=cluster[0].duplication_end
        )
    
    def generate_candidates(self, detection_results: Dict) -> List[ITDCandidate]:
        """Generate ITD candidates from detection results"""
        candidates = []
        
        # Process overlaps
        if 'overlaps' in detection_results and detection_results['overlaps']:
            overlap_candidates = self.extract_overlap_sequences(
                detection_results['overlaps'],
                detection_results['bam_file']
            )
            candidates.extend(overlap_candidates)
            logger.info(f"Generated {len(overlap_candidates)} candidates from overlaps")
        
        # Process soft clips
        if 'soft_clips' in detection_results and detection_results['soft_clips']:
            clip_candidates = self.process_soft_clips(detection_results['soft_clips'])
            candidates.extend(clip_candidates)
            logger.info(f"Generated {len(clip_candidates)} candidates from soft clips")
        
        # Correct ITD sizes for overlap-based candidates
        try:
            from itd_size_corrector import correct_itd_candidates_sizes
            logger.info("Applying ITD size correction to overlap-based candidates...")
            overlap_corrected = [c for c in candidates if c.support_type == 'overlap']
            if overlap_corrected:
                corrected_overlaps = correct_itd_candidates_sizes(
                    overlap_corrected, 
                    self.reference_sequence, 
                    detection_results['bam_file']
                )
                # Replace overlap candidates with corrected versions
                candidates = [c for c in candidates if c.support_type != 'overlap']
                candidates.extend(corrected_overlaps)
                logger.info(f"Size correction applied to {len(overlap_corrected)} overlap candidates")
        except Exception as e:
            logger.warning(f"ITD size correction failed: {e} - using original size estimates")
        
        # Merge similar candidates (both existing method and new biological event merging)
        merged_candidates = self.merge_candidates(candidates)
        
        # Apply additional biological event merging to prevent fragmentation
        biologically_merged = self._merge_similar_candidates(merged_candidates)
        
        # Filter by support
        filtered_candidates = [
            c for c in biologically_merged 
            if len(c.supporting_reads) >= self.min_support
        ]
        
        logger.info(f"ITD candidate pipeline: {len(candidates)} → {len(merged_candidates)} → {len(biologically_merged)} → {len(filtered_candidates)}")
        logger.info(f"Final ITD candidates: {len(filtered_candidates)}")
        
        logger.info("=== DEBUG: Checking ITD candidates ===")
        for i, cand in enumerate(filtered_candidates[:3]):
            logger.info(f"Candidate {i+1}: {cand.length}bp, pos={cand.position}, "
                       f"support={len(cand.supporting_reads)}, confidence={cand.confidence:.3f}")
            logger.info(f"  Sequence (first 50bp): '{cand.sequence[:50]}...'")
            if hasattr(cand, 'duplication_start') and cand.duplication_start:
                logger.info(f"  Duplication region: {cand.duplication_start}-{cand.duplication_end}")
        
        if filtered_candidates:
            logger.info(f"Top 3 ITD sequences found:")
            for i, cand in enumerate(filtered_candidates[:3]):
                logger.info(f"  {i+1}. {cand.length}bp at pos {cand.position}: '{cand.sequence[:50]}...'")
                logger.info(f"     Support: {len(cand.supporting_reads)} reads, Confidence: {cand.confidence:.3f}")
                logger.info(f"     Support type: {cand.support_type}")
        else:
            logger.info("No ITD candidates passed filtering criteria")

        return filtered_candidates

def generate_itd_candidates(detection_results: Dict, reference_sequence: str,
                           min_itd_length: int = 15, min_support: int = 2, 
                           max_itd_length: int = 500) -> List[ITDCandidate]:
    """Main function to generate ITD candidates"""
    generator = ITDGenerator(reference_sequence, min_itd_length, min_support, max_itd_length)
    return generator.generate_candidates(detection_results)
        

    
