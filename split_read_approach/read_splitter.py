#!/usr/bin/env python3
"""
Read Splitter Module
Splits long nanopore reads into artificial paired reads for ITD detection
"""

import logging
from typing import List, Dict, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class SplitConfig:
    """Configuration for read splitting"""
    split_method: str = "multiple"  # "center", "multiple", or "sliding"
    overlap: int = 20  # bp overlap between pairs
    min_split_size: int = 100  # Minimum size for each split
    num_splits: int = 3  # Number of split points for "multiple" method
    sliding_step: int = 50  # Step size for sliding window method

class ReadSplitter:
    """Split long reads into artificial paired reads"""
    
    def __init__(self, config: SplitConfig = None):
        self.config = config or SplitConfig()
        
    def split_read_center(self, read: Dict) -> List[Tuple[Dict, Dict]]:
        """Split read at center with overlap"""
        sequence = read['sequence']
        qualities = read['qualities']
        length = len(sequence)
        
        # Calculate split point
        mid_point = length // 2
        overlap_half = self.config.overlap // 2
        
        # Check if splits would be too small
        if mid_point - overlap_half < self.config.min_split_size:
            logger.debug(f"Read {read['name']} too short for center split")
            return []
        
        # Create R1 (first half + overlap)
        r1_end = mid_point + overlap_half
        r1 = {
            'name': f"{read['name']}_R1",
            'sequence': sequence[:r1_end],
            'qualities': qualities[:r1_end] if qualities else None,
            'original_name': read['name'],
            'split_type': 'R1',
            'split_method': 'center',
            'split_position': 0,
            'original_length': length
        }
        
        # Create R2 (second half + overlap)
        r2_start = mid_point - overlap_half
        r2 = {
            'name': f"{read['name']}_R2", 
            'sequence': sequence[r2_start:],
            'qualities': qualities[r2_start:] if qualities else None,
            'original_name': read['name'],
            'split_type': 'R2',
            'split_method': 'center',
            'split_position': r2_start,
            'original_length': length
        }
        
        return [(r1, r2)]
    
    def split_read_multiple(self, read: Dict) -> List[Tuple[Dict, Dict]]:
        """Split read at multiple points"""
        sequence = read['sequence']
        qualities = read['qualities']
        length = len(sequence)
        
        # Skip reads that are too short
        if length < 2 * self.config.min_split_size + self.config.overlap:
            return []
        
        pairs = []
        
        # Calculate split points based on read length
        if length < 400:  # Short reads - fewer splits
            split_points = [length // 2]  # Just center split
        elif self.config.num_splits == 3 and length >= 600:
            # 25%, 50%, 75% split points for longer reads
            split_points = [length // 4, length // 2, 3 * length // 4]
        else:
            # Adjusted split points to avoid too many pairs
            if length < 600:
                num_splits = min(2, self.config.num_splits)
            else:
                num_splits = self.config.num_splits
            
            step = length // (num_splits + 1)
            split_points = [i * step for i in range(1, num_splits + 1)]
        
        for i, split_point in enumerate(split_points):
            overlap_half = self.config.overlap // 2
            
            # Check minimum size constraints
            if (split_point - overlap_half < self.config.min_split_size or 
                length - split_point + overlap_half < self.config.min_split_size):
                continue
            
            # Create R1
            r1_end = min(split_point + overlap_half, length)
            r1 = {
                'name': f"{read['name']}_R1_{i}",
                'sequence': sequence[:r1_end],
                'qualities': qualities[:r1_end] if qualities else None,
                'original_name': read['name'],
                'split_type': 'R1',
                'split_method': 'multiple',
                'split_index': i,
                'split_position': 0,
                'original_length': length
            }
            
            # Create R2
            r2_start = max(0, split_point - overlap_half)
            r2 = {
                'name': f"{read['name']}_R2_{i}",
                'sequence': sequence[r2_start:],
                'qualities': qualities[r2_start:] if qualities else None,
                'original_name': read['name'],
                'split_type': 'R2',
                'split_method': 'multiple',
                'split_index': i,
                'split_position': r2_start,
                'original_length': length
            }
            
            pairs.append((r1, r2))
        
        return pairs
    
    def split_read_sliding(self, read: Dict) -> List[Tuple[Dict, Dict]]:
        """Split read using sliding window approach"""
        sequence = read['sequence']
        qualities = read['qualities']
        length = len(sequence)
        
        pairs = []
        
        # Slide through the read
        for start in range(0, length - 2 * self.config.min_split_size, self.config.sliding_step):
            # Calculate split point
            remaining = length - start
            split_point = start + remaining // 2
            overlap_half = self.config.overlap // 2
            
            # Check constraints
            if (split_point - start - overlap_half < self.config.min_split_size or
                length - split_point + overlap_half < self.config.min_split_size):
                continue
            
            # Create R1
            r1_start = start
            r1_end = min(split_point + overlap_half, length)
            r1 = {
                'name': f"{read['name']}_R1_s{start}",
                'sequence': sequence[r1_start:r1_end],
                'qualities': qualities[r1_start:r1_end] if qualities else None,
                'original_name': read['name'],
                'split_type': 'R1',
                'split_method': 'sliding',
                'split_position': r1_start,
                'window_start': start,
                'original_length': length
            }
            
            # Create R2
            r2_start = max(start, split_point - overlap_half)
            r2 = {
                'name': f"{read['name']}_R2_s{start}",
                'sequence': sequence[r2_start:],
                'qualities': qualities[r2_start:] if qualities else None,
                'original_name': read['name'],
                'split_type': 'R2',
                'split_method': 'sliding',
                'split_position': r2_start,
                'window_start': start,
                'original_length': length
            }
            
            pairs.append((r1, r2))
            
            # Limit number of windows per read
            if len(pairs) >= 10:
                break
        
        return pairs
    
    def split_reads(self, reads: List[Dict]) -> List[Tuple[Dict, Dict]]:
        """Split all reads using configured method"""
        all_pairs = []
        
        for read in reads:
            if self.config.split_method == "center":
                pairs = self.split_read_center(read)
            elif self.config.split_method == "multiple":
                pairs = self.split_read_multiple(read)
            elif self.config.split_method == "sliding":
                pairs = self.split_read_sliding(read)
            else:
                raise ValueError(f"Unknown split method: {self.config.split_method}")
            
            all_pairs.extend(pairs)
        
        logger.info(f"Generated {len(all_pairs)} read pairs from {len(reads)} reads using {self.config.split_method} method")
        
        # Log statistics
        if all_pairs:
            avg_r1_length = sum(len(r1['sequence']) for r1, _ in all_pairs) / len(all_pairs)
            avg_r2_length = sum(len(r2['sequence']) for _, r2 in all_pairs) / len(all_pairs)
            logger.info(f"Average lengths: R1={avg_r1_length:.1f}bp, R2={avg_r2_length:.1f}bp")
        
        return all_pairs

def split_reads(reads: List[Dict], method: str = "multiple", 
                overlap: int = 20, min_split_size: int = 100) -> List[Tuple[Dict, Dict]]:
    """Main function to split reads into artificial pairs"""
    config = SplitConfig(
        split_method=method,
        overlap=overlap,
        min_split_size=min_split_size
    )
    splitter = ReadSplitter(config)
    return splitter.split_reads(reads)