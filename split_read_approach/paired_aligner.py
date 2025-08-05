#!/usr/bin/env python3
"""
Paired Aligner Module
Aligns artificial paired reads and identifies ITD candidates from overlaps and soft clips
"""

import logging
import subprocess
import tempfile
import pysam
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
import re

logger = logging.getLogger(__name__)

class PairedAligner:
    """Align artificial paired reads and detect ITDs"""
    
    def __init__(self, reference_file: str, min_overlap: int = 15, min_softclip: int = 15):
        self.reference_file = reference_file
        self.min_overlap = min_overlap
        self.min_softclip = min_softclip
        self.temp_dir = tempfile.mkdtemp(prefix="flt3_itd_paired_")
        
    def write_pairs_to_fastq(self, pairs: List[Tuple[Dict, Dict]]) -> Tuple[str, str]:
        """Write paired reads to FASTQ files"""
        r1_file = Path(self.temp_dir) / "reads_R1.fastq"
        r2_file = Path(self.temp_dir) / "reads_R2.fastq"
        
        with open(r1_file, 'w') as f1, open(r2_file, 'w') as f2:
            for r1, r2 in pairs:
                # Write R1
                f1.write(f"@{r1['name']}\n")
                f1.write(f"{r1['sequence']}\n")
                f1.write("+\n")
                if r1['qualities']:
                    f1.write(''.join(chr(q + 33) for q in r1['qualities']) + "\n")
                else:
                    f1.write('I' * len(r1['sequence']) + "\n")  # Default quality
                
                # Write R2
                f2.write(f"@{r2['name']}\n")
                f2.write(f"{r2['sequence']}\n")
                f2.write("+\n")
                if r2['qualities']:
                    f2.write(''.join(chr(q + 33) for q in r2['qualities']) + "\n")
                else:
                    f2.write('I' * len(r2['sequence']) + "\n")
        
        return str(r1_file), str(r2_file)
    
    def align_pairs(self, r1_file: str, r2_file: str) -> str:
        """Align paired reads using minimap2"""
        output_sam = Path(self.temp_dir) / "aligned_pairs.sam"
        
        # Run minimap2 for paired-end alignment
        cmd = [
            'minimap2',
            '-ax', 'sr',  # Short read paired-end mode
            '-t', '4',    # Number of threads
            '--secondary=no',  # No secondary alignments
            '-k', '15',   # K-mer size
            self.reference_file,
            r1_file,
            r2_file
        ]
        
        with open(output_sam, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
            if result.returncode != 0:
                raise RuntimeError(f"Minimap2 failed: {result.stderr.decode()}")
        
        # Convert to sorted BAM
        output_bam = Path(self.temp_dir) / "aligned_pairs.bam"
        pysam.sort("-o", str(output_bam), str(output_sam))
        pysam.index(str(output_bam))
        
        return str(output_bam)
    
    def extract_soft_clips(self, read: pysam.AlignedRead) -> List[Dict]:
        """Extract soft-clipped sequences from read"""
        clips = []
        
        if read.cigartuples is None:
            return clips
        
        pos = 0
        ref_pos = read.reference_start
        
        for op, length in read.cigartuples:
            if op == 4:  # Soft clip
                if length >= self.min_softclip:
                    clip_seq = read.query_sequence[pos:pos+length]
                    clip_info = {
                        'sequence': clip_seq,
                        'length': length,
                        'position': 'start' if pos == 0 else 'end',
                        'ref_position': ref_pos,
                        'read_name': read.query_name
                    }
                    clips.append(clip_info)
            
            if op in [0, 2, 3, 7, 8]:  # Consumes reference
                ref_pos += length
            if op in [0, 1, 4, 7, 8]:  # Consumes query
                pos += length
        
        return clips
    
    def find_overlapping_pairs(self, bam_file: str) -> List[Dict]:
        """Find overlapping paired reads that might indicate ITDs"""
        overlaps = []
        
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Group reads by original read name
            read_groups = defaultdict(list)
            
            for read in bam:
                if read.is_unmapped:
                    continue
                    
                # Extract original read name
                original_name = read.query_name.rsplit('_R', 1)[0]
                read_groups[original_name].append(read)
            
            # Analyze each group
            for original_name, reads in read_groups.items():
                # Find R1 and R2 pairs
                r1_reads = [r for r in reads if '_R1' in r.query_name]
                r2_reads = [r for r in reads if '_R2' in r.query_name]
                
                for r1 in r1_reads:
                    for r2 in r2_reads:
                        # Check if they're from the same split
                        if not self._same_split(r1.query_name, r2.query_name):
                            continue
                        
                        # Check for overlap
                        overlap = self._calculate_overlap(r1, r2)
                        if overlap and overlap['length'] >= self.min_overlap:
                            overlap['original_read'] = original_name
                            overlaps.append(overlap)
        
        return overlaps
    
    def _same_split(self, name1: str, name2: str) -> bool:
        """Check if two reads are from the same split"""
        # Extract split indices if present
        match1 = re.search(r'_R[12]_(\d+)', name1)
        match2 = re.search(r'_R[12]_(\d+)', name2)
        
        if match1 and match2:
            return match1.group(1) == match2.group(1)
        
        # For simple splits without index
        base1 = name1.rsplit('_R', 1)[0]
        base2 = name2.rsplit('_R', 1)[0]
        return base1 == base2
    
    def _calculate_overlap(self, r1: pysam.AlignedRead, r2: pysam.AlignedRead) -> Optional[Dict]:
        """Calculate overlap between paired reads"""
        if r1.reference_start is None or r1.reference_end is None:
            return None
        if r2.reference_start is None or r2.reference_end is None:
            return None
        
        # Check if reads overlap
        overlap_start = max(r1.reference_start, r2.reference_start)
        overlap_end = min(r1.reference_end, r2.reference_end)
        
        if overlap_start < overlap_end:
            # Extract overlapping sequence
            # This is simplified - would need proper CIGAR parsing for accuracy
            overlap_length = overlap_end - overlap_start
            
            # Filter out small overlaps that are likely from normal alignment
            if overlap_length < self.min_overlap:
                return None
            
            # For split reads from the same original read, we expect specific patterns
            # Filter overlaps that don't match expected ITD patterns
            r1_aligned_length = r1.reference_end - r1.reference_start
            r2_aligned_length = r2.reference_end - r2.reference_start
            
            # If both reads align to similar lengths, it's likely not an ITD
            if abs(r1_aligned_length - r2_aligned_length) < 10:
                return None
            
            return {
                'r1_name': r1.query_name,
                'r2_name': r2.query_name,
                'overlap_start': overlap_start,
                'overlap_end': overlap_end,
                'length': overlap_length,
                'r1_start': r1.reference_start,
                'r1_end': r1.reference_end,
                'r2_start': r2.reference_start,
                'r2_end': r2.reference_end,
                'r1_aligned_length': r1_aligned_length,
                'r2_aligned_length': r2_aligned_length
            }
        
        return None
    
    def detect_itd_candidates(self, pairs: List[Tuple[Dict, Dict]]) -> Dict:
        """Main function to detect ITD candidates from paired reads"""
        logger.info("Starting paired alignment for ITD detection")
        
        # Write pairs to FASTQ
        r1_file, r2_file = self.write_pairs_to_fastq(pairs)
        
        # Align pairs
        bam_file = self.align_pairs(r1_file, r2_file)
        
        # Find overlaps
        overlaps = self.find_overlapping_pairs(bam_file)
        logger.info(f"Found {len(overlaps)} overlapping pairs")
        
        # Extract soft clips
        soft_clips = []
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                clips = self.extract_soft_clips(read)
                soft_clips.extend(clips)
        
        logger.info(f"Found {len(soft_clips)} soft-clipped sequences")
        
        # Combine results
        results = {
            'overlaps': overlaps,
            'soft_clips': soft_clips,
            'bam_file': bam_file,
            'temp_dir': self.temp_dir  # Add temp directory for cleanup tracking
        }
        
        return results
    
    def cleanup(self):
        """Clean up temporary files"""
        import shutil
        if Path(self.temp_dir).exists():
            shutil.rmtree(self.temp_dir)

def detect_itds_from_pairs(pairs: List[Tuple[Dict, Dict]], reference_file: str,
                          min_overlap: int = 15, min_softclip: int = 15) -> Dict:
    """Main function to detect ITDs from paired reads"""
    aligner = PairedAligner(reference_file, min_overlap, min_softclip)
    results = aligner.detect_itd_candidates(pairs)
    # Do not clean up here; cleanup should be handled after downstream processing is complete
    return results