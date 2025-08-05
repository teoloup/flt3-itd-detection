#!/usr/bin/env python3
"""
BAM Read Extractor Module
Extracts reads mapping to FLT3 region and performs primer trimming
"""

import logging
import pysam
from pathlib import Path
from typing import List, Dict, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)


class FLT3ReadExtractor:
    """Extract and process reads from BAM file"""
    
    def __init__(self, bam_file: str, genome_build: str = "hg38", min_mapping_quality: int = 20, config=None):
        self.bam_file = bam_file
        self.genome_build = genome_build
        self.min_mapping_quality = min_mapping_quality
        if config is None:
            logger.warning("No config object provided to FLT3ReadExtractor. Falling back to minimal default settings. This may cause resource or parameter mismatches.")
            # Minimal fallback config with required attributes
            class MinimalConfig:
                threads = 4
                amplicon_length = 336
                forward_primer = "AGCAATTTAGGTATGAAAGCCAGC"
                reverse_primer = "CTGTACCTTTCAGCATTTTGACG" 
                start_hg38 = 28033301
                end_hg38 = 28034800
                start_hg19 = 28607438
                end_hg19 = 28608937
            self.config = MinimalConfig()
        else:
            self.config = config
        self.setup_coordinates()
        
    def setup_coordinates(self):
        """Set up genomic coordinates based on genome build"""
        if self.genome_build == "hg38":
            self.start = self.config.start_hg38
            self.end = self.config.end_hg38
        else:  # hg19
            self.start = self.config.start_hg19
            self.end = self.config.end_hg19
            
        # Check chromosome naming in BAM
        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            references = bam.references
            if "chr13" in references:
                self.chromosome = "chr13"
            elif "13" in references:
                self.chromosome = "13"
            else:
                raise ValueError("Cannot find chromosome 13 in BAM file")
                
        self.region = f"{self.chromosome}:{self.start}-{self.end}"
        logger.info(f"Using region: {self.region}")
    
    def extract_reads(self) -> List[Dict]:
        """Extract reads mapping to FLT3 region using samtools view + fastq, then parse to dicts"""
        import tempfile
        import subprocess
        import os
        import shutil

        temp_dir = tempfile.mkdtemp(prefix="flt3_samtools_")
        fastq_out = os.path.join(temp_dir, "region_reads.fastq")
        region = f"{self.chromosome}:{self.start}-{self.end}"
        # Compose samtools view (region + quality filter) piped to samtools fastq
        samtools_view_cmd = [
            "samtools", "view",
            "-h",  # include header so samtools fastq works
            "-F", "4",  # skip unmapped, secondary, supplementary
            "-q", str(self.min_mapping_quality),
            self.bam_file,
            region
        ]
        # Use shell redirection for samtools fastq output
        fastq_cmd_str = f"samtools fastq - > '{fastq_out}'"
        logger.info(f"Running samtools view + fastq for region: {' '.join(samtools_view_cmd)} | {fastq_cmd_str}")
        view_proc = subprocess.Popen(samtools_view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        fastq_proc = subprocess.Popen(fastq_cmd_str, shell=True, stdin=view_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, fastq_stderr = fastq_proc.communicate()
        view_proc.stdout.close()
        view_stderr = view_proc.stderr.read()
        view_proc.stderr.close()
        view_proc.wait()
        if fastq_proc.returncode != 0 or view_proc.returncode != 0:
            logger.error(f"Samtools view/fastq failed: {view_stderr.decode()} {fastq_stderr.decode()}")
            shutil.rmtree(temp_dir)
            raise RuntimeError("Samtools view/fastq failed for region extraction.")

        # Parse FASTQ to read dicts
        reads = []
        with open(fastq_out, 'r') as f:
            while True:
                name = f.readline().strip()
                if not name:
                    break
                seq = f.readline().strip()
                plus = f.readline()
                qual = f.readline().strip()
                if name.startswith('@'):
                    read_name = name[1:]
                else:
                    read_name = name
                reads.append({
                    'name': read_name,
                    'sequence': seq,
                    'qualities': [ord(q) - 33 for q in qual] if qual else None,
                    'length': len(seq)
                })
        logger.info(f"Extracted {len(reads)} reads from FLT3 region using samtools view + fastq")
        shutil.rmtree(temp_dir)
        return reads
    
    def find_primer_positions(self, sequence: str) -> Tuple[int, int]:
        """Find primer positions in sequence allowing for mismatches"""
        forward_pos = -1
        reverse_pos = -1
        
        # Forward primer search (allow up to 2 mismatches)
        for i in range(len(sequence) - len(self.config.forward_primer) + 1):
            mismatches = 0
            for j in range(len(self.config.forward_primer)):
                if sequence[i+j] != self.config.forward_primer[j]:
                    mismatches += 1
                    if mismatches > 2:
                        break
            if mismatches <= 2:
                forward_pos = i
                break
        
        # Reverse primer search (search for reverse complement)
        rev_primer_rc = self.reverse_complement(self.config.reverse_primer)
        for i in range(len(sequence) - len(rev_primer_rc) + 1):
            mismatches = 0
            for j in range(len(rev_primer_rc)):
                if sequence[i+j] != rev_primer_rc[j]:
                    mismatches += 1
                    if mismatches > 2:
                        break
            if mismatches <= 2:
                reverse_pos = i + len(rev_primer_rc)
                break
                
        return forward_pos, reverse_pos
    
    def reverse_complement(self, seq: str) -> str:
        """Get reverse complement of sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                     'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                     'N': 'N', 'n': 'n'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    

    def trim_primers(self, reads: List[Dict], error_rate: float = None) -> Tuple[List[Dict], str]:
        """Trim primers using cutadapt with linked adapters and Nanopore error tolerance. Always keep trimmed FASTQ file for validation."""
        import tempfile
        import subprocess
        import os
        import shutil

        # Use error_rate from config if available, else default to 0.20
        if error_rate is None:
            error_rate = getattr(self.config, 'cutadapt_error_rate', 0.20)
        threads = self.config.threads
        amplicon_length = self.config.amplicon_length
        fwd = self.config.forward_primer
        rev = self.reverse_complement(self.config.reverse_primer)
        logger.debug(f"Using {threads} threads for cutadapt (from config)")

        # Prepare temporary files
        temp_dir = tempfile.mkdtemp(prefix="flt3_cutadapt_")
        fastq_in = os.path.join(temp_dir, "input.fastq")
        fastq_out = os.path.join(temp_dir, "trimmed.fastq")

        # Write reads to FASTQ
        with open(fastq_in, 'w') as f:
            for read in reads:
                f.write(f"@{read['name']}\n{read['sequence']}\n+\n")
                if read['qualities']:
                    f.write(''.join(chr(q + 33) for q in read['qualities']) + "\n")
                else:
                    f.write('I' * len(read['sequence']) + "\n")

        # Prepare cutadapt command
        # Using linked adapters to find and trim everything OUTSIDE the amplicon
        # This keeps the amplicon (including primers) and removes flanking sequences
        cutadapt_cmd = [
            "cutadapt",
            "-g", f"{fwd}...{rev}",  # linked adapters: find forward then reverse
            "-e", str(error_rate),
            "-j", str(threads),  # number of parallel jobs
            "-m", str(amplicon_length - 30),  # min length after trimming
            "--rc",  # reverse complement search also
            "--rename", "{header}",    # keep original read names
            "--discard-untrimmed",  # discard reads without both primers
            "--action=retain",  # keep primer sequences
            "-o", fastq_out,
            fastq_in
        ]

        logger.info(f"Running cutadapt for primer trimming: {' '.join(cutadapt_cmd)}")
        logger.debug(f"Full cutadapt command: {cutadapt_cmd}")  # Debug the actual command list
        result = subprocess.run(cutadapt_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            logger.error(f"Cutadapt failed: {result.stderr.decode()}")
            shutil.rmtree(temp_dir)
            raise RuntimeError("Cutadapt failed for primer trimming.")

        # Parse trimmed FASTQ
        trimmed_reads = []
        with open(fastq_out, 'r') as f:
            while True:
                name = f.readline().strip()
                if not name:
                    break
                seq = f.readline().strip()
                plus = f.readline()
                qual = f.readline().strip()
                if name.startswith('@'):
                    read_name = name[1:]
                else:
                    read_name = name
                trimmed_reads.append({
                    'name': read_name,
                    'sequence': seq,
                    'qualities': [ord(q) - 33 for q in qual] if qual else None,
                    'length': len(seq),
                    'trimmed': True
                })

        logger.info(f"Cutadapt trimmed {len(trimmed_reads)} reads (input: {len(reads)})")
        # Do not delete temp_dir here; FASTQ is needed for validation
        return trimmed_reads, fastq_out
    
    def process_reads(self, min_length: int = 250) -> Tuple[List[Dict], str]:
        """Extract reads and trim primers using cutadapt. Always keep trimmed FASTQ file for validation."""
        reads = self.extract_reads()
        trimmed_reads, fastq_out = self.trim_primers(reads)
        # Filter by minimum length after trimming
        processed_reads = [r for r in trimmed_reads if r['length'] >= min_length]
        logger.info(f"Processed {len(processed_reads)} reads after cutadapt trimming and length filtering")
        # Log statistics
        trimmed_count = sum(1 for r in processed_reads if r.get('trimmed', False))
        avg_length = sum(r['length'] for r in processed_reads) / len(processed_reads) if processed_reads else 0
        logger.info(f"Trimming statistics: {trimmed_count}/{len(processed_reads)} reads trimmed")
        logger.info(f"Average read length after processing: {avg_length:.1f} bp")
        return processed_reads, fastq_out

def extract_flt3_reads(bam_file: str, genome_build: str = "hg38", 
                      min_mapping_quality: int = 20, min_length: int = 250, config=None) -> Tuple[List[Dict], str]:
    """Main function to extract and process FLT3 reads, using main config if provided.
    Returns processed_reads and the trimmed FASTQ file path for centralized cleanup."""
    extractor = FLT3ReadExtractor(bam_file, genome_build, min_mapping_quality, config=config)
    processed_reads, fastq_file = extractor.process_reads(min_length)
    return processed_reads, fastq_file