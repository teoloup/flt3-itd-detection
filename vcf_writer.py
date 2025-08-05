#!/usr/bin/env python3
"""
VCF Writer Module
Writes ITD detection results to VCF format
"""

import logging
from datetime import datetime
from typing import List, Dict
from pathlib import Path

logger = logging.getLogger(__name__)

class VCFWriter:
    """Write ITD results to VCF format"""
    
    def __init__(self, reference_name: str = "FLT3", sample_name: str = "Sample",
                 genome_chr: str = "chr13", genome_start: int = 28033881):
        self.reference_name = reference_name
        self.sample_name = sample_name
        self.genome_chr = genome_chr
        self.genome_start = genome_start  # Genome coordinate offset
        
    def write_header(self, f, reference_length: int, command_line: str = ""):
        """Write VCF header"""
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        f.write(f"##source=FLT3_ITD_Detector_v1.0\n")
        f.write(f"##reference={self.genome_chr}\n")
        f.write(f"##contig=<ID={self.genome_chr},length=248956422>\n")  # chr13 length
        
        # INFO fields
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        f.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the ITD">\n')
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
        f.write('##INFO=<ID=DUPSEQ,Number=1,Type=String,Description="Duplicated sequence">\n')
        f.write('##INFO=<ID=DUPSTART,Number=1,Type=Integer,Description="Start position of duplicated region">\n')
        f.write('##INFO=<ID=DUPEND,Number=1,Type=Integer,Description="End position of duplicated region">\n')
        f.write('##INFO=<ID=CONFIDENCE,Number=1,Type=Float,Description="Confidence score">\n')
        f.write('##INFO=<ID=PRIMARY,Number=0,Type=Flag,Description="Primary ITD (matches size analysis)">\n')
        
        # FORMAT fields
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for ref and alt alleles">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">\n')
        f.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
        f.write('##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of supporting reads">\n')
        
        # Column headers
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}\n")
    
    def format_itd_variant(self, validation_result: 'ValidationResult', 
                          reference_sequence: str, idx: int) -> str:
        """Format a single ITD variant for VCF"""
        # For VCF insertion format, we need the position BEFORE the insertion
        local_insertion_pos = validation_result.insertion_position
        
        # VCF position is the base before the insertion (0-based to 1-based conversion)
        vcf_local_pos = max(0, local_insertion_pos - 1)
        genome_pos = self.genome_start + vcf_local_pos
        
        # Get reference base at the position before insertion
        if vcf_local_pos < len(reference_sequence):
            ref_base = reference_sequence[vcf_local_pos]
        else:
            # Fallback if position is at the end
            ref_base = reference_sequence[-1]
        
        # Create ALT allele (ref base + inserted sequence)
        alt_allele = ref_base + validation_result.itd_candidate.sequence
        
        # Calculate actual insertion length (length of inserted sequence)
        # This ensures SVLEN matches the actual VCF ALT allele length
        actual_insertion_length = len(validation_result.itd_candidate.sequence)
        
        # Calculate quality score (Phred-scaled confidence)
        qual_score = min(99, int(-10 * (1 - validation_result.validation_confidence) 
                                 if validation_result.validation_confidence > 0 else 0))
        
        # Build INFO field
        info_fields = [
            f"SVTYPE=DUP",
            f"SVLEN={actual_insertion_length}",
            f"END={genome_pos + actual_insertion_length}",
            f"DUPSEQ={validation_result.itd_candidate.sequence}",
            f"CONFIDENCE={validation_result.validation_confidence:.3f}"
        ]
        
        # Add PRIMARY flag if candidate is marked as primary
        if hasattr(validation_result.itd_candidate, 'is_primary') and validation_result.itd_candidate.is_primary:
            info_fields.append("PRIMARY")
        
        if validation_result.itd_candidate.duplication_start is not None:
            # Convert duplication coordinates to genome coordinates too
            dup_start_genome = self.genome_start + validation_result.itd_candidate.duplication_start
            dup_end_genome = self.genome_start + validation_result.itd_candidate.duplication_end
            info_fields.append(f"DUPSTART={dup_start_genome}")
            info_fields.append(f"DUPEND={dup_end_genome}")
        
        info = ";".join(info_fields)
        
        # Build FORMAT fields
        format_field = "GT:AD:DP:AF:SR"
        
        # Calculate allelic depths
        ref_depth = validation_result.total_coverage - validation_result.itd_coverage
        alt_depth = validation_result.itd_coverage
        
        # Determine genotype
        if validation_result.allele_frequency > 0.9:
            genotype = "1/1"  # Homozygous
        elif validation_result.allele_frequency > 0.1:
            genotype = "0/1"  # Heterozygous
        else:
            genotype = "0/1"  # Low frequency het
        
        sample_data = f"{genotype}:{ref_depth},{alt_depth}:{validation_result.total_coverage}:" \
                     f"{validation_result.allele_frequency:.3f}:{validation_result.itd_coverage}"
        
        # Build VCF line with genome coordinates
        vcf_line = f"{self.genome_chr}\t{genome_pos}\tITD_{idx}\t{ref_base}\t{alt_allele}\t" \
                  f"{qual_score}\tPASS\t{info}\t{format_field}\t{sample_data}"
        
        return vcf_line
    
    def write_vcf(self, validation_results: List['ValidationResult'], 
                  reference_sequence: str, output_file: str, 
                  command_line: str = ""):
        """Write complete VCF file"""
        with open(output_file, 'w') as f:
            # Write header
            self.write_header(f, len(reference_sequence), command_line)
            
            # Write variants
            for idx, result in enumerate(validation_results, 1):
                vcf_line = self.format_itd_variant(result, reference_sequence, idx)
                f.write(vcf_line + "\n")
        
        logger.info(f"Wrote {len(validation_results)} variants to {output_file}")

def write_vcf_results(validation_results: List['ValidationResult'],
                     reference_sequence: str,
                     output_file: str,
                     sample_name: str = "Sample",
                     genome_chr: str = "chr13", 
                     genome_start: int = 28033881) -> None:
    """Main function to write VCF results with genome coordinates"""
    writer = VCFWriter(sample_name=sample_name, genome_chr=genome_chr, genome_start=genome_start)
    writer.write_vcf(validation_results, reference_sequence, output_file)