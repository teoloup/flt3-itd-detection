#!/usr/bin/env python3
"""
Configuration Module for FLT3 ITD Detection
Clean CIGAR-based architecture - no legacy code
"""

import argparse
from dataclasses import dataclass
from typing import Optional

@dataclass
class FLT3Config:
    """Configuration for FLT3 ITD detection using CIGAR-based approach"""
    
    # Required fields
    bam_file: str
    output_dir: str
    
    # Core parameters
    sample_name: str = "Sample"
    genome_build: str = "hg38"
    reference_sequence: Optional[str] = None
    
    # ITD detection parameters
    min_itd_length: int = 15
    max_itd_length: int = 500
    min_supporting_reads: int = 3
    min_allele_frequency: float = 0.01
    max_candidates: int = 50
    
    # Quality filters
    min_mapping_quality: int = 20
    min_read_length: int = 250
    wt_tolerance: int = 20  # Wild-type size tolerance (±bp)
    
    # CIGAR detection parameters
    min_softclip_length: int = 50
    enable_softclip_fallback: bool = True
    
    # Output options
    write_html: bool = True
    write_vcf: bool = True
    bamout: bool = False  # Export BAM/FASTA for IGV viewing
    debug: bool = False
    
    # Processing
    threads: int = 4
    
    # Genomic coordinates for VCF output
    chromosome: str = "chr13"
    start_hg38: int = 28033301
    end_hg38: int = 28034800
    start_hg19: int = 28607438
    end_hg19: int = 28608937
    
    # BLAT coordinates for VCF annotation (actual reference coordinates)
    # These are the coordinates where our amplicon reference maps to the genome
    vcf_chr: str = "chr13"
    vcf_start_hg38: int = 28033881  # Start coordinate for VCF output (matches genome browsers)
    vcf_end_hg38: int = 28034216    # End coordinate for VCF output  
    vcf_start_hg19: int = 28608018  # Start coordinate for VCF output (hg19 equivalent)
    vcf_end_hg19: int = 28608353    # End coordinate for VCF output
    
    # Amplicon settings for primer trimming
    forward_primer: str = "AGCAATTTAGGTATGAAAGCCAGC"
    reverse_primer: str = "CTGTACCTTTCAGCATTTTGACG"
    amplicon_length: int = 336
    
    def __post_init__(self):
        """Set reference sequence based on genome build"""
        if self.reference_sequence is None:
            # FLT3 exon 14-15 reference sequence (reverse complement to match genome browsers)
            # hg38: chr13:28033881-28034216 (336 bp)
            self.reference_sequence = (
                "CTGTACCTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTG"
                "CTAATTCCATAAGCTGTTGCGTTCATCACTTTTCCAAAAGCACCTGATCC"
                "TAGTACCTTCCCTGCAAAGACAAATGGTGAGTACGTGCATTTTAAAGATT"
                "TTCCAATGGAAAAGAAATGCTGCAGAAACATTTGGCACATTCCATTCTTA"
                "CCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATA"
                "TTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCT"
                "GTACCATCTGTAGCTGGCTTTCATACCTAAATTGCT"
            )
    
    @property
    def vcf_start_pos(self) -> int:
        """Get VCF start position based on genome build"""
        return self.vcf_start_hg38 if self.genome_build == "hg38" else self.vcf_start_hg19
    
    @property  
    def vcf_end_pos(self) -> int:
        """Get VCF end position based on genome build"""
        return self.vcf_end_hg38 if self.genome_build == "hg38" else self.vcf_end_hg19

def create_argument_parser():
    """Create command line argument parser for CIGAR-based ITD detection"""
    parser = argparse.ArgumentParser(
        description="FLT3 ITD Detection using CIGAR-based approach",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -b sample.bam -o results/ -n SampleName
  %(prog)s -b sample.bam -o results/ --min-supporting-reads 5 --min-allele-frequency 0.05
  %(prog)s -b sample.bam -o results/ --bamout --debug
  %(prog)s -b sample.bam -o results/ --max-itd-length 200 --max-candidates 20
        """
    )
    
    # Required arguments
    parser.add_argument('-b', '--bam', required=True,
                        help='Input BAM file with reads aligned to FLT3 region')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory')
    
    # Core arguments
    parser.add_argument('-n', '--name', default='Sample',
                        help='Sample name (default: Sample)')
    parser.add_argument('-g', '--genome', choices=['hg19', 'hg38'], default='hg38',
                        help='Reference genome build (default: hg38)')
    
    # ITD detection parameters
    parser.add_argument('--min-itd-length', type=int, default=15,
                        help='Minimum ITD length (default: 15)')
    parser.add_argument('--max-itd-length', type=int, default=500,
                        help='Maximum ITD length (default: 500)')
    parser.add_argument('--min-supporting-reads', type=int, default=3,
                        help='Minimum supporting reads (default: 3)')
    parser.add_argument('--min-allele-frequency', type=float, default=0.01,
                        help='Minimum allele frequency (default: 0.01)')
    parser.add_argument('--max-candidates', type=int, default=50,
                        help='Maximum candidates to validate (default: 50)')
    parser.add_argument('--wt-tolerance', type=int, default=20,
                        help='Maximum wild-type tolerance (default: 20)')
    
    # Quality filters
    parser.add_argument('--min-mapping-quality', type=int, default=20,
                        help='Minimum mapping quality (default: 20)')
    parser.add_argument('--min-read-length', type=int, default=250,
                        help='Minimum read length after trimming (default: 250)')
    
    # CIGAR detection options
    parser.add_argument('--min-softclip-length', type=int, default=50,
                        help='Minimum soft-clip length for secondary detection (default: 50)')
    parser.add_argument('--no-softclip-fallback', action='store_true',
                        help='Disable soft-clip fallback analysis')
    
    # Output options
    parser.add_argument('--no-html', action='store_true',
                        help='Do not generate HTML report')
    parser.add_argument('--no-vcf', action='store_true',
                        help='Do not generate VCF file')
    parser.add_argument('--bamout', action='store_true',
                        help='Export reference FASTA and BAM files for IGV viewing')
    
    # Processing options
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of threads (default: 4)')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug logging')
    
    return parser

def parse_arguments():
    """Parse command line arguments and create config"""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Create config object with clean parameters
    config = FLT3Config(
        bam_file=args.bam,
        output_dir=args.output,
        sample_name=args.name,
        genome_build=args.genome,
        min_itd_length=args.min_itd_length,
        max_itd_length=args.max_itd_length,
        min_supporting_reads=args.min_supporting_reads,
        min_allele_frequency=args.min_allele_frequency,
        max_candidates=args.max_candidates,
        min_mapping_quality=args.min_mapping_quality,
        min_read_length=args.min_read_length,
        min_softclip_length=args.min_softclip_length,
        wt_tolerance=args.wt_tolerance,
        enable_softclip_fallback=not args.no_softclip_fallback,
        write_html=not args.no_html,
        write_vcf=not args.no_vcf,
        bamout=args.bamout,
        threads=args.threads,
        debug=args.debug
    )
    
    return config