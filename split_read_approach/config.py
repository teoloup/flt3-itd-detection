#!/usr/bin/env python3
"""
Configuration Module
Manages all configuration parameters for FLT3 ITD detection
"""

import argparse
from dataclasses import dataclass
from typing import Optional

@dataclass
class FLT3Config:
    # Required fields (no defaults)
    bam_file: str
    output_dir: str

    # IGV/validation output
    bamout: bool = False

    # Optional fields (with defaults)
    chromosome: str = "chr13"
    # Default coordinates for hg38
    start_hg38: int = 28033301
    end_hg38: int = 28034800
    # Default coordinates for hg19
    start_hg19: int = 28607438
    end_hg19: int = 28608937
    """Configuration for FLT3 ITD detection"""
    # Input/Output
    sample_name: str = "Sample"

    # Reference settings
    genome_build: str = "hg38"
    reference_sequence: Optional[str] = None

    # Amplicon settings
    forward_primer: str = "AGCAATTTAGGTATGAAAGCCAGC"
    reverse_primer: str = "CTGTACCTTTCAGCATTTTGACG"
    amplicon_length: int = 336

    # Quality filters
    min_mapping_quality: int = 20
    min_read_length: int = 250
    min_base_quality: int = 10
    # ITD detection parameters
    min_itd_length: int = 15
    max_itd_length: int = 200
    min_supporting_reads: int = 3
    min_allele_frequency: float = 0.01
    max_candidates_per_position: int = 5  # Limit candidates at similar positions

    # Split read parameters
    split_method: str = "multiple"  # "center", "multiple", or "sliding"
    split_overlap: int = 20
    min_split_size: int = 100
    num_splits: int = 3

    # Alignment parameters
    min_overlap_length: int = 15
    min_softclip_length: int = 15
    
    # Validation parameters
    validation_min_reads: int = 3
    confidence_threshold: float = 0.5
    min_high_confidence_ratio: float = 0.5  # Ratio of high-confidence ITD reads
    
    # Output options
    write_html: bool = True
    write_vcf: bool = True
    keep_temp: bool = False  # Use keep_temp for consistency
    debug: bool = False
    fastq_file: Optional[str] = None  # Path to trimmed FASTQ for validator
    
    # Processing options
    threads: int = 24
    
    def __post_init__(self):
        """Set reference sequence based on genome build"""
        if self.reference_sequence is None:
            # This is the FLT3 exon 14-15 reference sequence
            self.reference_sequence = (
                "AGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCC"
                "TCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCT"
                "CAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGT"
                "GCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGC"
                "ACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTTT"
                "TGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCT"
                "CAATCCAGGTTGCCGTCAAAATGCTGAAAGGTACAG"
            )

def create_argument_parser():

    """Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description="FLT3 ITD Detection using split-read approach for Nanopore data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -b sample.bam -o results/ -n SampleName
  %(prog)s -b sample.bam -o results/ --min-supporting-reads 5 --min-allele-frequency 0.05
  %(prog)s -b sample.bam -o results/ --split-method sliding --debug
        """
    )
    
    # Required arguments
    parser.add_argument('-b', '--bam', required=True,
                        help='Input BAM file with reads aligned to FLT3 region')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory')
    
    # Optional arguments
    parser.add_argument('-n', '--name', default='Sample',
                        help='Sample name (default: Sample)')
    parser.add_argument('-g', '--genome', choices=['hg19', 'hg38'], default='hg38',
                        help='Reference genome build (default: hg38)')
    parser.add_argument('-R','--reference-fasta', type=str, default=None,
                        help='Reference FASTA file for samtools view -T (required if BAM header lacks @SQ lines)')
    # Quality filters
    parser.add_argument('--min-mapping-quality', type=int, default=20,
                        help='Minimum mapping quality (default: 20)')
    parser.add_argument('--min-read-length', type=int, default=250,
                        help='Minimum read length after trimming (default: 250)')
    
    # ITD detection parameters
    parser.add_argument('--min-itd-length', type=int, default=15,
                        help='Minimum ITD length (default: 15)')
    parser.add_argument('--max-itd-length', type=int, default=500,
                        help='Maximum ITD length (default: 500)')
    parser.add_argument('--min-supporting-reads', type=int, default=3,
                        help='Minimum supporting reads (default: 3)')
    parser.add_argument('--min-allele-frequency', type=float, default=0.01,
                        help='Minimum allele frequency (default: 0.01)')
    
    # Split read parameters
    parser.add_argument('--split-method', choices=['center', 'multiple', 'sliding'],
                        default='multiple',
                        help='Read splitting method (default: multiple)')
    parser.add_argument('--split-overlap', type=int, default=20,
                        help='Overlap between split pairs (default: 20)')
    parser.add_argument('--num-splits', type=int, default=3,
                        help='Number of split points for multiple method (default: 3)')
    
    # Output options
    parser.add_argument('--no-html', action='store_true',
                        help='Do not generate HTML report')
    parser.add_argument('--no-vcf', action='store_true',
                        help='Do not generate VCF file')
    parser.add_argument('--bamout', action='store_true',
                        help='For each validated ITD, output a reference FASTA and BAM of supporting reads for IGV viewing (in results directory)')

    
    # Processing options
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of threads (default: 4)')
    parser.add_argument('--keep-temp', action='store_true',
                        help='Keep temporary files')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug logging')
    
    parser.add_argument('--conservative', action='store_true',
                        help='Use conservative settings (fewer splits, stricter filtering)')
    parser.add_argument('--max-candidates', type=int, default=50,
                        help='Maximum candidates to validate (default: 50)')
    
    return parser

def parse_arguments():
    """Parse command line arguments"""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Apply conservative settings if requested
    if args.conservative:
        num_splits = 2  # Fewer splits
        min_supporting_reads = 5  # Higher support requirement
        min_confidence = 0.7  # Higher confidence threshold
    else:
        num_splits = args.num_splits
        min_supporting_reads = args.min_supporting_reads
        min_confidence = 0.5
    
    # Create config object
    config = FLT3Config(
        bam_file=args.bam,
        output_dir=args.output,
        sample_name=args.name,
        genome_build=args.genome,
        min_mapping_quality=args.min_mapping_quality,
        min_read_length=args.min_read_length,
        min_itd_length=args.min_itd_length,
        max_itd_length=args.max_itd_length,
        min_supporting_reads=min_supporting_reads,
        min_allele_frequency=args.min_allele_frequency,
        split_method=args.split_method,
        split_overlap=args.split_overlap,
        num_splits=num_splits,
        write_html=not args.no_html,
        write_vcf=not args.no_vcf,
        threads=args.threads,
        keep_temp=args.keep_temp,
        debug=args.debug,
        max_candidates_per_position=5 if args.conservative else 10,
        confidence_threshold=min_confidence,
        min_high_confidence_ratio=0.7 if args.conservative else 0.5,
        bamout=args.bamout
    )
    
    return config