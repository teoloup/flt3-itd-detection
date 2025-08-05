#!/usr/bin/env python3
"""
Integration Adapter Module
Bridges between old and new ITD detection architectures
"""

import logging
from typing import List, Dict, Any, Union
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class LegacyValidationResult:
    """Validation result in old format for compatibility with VCF/HTML writers"""
    itd_candidate: Any  # Will contain the candidate object
    is_valid: bool
    validation_confidence: float
    allele_frequency: float
    total_coverage: int
    itd_coverage: int
    insertion_position: int
    duplication_length: int
    duplication_start: int = None
    duplication_end: int = None

class IntegrationAdapter:
    """Adapter to bridge old and new ITD detection systems"""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def convert_new_to_legacy_results(self, new_candidates: List, validation_results: List) -> List[LegacyValidationResult]:
        """
        Convert new architecture results to legacy format for VCF/HTML writers
        
        Args:
            new_candidates: List of candidates from new pipeline
            validation_results: List of validation results from new pipeline
            
        Returns:
            List of LegacyValidationResult objects
        """
        legacy_results = []
        
        for i, candidate in enumerate(new_candidates):
            # Get corresponding validation result
            validation = validation_results[i] if i < len(validation_results) else None
            
            # Create legacy ITD candidate object
            legacy_candidate = self._create_legacy_candidate(candidate)
            
            # Create legacy validation result
            legacy_result = LegacyValidationResult(
                itd_candidate=legacy_candidate,
                is_valid=validation.is_valid if validation else True,
                validation_confidence=validation.confidence if validation else candidate.confidence,
                allele_frequency=validation.allele_frequency if validation else 0.1,
                total_coverage=validation.total_reads if validation else len(candidate.supporting_reads) * 2,
                itd_coverage=validation.itd_reads if validation else len(candidate.supporting_reads),
                insertion_position=candidate.position,
                duplication_length=candidate.length,
                duplication_start=getattr(candidate, 'duplication_start', None),
                duplication_end=getattr(candidate, 'duplication_end', None)
            )
            
            legacy_results.append(legacy_result)
            
        return legacy_results
    
    def _create_legacy_candidate(self, new_candidate) -> Any:
        """Create legacy candidate object from new candidate"""
        
        class LegacyITDCandidate:
            def __init__(self, new_cand):
                self.sequence = new_cand.sequence
                self.length = new_cand.length
                self.position = new_cand.position
                self.supporting_reads = new_cand.supporting_reads
                self.confidence = new_cand.confidence
                self.support_type = getattr(new_cand, 'support_type', 'cigar_insertion')
                self.insertion_site = getattr(new_cand, 'insertion_site', new_cand.position)
                self.duplication_start = getattr(new_cand, 'duplication_start', None)
                self.duplication_end = getattr(new_cand, 'duplication_end', None)
        
        return LegacyITDCandidate(new_candidate)
    
    def extract_reference_file_from_config(self, config) -> str:
        """Extract or create reference FASTA file path from config"""
        
        # Check if reference FASTA already specified
        if hasattr(config, 'reference_fasta') and config.reference_fasta:
            return config.reference_fasta
        
        # Create temporary reference file
        import tempfile
        from pathlib import Path
        
        temp_ref = Path(config.output_dir) / f"{config.sample_name}_reference.fasta"
        
        with open(temp_ref, 'w') as f:
            f.write(f">FLT3_reference\n")
            f.write(f"{config.reference_sequence}\n")
        
        self.logger.info(f"Created temporary reference file: {temp_ref}")
        return str(temp_ref)
    
    def create_trimmed_bam_path(self, config) -> str:
        """Create path for trimmed BAM file that new pipeline expects"""
        
        # The new pipeline expects a BAM file, but we need to point to the extracted/trimmed reads
        # For now, we'll use the original BAM file and let the new pipeline handle extraction
        # In future, we could create a trimmed BAM, but that's more complex
        
        return config.bam_file  # Use original BAM for now
    
    def prepare_new_pipeline_args(self, config) -> Dict[str, Any]:
        """Prepare arguments for the new pipeline from legacy config"""
        
        # Extract reference file path
        reference_file = self.extract_reference_file_from_config(config)
        
        # Create arguments dictionary
        args = {
            'bam_file': self.create_trimmed_bam_path(config),
            'reference_sequence': config.reference_sequence,
            'reference_file': reference_file,
            'min_itd_length': config.min_itd_length,
            'max_itd_length': config.max_itd_length,
            'min_support': config.min_supporting_reads,
            'min_softclip_length': getattr(config, 'min_softclip_length', 50),
            'enable_softclip_fallback': getattr(config, 'enable_softclip_fallback', True),
            'output_prefix': f"{config.output_dir}/{config.sample_name}"
        }
        
        return args
    
    def validate_new_pipeline_dependencies(self) -> List[str]:
        """Check if new pipeline dependencies are available"""
        missing_deps = []
        
        # Check for required Python packages
        try:
            import pysam
        except ImportError:
            missing_deps.append("pysam")
        
        try:
            import numpy
        except ImportError:
            missing_deps.append("numpy")
        
        # Check for external tools
        import shutil
        if not shutil.which('minimap2'):
            missing_deps.append("minimap2 (not in PATH)")
        
        if not shutil.which('samtools'):
            missing_deps.append("samtools (not in PATH)")
        
        return missing_deps
    
    def log_pipeline_transition(self, config):
        """Log the transition from old to new pipeline"""
        self.logger.info("=" * 60)
        self.logger.info("TRANSITIONING TO NEW CIGAR-BASED ITD DETECTION")
        self.logger.info("=" * 60)
        self.logger.info("Old approach: Split-read overlap analysis")
        self.logger.info("New approach: Direct CIGAR + Soft-clip + Dual-reference validation")
        self.logger.info(f"Expected size accuracy improvement for {config.sample_name}")
        self.logger.info("=" * 60)


def create_integration_adapter() -> IntegrationAdapter:
    """Factory function to create integration adapter"""
    return IntegrationAdapter()


# Compatibility functions for easy migration
def run_new_pipeline_with_legacy_config(config) -> tuple:
    """
    Run new pipeline with legacy config and return legacy-compatible results
    
    Args:
        config: Legacy FLT3Config object
        
    Returns:
        Tuple of (legacy_validation_results, detection_summary)
    """
    adapter = IntegrationAdapter()
    
    # Check dependencies
    missing_deps = adapter.validate_new_pipeline_dependencies()
    if missing_deps:
        raise RuntimeError(f"Missing dependencies for new pipeline: {', '.join(missing_deps)}")
    
    # Log transition
    adapter.log_pipeline_transition(config)
    
    # Prepare arguments for new pipeline
    pipeline_args = adapter.prepare_new_pipeline_args(config)
    
    # Import and run new pipeline
    try:
        from flt3_itd_pipeline import run_flt3_itd_pipeline
        
        result = run_flt3_itd_pipeline(**pipeline_args)
        
        # Convert results to legacy format
        legacy_results = adapter.convert_new_to_legacy_results(
            result.validated_candidates,
            result.validation_results
        )
        
        return legacy_results, result.detection_summary
        
    except ImportError as e:
        logger.error(f"Failed to import new pipeline modules: {e}")
        logger.error("Falling back to old pipeline...")
        raise
    
    except Exception as e:
        logger.error(f"New pipeline failed: {e}")
        logger.error("Consider falling back to old pipeline")
        raise


if __name__ == "__main__":
    # Test the adapter
    import argparse
    
    parser = argparse.ArgumentParser(description="Test integration adapter")
    parser.add_argument("--test-deps", action="store_true", help="Test dependencies")
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    adapter = IntegrationAdapter()
    
    if args.test_deps:
        missing = adapter.validate_new_pipeline_dependencies()
        if missing:
            print(f"Missing dependencies: {missing}")
        else:
            print("All dependencies available!")
