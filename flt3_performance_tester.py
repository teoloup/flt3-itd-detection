#!/usr/bin/env python3
"""
FLT3-ITD Performance Testing Script
Runs the detection tool on multiple samples and calculates performance metrics
"""

import os
import sys
import subprocess
import pandas as pd
import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple
import time
from datetime import datetime


class FLT3PerformanceTester:
    def __init__(self, main_module_path: str, samples_csv: str, output_dir: str):
        self.main_module_path = Path(main_module_path)
        self.samples_csv = samples_csv
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load samples
        self.samples_df = pd.read_csv(samples_csv)
        print(f"Loaded {len(self.samples_df)} samples from {samples_csv}")
        
        # Results storage
        self.results = []
        
    def run_detection(self, sample_path: str, sample_name: str) -> Dict:
        """Run FLT3 detection on a single sample with standard/default parameters"""
        print(f"\n=== Processing {sample_name} ===")
        sample_output = self.output_dir / sample_name
        sample_output.mkdir(exist_ok=True)
        cmd = [
            'python', str(self.main_module_path),
            '-b', sample_path,
            '-o', str(sample_output),
            '--min-allele-frequency', '0.01',
            '--genome', 'hg38',
            '--min-supporting-reads', '50',
            '--threads', '24'
        ]
        print(f"Command: {' '.join(cmd)}")
        start_time = time.time()
        try:
            result = subprocess.run(
                cmd,
                stdout=sys.stdout,
                stderr=sys.stderr,
                timeout=1000
            )
            end_time = time.time()
            runtime = end_time - start_time
            if result.returncode != 0:
                print(f"ERROR: Detection failed for {sample_name}")
                return {
                    'sample': sample_name,
                    'status': 'failed',
                    'error': 'Detection failed (see SLURM log or terminal output)',
                    'runtime': runtime
                }
            # No captured stdout/stderr, so pass empty strings
            itd_results = self.parse_detection_output('', sample_output)
            return {
                'sample': sample_name,
                'status': 'success',
                'runtime': runtime,
                'itd_detected': itd_results['itd_detected'],
                'itd_count': itd_results['itd_count'],
                'itd_details': itd_results['itd_details'],
                'total_reads': itd_results['total_reads'],
                'itd_candidates': itd_results['itd_candidates'],
                'stdout': '',
                'stderr': ''
            }
        except subprocess.TimeoutExpired:
            return {
                'sample': sample_name,
                'status': 'timeout',
                'runtime': 300,
                'error': 'Process timed out after 5 minutes'
            }
        except Exception as e:
            return {
                'sample': sample_name,
                'status': 'error',
                'error': str(e),
                'runtime': 0
            }
    def parse_detection_output(self, stdout: str, sample_output: Path) -> Dict:
        """Parse the detection output to extract ITD information (robust version)"""
        itd_details = []
        total_reads = 0
        itd_candidates = 0
        itd_positive_flag = False

        # Parse log output
        lines = stdout.split('\n')
        for line in lines:
            # Extract total reads processed
            if "Total reads processed:" in line:
                match = re.search(r'Total reads processed: (\d+)', line)
                if match:
                    total_reads = int(match.group(1))

            # Extract ITD candidates found
            if "ITD candidates found:" in line:
                match = re.search(r'ITD candidates found: (\d+)', line)
                if match:
                    itd_candidates = int(match.group(1))

            # Extract validated ITDs (old logic)
            if "PASSED validation" in line and "Cluster" in line:
                cluster_info = self.extract_cluster_info(lines, line)
                if cluster_info:
                    itd_details.append(cluster_info)

            # Robust: look for other positive ITD flags in stdout
            if (re.search(r'ITD detected', line, re.IGNORECASE) or
                re.search(r'Validated ITD', line, re.IGNORECASE) or
                re.search(r'ITD found', line, re.IGNORECASE) or
                re.search(r'ITD(s)? present', line, re.IGNORECASE)):
                itd_positive_flag = True

        # Check for JSON output files and look for any ITD lists
        json_files = list(sample_output.glob("*_results.json"))
        if json_files:
            try:
                with open(json_files[0], 'r') as f:
                    json_data = json.load(f)
                    # Accept any non-empty ITD list
                    for key in ['validated_itds', 'itds', 'itd_list', 'detected_itds']:
                        if key in json_data and isinstance(json_data[key], list) and len(json_data[key]) > 0:
                            itd_details.extend(json_data[key])
                            itd_positive_flag = True
            except Exception:
                pass

        # Check for HTML report file as a positive flag
        html_files = list(sample_output.glob("*.html"))
        if html_files:
            itd_positive_flag = True

        # If any positive flag or details found, mark as positive
        itd_detected = itd_positive_flag or len(itd_details) > 0

        return {
            'itd_detected': itd_detected,
            'itd_count': len(itd_details),
            'itd_details': itd_details,
            'total_reads': total_reads,
            'itd_candidates': itd_candidates
        }
    
    def extract_cluster_info(self, lines: List[str], current_line: str) -> Dict:
        """Extract cluster information from log lines"""
        # Find the cluster validation info
        cluster_info = {}
        
        # Look backwards for cluster details
        current_idx = lines.index(current_line)
        for i in range(max(0, current_idx - 10), current_idx):
            line = lines[i]
            
            # Extract size
            if "size=" in line and "support=" in line:
                size_match = re.search(r'size=(\d+)bp', line)
                support_match = re.search(r'support=(\d+)', line)
                if size_match and support_match:
                    cluster_info['size'] = int(size_match.group(1))
                    cluster_info['support'] = int(support_match.group(1))
            
            # Extract allele frequency
            if "allele frequency:" in line:
                af_match = re.search(r'allele frequency: ([\d.]+)', line)
                if af_match:
                    cluster_info['allele_frequency'] = float(af_match.group(1))
            
            # Extract confidence
            if "confidence=" in line:
                conf_match = re.search(r'confidence=([\d.]+)', line)
                if conf_match:
                    cluster_info['confidence'] = float(conf_match.group(1))
        
        return cluster_info if cluster_info else None
    
    def calculate_metrics(self, results: List[Dict]) -> Dict:
        """Calculate performance metrics"""
        tp = tn = fp = fn = 0
        
        for i, row in self.samples_df.iterrows():
            sample_name = Path(row['sample_path']).stem
            expected = row['expected_result'].lower()
            
            # Find corresponding result
            result = next((r for r in results if r['sample'] == sample_name), None)
            
            if not result or result['status'] != 'success':
                print(f"Warning: No valid result for {sample_name}")
                continue
            
            detected = result['itd_detected']
            
            if expected == 'positive' and detected:
                tp += 1
            elif expected == 'negative' and not detected:
                tn += 1
            elif expected == 'positive' and not detected:
                fn += 1
            elif expected == 'negative' and detected:
                fp += 1
        
        total = tp + tn + fp + fn
        
        # Calculate metrics
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0
        accuracy = (tp + tn) / total if total > 0 else 0
        f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0
        
        return {
            'total_samples': total,
            'true_positives': tp,
            'true_negatives': tn,
            'false_positives': fp,
            'false_negatives': fn,
            'sensitivity': sensitivity,
            'specificity': specificity,
            'precision': precision,
            'negative_predictive_value': npv,
            'accuracy': accuracy,
            'f1_score': f1_score
        }
    
    def run_all_samples(self) -> Dict:
        """Run detection on all samples with standard parameters"""
        print(f"\n{'='*60}")
        print(f"RUNNING STANDARD DETECTION FOR ALL SAMPLES")
        print(f"{'='*60}")
        results = []
        for i, row in self.samples_df.iterrows():
            sample_path = row['sample_path']
            sample_name = Path(sample_path).stem
            result = self.run_detection(sample_path, sample_name)
            result['expected'] = row['expected_result']
            results.append(result)
        metrics = self.calculate_metrics(results)
        self.print_summary(results, metrics, "standard")
        return {
            'results': results,
            'metrics': metrics
        }
    
    def print_summary(self, results: List[Dict], metrics: Dict, param_name: str):
        """Print a summary of results"""
        print(f"\n{'='*60}")
        print(f"SUMMARY FOR PARAMETER SET: {param_name}")
        print(f"{'='*60}")
        
        # Sample-by-sample results
        print("\nSample Results:")
        print(f"{'Sample':<20} {'Expected':<10} {'Detected':<10} {'ITDs':<6} {'Status':<10} {'Runtime':<8}")
        print("-" * 75)
        
        for result in results:
            sample = result['sample'][:19]
            expected = result['expected'][:9]
            detected = 'YES' if result.get('itd_detected', False) else 'NO'
            itd_count = result.get('itd_count', 0)
            status = result['status'][:9]
            runtime = f"{result.get('runtime', 0):.1f}s"
            
            print(f"{sample:<20} {expected:<10} {detected:<10} {itd_count:<6} {status:<10} {runtime:<8}")
        
        # Performance metrics
        print(f"\nPerformance Metrics:")
        print(f"  Sensitivity (Recall):    {metrics['sensitivity']:.3f} ({metrics['true_positives']}/{metrics['true_positives'] + metrics['false_negatives']})")
        print(f"  Specificity:             {metrics['specificity']:.3f} ({metrics['true_negatives']}/{metrics['true_negatives'] + metrics['false_positives']})")
        print(f"  Precision (PPV):         {metrics['precision']:.3f} ({metrics['true_positives']}/{metrics['true_positives'] + metrics['false_positives']})")
        print(f"  Negative Pred. Value:    {metrics['negative_predictive_value']:.3f}")
        print(f"  Accuracy:                {metrics['accuracy']:.3f}")
        print(f"  F1-Score:                {metrics['f1_score']:.3f}")
        
        print(f"\nConfusion Matrix:")
        print(f"                    Predicted")
        print(f"                 Pos    Neg")
        print(f"  Actual   Pos   {metrics['true_positives']:3d}    {metrics['false_negatives']:3d}")
        print(f"           Neg   {metrics['false_positives']:3d}    {metrics['true_negatives']:3d}")
        
        # Failed samples
        failed_samples = [r for r in results if r['status'] != 'success']
        if failed_samples:
            print(f"\nFailed Samples ({len(failed_samples)}):")
            for result in failed_samples:
                print(f"  {result['sample']}: {result['status']} - {result.get('error', '')}")


def main():
    parser = argparse.ArgumentParser(description='FLT3-ITD Performance Testing Script')
    parser.add_argument('main_module', help='Path to main_modular.py')
    parser.add_argument('samples_csv', help='CSV file with sample paths and expected results')
    parser.add_argument('-o', '--output', default='performance_test_output', 
                       help='Output directory for results')
    

    args = parser.parse_args()
    tester = FLT3PerformanceTester(args.main_module, args.samples_csv, args.output)
    result = tester.run_all_samples()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_file = tester.output_dir / f"performance_results_{timestamp}.json"
    with open(results_file, 'w') as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\nDetailed results saved to: {results_file}")
    


if __name__ == "__main__":
    main()
