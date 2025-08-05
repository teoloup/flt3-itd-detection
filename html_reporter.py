#!/usr/bin/env python3
"""
HTML Report Generator Module
Creates visual HTML reports for ITD detection results with IGV integration
"""

import logging
import base64
from typing import List, Dict, TYPE_CHECKING
from pathlib import Path
from datetime import datetime

if TYPE_CHECKING:
    from reference_validator import ValidationResult

logger = logging.getLogger(__name__)

class HTMLReporter:
    """Generate HTML reports for ITD results"""
    
    def __init__(self):
        self.css_style = """
        <style>
            body {
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f5f5f5;
            }
            .container {
                max-width: 1200px;
                margin: 0 auto;
                background-color: white;
                padding: 20px;
                box-shadow: 0 0 10px rgba(0,0,0,0.1);
            }
            h1, h2, h3 {
                color: #333;
            }
            .summary-box {
                background-color: #e8f4f8;
                padding: 15px;
                border-radius: 5px;
                margin: 20px 0;
            }
            .itd-card {
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 5px;
                padding: 15px;
                margin: 15px 0;
            }
            .sequence-display {
                font-family: 'Courier New', monospace;
                background-color: #f0f0f0;
                padding: 10px;
                border-radius: 3px;
                word-break: break-word;
                margin: 10px 0;
            }
            .itd-highlight {
                background-color: #ffeb3b;
                font-weight: bold;
            }
            .metric {
                display: inline-block;
                margin: 10px 20px 10px 0;
            }
            .metric-label {
                font-weight: bold;
                color: #666;
            }
            .metric-value {
                font-size: 1.2em;
                color: #333;
            }
            table {
                border-collapse: collapse;
                width: 100%;
                margin: 15px 0;
            }
            th, td {
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }
            th {
                background-color: #4CAF50;
                color: white;
            }
            tr:nth-child(even) {
                background-color: #f2f2f2;
            }
            .no-itd {
                color: #4CAF50;
                font-size: 1.2em;
                text-align: center;
                padding: 20px;
            }
            .warning {
                background-color: #fff3cd;
                border: 1px solid #ffeaa7;
                color: #856404;
                padding: 10px;
                border-radius: 5px;
                margin: 10px 0;
            }
            .confidence-high { color: #28a745; }
            .confidence-medium { color: #ffc107; }
            .confidence-low { color: #dc3545; }
            
            /* Size Analysis Styles */
            .analysis-section {
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 8px;
                padding: 20px;
                margin: 20px 0;
            }
            .size-stats {
                background-color: #e3f2fd;
                padding: 15px;
                border-radius: 5px;
                margin: 15px 0;
            }
            .stats-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 10px;
                margin-top: 10px;
            }
            .stat-item {
                background-color: white;
                padding: 8px 12px;
                border-radius: 4px;
                border-left: 4px solid #2196F3;
            }
            .stat-label {
                font-weight: bold;
                color: #666;
                display: block;
            }
            .stat-value {
                font-size: 1.1em;
                color: #333;
                font-weight: bold;
            }
            .itd-signatures {
                background-color: #fff3e0;
                padding: 15px;
                border-radius: 5px;
                margin: 15px 0;
            }
            .signature-stats {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 10px;
                margin: 10px 0;
            }
            .sig-item {
                background-color: white;
                padding: 8px 12px;
                border-radius: 4px;
                border-left: 4px solid #FF9800;
            }
            .sig-label {
                font-weight: bold;
                color: #666;
                display: block;
            }
            .sig-value {
                font-size: 1.1em;
                color: #333;
                font-weight: bold;
            }
            .size-classes {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                gap: 8px;
                margin: 10px 0;
            }
            .class-item {
                background-color: white;
                padding: 6px 10px;
                border-radius: 3px;
                border-left: 3px solid #4CAF50;
                font-size: 0.9em;
            }
            .class-label {
                font-weight: bold;
                color: #666;
            }
            .class-value {
                color: #333;
                font-weight: bold;
            }
            .size-plots {
                margin: 20px 0;
            }
            .plot-container {
                margin: 20px 0;
                text-align: center;
            }
            .analysis-plot {
                max-width: 100%;
                height: auto;
                border: 1px solid #ddd;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            .detailed-report {
                background-color: #f5f5f5;
                padding: 15px;
                border-radius: 5px;
                margin: 15px 0;
            }
            .report-content {
                background-color: white;
                padding: 15px;
                border-radius: 3px;
                border: 1px solid #ddd;
                max-height: 400px;
                overflow-y: auto;
            }
            .report-content pre {
                margin: 0;
                font-family: 'Courier New', monospace;
                font-size: 0.9em;
                line-height: 1.4;
                white-space: pre-wrap;
            }
            
            /* Size Analysis Styles */
            .analysis-section {
        </style>
        """
    
    def generate_report(self, validation_results: List['ValidationResult'],
                       sample_name: str, reference_sequence: str,
                       total_reads: int, output_file: str,
                       size_analysis_results: Dict = None,
                       size_analysis_dir: str = None):
        """Generate complete HTML report"""
        html_content = self._create_html_header(sample_name)
        html_content += self._create_summary_section(validation_results, total_reads)
        
        # Add size analysis section if available
        if size_analysis_results and size_analysis_dir:
            html_content += self._create_size_analysis_section(size_analysis_results, size_analysis_dir, sample_name)
        
        if validation_results:
            html_content += self._create_itd_details(validation_results, reference_sequence)
            html_content += self._create_visualization(validation_results, reference_sequence)
        else:
            html_content += '<div class="no-itd">No ITDs detected</div>'
        
        html_content += self._create_footer()
        
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"HTML report written to {output_file}")
    
    def _create_html_header(self, sample_name: str) -> str:
        """Create HTML header"""
        return f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>FLT3 ITD Report - {sample_name}</title>
            {self.css_style}
        </head>
        <body>
            <div class="container">
                <h1>FLT3 ITD Detection Report</h1>
                <p><strong>Sample:</strong> {sample_name}</p>
                <p><strong>Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        """
    
    def _create_summary_section(self, results: List['ValidationResult'], 
                               total_reads: int) -> str:
        """Create summary statistics section"""
        num_itds = len(results)
        
        if num_itds > 0:
            max_af = max(r.allele_frequency for r in results)
            # Use the highest ITD coverage instead of summing to avoid double counting
            # since merged candidates may include the same supporting reads
            max_itd_reads = max(r.itd_coverage for r in results)
            avg_confidence = sum(r.validation_confidence for r in results) / num_itds
        else:
            max_af = 0
            max_itd_reads = 0
            avg_confidence = 0
        
        summary = f"""
        <div class="summary-box">
            <h2>Summary</h2>
            <div class="metric">
                <span class="metric-label">ITDs Detected:</span>
                <span class="metric-value">{num_itds}</span>
            </div>
            <div class="metric">
                <span class="metric-label">Total Reads:</span>
                <span class="metric-value">{total_reads}</span>
            </div>
            <div class="metric">
                <span class="metric-label">Max ITD Supporting Reads:</span>
                <span class="metric-value">{max_itd_reads}</span>
            </div>
        """
        
        if num_itds > 0:
            summary += f"""
            <div class="metric">
                <span class="metric-label">Highest Allele Frequency:</span>
                <span class="metric-value">{max_af:.1%}</span>
            </div>
            <div class="metric">
                <span class="metric-label">Average Confidence:</span>
                <span class="metric-value">{avg_confidence:.1%}</span>
            </div>
            """
        
        summary += "</div>"
        return summary
    
    def _create_size_analysis_section(self, size_analysis_results: Dict, 
                                    size_analysis_dir: str, sample_name: str) -> str:
        """Create size analysis section with plots and statistics"""
        size_dir = Path(size_analysis_dir)
        
        # Look for plot files
        plot_files = {
            'distribution': None,
            'itd_analysis': None
        }
        
        # Find the plot files
        for plot_file in size_dir.glob("*.png"):
            if "size_distribution" in plot_file.name:
                plot_files['distribution'] = plot_file
            elif "itd_analysis" in plot_file.name:
                plot_files['itd_analysis'] = plot_file
        
        # Look for the text report
        report_file = None
        for report in size_dir.glob("*_size_analysis_report.txt"):
            report_file = report
            break
        
        section = """
        <div class="analysis-section">
            <h2>Read Size Distribution Analysis</h2>
            <p>This analysis examines the distribution of read lengths to detect potential ITD signatures.</p>
        """
        
        # Add statistics if available
        if 'basic_stats' in size_analysis_results:
            stats = size_analysis_results['basic_stats']
            section += f"""
            <div class="size-stats">
                <h3>Size Statistics</h3>
                <div class="stats-grid">
                    <div class="stat-item">
                        <span class="stat-label">Total Reads:</span>
                        <span class="stat-value">{stats.get('total_reads', 0):,}</span>
                    </div>
                    <div class="stat-item">
                        <span class="stat-label">Mean Length:</span>
                        <span class="stat-value">{stats.get('mean_length', 0):.1f} bp</span>
                    </div>
                    <div class="stat-item">
                        <span class="stat-label">Median Length:</span>
                        <span class="stat-value">{stats.get('median_length', 0):.1f} bp</span>
                    </div>
                    <div class="stat-item">
                        <span class="stat-label">Length Range:</span>
                        <span class="stat-value">{stats.get('min_length', 0):.0f} - {stats.get('max_length', 0):.0f} bp</span>
                    </div>
                </div>
            </div>
            """
        
        # Add ITD signatures if available
        if 'itd_signatures' in size_analysis_results:
            itd_sigs = size_analysis_results['itd_signatures']
            section += f"""
            <div class="itd-signatures">
                <h3>ITD Signature Analysis</h3>
                <div class="signature-stats">
                    <div class="sig-item">
                        <span class="sig-label">Estimated ITD Frequency:</span>
                        <span class="sig-value">{itd_sigs.get('estimated_itd_frequency', 0) * 100:.1f}%</span>
                    </div>
                    <div class="sig-item">
                        <span class="sig-label">Potential ITD Reads:</span>
                        <span class="sig-value">{itd_sigs.get('potential_itd_reads', 0):,}</span>
                    </div>
                    <div class="sig-item">
                        <span class="sig-label">Wild-type Reads:</span>
                        <span class="sig-value">{itd_sigs.get('wild_type_reads', 0):,}</span>
                    </div>
                </div>
            """
            
            # Add size classes breakdown
            if 'size_classes' in itd_sigs:
                size_classes = itd_sigs['size_classes']
                section += """
                <h4>Size Class Distribution</h4>
                <div class="size-classes">
                """
                for class_name, count in size_classes.items():
                    if count > 0:
                        section += f"""
                        <div class="class-item">
                            <span class="class-label">{class_name.replace('_', ' ').title()}:</span>
                            <span class="class-value">{count:,}</span>
                        </div>
                        """
                section += "</div>"
            section += "</div>"
        
        # Add plots
        if plot_files['distribution'] or plot_files['itd_analysis']:
            section += '<div class="size-plots">'
            
            if plot_files['distribution']:
                section += f"""
                <div class="plot-container">
                    <h3>Size Distribution Plot</h3>
                    <img src="data:image/png;base64,{self._encode_image(plot_files['distribution'])}" 
                         alt="Size Distribution" class="analysis-plot">
                </div>
                """
            
            if plot_files['itd_analysis']:
                section += f"""
                <div class="plot-container">
                    <h3>ITD Analysis Plot</h3>
                    <img src="data:image/png;base64,{self._encode_image(plot_files['itd_analysis'])}" 
                         alt="ITD Analysis" class="analysis-plot">
                </div>
                """
            
            section += '</div>'
        
        # Add detailed report if available
        if report_file and report_file.exists():
            try:
                with open(report_file, 'r') as f:
                    report_content = f.read()
                
                section += f"""
                <div class="detailed-report">
                    <h3>Detailed Size Analysis Report</h3>
                    <div class="report-content">
                        <pre>{report_content}</pre>
                    </div>
                </div>
                """
            except Exception as e:
                logger.warning(f"Could not read size analysis report: {e}")
        
        section += "</div>"
        return section
    
    def _encode_image(self, image_path: Path) -> str:
        """Encode image as base64 for embedding in HTML"""
        try:
            with open(image_path, 'rb') as img_file:
                return base64.b64encode(img_file.read()).decode('utf-8')
        except Exception as e:
            logger.error(f"Error encoding image {image_path}: {e}")
            return ""
    
    def _create_itd_details(self, results: List['ValidationResult'], 
                           reference_sequence: str) -> str:
        """Create detailed ITD information"""
        details = "<h2>ITD Details</h2>"
        
        for i, result in enumerate(results, 1):
            confidence_class = self._get_confidence_class(result.validation_confidence)
            
            details += f"""
            <div class="itd-card">
                <h3>ITD {i} {'🔺 PRIMARY' if hasattr(result.itd_candidate, 'is_primary') and result.itd_candidate.is_primary else ''}</h3>
                <table>
                    <tr>
                        <th>Property</th>
                        <th>Value</th>
                    </tr>
                    <tr>
                        <td>Length</td>
                        <td>{len(result.itd_candidate.sequence)} bp (sequence) / {result.duplication_length} bp (calculated)</td>
                    </tr>
                    <tr>
                        <td>Position</td>
                        <td>{result.insertion_position}</td>
                    </tr>
                    <tr>
                        <td>Allele Frequency</td>
                        <td>{result.allele_frequency:.1%}</td>
                    </tr>
                    <tr>
                        <td>Supporting Reads</td>
                        <td>{result.itd_coverage} / {result.total_coverage}</td>
                    </tr>
                    <tr>
                        <td>Confidence</td>
                        <td class="{confidence_class}">{result.validation_confidence:.1%}</td>
                    </tr>
                    {'<tr><td>Primary ITD</td><td><span style="color: #d63384; font-weight: bold;">Yes (matches size analysis)</span></td></tr>' if hasattr(result.itd_candidate, 'is_primary') and result.itd_candidate.is_primary else ''}
                </table>
                
                <h4>ITD Sequence (from reads):</h4>
                <div class="sequence-display">
                    {self._format_sequence(result.itd_candidate.sequence)}
                </div>
                
                <h4>Context:</h4>
                {self._create_context_display(result, reference_sequence)}
            </div>
            """
        
        return details

    def _create_context_display(self, result: 'ValidationResult', reference_sequence: str) -> str:
        """Show ITD in sequence context"""
        pos = result.insertion_position
        context_size = 50

        # Get context around insertion point
        start = max(0, pos - context_size)
        end = min(len(reference_sequence), pos + context_size)

        left_context = reference_sequence[start:pos]
        right_context = reference_sequence[pos:end]
        itd_sequence = result.itd_candidate.sequence  # Use actual duplication

        context_html = f"""
        <div class="sequence-display">
            <span>{self._format_sequence(left_context)}</span>
            <span class="itd-highlight">{self._format_sequence(itd_sequence)}</span>
            <span>{self._format_sequence(right_context)}</span>
        </div>
        """
        return context_html

    def _create_visualization(self, results: List['ValidationResult'], 
                             reference_sequence: str) -> str:
        """Create visual representation of ITDs"""
        viz = """
        <h2>ITD Visualization</h2>
        <div class="itd-card">
            <p>Reference sequence with ITD positions marked:</p>
            <div style="position: relative; margin: 20px 0;">
        """
        
        ref_length = len(reference_sequence)
        scale_factor = 800 / ref_length
        
        viz += """
            <div style="background-color: #ccc; height: 30px; position: relative;">
                <span style="position: absolute; left: 5px; top: 5px;">FLT3 Reference</span>
            </div>
        """
        
        for i, result in enumerate(results, 1):
            pos = result.insertion_position
            left_px = int(pos * scale_factor)
            width_px = max(10, int(result.duplication_length * scale_factor))
            
            viz += f"""
            <div style="position: absolute; left: {left_px}px; top: 40px; 
                        width: {width_px}px; height: 20px; 
                        background-color: #ff6b6b; opacity: 0.7;
                        border: 1px solid #c92a2a;" 
                 title="ITD {i}: {result.duplication_length}bp at position {pos}">
            </div>
            """
        
        viz += """
            </div>
            <p style="margin-top: 80px;"><em>Red bars indicate ITD positions and relative sizes</em></p>
        </div>
        """
        
        return viz
    
    def _format_sequence(self, sequence: str, chunk_size: int = 10) -> str:
        chunks = [sequence[i:i+chunk_size] for i in range(0, len(sequence), chunk_size)]
        return ' '.join(chunks)
    
    def _get_confidence_class(self, confidence: float) -> str:
        if confidence >= 0.8:
            return "confidence-high"
        elif confidence >= 0.5:
            return "confidence-medium"
        else:
            return "confidence-low"
    
    def _create_footer(self) -> str:
        return """
            </div>
        </body>
        </html>
        """

def generate_html_report(validation_results: List['ValidationResult'],
                        sample_name: str,
                        reference_sequence: str,
                        total_reads: int,
                        output_file: str,
                        size_analysis_results: Dict = None,
                        size_analysis_dir: str = None) -> None:
    reporter = HTMLReporter()
    reporter.generate_report(validation_results, sample_name, 
                           reference_sequence, total_reads, output_file,
                           size_analysis_results, size_analysis_dir)
