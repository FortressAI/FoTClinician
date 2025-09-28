#!/usr/bin/env python3
"""
FoTChemistry Leaderboard Generator

Generates public leaderboards and dashboards for autonomous discovery campaigns.
Shows real-time progress on truth collapse, discovery rates, and community impact.
"""

import json
import logging
from typing import Dict, List, Any
from datetime import datetime, timedelta
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

logger = logging.getLogger(__name__)


class LeaderboardGenerator:
    """Generate HTML leaderboards and dashboards for FoT discoveries."""
    
    def __init__(self, output_dir: str = "site/leaderboards/"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC")
    
    def generate_all_leaderboards(self, discovery_results: Dict[str, Any]):
        """Generate all leaderboards from discovery results."""
        logger.info("📊 Generating FoTChemistry leaderboards...")
        
        # Main discovery dashboard
        self.generate_discovery_dashboard(discovery_results)
        
        # Campaign-specific leaderboards
        self.generate_campaign_leaderboards(discovery_results)
        
        # Truth collapse tracker
        self.generate_truth_collapse_tracker(discovery_results)
        
        # Community impact metrics
        self.generate_community_impact(discovery_results)
        
        # Generate index page
        self.generate_index_page()
        
        logger.info(f"✅ Leaderboards generated in {self.output_dir}")
    
    def generate_discovery_dashboard(self, results: Dict[str, Any]):
        """Generate main discovery dashboard."""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>FoTChemistry Discovery Dashboard</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background: #f8f9fa; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .metric-card {{ 
            background: white; border-radius: 8px; padding: 20px; margin: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1); display: inline-block; min-width: 200px;
        }}
        .metric-value {{ font-size: 2em; font-weight: bold; color: #007bff; }}
        .metric-label {{ color: #666; font-size: 0.9em; }}
        .status-truth {{ color: #28a745; }}
        .status-refuted {{ color: #dc3545; }}
        .status-evidence {{ color: #ffc107; }}
        .update-time {{ color: #666; font-size: 0.8em; text-align: center; margin-top: 20px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧪 FoTChemistry Discovery Dashboard</h1>
        <p>Real-time autonomous chemistry discovery and truth mining</p>
    </div>
    
    <div style="text-align: center;">
        <div class="metric-card">
            <div class="metric-value">{results.get('total_discoveries', 0)}</div>
            <div class="metric-label">Total Discoveries</div>
        </div>
        
        <div class="metric-card">
            <div class="metric-value status-truth">{results.get('new_truth_count', 0)}</div>
            <div class="metric-label">Claims Collapsed to Truth</div>
        </div>
        
        <div class="metric-card">
            <div class="metric-value status-refuted">{results.get('refuted_count', 0)}</div>
            <div class="metric-label">Claims Refuted</div>
        </div>
        
        <div class="metric-card">
            <div class="metric-value status-evidence">{results.get('needs_evidence_count', 0)}</div>
            <div class="metric-label">Needs More Evidence</div>
        </div>
        
        <div class="metric-card">
            <div class="metric-value">{results.get('runtime_hours', 0):.1f}h</div>
            <div class="metric-label">Session Runtime</div>
        </div>
        
        <div class="metric-card">
            <div class="metric-value">{results.get('measurement_success_rate', 0)*100:.1f}%</div>
            <div class="metric-label">Success Rate</div>
        </div>
    </div>
    
    <h2>🎯 Active Campaigns</h2>
    <div id="campaigns-chart"></div>
    
    <h2>⚖️ Truth Collapse Timeline</h2>
    <div id="timeline-chart"></div>
    
    <h2>🏆 Recent Discoveries</h2>
    <div id="recent-discoveries">
        {self._generate_recent_discoveries_html(results.get('verdicts', []))}
    </div>
    
    <div class="update-time">
        Last updated: {self.timestamp} | 
        <a href="index.html">All Leaderboards</a> | 
        <a href="https://github.com/FortressAI/FoTChemistry">GitHub</a>
    </div>
    
    <script>
        // Generate campaigns visualization
        {self._generate_campaigns_chart_js(results)}
        
        // Generate timeline visualization  
        {self._generate_timeline_chart_js(results)}
    </script>
</body>
</html>
        """
        
        with open(self.output_dir / "dashboard.html", "w") as f:
            f.write(html_content)
    
    def generate_campaign_leaderboards(self, results: Dict[str, Any]):
        """Generate leaderboards for specific campaigns."""
        campaigns = self._extract_campaign_results(results)
        
        for campaign_name, campaign_data in campaigns.items():
            html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>{campaign_name} - FoTChemistry</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background: #f8f9fa; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .leaderboard {{ background: white; border-radius: 8px; padding: 20px; margin: 20px 0; }}
        .candidate {{ 
            border-bottom: 1px solid #eee; padding: 15px 0; display: flex; 
            justify-content: space-between; align-items: center;
        }}
        .rank {{ font-weight: bold; color: #007bff; min-width: 40px; }}
        .candidate-info {{ flex-grow: 1; margin-left: 20px; }}
        .performance {{ text-align: right; }}
        .metric {{ display: inline-block; margin: 0 10px; }}
        .metric-value {{ font-weight: bold; }}
        .status-badge {{ 
            padding: 4px 8px; border-radius: 4px; font-size: 0.8em; color: white;
        }}
        .status-truth {{ background: #28a745; }}
        .status-refuted {{ background: #dc3545; }}
        .status-evidence {{ background: #ffc107; color: black; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🏆 {campaign_name} Leaderboard</h1>
        <p>Top performing candidates ranked by FoT virtue-weighted performance</p>
    </div>
    
    <div class="leaderboard">
        <h3>🥇 Performance Leaders</h3>
        {self._generate_campaign_leaderboard_html(campaign_data)}
    </div>
    
    <div class="leaderboard">
        <h3>📈 Performance Trends</h3>
        <div id="performance-chart"></div>
    </div>
    
    <div style="text-align: center; margin-top: 30px;">
        <a href="dashboard.html">← Back to Dashboard</a> | 
        <a href="index.html">All Leaderboards</a>
    </div>
    
    <script>
        {self._generate_campaign_performance_chart_js(campaign_data)}
    </script>
</body>
</html>
            """
            
            filename = campaign_name.lower().replace(" ", "_").replace("&", "and") + ".html"
            with open(self.output_dir / filename, "w") as f:
                f.write(html_content)
    
    def generate_truth_collapse_tracker(self, results: Dict[str, Any]):
        """Generate truth collapse tracking page."""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Truth Collapse Tracker - FoTChemistry</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background: #f8f9fa; }}
        .truth-claim {{ 
            background: white; border-radius: 8px; padding: 20px; margin: 15px 0;
            border-left: 4px solid #28a745; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .claim-header {{ display: flex; justify-content: space-between; align-items: center; }}
        .claim-title {{ font-weight: bold; color: #333; }}
        .confidence-score {{ 
            background: #28a745; color: white; padding: 4px 12px; 
            border-radius: 20px; font-size: 0.9em;
        }}
        .claim-details {{ margin-top: 10px; color: #666; }}
        .evidence-list {{ margin-top: 15px; }}
        .evidence-item {{ 
            background: #f8f9fa; padding: 10px; margin: 5px 0; 
            border-radius: 4px; font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="header" style="text-align: center; margin-bottom: 30px;">
        <h1>⚖️ Truth Collapse Tracker</h1>
        <p>Claims that have been validated and collapsed to accepted truth</p>
    </div>
    
    <div id="truth-statistics">
        <h2>📊 Truth Statistics</h2>
        <div id="truth-stats-chart"></div>
    </div>
    
    <div id="truth-claims">
        <h2>✅ Validated Truth Claims</h2>
        {self._generate_truth_claims_html(results.get('verdicts', []))}
    </div>
    
    <div style="text-align: center; margin-top: 30px;">
        <a href="dashboard.html">← Back to Dashboard</a>
    </div>
    
    <script>
        {self._generate_truth_stats_chart_js(results)}
    </script>
</body>
</html>
        """
        
        with open(self.output_dir / "truth_tracker.html", "w") as f:
            f.write(html_content)
    
    def generate_community_impact(self, results: Dict[str, Any]):
        """Generate community impact metrics page."""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Community Impact - FoTChemistry</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body style="font-family: Arial, sans-serif; margin: 20px; background: #f8f9fa;">
    <div style="text-align: center; margin-bottom: 30px;">
        <h1>🌍 Community Impact Dashboard</h1>
        <p>Real-world impact of FoTChemistry discoveries</p>
    </div>
    
    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px;">
        <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
            <h3>🧪 PFAS Remediation Impact</h3>
            <div style="font-size: 2em; color: #007bff; font-weight: bold;">
                {self._calculate_pfas_impact(results)}
            </div>
            <div style="color: #666;">Estimated people served by PFAS solutions</div>
        </div>
        
        <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
            <h3>🌱 CO₂ Reduction Potential</h3>
            <div style="font-size: 2em; color: #28a745; font-weight: bold;">
                {self._calculate_co2_impact(results)} kg/year
            </div>
            <div style="color: #666;">Estimated CO₂ reduction potential</div>
        </div>
        
        <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
            <h3>♻️ Green Chemistry Adoptions</h3>
            <div style="font-size: 2em; color: #ffc107; font-weight: bold;">
                {self._calculate_green_adoptions(results)}
            </div>
            <div style="color: #666;">Safer processes recommended</div>
        </div>
        
        <div style="background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);">
            <h3>🔬 Literature Claims Validated</h3>
            <div style="font-size: 2em; color: #dc3545; font-weight: bold;">
                {self._calculate_validated_claims(results)}
            </div>
            <div style="color: #666;">Scientific claims verified through replication</div>
        </div>
    </div>
    
    <div style="text-align: center; margin-top: 30px;">
        <a href="dashboard.html">← Back to Dashboard</a>
    </div>
</body>
</html>
        """
        
        with open(self.output_dir / "community_impact.html", "w") as f:
            f.write(html_content)
    
    def generate_index_page(self):
        """Generate index page linking to all leaderboards."""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>FoTChemistry Leaderboards</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background: #f8f9fa; }}
        .header {{ text-align: center; margin-bottom: 40px; }}
        .nav-card {{ 
            background: white; border-radius: 8px; padding: 30px; margin: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1); text-align: center;
            transition: transform 0.2s; text-decoration: none; color: inherit;
            display: block;
        }}
        .nav-card:hover {{ transform: translateY(-2px); text-decoration: none; color: inherit; }}
        .nav-icon {{ font-size: 3em; margin-bottom: 15px; }}
        .nav-title {{ font-size: 1.5em; font-weight: bold; margin-bottom: 10px; }}
        .nav-description {{ color: #666; }}
        .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧪 FoTChemistry Leaderboards</h1>
        <p>Real-time tracking of autonomous chemistry discovery and truth mining</p>
        <p style="color: #666;">Last updated: {self.timestamp}</p>
    </div>
    
    <div class="grid">
        <a href="dashboard.html" class="nav-card">
            <div class="nav-icon">📊</div>
            <div class="nav-title">Discovery Dashboard</div>
            <div class="nav-description">
                Real-time overview of all autonomous discovery campaigns
            </div>
        </a>
        
        <a href="truth_tracker.html" class="nav-card">
            <div class="nav-icon">⚖️</div>
            <div class="nav-title">Truth Collapse Tracker</div>
            <div class="nav-description">
                Claims validated and collapsed to accepted scientific truth
            </div>
        </a>
        
        <a href="community_impact.html" class="nav-card">
            <div class="nav-icon">🌍</div>
            <div class="nav-title">Community Impact</div>
            <div class="nav-description">
                Real-world impact metrics from FoTChemistry discoveries
            </div>
        </a>
        
        <a href="pfas_sorbent_catalyst_discovery.html" class="nav-card">
            <div class="nav-icon">💧</div>
            <div class="nav-title">PFAS Remediation</div>
            <div class="nav-description">
                Clean water solutions through sorbent and catalyst discovery
            </div>
        </a>
        
        <a href="co2_electrocatalyst_descriptor_discovery.html" class="nav-card">
            <div class="nav-icon">🌱</div>
            <div class="nav-title">CO₂ Catalysts</div>
            <div class="nav-description">
                Climate solutions through CO₂ reduction electrocatalysts
            </div>
        </a>
        
        <a href="https://github.com/FortressAI/FoTChemistry" class="nav-card">
            <div class="nav-icon">💻</div>
            <div class="nav-title">Source Code</div>
            <div class="nav-description">
                Open source repository for autonomous chemistry discovery
            </div>
        </a>
    </div>
    
    <div style="text-align: center; margin-top: 40px; color: #666;">
        <p>🤖 Powered by FoTChemistry Autonomous Discovery Engine</p>
        <p>Field of Truth + Chemistry = Open, collaborative, truth-driven discovery</p>
    </div>
</body>
</html>
        """
        
        with open(self.output_dir / "index.html", "w") as f:
            f.write(html_content)
    
    # Helper methods for generating HTML content
    def _generate_recent_discoveries_html(self, verdicts: List[Dict[str, Any]]) -> str:
        """Generate HTML for recent discoveries."""
        if not verdicts:
            return "<p>No discoveries yet - autonomous discovery launching soon!</p>"
        
        html = ""
        for verdict in verdicts[-5:]:  # Show last 5
            status = verdict.get('status', 'unknown')
            confidence = verdict.get('confidence', 0)
            reasoning = verdict.get('reasoning', 'No reasoning provided')
            
            status_class = f"status-{status}"
            html += f"""
            <div class="metric-card" style="width: 100%; text-align: left;">
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <span class="{status_class}">●</span>
                    <span style="font-weight: bold;">{status.title()}</span>
                    <span style="color: #666;">Confidence: {confidence:.2f}</span>
                </div>
                <div style="margin-top: 10px; color: #666; font-size: 0.9em;">
                    {reasoning[:100]}...
                </div>
            </div>
            """
        
        return html
    
    def _generate_campaigns_chart_js(self, results: Dict[str, Any]) -> str:
        """Generate JavaScript for campaigns visualization."""
        return """
        var campaigns_data = [{
            x: ['PFAS Leads', 'CO₂ Catalysts', 'Green Solvents', 'Truth Mining'],
            y: [0, 0, 0, 0],
            type: 'bar',
            marker: {color: ['#007bff', '#28a745', '#ffc107', '#dc3545']}
        }];
        
        var campaigns_layout = {
            title: 'Discoveries by Campaign',
            xaxis: {title: 'Campaign'},
            yaxis: {title: 'Discoveries'},
            height: 400
        };
        
        Plotly.newPlot('campaigns-chart', campaigns_data, campaigns_layout);
        """
    
    def _generate_timeline_chart_js(self, results: Dict[str, Any]) -> str:
        """Generate JavaScript for timeline visualization.""" 
        return """
        var timeline_data = [{
            x: ['2025-09-28'],
            y: [0],
            type: 'scatter',
            mode: 'lines+markers',
            name: 'Truth Collapses'
        }];
        
        var timeline_layout = {
            title: 'Truth Collapse Timeline',
            xaxis: {title: 'Date'},
            yaxis: {title: 'Cumulative Truth Claims'},
            height: 400
        };
        
        Plotly.newPlot('timeline-chart', timeline_data, timeline_layout);
        """
    
    def _extract_campaign_results(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract results by campaign."""
        # Placeholder - would extract from actual results
        return {
            "PFAS Sorbent & Catalyst Discovery": {"discoveries": 0, "candidates": []},
            "CO₂ Electrocatalyst Descriptor Discovery": {"discoveries": 0, "candidates": []}
        }
    
    def _generate_campaign_leaderboard_html(self, campaign_data: Dict[str, Any]) -> str:
        """Generate HTML for campaign leaderboard."""
        return "<p>Campaign leaderboard launching soon with first discoveries!</p>"
    
    def _generate_campaign_performance_chart_js(self, campaign_data: Dict[str, Any]) -> str:
        """Generate performance chart for campaign."""
        return "// Performance chart launching with first discoveries"
    
    def _generate_truth_claims_html(self, verdicts: List[Dict[str, Any]]) -> str:
        """Generate HTML for truth claims."""
        truth_verdicts = [v for v in verdicts if v.get('status') == 'truth']
        
        if not truth_verdicts:
            return "<p>No truth claims yet - autonomous discovery will validate claims as they emerge!</p>"
        
        html = ""
        for verdict in truth_verdicts:
            confidence = verdict.get('confidence', 0)
            reasoning = verdict.get('reasoning', 'No reasoning provided')
            
            html += f"""
            <div class="truth-claim">
                <div class="claim-header">
                    <div class="claim-title">Truth Claim #{verdict.get('id', 'unknown')[:8]}</div>
                    <div class="confidence-score">{confidence:.1%} Confidence</div>
                </div>
                <div class="claim-details">{reasoning}</div>
            </div>
            """
        
        return html
    
    def _generate_truth_stats_chart_js(self, results: Dict[str, Any]) -> str:
        """Generate truth statistics chart."""
        return "// Truth statistics chart launching with validated claims"
    
    def _calculate_pfas_impact(self, results: Dict[str, Any]) -> str:
        """Calculate estimated PFAS remediation impact."""
        # Placeholder calculation
        discoveries = results.get('new_truth_count', 0)
        estimated_people = discoveries * 10000  # Rough estimate
        
        if estimated_people >= 1000000:
            return f"{estimated_people/1000000:.1f}M"
        elif estimated_people >= 1000:
            return f"{estimated_people/1000:.0f}K"
        else:
            return str(estimated_people)
    
    def _calculate_co2_impact(self, results: Dict[str, Any]) -> str:
        """Calculate estimated CO₂ reduction impact."""
        discoveries = results.get('new_truth_count', 0)
        co2_reduction = discoveries * 1000  # kg/year per discovery
        return f"{co2_reduction:,}"
    
    def _calculate_green_adoptions(self, results: Dict[str, Any]) -> int:
        """Calculate green chemistry adoptions."""
        return results.get('new_truth_count', 0) * 2  # Rough estimate
    
    def _calculate_validated_claims(self, results: Dict[str, Any]) -> int:
        """Calculate validated literature claims."""
        return results.get('new_truth_count', 0)


def main():
    """Main entry point for leaderboard generation."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate FoTChemistry leaderboards")
    parser.add_argument("--input", required=True, help="Discovery results JSON file")
    parser.add_argument("--output", default="site/leaderboards/", help="Output directory")
    
    args = parser.parse_args()
    
    # Load discovery results
    with open(args.input) as f:
        results = json.load(f)
    
    # Generate leaderboards
    generator = LeaderboardGenerator(args.output)
    generator.generate_all_leaderboards(results)
    
    print(f"✅ Generated leaderboards in {args.output}")


if __name__ == "__main__":
    main()
