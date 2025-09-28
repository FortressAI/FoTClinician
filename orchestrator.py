#!/usr/bin/env python3
"""
FoTChemistry Autonomous Discovery Orchestrator

The brain of the perpetual discovery machine. Coordinates agents to:
1. Scout for signals in data streams
2. Hypothesize testable claims  
3. Generate candidate solutions
4. Measure via computation/experiment
5. Judge results via FoT collapse rules
6. Publish validated discoveries to AKG

This runs continuously, hunting for truth without human prompts.
"""

import json
import time
import pathlib
import random
import logging
import argparse
import glob
import yaml
from typing import Dict, Any, List, Optional
from datetime import datetime, timedelta
from dataclasses import dataclass, asdict

# Import FoT agent modules
from agents.scout.sources import SourceRegistry
from agents.architect.claims import ClaimFactory
from agents.alchemist.generator import CandidateGenerator
from agents.physicist.measure import MeasureQueue
from agents.statistician.evaluate import CollapseRules
from agents.ethics.guard import EthicsGate
from akg.client import AKG

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'logs/discovery_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


@dataclass
class DiscoveryMetrics:
    """Track autonomous discovery performance."""
    discoveries_attempted: int = 0
    discoveries_completed: int = 0
    claims_generated: int = 0
    claims_collapsed: int = 0
    truth_collapsed: int = 0
    refuted_claims: int = 0
    needs_evidence: int = 0
    ethics_blocks: int = 0
    measurement_failures: int = 0
    runtime_seconds: float = 0.0


class AutonomousOrchestrator:
    """
    The autonomous discovery orchestrator - the brain of FoTChemistry.
    
    This class coordinates all agents to run continuous discovery campaigns
    without human intervention, following Field of Truth principles.
    """
    
    def __init__(self, campaigns: List[pathlib.Path], budget_seconds: int, mode: str = "autonomous"):
        """
        Initialize the autonomous discovery system.
        
        Args:
            campaigns: List of campaign configuration files
            budget_seconds: Total time budget for discovery session
            mode: Discovery mode ("autonomous", "interactive", "validation")
        """
        self.mode = mode
        self.budget_seconds = budget_seconds
        self.start_time = time.time()
        self.deadline = self.start_time + budget_seconds
        self.metrics = DiscoveryMetrics()
        
        # Initialize FoT agent ecosystem
        logger.info("🧠 Initializing FoTChemistry autonomous discovery agents...")
        
        try:
            self.akg = AKG()
            self.sources = SourceRegistry()
            self.factory = ClaimFactory()
            self.generator = CandidateGenerator()
            self.measure_queue = MeasureQueue(parallel=4)
            self.collapse_rules = CollapseRules()
            self.ethics_gate = EthicsGate()
            
            logger.info("✅ All agents initialized successfully")
            
        except Exception as e:
            logger.error(f"❌ Failed to initialize agents: {e}")
            raise
        
        # Load and validate campaigns
        self.campaigns = self._load_campaigns(campaigns)
        logger.info(f"📋 Loaded {len(self.campaigns)} discovery campaigns")
        
    def _load_campaigns(self, campaign_paths: List[pathlib.Path]) -> List[Dict[str, Any]]:
        """Load and validate campaign configurations."""
        campaigns = []
        
        for path_pattern in campaign_paths:
            for path in glob.glob(str(path_pattern)):
                try:
                    with open(path, 'r') as f:
                        campaign = yaml.safe_load(f)
                    
                    # Validate required fields
                    required = ['name', 'batch_size', 'objective', 'sources', 'generator', 'measure']
                    missing = [field for field in required if field not in campaign]
                    if missing:
                        logger.warning(f"⚠️ Campaign {path} missing fields: {missing}")
                        continue
                    
                    campaign['config_path'] = path
                    campaigns.append(campaign)
                    logger.info(f"📄 Loaded campaign: {campaign['name']}")
                    
                except Exception as e:
                    logger.error(f"❌ Failed to load campaign {path}: {e}")
        
        return campaigns
    
    def calculate_curiosity_score(self, claim: Dict[str, Any]) -> float:
        """
        Calculate curiosity score for claim prioritization.
        
        Uses expected information gain weighted by virtue and cost:
        curiosity = (uncertainty * virtue_weight) / computational_cost
        
        Args:
            claim: Claim dictionary with uncertainty, virtue_weight, estimated_cost
            
        Returns:
            Curiosity score (higher = more interesting to pursue)
        """
        uncertainty = claim.get("uncertainty", 0.5)
        virtue_weight = claim.get("virtue_weight", 0.5)
        estimated_cost = max(claim.get("estimated_cost", 1.0), 0.1)  # Avoid division by zero
        
        # FoT virtue weighting combines epistemic and moral considerations
        virtue_boost = 0.5 + virtue_weight  # Boost for high-virtue discoveries
        
        curiosity = (uncertainty * virtue_boost) / estimated_cost
        
        logger.debug(f"🤔 Curiosity score: {curiosity:.3f} "
                    f"(u={uncertainty:.2f}, v={virtue_weight:.2f}, c={estimated_cost:.2f})")
        
        return curiosity
    
    def run_discovery_cycle(self, campaign: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Run one complete discovery cycle for a campaign.
        
        Steps:
        1. Scout for new signals
        2. Generate hypotheses/claims
        3. Rank by curiosity
        4. Generate test candidates
        5. Measure candidates
        6. Evaluate and collapse claims
        7. Archive results
        
        Args:
            campaign: Campaign configuration
            
        Returns:
            List of discovery verdicts
        """
        campaign_name = campaign['name']
        logger.info(f"🔬 Starting discovery cycle: {campaign_name}")
        
        try:
            # Step 1: Scout for signals
            logger.info("👁️ Scouting for new signals...")
            signals = self.sources.poll(campaign['sources'])
            logger.info(f"📡 Found {len(signals)} signals")
            
            if not signals:
                logger.info("🤷 No new signals found, skipping cycle")
                return []
            
            # Step 2: Generate claims from signals
            logger.info("💡 Generating hypotheses...")
            claims = []
            for signal in signals:
                try:
                    claim = self.factory.make_claim(signal, campaign)
                    if claim:
                        claims.append(claim)
                        self.metrics.claims_generated += 1
                except Exception as e:
                    logger.warning(f"⚠️ Failed to generate claim from signal: {e}")
            
            logger.info(f"🧠 Generated {len(claims)} testable claims")
            
            if not claims:
                logger.info("🤷 No valid claims generated, skipping cycle")
                return []
            
            # Step 3: Rank by curiosity (FoT-weighted exploration)
            logger.info("🎯 Ranking claims by curiosity...")
            claims_with_scores = [(claim, self.calculate_curiosity_score(claim)) for claim in claims]
            claims_with_scores.sort(key=lambda x: x[1], reverse=True)
            
            # Select top claims within batch size
            batch_size = campaign.get('batch_size', 10)
            selected_claims = [claim for claim, score in claims_with_scores[:batch_size]]
            
            logger.info(f"🎲 Selected {len(selected_claims)} high-curiosity claims for testing")
            
            # Step 4: Generate test candidates
            logger.info("🧪 Generating test candidates...")
            tasks = []
            
            for claim in selected_claims:
                try:
                    candidates = self.generator.propose_candidates(claim, campaign)
                    
                    for candidate in candidates:
                        # Ethics screening
                        if self.ethics_gate.validate_candidate(candidate, claim):
                            tasks.append({
                                "claim": claim,
                                "candidate": candidate,
                                "campaign": campaign_name
                            })
                        else:
                            self.metrics.ethics_blocks += 1
                            logger.warning(f"🚫 Ethics gate blocked candidate")
                            
                except Exception as e:
                    logger.warning(f"⚠️ Failed to generate candidates for claim: {e}")
            
            logger.info(f"⚗️ Generated {len(tasks)} measurement tasks")
            
            if not tasks:
                logger.info("🤷 No valid measurement tasks, skipping cycle")
                return []
            
            # Step 5: Measure candidates
            logger.info("📏 Running measurements...")
            timeout = campaign.get('per_task_timeout', 600)  # 10 min default
            
            try:
                results = self.measure_queue.execute_tasks(tasks, timeout=timeout)
                successful_results = [r for r in results if r.get('success', False)]
                
                logger.info(f"📊 Completed {len(successful_results)}/{len(tasks)} measurements successfully")
                self.metrics.discoveries_completed += len(successful_results)
                self.metrics.measurement_failures += len(tasks) - len(successful_results)
                
            except Exception as e:
                logger.error(f"❌ Measurement execution failed: {e}")
                return []
            
            # Step 6: Evaluate and collapse claims
            logger.info("⚖️ Evaluating results and collapsing claims...")
            verdicts = []
            
            for result in successful_results:
                try:
                    verdict = self.collapse_rules.judge_claim(
                        result['claim'], 
                        result['candidate'], 
                        result['measurement_result']
                    )
                    
                    # Update metrics
                    if verdict['status'] == 'truth':
                        self.metrics.truth_collapsed += 1
                    elif verdict['status'] == 'refuted':
                        self.metrics.refuted_claims += 1
                    else:
                        self.metrics.needs_evidence += 1
                    
                    self.metrics.claims_collapsed += 1
                    verdicts.append(verdict)
                    
                    # Write to AKG with full provenance
                    self.akg.record_discovery_verdict(verdict)
                    
                    logger.info(f"⚖️ Claim verdict: {verdict['status']} "
                              f"(confidence: {verdict.get('confidence', 0):.2f})")
                    
                except Exception as e:
                    logger.warning(f"⚠️ Failed to evaluate result: {e}")
            
            logger.info(f"✅ Discovery cycle complete: {len(verdicts)} verdicts generated")
            return verdicts
            
        except Exception as e:
            logger.error(f"❌ Discovery cycle failed for {campaign_name}: {e}")
            return []
    
    def run_autonomous_discovery(self) -> Dict[str, Any]:
        """
        Run autonomous discovery across all campaigns until budget exhausted.
        
        Returns:
            Discovery session summary
        """
        logger.info(f"🚀 Starting autonomous discovery session")
        logger.info(f"⏰ Budget: {self.budget_seconds/3600:.1f} hours")
        logger.info(f"📋 Campaigns: {len(self.campaigns)}")
        
        session_verdicts = []
        cycles_completed = 0
        
        try:
            while time.time() < self.deadline and self.campaigns:
                remaining_time = self.deadline - time.time()
                logger.info(f"⏱️ Time remaining: {remaining_time/60:.1f} minutes")
                
                # Fair round-robin through campaigns
                random.shuffle(self.campaigns)
                
                for campaign in self.campaigns:
                    if time.time() >= self.deadline:
                        logger.info("⏰ Time budget exhausted")
                        break
                    
                    cycle_verdicts = self.run_discovery_cycle(campaign)
                    session_verdicts.extend(cycle_verdicts)
                    cycles_completed += 1
                    
                    # Brief pause between campaigns for resource management
                    time.sleep(5)
            
            # Calculate final metrics
            self.metrics.runtime_seconds = time.time() - self.start_time
            
            summary = {
                "session_type": "autonomous_discovery",
                "start_time": datetime.fromtimestamp(self.start_time).isoformat(),
                "end_time": datetime.now().isoformat(),
                "runtime_hours": self.metrics.runtime_seconds / 3600,
                "campaigns_run": len(self.campaigns),
                "cycles_completed": cycles_completed,
                "total_discoveries": len(session_verdicts),
                "collapsed_claims": self.metrics.claims_collapsed,
                "new_truth_count": self.metrics.truth_collapsed,
                "refuted_count": self.metrics.refuted_claims,
                "needs_evidence_count": self.metrics.needs_evidence,
                "ethics_blocks": self.metrics.ethics_blocks,
                "measurement_success_rate": (
                    self.metrics.discoveries_completed / 
                    max(self.metrics.discoveries_attempted, 1)
                ),
                "verdicts": session_verdicts,
                "metrics": asdict(self.metrics)
            }
            
            logger.info("🎉 Autonomous discovery session complete!")
            logger.info(f"📊 Discoveries: {len(session_verdicts)}")
            logger.info(f"✅ Truth collapsed: {self.metrics.truth_collapsed}")
            logger.info(f"❌ Claims refuted: {self.metrics.refuted_claims}")
            logger.info(f"🔬 Needs evidence: {self.metrics.needs_evidence}")
            
            return summary
            
        except Exception as e:
            logger.error(f"❌ Autonomous discovery session failed: {e}")
            raise


def parse_time_budget(budget_str: str) -> int:
    """Parse time budget string to seconds."""
    budget_map = {
        "30m": 1800,
        "1h": 3600, 
        "2h": 7200,
        "4h": 14400,
        "8h": 28800
    }
    return budget_map.get(budget_str, 3600)


def main():
    """Main entry point for autonomous discovery orchestrator."""
    parser = argparse.ArgumentParser(
        description="FoTChemistry Autonomous Discovery Orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "--campaigns",
        nargs="+",
        required=True,
        help="Campaign configuration files (glob patterns supported)"
    )
    
    parser.add_argument(
        "--budget", 
        default="2h",
        choices=["30m", "1h", "2h", "4h", "8h"],
        help="Time budget for discovery session"
    )
    
    parser.add_argument(
        "--mode",
        default="autonomous",
        choices=["autonomous", "interactive", "validation"],
        help="Discovery mode"
    )
    
    parser.add_argument(
        "--output",
        default="discovery_results.json",
        help="Output file for discovery results"
    )
    
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )
    
    args = parser.parse_args()
    
    # Configure logging level
    logging.getLogger().setLevel(getattr(logging, args.log_level))
    
    # Create logs directory
    pathlib.Path("logs").mkdir(exist_ok=True)
    
    try:
        # Parse time budget
        budget_seconds = parse_time_budget(args.budget)
        
        # Initialize orchestrator
        orchestrator = AutonomousOrchestrator(
            campaigns=args.campaigns,
            budget_seconds=budget_seconds,
            mode=args.mode
        )
        
        # Run autonomous discovery
        results = orchestrator.run_autonomous_discovery()
        
        # Save results
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        logger.info(f"💾 Results saved to {args.output}")
        
        print(f"\n🎉 Autonomous Discovery Complete!")
        print(f"📊 Total discoveries: {results['total_discoveries']}")
        print(f"✅ Claims collapsed to truth: {results['new_truth_count']}")
        print(f"⏰ Runtime: {results['runtime_hours']:.2f} hours")
        
    except KeyboardInterrupt:
        logger.info("🛑 Discovery interrupted by user")
    except Exception as e:
        logger.error(f"💥 Discovery failed: {e}")
        raise


if __name__ == "__main__":
    main()
