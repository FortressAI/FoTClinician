#!/usr/bin/env python3
"""
COMPLETE DISCOVERY VALIDATION PIPELINE
End-to-end validation of molecular discoveries for public benefit

🎯 VALIDATION STAGES:
1. Extract discoveries from logs
2. Canonicalize and deduplicate
3. Novelty validation against databases
4. Reality filters (SA, PAINS, ADMET, safety)
5. Public benefit assessment
6. Generate final discovery report

Author: FoT Research Team
Purpose: Transform generated candidates into validated discoveries
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime

from validate_overnight_discoveries import OvernightDiscoveryValidator
from reality_filters import filter_discovery_batch
from novelty_validation_engine import validate_discovery_batch

logger = logging.getLogger(__name__)

class CompleteValidationPipeline:
    """
    Complete validation pipeline for molecular discoveries
    """
    
    def __init__(self, log_file: Path = Path("continuous_chemistry_discovery.log")):
        """Initialize validation pipeline"""
        self.log_file = log_file
        self.validator = OvernightDiscoveryValidator(log_file)
        
        logger.info("✅ Complete validation pipeline initialized")
    
    def run_full_validation(self) -> Dict[str, Any]:
        """Run complete validation pipeline"""
        logger.info("🚀 Starting complete discovery validation pipeline...")
        
        # Stage 1: Extract and canonicalize discoveries
        logger.info("📖 Stage 1: Extracting discoveries from logs...")
        raw_discoveries = self.validator.extract_discoveries_from_log()
        if not raw_discoveries:
            return {"error": "No discoveries found in log file"}
        
        unique_discoveries = self.validator.canonicalize_and_deduplicate(raw_discoveries)
        logger.info(f"✅ Stage 1 complete: {len(unique_discoveries)} unique candidates")
        
        # Stage 2: Novelty validation
        logger.info("🔍 Stage 2: Validating novelty against databases...")
        smiles_list = [d["smiles"] for d in unique_discoveries]
        novelty_results = validate_discovery_batch(smiles_list)
        
        # Filter for novel compounds only
        novel_candidates = []
        for i, discovery in enumerate(unique_discoveries):
            novelty_result = novelty_results[i]
            if novelty_result.is_novel:
                novel_candidates.append({
                    **discovery,
                    "novelty_result": novelty_result
                })
        
        logger.info(f"✅ Stage 2 complete: {len(novel_candidates)} novel candidates")
        
        if not novel_candidates:
            logger.warning("⚠️ No novel candidates found - all compounds are known")
            return self._generate_no_novel_report(unique_discoveries, novelty_results)
        
        # Stage 3: Reality filters
        logger.info("🔍 Stage 3: Applying reality filters...")
        novel_smiles = [c["smiles"] for c in novel_candidates]
        reality_results = filter_discovery_batch(novel_smiles)
        
        # Filter for compounds that pass reality checks
        validated_candidates = []
        for i, candidate in enumerate(novel_candidates):
            reality_result = reality_results[i]
            if reality_result.passes_filters:
                validated_candidates.append({
                    **candidate,
                    "reality_result": reality_result
                })
        
        logger.info(f"✅ Stage 3 complete: {len(validated_candidates)} validated candidates")
        
        # Stage 4: Public benefit assessment and ranking
        logger.info("🌟 Stage 4: Assessing public benefit...")
        final_discoveries = self._assess_and_rank_discoveries(validated_candidates)
        
        # Stage 5: Generate comprehensive report
        logger.info("📊 Stage 5: Generating discovery report...")
        report = self._generate_discovery_report(
            raw_discoveries, unique_discoveries, novel_candidates, 
            validated_candidates, final_discoveries
        )
        
        # Stage 6: Save results
        self._save_validation_results(report)
        
        logger.info("✅ Complete validation pipeline finished!")
        return report
    
    def _assess_and_rank_discoveries(self, candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Assess and rank discoveries by public benefit potential"""
        
        for candidate in candidates:
            # Calculate composite benefit score
            novelty_score = candidate["novelty_result"].novelty_score
            public_benefit = candidate["novelty_result"].public_benefit_score
            safety_score = candidate["reality_result"].safety_score
            lead_likeness = candidate["reality_result"].lead_likeness_score
            sa_score = candidate["reality_result"].synthetic_accessibility_score
            
            # Synthetic accessibility (inverted - lower is better)
            sa_benefit = max(0, (10 - sa_score) / 10)
            
            # Composite public benefit score
            composite_score = (
                novelty_score * 0.25 +          # 25% novelty
                public_benefit * 0.30 +         # 30% public benefit
                safety_score * 0.20 +           # 20% safety
                lead_likeness * 0.15 +          # 15% lead-likeness
                sa_benefit * 0.10               # 10% synthetic accessibility
            )
            
            candidate["composite_benefit_score"] = composite_score
            
            # Assign discovery class
            if composite_score >= 0.8:
                candidate["discovery_class"] = "Exceptional"
            elif composite_score >= 0.7:
                candidate["discovery_class"] = "High Priority"
            elif composite_score >= 0.6:
                candidate["discovery_class"] = "Promising"
            else:
                candidate["discovery_class"] = "Marginal"
        
        # Sort by composite benefit score
        candidates.sort(key=lambda x: x["composite_benefit_score"], reverse=True)
        
        return candidates
    
    def _generate_discovery_report(self, raw_discoveries: List[Dict], 
                                 unique_discoveries: List[Dict],
                                 novel_candidates: List[Dict],
                                 validated_candidates: List[Dict],
                                 final_discoveries: List[Dict]) -> Dict[str, Any]:
        """Generate comprehensive discovery report"""
        
        timestamp = datetime.now()
        
        # Pipeline statistics
        pipeline_stats = {
            "raw_discoveries": len(raw_discoveries),
            "unique_canonical": len(unique_discoveries),
            "novel_candidates": len(novel_candidates),
            "reality_validated": len(validated_candidates),
            "final_discoveries": len(final_discoveries),
            "novelty_rate": len(novel_candidates) / len(unique_discoveries) if unique_discoveries else 0,
            "validation_rate": len(validated_candidates) / len(novel_candidates) if novel_candidates else 0,
            "overall_success_rate": len(final_discoveries) / len(raw_discoveries) if raw_discoveries else 0
        }
        
        # Discovery classification
        discovery_classes = {}
        for discovery in final_discoveries:
            class_name = discovery["discovery_class"]
            if class_name not in discovery_classes:
                discovery_classes[class_name] = []
            discovery_classes[class_name].append(discovery)
        
        # Top discoveries (up to 10)
        top_discoveries = final_discoveries[:10]
        
        # Failure analysis
        failure_reasons = {
            "not_novel": len(unique_discoveries) - len(novel_candidates),
            "failed_reality_filters": len(novel_candidates) - len(validated_candidates),
            "low_benefit_score": 0  # All validated candidates are included
        }
        
        # Generate actionable recommendations
        recommendations = self._generate_recommendations(
            pipeline_stats, discovery_classes, failure_reasons
        )
        
        report = {
            "validation_metadata": {
                "timestamp": timestamp.isoformat(),
                "log_file": str(self.log_file),
                "pipeline_version": "1.0.0",
                "validation_stages": [
                    "extraction_canonicalization",
                    "novelty_validation", 
                    "reality_filters",
                    "benefit_assessment"
                ]
            },
            "pipeline_statistics": pipeline_stats,
            "discovery_classifications": {
                class_name: len(discoveries) 
                for class_name, discoveries in discovery_classes.items()
            },
            "top_discoveries": [
                {
                    "rank": i + 1,
                    "smiles": d["smiles"],
                    "discovery_class": d["discovery_class"],
                    "composite_benefit_score": d["composite_benefit_score"],
                    "novelty_score": d["novelty_result"].novelty_score,
                    "public_benefit_score": d["novelty_result"].public_benefit_score,
                    "safety_score": d["reality_result"].safety_score,
                    "synthetic_accessibility": d["reality_result"].synthetic_accessibility_score,
                    "inchi_key": d["novelty_result"].inchi_key
                }
                for i, d in enumerate(top_discoveries)
            ],
            "failure_analysis": failure_reasons,
            "recommendations": recommendations,
            "detailed_discoveries": final_discoveries,
            "validation_summary": {
                "total_processed": len(raw_discoveries),
                "real_discoveries": len(final_discoveries),
                "success_message": self._generate_success_message(len(final_discoveries), pipeline_stats)
            }
        }
        
        return report
    
    def _generate_no_novel_report(self, unique_discoveries: List[Dict], 
                                novelty_results: List) -> Dict[str, Any]:
        """Generate report when no novel compounds are found"""
        
        # Analyze why no novel compounds were found
        database_matches = {}
        for result in novelty_results:
            for match in result.database_matches:
                db = match["database"]
                if db not in database_matches:
                    database_matches[db] = 0
                database_matches[db] += 1
        
        return {
            "validation_metadata": {
                "timestamp": datetime.now().isoformat(),
                "log_file": str(self.log_file),
                "result": "no_novel_discoveries"
            },
            "pipeline_statistics": {
                "raw_discoveries": len(unique_discoveries),
                "unique_canonical": len(unique_discoveries),
                "novel_candidates": 0,
                "novelty_rate": 0.0
            },
            "database_coverage": database_matches,
            "recommendations": [
                "Adjust molecular generation parameters to explore less-known chemical space",
                "Implement more sophisticated molecular transformations",
                "Consider targeting specific therapeutic areas or material properties",
                "Increase diversity in seed molecules and generation strategies",
                "Review scoring function to reward structural novelty"
            ],
            "calibration_assessment": {
                "engine_status": "functional",
                "convergence_behavior": "converging to known compounds",
                "next_steps": "recalibrate for novelty exploration"
            }
        }
    
    def _generate_recommendations(self, stats: Dict, classes: Dict, 
                                failures: Dict) -> List[str]:
        """Generate actionable recommendations based on results"""
        recommendations = []
        
        if stats["final_discoveries"] == 0:
            recommendations.extend([
                "No validated discoveries found - system needs recalibration",
                "Consider relaxing reality filter thresholds for initial exploration",
                "Implement more diverse molecular generation strategies"
            ])
        elif stats["final_discoveries"] < 5:
            recommendations.extend([
                "Low discovery yield - optimize generation parameters",
                "Focus on specific therapeutic targets or material properties",
                "Consider expanding chemical space exploration"
            ])
        else:
            recommendations.extend([
                f"Excellent yield: {stats['final_discoveries']} validated discoveries",
                "Prioritize top-ranked candidates for experimental validation",
                "Consider scaling up discovery campaigns"
            ])
        
        if stats["novelty_rate"] < 0.05:  # Less than 5% novel
            recommendations.append("Very low novelty rate - adjust generation to explore unknown chemical space")
        
        if failures["failed_reality_filters"] > failures["not_novel"]:
            recommendations.append("Many novel candidates failed reality filters - consider adjusting filter thresholds")
        
        return recommendations
    
    def _generate_success_message(self, discovery_count: int, stats: Dict) -> str:
        """Generate appropriate success message"""
        if discovery_count == 0:
            return "No validated discoveries found. System calibration successful - converged to known compounds."
        elif discovery_count < 5:
            return f"Found {discovery_count} validated discoveries. Promising start - optimize for higher yield."
        else:
            return f"Excellent results: {discovery_count} validated discoveries ready for experimental validation!"
    
    def _save_validation_results(self, report: Dict[str, Any]):
        """Save validation results to files"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save complete report
        report_file = Path(f"complete_validation_report_{timestamp}.json")
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        
        # Save top discoveries for easy access
        if report["top_discoveries"]:
            top_file = Path(f"top_discoveries_{timestamp}.json")
            with open(top_file, 'w') as f:
                json.dump(report["top_discoveries"], f, indent=2)
        
        # Save summary for quick review
        summary = {
            "timestamp": report["validation_metadata"]["timestamp"],
            "total_discoveries": report["pipeline_statistics"]["final_discoveries"],
            "success_rate": report["pipeline_statistics"]["overall_success_rate"],
            "top_3_discoveries": report["top_discoveries"][:3],
            "recommendations": report["recommendations"]
        }
        
        summary_file = Path(f"validation_summary_{timestamp}.json")
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"💾 Validation results saved:")
        logger.info(f"   📄 Complete report: {report_file}")
        logger.info(f"   🏆 Top discoveries: {top_file if report['top_discoveries'] else 'None'}")
        logger.info(f"   📊 Summary: {summary_file}")

def main():
    """Main validation workflow"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    pipeline = CompleteValidationPipeline()
    report = pipeline.run_full_validation()
    
    if "error" in report:
        logger.error(f"❌ Validation failed: {report['error']}")
        return
    
    # Print executive summary
    print("\n" + "="*80)
    print("🧬 FOTCHEMISTRY DISCOVERY VALIDATION REPORT")
    print("="*80)
    
    stats = report["pipeline_statistics"]
    print(f"📊 PIPELINE STATISTICS:")
    print(f"   Raw discoveries extracted: {stats['raw_discoveries']}")
    print(f"   Unique canonical compounds: {stats['unique_canonical']}")
    print(f"   Novel candidates: {stats['novel_candidates']} ({stats['novelty_rate']*100:.1f}%)")
    print(f"   Reality-validated: {stats['reality_validated']}")
    print(f"   Final discoveries: {stats['final_discoveries']}")
    print(f"   Overall success rate: {stats['overall_success_rate']*100:.1f}%")
    
    if report["top_discoveries"]:
        print(f"\n🏆 TOP DISCOVERIES:")
        for discovery in report["top_discoveries"][:5]:
            print(f"   {discovery['rank']}. {discovery['smiles']}")
            print(f"      Class: {discovery['discovery_class']}")
            print(f"      Benefit Score: {discovery['composite_benefit_score']:.3f}")
            print(f"      Novelty: {discovery['novelty_score']:.3f} | Safety: {discovery['safety_score']:.3f}")
    
    print(f"\n🎯 RECOMMENDATIONS:")
    for rec in report["recommendations"]:
        print(f"   • {rec}")
    
    print(f"\n✅ {report['validation_summary']['success_message']}")
    print("="*80)

if __name__ == "__main__":
    main()
