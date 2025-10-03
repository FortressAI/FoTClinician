#!/usr/bin/env python3
"""
VALIDATE OVERNIGHT DISCOVERIES
Extract and validate all molecular candidates from overnight discovery logs

🎯 FUNCTIONS:
- Extract SMILES from discovery logs
- Canonicalize and deduplicate compounds
- Run comprehensive novelty validation
- Generate discovery quality report
- Identify truly novel candidates for further validation

Author: FoT Research Team
Purpose: Convert generated candidates into validated discoveries
"""

import re
import json
import logging
from pathlib import Path
from typing import List, Dict, Any, Set
from collections import Counter
from datetime import datetime

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

from novelty_validation_engine import NoveltyValidationEngine, validate_discovery_batch

logger = logging.getLogger(__name__)

class OvernightDiscoveryValidator:
    """
    Extract and validate overnight molecular discoveries
    """
    
    def __init__(self, log_file: Path = Path("continuous_chemistry_discovery.log")):
        """Initialize validator"""
        self.log_file = log_file
        self.novelty_engine = NoveltyValidationEngine()
        
        # Discovery extraction patterns
        self.discovery_patterns = [
            r"✅ Discovery \d+: ([A-Za-z0-9\(\)\[\]=\-\+\#@\.\\\/:]+) \(score: ([\d\.]+)\)",
            r"Discovery.*?([A-Za-z0-9\(\)\[\]=\-\+\#@\.\\\/:]+).*?score.*?([\d\.]+)",
        ]
        
        logger.info("✅ Overnight discovery validator initialized")
    
    def extract_discoveries_from_log(self) -> List[Dict[str, Any]]:
        """Extract all discoveries from the log file"""
        discoveries = []
        
        if not self.log_file.exists():
            logger.error(f"❌ Log file not found: {self.log_file}")
            return discoveries
        
        logger.info(f"📖 Reading discovery log: {self.log_file}")
        
        try:
            with open(self.log_file, 'r') as f:
                log_content = f.read()
            
            # Extract discoveries using regex patterns
            for pattern in self.discovery_patterns:
                matches = re.findall(pattern, log_content)
                
                for match in matches:
                    smiles = match[0].strip()
                    score = float(match[1])
                    
                    # Basic SMILES validation
                    if self._is_valid_smiles(smiles):
                        discoveries.append({
                            "smiles": smiles,
                            "score": score,
                            "extraction_method": "log_regex"
                        })
            
            logger.info(f"📊 Extracted {len(discoveries)} raw discoveries from log")
            
        except Exception as e:
            logger.error(f"❌ Failed to read log file: {e}")
        
        return discoveries
    
    def canonicalize_and_deduplicate(self, discoveries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Canonicalize SMILES and remove duplicates"""
        if not HAS_RDKIT:
            logger.error("❌ RDKit not available - cannot canonicalize SMILES")
            return discoveries
        
        canonical_discoveries = {}
        invalid_smiles = []
        
        logger.info("🔧 Canonicalizing and deduplicating SMILES...")
        
        for discovery in discoveries:
            smiles = discovery["smiles"]
            
            try:
                # Parse and canonicalize
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    canonical_smiles = Chem.MolToSmiles(mol)
                    
                    # Keep the highest scoring instance of each unique molecule
                    if canonical_smiles not in canonical_discoveries:
                        canonical_discoveries[canonical_smiles] = {
                            "smiles": canonical_smiles,
                            "original_smiles": smiles,
                            "score": discovery["score"],
                            "extraction_method": discovery["extraction_method"],
                            "duplicate_count": 1
                        }
                    else:
                        # Update if this instance has a higher score
                        existing = canonical_discoveries[canonical_smiles]
                        existing["duplicate_count"] += 1
                        if discovery["score"] > existing["score"]:
                            existing["score"] = discovery["score"]
                            existing["original_smiles"] = smiles
                else:
                    invalid_smiles.append(smiles)
            
            except Exception as e:
                logger.warning(f"⚠️ Failed to canonicalize {smiles}: {e}")
                invalid_smiles.append(smiles)
        
        unique_discoveries = list(canonical_discoveries.values())
        
        logger.info(f"✅ Canonicalization complete:")
        logger.info(f"   📊 Raw discoveries: {len(discoveries)}")
        logger.info(f"   🔧 Unique canonical: {len(unique_discoveries)}")
        logger.info(f"   ❌ Invalid SMILES: {len(invalid_smiles)}")
        
        if invalid_smiles:
            logger.warning(f"⚠️ Invalid SMILES found: {invalid_smiles[:5]}...")
        
        return unique_discoveries
    
    def validate_overnight_batch(self) -> Dict[str, Any]:
        """Complete validation of overnight discoveries"""
        logger.info("🚀 Starting overnight discovery validation...")
        
        # Step 1: Extract discoveries from log
        raw_discoveries = self.extract_discoveries_from_log()
        if not raw_discoveries:
            return {"error": "No discoveries found in log file"}
        
        # Step 2: Canonicalize and deduplicate
        unique_discoveries = self.canonicalize_and_deduplicate(raw_discoveries)
        
        # Step 3: Novelty validation
        smiles_list = [d["smiles"] for d in unique_discoveries]
        novelty_results = validate_discovery_batch(smiles_list)
        
        # Step 4: Combine results
        validated_discoveries = []
        for i, discovery in enumerate(unique_discoveries):
            novelty_result = novelty_results[i]
            
            validated_discovery = {
                **discovery,
                "inchi_key": novelty_result.inchi_key,
                "is_novel": novelty_result.is_novel,
                "novelty_score": novelty_result.novelty_score,
                "public_benefit_score": novelty_result.public_benefit_score,
                "database_matches": novelty_result.database_matches,
                "validation_timestamp": novelty_result.validation_timestamp.isoformat()
            }
            validated_discoveries.append(validated_discovery)
        
        # Step 5: Generate summary statistics
        summary = self._generate_validation_summary(validated_discoveries)
        
        # Step 6: Save results
        self._save_validation_results(validated_discoveries, summary)
        
        return {
            "summary": summary,
            "validated_discoveries": validated_discoveries
        }
    
    def _generate_validation_summary(self, discoveries: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Generate comprehensive validation summary"""
        total_count = len(discoveries)
        novel_discoveries = [d for d in discoveries if d["is_novel"]]
        known_discoveries = [d for d in discoveries if not d["is_novel"]]
        high_benefit = [d for d in discoveries if d["public_benefit_score"] > 0.7]
        
        # Score distributions
        scores = [d["score"] for d in discoveries]
        novelty_scores = [d["novelty_score"] for d in discoveries]
        benefit_scores = [d["public_benefit_score"] for d in discoveries]
        
        # Database match analysis
        database_coverage = Counter()
        for discovery in discoveries:
            for match in discovery["database_matches"]:
                database_coverage[match["database"]] += 1
        
        # Top novel candidates
        novel_sorted = sorted(novel_discoveries, 
                            key=lambda x: (x["novelty_score"], x["public_benefit_score"], x["score"]), 
                            reverse=True)
        top_novel = novel_sorted[:10]  # Top 10 novel candidates
        
        summary = {
            "validation_timestamp": datetime.now().isoformat(),
            "total_candidates": total_count,
            "novel_candidates": len(novel_discoveries),
            "known_candidates": len(known_discoveries),
            "novelty_rate": len(novel_discoveries) / total_count if total_count > 0 else 0,
            "high_benefit_candidates": len(high_benefit),
            "score_statistics": {
                "min": min(scores) if scores else 0,
                "max": max(scores) if scores else 0,
                "mean": sum(scores) / len(scores) if scores else 0
            },
            "novelty_statistics": {
                "min": min(novelty_scores) if novelty_scores else 0,
                "max": max(novelty_scores) if novelty_scores else 0,
                "mean": sum(novelty_scores) / len(novelty_scores) if novelty_scores else 0
            },
            "benefit_statistics": {
                "min": min(benefit_scores) if benefit_scores else 0,
                "max": max(benefit_scores) if benefit_scores else 0,
                "mean": sum(benefit_scores) / len(benefit_scores) if benefit_scores else 0
            },
            "database_coverage": dict(database_coverage),
            "top_novel_candidates": [
                {
                    "smiles": d["smiles"],
                    "score": d["score"],
                    "novelty_score": d["novelty_score"],
                    "public_benefit_score": d["public_benefit_score"]
                }
                for d in top_novel
            ]
        }
        
        return summary
    
    def _save_validation_results(self, discoveries: List[Dict[str, Any]], summary: Dict[str, Any]):
        """Save validation results to files"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save full results
        results_file = Path(f"overnight_validation_results_{timestamp}.json")
        with open(results_file, 'w') as f:
            json.dump({
                "summary": summary,
                "discoveries": discoveries
            }, f, indent=2)
        
        # Save summary report
        summary_file = Path(f"overnight_validation_summary_{timestamp}.json")
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Save novel candidates only
        novel_discoveries = [d for d in discoveries if d["is_novel"]]
        if novel_discoveries:
            novel_file = Path(f"novel_candidates_{timestamp}.json")
            with open(novel_file, 'w') as f:
                json.dump(novel_discoveries, f, indent=2)
        
        logger.info(f"💾 Validation results saved:")
        logger.info(f"   📄 Full results: {results_file}")
        logger.info(f"   📊 Summary: {summary_file}")
        if novel_discoveries:
            logger.info(f"   🆕 Novel candidates: {novel_file}")
    
    def _is_valid_smiles(self, smiles: str) -> bool:
        """Basic SMILES validation"""
        if not smiles or len(smiles) < 1:
            return False
        
        # Check for reasonable characters
        valid_chars = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[]=-+#@.\\/:*")
        if not all(c in valid_chars for c in smiles):
            return False
        
        # Try RDKit parsing if available
        if HAS_RDKIT:
            try:
                mol = Chem.MolFromSmiles(smiles)
                return mol is not None
            except:
                return False
        
        return True

def main():
    """Main validation workflow"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    validator = OvernightDiscoveryValidator()
    results = validator.validate_overnight_batch()
    
    if "error" in results:
        logger.error(f"❌ Validation failed: {results['error']}")
        return
    
    summary = results["summary"]
    
    # Print summary report
    print("\n" + "="*60)
    print("🧬 OVERNIGHT DISCOVERY VALIDATION REPORT")
    print("="*60)
    print(f"📊 Total Candidates: {summary['total_candidates']}")
    print(f"🆕 Novel Candidates: {summary['novel_candidates']} ({summary['novelty_rate']*100:.1f}%)")
    print(f"📚 Known Compounds: {summary['known_candidates']}")
    print(f"🌟 High Benefit: {summary['high_benefit_candidates']}")
    print()
    print("📈 Score Statistics:")
    print(f"   Discovery Score: {summary['score_statistics']['min']:.3f} - {summary['score_statistics']['max']:.3f} (avg: {summary['score_statistics']['mean']:.3f})")
    print(f"   Novelty Score: {summary['novelty_statistics']['min']:.3f} - {summary['novelty_statistics']['max']:.3f} (avg: {summary['novelty_statistics']['mean']:.3f})")
    print(f"   Benefit Score: {summary['benefit_statistics']['min']:.3f} - {summary['benefit_statistics']['max']:.3f} (avg: {summary['benefit_statistics']['mean']:.3f})")
    print()
    print("🏆 Top Novel Candidates:")
    for i, candidate in enumerate(summary['top_novel_candidates'][:5], 1):
        print(f"   {i}. {candidate['smiles']} (Novelty: {candidate['novelty_score']:.3f}, Benefit: {candidate['public_benefit_score']:.3f})")
    print()
    print("🎯 NEXT STEPS:")
    if summary['novel_candidates'] > 0:
        print("   1. Review novel candidates for synthetic accessibility")
        print("   2. Run orthogonal computational validation (docking, QSAR)")
        print("   3. Design experimental validation protocols")
        print("   4. Select top 3-5 candidates for wet-lab synthesis")
    else:
        print("   1. No novel candidates found - system needs recalibration")
        print("   2. Adjust generation parameters to explore less-known chemical space")
        print("   3. Implement more sophisticated molecular transformations")
    print("="*60)

if __name__ == "__main__":
    main()
