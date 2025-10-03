#!/usr/bin/env python3
"""
REALITY FILTERS FOR MOLECULAR DISCOVERIES
Comprehensive screening for synthetic accessibility, safety, and practical utility

🎯 FILTER CATEGORIES:
- Synthetic Accessibility Scoring (SA_Score)
- PAINS (Pan Assay Interference Compounds) alerts
- Structural alerts for toxicity
- ADMET property screening
- Reactive functional group detection
- Lead-likeness assessment

Author: FoT Research Team
Purpose: Filter out impractical or dangerous molecular candidates
"""

import logging
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Fragments
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem import rdRGroupDecomposition as rdRGD
    from rdkit.Chem.Crippen import MolLogP
    from rdkit.Chem.Lipinski import NumHDonors, NumHAcceptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    # Try to import SA_Score (synthetic accessibility)
    from rdkit.Contrib.SA_Score import sascorer
    HAS_SA_SCORE = True
except ImportError:
    HAS_SA_SCORE = False

logger = logging.getLogger(__name__)

@dataclass
class RealityFilterResult:
    """Results of reality filtering"""
    smiles: str
    passes_filters: bool
    synthetic_accessibility_score: float  # 1 = easy, 10 = very difficult
    pains_alerts: List[str]
    structural_alerts: List[str]
    admet_properties: Dict[str, float]
    lead_likeness_score: float
    safety_score: float
    filter_timestamp: datetime
    filter_details: Dict[str, Any]

class RealityFilterEngine:
    """
    Comprehensive reality filtering for molecular discoveries
    """
    
    def __init__(self):
        """Initialize reality filter engine"""
        
        # Synthetic accessibility thresholds
        self.sa_score_threshold = 6.0  # Compounds with SA > 6 are very difficult to synthesize
        
        # ADMET thresholds (drug-like ranges)
        self.admet_thresholds = {
            "molecular_weight": (150, 500),      # Daltons
            "logp": (-2, 5),                     # Partition coefficient
            "hbd": (0, 5),                       # Hydrogen bond donors
            "hba": (0, 10),                      # Hydrogen bond acceptors
            "tpsa": (20, 140),                   # Topological polar surface area
            "rotatable_bonds": (0, 10),          # Rotatable bonds
            "aromatic_rings": (0, 4),            # Aromatic rings
            "aliphatic_rings": (0, 4),           # Aliphatic rings
        }
        
        # Initialize RDKit filter catalogs
        if HAS_RDKIT:
            self._initialize_filter_catalogs()
        
        # Reactive/problematic functional groups (SMARTS patterns)
        self.reactive_groups = {
            "aldehyde": "[CX3H1](=O)[#6]",
            "acyl_chloride": "[CX3](=[OX1])[Cl]",
            "anhydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",
            "isocyanate": "[NX2]=[CX2]=[OX1]",
            "epoxide": "[OX2r3]1[#6r3][#6r3]1",
            "aziridine": "[NX3r3]1[#6r3][#6r3]1",
            "peroxide": "[OX2][OX2]",
            "nitro_aromatic": "[cX3]([nX2+])([nX1-])",
            "quinone": "[#6]1=[#6][#6](=[OX1])[#6]=[#6][#6]1=[OX1]",
            "michael_acceptor": "[CX3]=[CX3][CX3]=[OX1]",
            "alkyl_halide_primary": "[CX4H2][F,Cl,Br,I]",
            "sulfonyl_chloride": "[SX4](=[OX1])(=[OX1])[Cl]",
            "phosphorus_halide": "[PX4]([F,Cl,Br,I])([F,Cl,Br,I])([F,Cl,Br,I])",
        }
        
        # Toxicophore patterns (simplified)
        self.toxicophores = {
            "aromatic_amine": "[cX3][NX3H2]",
            "nitroaromatic": "[cX3][NX3+](=[OX1])[OX1-]",
            "aromatic_nitro": "[cX3][NX3+](=[OX1])[OX1-]",
            "hydrazine": "[NX3][NX3]",
            "azo": "[NX2]=[NX2]",
            "nitroso": "[NX2]=[OX1]",
            "thiourea": "[NX3][CX3](=[SX1])[NX3]",
            "polyhalogen": "[CX4]([F,Cl,Br,I])([F,Cl,Br,I])([F,Cl,Br,I])",
        }
        
        logger.info("✅ Reality filter engine initialized")
    
    def _initialize_filter_catalogs(self):
        """Initialize RDKit filter catalogs"""
        try:
            # PAINS filters
            self.pains_params = FilterCatalogParams()
            self.pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
            self.pains_catalog = FilterCatalog(self.pains_params)
            
            # Brenk filters (for lead-likeness)
            self.brenk_params = FilterCatalogParams()
            self.brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
            self.brenk_catalog = FilterCatalog(self.brenk_params)
            
            logger.info("✅ RDKit filter catalogs initialized")
        except Exception as e:
            logger.warning(f"⚠️ Failed to initialize filter catalogs: {e}")
            self.pains_catalog = None
            self.brenk_catalog = None
    
    def apply_reality_filters(self, smiles: str) -> RealityFilterResult:
        """Apply comprehensive reality filters to a SMILES string"""
        
        if not HAS_RDKIT:
            logger.error("❌ RDKit not available - cannot apply reality filters")
            return self._create_failed_result(smiles, "RDKit not available")
        
        # Parse molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return self._create_failed_result(smiles, "Invalid SMILES")
        
        # Sanitize molecule
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            return self._create_failed_result(smiles, f"Sanitization failed: {e}")
        
        # Apply individual filters
        sa_score = self._calculate_synthetic_accessibility(mol)
        pains_alerts = self._check_pains_alerts(mol)
        structural_alerts = self._check_structural_alerts(mol)
        admet_properties = self._calculate_admet_properties(mol)
        lead_likeness = self._calculate_lead_likeness(mol)
        safety_score = self._calculate_safety_score(mol)
        
        # Determine if molecule passes all filters
        passes_filters = self._evaluate_overall_pass(
            sa_score, pains_alerts, structural_alerts, 
            admet_properties, lead_likeness, safety_score
        )
        
        # Create detailed filter information
        filter_details = {
            "sa_score_pass": sa_score <= self.sa_score_threshold,
            "pains_pass": len(pains_alerts) == 0,
            "structural_alerts_pass": len(structural_alerts) == 0,
            "admet_pass": self._check_admet_pass(admet_properties),
            "lead_likeness_pass": lead_likeness >= 0.5,
            "safety_pass": safety_score >= 0.6,
            "filter_criteria": {
                "max_sa_score": self.sa_score_threshold,
                "admet_thresholds": self.admet_thresholds,
                "min_lead_likeness": 0.5,
                "min_safety_score": 0.6
            }
        }
        
        return RealityFilterResult(
            smiles=smiles,
            passes_filters=passes_filters,
            synthetic_accessibility_score=sa_score,
            pains_alerts=pains_alerts,
            structural_alerts=structural_alerts,
            admet_properties=admet_properties,
            lead_likeness_score=lead_likeness,
            safety_score=safety_score,
            filter_timestamp=datetime.now(),
            filter_details=filter_details
        )
    
    def _calculate_synthetic_accessibility(self, mol) -> float:
        """Calculate synthetic accessibility score"""
        if HAS_SA_SCORE:
            try:
                return sascorer.calculateScore(mol)
            except Exception as e:
                logger.warning(f"⚠️ SA_Score calculation failed: {e}")
        
        # Fallback: simple heuristic based on complexity
        num_atoms = mol.GetNumAtoms()
        num_rings = Descriptors.RingCount(mol)
        num_heteroatoms = Descriptors.NumHeteroatoms(mol)
        num_rotatable = Descriptors.NumRotatableBonds(mol)
        
        # Simple complexity score (higher = more difficult)
        complexity = (
            num_atoms * 0.1 +
            num_rings * 0.5 +
            num_heteroatoms * 0.2 +
            num_rotatable * 0.1
        )
        
        # Convert to SA_Score-like scale (1-10)
        sa_score = min(10.0, max(1.0, complexity))
        return sa_score
    
    def _check_pains_alerts(self, mol) -> List[str]:
        """Check for PAINS (Pan Assay Interference Compounds) alerts"""
        alerts = []
        
        if self.pains_catalog:
            try:
                entries = self.pains_catalog.GetMatches(mol)
                for entry in entries:
                    alerts.append(entry.GetDescription())
            except Exception as e:
                logger.warning(f"⚠️ PAINS check failed: {e}")
        
        return alerts
    
    def _check_structural_alerts(self, mol) -> List[str]:
        """Check for problematic structural features"""
        alerts = []
        
        # Check reactive groups
        for group_name, smarts in self.reactive_groups.items():
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    alerts.append(f"Reactive group: {group_name}")
            except Exception as e:
                logger.warning(f"⚠️ Reactive group check failed for {group_name}: {e}")
        
        # Check toxicophores
        for tox_name, smarts in self.toxicophores.items():
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    alerts.append(f"Toxicophore: {tox_name}")
            except Exception as e:
                logger.warning(f"⚠️ Toxicophore check failed for {tox_name}: {e}")
        
        # Additional structural checks
        if Descriptors.NumRadicalElectrons(mol) > 0:
            alerts.append("Contains radical electrons")
        
        if Descriptors.NumValenceElectrons(mol) % 2 != 0:
            alerts.append("Odd number of valence electrons")
        
        return alerts
    
    def _calculate_admet_properties(self, mol) -> Dict[str, float]:
        """Calculate ADMET-relevant molecular properties"""
        properties = {}
        
        try:
            properties["molecular_weight"] = Descriptors.MolWt(mol)
            properties["logp"] = Descriptors.MolLogP(mol)
            properties["hbd"] = Descriptors.NumHDonors(mol)
            properties["hba"] = Descriptors.NumHAcceptors(mol)
            properties["tpsa"] = Descriptors.TPSA(mol)
            properties["rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
            properties["aromatic_rings"] = Descriptors.NumAromaticRings(mol)
            properties["aliphatic_rings"] = Descriptors.NumAliphaticRings(mol)
            properties["formal_charge"] = Chem.rdmolops.GetFormalCharge(mol)
            properties["num_atoms"] = mol.GetNumAtoms()
            properties["num_heavy_atoms"] = mol.GetNumHeavyAtoms()
            
            # Additional descriptors
            properties["fraction_csp3"] = Descriptors.FractionCsp3(mol)
            properties["num_heterocycles"] = Descriptors.NumHeterocycles(mol)
            properties["num_saturated_rings"] = Descriptors.NumSaturatedRings(mol)
            
        except Exception as e:
            logger.warning(f"⚠️ ADMET property calculation failed: {e}")
        
        return properties
    
    def _calculate_lead_likeness(self, mol) -> float:
        """Calculate lead-likeness score"""
        try:
            # Lead-like criteria (more permissive than drug-like)
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            rotatable = Descriptors.NumRotatableBonds(mol)
            
            score = 1.0
            
            # Molecular weight: 200-450 Da (optimal for leads)
            if not (200 <= mw <= 450):
                score -= 0.2
            
            # LogP: -1 to 4 (slightly more permissive)
            if not (-1 <= logp <= 4):
                score -= 0.2
            
            # Hydrogen bond donors: 0-4
            if hbd > 4:
                score -= 0.2
            
            # Hydrogen bond acceptors: 0-8
            if hba > 8:
                score -= 0.2
            
            # Rotatable bonds: 0-8
            if rotatable > 8:
                score -= 0.2
            
            return max(0.0, score)
            
        except Exception as e:
            logger.warning(f"⚠️ Lead-likeness calculation failed: {e}")
            return 0.5  # Neutral score
    
    def _calculate_safety_score(self, mol) -> float:
        """Calculate safety score based on structural features"""
        safety_score = 1.0
        
        try:
            # Penalize for reactive/toxic groups
            for group_name, smarts in {**self.reactive_groups, **self.toxicophores}.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    if "toxic" in group_name or "reactive" in group_name:
                        safety_score -= 0.3  # Major penalty
                    else:
                        safety_score -= 0.1  # Minor penalty
            
            # Penalize for unusual formal charges
            formal_charge = abs(Chem.rdmolops.GetFormalCharge(mol))
            if formal_charge > 2:
                safety_score -= 0.2
            
            # Penalize for too many heteroatoms (potential for reactivity)
            num_heteroatoms = Descriptors.NumHeteroatoms(mol)
            num_heavy_atoms = mol.GetNumHeavyAtoms()
            if num_heavy_atoms > 0:
                heteroatom_ratio = num_heteroatoms / num_heavy_atoms
                if heteroatom_ratio > 0.5:  # More than 50% heteroatoms
                    safety_score -= 0.2
            
            # Bonus for simple, stable structures
            num_rings = Descriptors.RingCount(mol)
            if num_rings <= 2 and num_heavy_atoms <= 20:
                safety_score += 0.1
            
            return max(0.0, min(1.0, safety_score))
            
        except Exception as e:
            logger.warning(f"⚠️ Safety score calculation failed: {e}")
            return 0.5  # Neutral score
    
    def _check_admet_pass(self, properties: Dict[str, float]) -> bool:
        """Check if ADMET properties are within acceptable ranges"""
        for prop, value in properties.items():
            if prop in self.admet_thresholds:
                min_val, max_val = self.admet_thresholds[prop]
                if not (min_val <= value <= max_val):
                    return False
        return True
    
    def _evaluate_overall_pass(self, sa_score: float, pains_alerts: List[str], 
                              structural_alerts: List[str], admet_properties: Dict[str, float],
                              lead_likeness: float, safety_score: float) -> bool:
        """Evaluate if molecule passes all reality filters"""
        
        # Must pass all critical filters
        if sa_score > self.sa_score_threshold:
            return False
        
        if len(pains_alerts) > 0:
            return False
        
        if len(structural_alerts) > 0:
            return False
        
        if not self._check_admet_pass(admet_properties):
            return False
        
        if lead_likeness < 0.5:
            return False
        
        if safety_score < 0.6:
            return False
        
        return True
    
    def _create_failed_result(self, smiles: str, reason: str) -> RealityFilterResult:
        """Create a failed filter result"""
        return RealityFilterResult(
            smiles=smiles,
            passes_filters=False,
            synthetic_accessibility_score=10.0,  # Maximum difficulty
            pains_alerts=[f"Filter error: {reason}"],
            structural_alerts=[],
            admet_properties={},
            lead_likeness_score=0.0,
            safety_score=0.0,
            filter_timestamp=datetime.now(),
            filter_details={"error": reason}
        )

def filter_discovery_batch(smiles_list: List[str]) -> List[RealityFilterResult]:
    """Apply reality filters to a batch of molecular discoveries"""
    engine = RealityFilterEngine()
    results = []
    
    logger.info(f"🔍 Applying reality filters to {len(smiles_list)} candidates")
    
    for i, smiles in enumerate(smiles_list):
        if (i + 1) % 10 == 0:
            logger.info(f"Progress: {i+1}/{len(smiles_list)}")
        
        result = engine.apply_reality_filters(smiles)
        results.append(result)
    
    # Summary statistics
    passed_count = sum(1 for r in results if r.passes_filters)
    high_sa_count = sum(1 for r in results if r.synthetic_accessibility_score <= 4.0)
    safe_count = sum(1 for r in results if r.safety_score >= 0.8)
    
    logger.info(f"✅ Reality filtering complete:")
    logger.info(f"   📊 Total candidates: {len(results)}")
    logger.info(f"   ✅ Passed all filters: {passed_count} ({passed_count/len(results)*100:.1f}%)")
    logger.info(f"   🧪 Easy to synthesize: {high_sa_count} ({high_sa_count/len(results)*100:.1f}%)")
    logger.info(f"   🛡️ High safety score: {safe_count} ({safe_count/len(results)*100:.1f}%)")
    
    return results

if __name__ == "__main__":
    # Test with known compounds
    test_smiles = [
        "C",                    # Methane (should pass)
        "CCO",                  # Ethanol (should pass)
        "c1ccccc1",            # Benzene (might have alerts)
        "CC(=O)Cl",            # Acetyl chloride (should fail - reactive)
        "c1ccc(N)cc1",         # Aniline (might have toxicophore alert)
        "CCCCCCCCCCCCCCCCCCCC", # Very long chain (might fail SA)
    ]
    
    logging.basicConfig(level=logging.INFO)
    results = filter_discovery_batch(test_smiles)
    
    print("\n" + "="*60)
    print("🔍 REALITY FILTER TEST RESULTS")
    print("="*60)
    
    for result in results:
        print(f"\n{result.smiles}:")
        print(f"  ✅ Passes Filters: {result.passes_filters}")
        print(f"  🧪 SA Score: {result.synthetic_accessibility_score:.2f}")
        print(f"  🛡️ Safety Score: {result.safety_score:.2f}")
        print(f"  📊 Lead-likeness: {result.lead_likeness_score:.2f}")
        if result.pains_alerts:
            print(f"  ⚠️ PAINS Alerts: {result.pains_alerts}")
        if result.structural_alerts:
            print(f"  🚨 Structural Alerts: {result.structural_alerts}")
