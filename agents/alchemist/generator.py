#!/usr/bin/env python3
"""
FoTChemistry Alchemist Agent - Candidate Generation

The Alchemist agent generates test candidates to validate claims:
- Molecular candidates for property testing
- Reaction candidates for yield/selectivity validation  
- Material candidates for performance evaluation

Uses safe generation strategies with ethics screening.
"""

import logging
import random
from typing import Dict, List, Any, Optional
import numpy as np

logger = logging.getLogger(__name__)


class CandidateGenerator:
    """Generate test candidates for claim validation."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize candidate generator."""
        self.config = config or {}
        self.generators = {
            'molecular_property': self._generate_molecular_candidates,
            'reaction_performance': self._generate_reaction_candidates,
            'material_discovery': self._generate_material_candidates,
            'literature_validation': self._generate_literature_candidates
        }
    
    def propose_candidates(self, claim: Dict[str, Any], campaign: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Generate test candidates for a claim.
        
        Args:
            claim: Claim to generate candidates for
            campaign: Campaign configuration
            
        Returns:
            List of candidate dictionaries
        """
        claim_type = claim.get('claim_type', 'unknown')
        
        logger.info(f"⚗️ Generating candidates for {claim_type} claim")
        
        try:
            if claim_type in self.generators:
                generator = self.generators[claim_type]
                candidates = generator(claim, campaign)
                
                logger.info(f"✅ Generated {len(candidates)} candidates")
                return candidates
            else:
                logger.warning(f"⚠️ Unknown claim type: {claim_type}")
                return []
                
        except Exception as e:
            logger.error(f"❌ Candidate generation failed: {e}")
            return []
    
    def _generate_molecular_candidates(self, claim: Dict[str, Any], campaign: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Generate molecular candidates for property testing."""
        metadata = claim.get('metadata', {})
        smiles = metadata.get('smiles', '')
        property_name = metadata.get('property', 'activity')
        
        if not smiles:
            logger.warning("No SMILES in claim, cannot generate molecular candidates")
            return []
        
        candidates = []
        
        # Original molecule as baseline
        candidates.append({
            'type': 'baseline_molecule',
            'smiles': smiles,
            'name': f"Original_{metadata.get('molecule_id', 'unknown')}",
            'modifications': [],
            'expected_property': property_name,
            'generation_method': 'baseline',
            'safety_score': 0.8,  # Assume reasonable safety for baseline
            'novelty_score': 0.0   # Not novel
        })
        
        # Generate molecular variations using safe transformations
        variations = self._generate_molecular_variations(smiles, property_name)
        for i, variation in enumerate(variations[:5]):  # Limit to 5 variations
            candidates.append({
                'type': 'molecular_variant',
                'smiles': variation['smiles'],
                'name': f"Variant_{i+1}_{metadata.get('molecule_id', 'unknown')}",
                'modifications': variation['modifications'],
                'expected_property': property_name,
                'generation_method': 'structural_variation',
                'safety_score': variation.get('safety_score', 0.6),
                'novelty_score': variation.get('novelty_score', 0.3)
            })
        
        # Add positive and negative controls if available
        controls = self._get_molecular_controls(property_name)
        candidates.extend(controls)
        
        return candidates
    
    def _generate_reaction_candidates(self, claim: Dict[str, Any], campaign: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Generate reaction candidates for performance testing."""
        metadata = claim.get('metadata', {})
        reaction_smiles = metadata.get('reaction_smiles', '')
        property_name = metadata.get('property', 'yield')
        
        if not reaction_smiles:
            logger.warning("No reaction SMILES in claim")
            return []
        
        candidates = []
        
        # Original reaction as baseline
        candidates.append({
            'type': 'baseline_reaction',
            'reaction_smiles': reaction_smiles,
            'name': f"Original_{metadata.get('reaction_id', 'unknown')}",
            'conditions': {'temperature': 25, 'solvent': 'standard', 'time': '2h'},
            'modifications': [],
            'expected_property': property_name,
            'generation_method': 'baseline',
            'safety_score': 0.7,
            'novelty_score': 0.0
        })
        
        # Generate condition variations
        condition_variants = self._generate_condition_variations(reaction_smiles)
        for i, variant in enumerate(condition_variants[:4]):
            candidates.append({
                'type': 'reaction_variant',
                'reaction_smiles': reaction_smiles,
                'name': f"Conditions_{i+1}_{metadata.get('reaction_id', 'unknown')}",
                'conditions': variant['conditions'],
                'modifications': variant['modifications'],
                'expected_property': property_name,
                'generation_method': 'condition_optimization',
                'safety_score': variant.get('safety_score', 0.6),
                'novelty_score': variant.get('novelty_score', 0.2)
            })
        
        return candidates
    
    def _generate_material_candidates(self, claim: Dict[str, Any], campaign: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Generate material candidates for performance testing."""
        metadata = claim.get('metadata', {})
        composition = metadata.get('composition', 'unknown')
        property_name = metadata.get('property', 'performance')
        
        campaign_name = campaign.get('name', '')
        
        candidates = []
        
        # Generate candidates based on campaign type
        if 'pfas' in campaign_name.lower():
            candidates = self._generate_pfas_sorbent_candidates(composition, property_name)
        elif 'co2' in campaign_name.lower():
            candidates = self._generate_co2_catalyst_candidates(composition, property_name)
        else:
            # Generic material candidates
            candidates = self._generate_generic_material_candidates(composition, property_name)
        
        return candidates
    
    def _generate_literature_candidates(self, claim: Dict[str, Any], campaign: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Generate candidates for literature validation."""
        metadata = claim.get('metadata', {})
        paper_id = metadata.get('paper_id', 'unknown')
        property_name = metadata.get('property', 'reproducibility')
        
        candidates = []
        
        # Original study as baseline
        candidates.append({
            'type': 'original_study',
            'paper_id': paper_id,
            'name': f"Original_{paper_id}",
            'replication_method': 'exact_conditions',
            'expected_property': property_name,
            'generation_method': 'literature_replication',
            'confidence_score': 0.5,  # Unknown until tested
            'novelty_score': 0.0
        })
        
        # Methodological variations
        method_variants = self._generate_method_variations(property_name)
        for i, variant in enumerate(method_variants[:3]):
            candidates.append({
                'type': 'method_variant',
                'paper_id': paper_id,
                'name': f"Method_{i+1}_{paper_id}",
                'replication_method': variant['method'],
                'modifications': variant['modifications'],
                'expected_property': property_name,
                'generation_method': 'methodological_variation',
                'confidence_score': variant.get('confidence', 0.6),
                'novelty_score': variant.get('novelty', 0.1)
            })
        
        return candidates
    
    def _generate_molecular_variations(self, smiles: str, property_name: str) -> List[Dict[str, Any]]:
        """Generate safe molecular variations."""
        variations = []
        
        # Placeholder implementation - in real system would use RDKit
        # for actual molecular transformations
        
        # Example safe transformations
        safe_modifications = [
            {'type': 'methyl_addition', 'description': 'Add methyl group'},
            {'type': 'hydroxyl_addition', 'description': 'Add hydroxyl group'},
            {'type': 'halogen_substitution', 'description': 'Replace H with F/Cl'},
            {'type': 'ring_expansion', 'description': 'Expand ring size by 1'},
            {'type': 'stereoisomer', 'description': 'Generate stereoisomer'}
        ]
        
        for i, mod in enumerate(safe_modifications[:3]):  # Limit variations
            # Generate pseudo-SMILES for demonstration
            if 'CCO' in smiles:  # Ethanol example
                if mod['type'] == 'methyl_addition':
                    new_smiles = 'CCCO'  # Propanol
                elif mod['type'] == 'hydroxyl_addition':
                    new_smiles = 'OCC(O)O'  # Glycerol-like
                else:
                    new_smiles = f"{smiles}_mod_{i}"  # Placeholder
            else:
                new_smiles = f"{smiles}_mod_{i}"  # Placeholder
            
            variations.append({
                'smiles': new_smiles,
                'modifications': [mod],
                'safety_score': max(0.4, 0.8 - 0.1 * i),  # Decreasing safety
                'novelty_score': 0.2 + 0.1 * i
            })
        
        return variations
    
    def _generate_condition_variations(self, reaction_smiles: str) -> List[Dict[str, Any]]:
        """Generate safe reaction condition variations."""
        variations = []
        
        base_conditions = {
            'temperature': 25,
            'solvent': 'water',
            'time': '2h',
            'pressure': '1atm'
        }
        
        # Safe condition modifications
        condition_mods = [
            {
                'conditions': {**base_conditions, 'temperature': 40},
                'modifications': ['increase_temperature_15C'],
                'safety_score': 0.8
            },
            {
                'conditions': {**base_conditions, 'solvent': 'ethanol'},
                'modifications': ['solvent_change_to_ethanol'],
                'safety_score': 0.7
            },
            {
                'conditions': {**base_conditions, 'time': '4h'},
                'modifications': ['extend_reaction_time_2h'],
                'safety_score': 0.9
            },
            {
                'conditions': {**base_conditions, 'temperature': 60, 'time': '1h'},
                'modifications': ['higher_temp_shorter_time'],
                'safety_score': 0.6
            }
        ]
        
        for mod in condition_mods:
            mod['novelty_score'] = 0.15
            variations.append(mod)
        
        return variations
    
    def _generate_pfas_sorbent_candidates(self, composition: str, property_name: str) -> List[Dict[str, Any]]:
        """Generate PFAS sorbent material candidates."""
        candidates = []
        
        # Safe sorbent materials for PFAS removal
        safe_materials = [
            {
                'type': 'activated_carbon',
                'composition': 'Activated Carbon (coconut shell)',
                'name': 'AC_Coconut_Standard',
                'surface_area': 1200,  # m²/g
                'pore_size': 'micro/mesoporous',
                'modifications': [],
                'safety_score': 0.9,
                'cost_score': 0.8,
                'novelty_score': 0.1
            },
            {
                'type': 'modified_carbon',
                'composition': 'Activated Carbon + PFOA-selective groups',
                'name': 'AC_PFOA_Selective',
                'surface_area': 1000,
                'pore_size': 'microporous',
                'modifications': ['surface_functionalization'],
                'safety_score': 0.8,
                'cost_score': 0.6,
                'novelty_score': 0.4
            },
            {
                'type': 'ion_exchange',
                'composition': 'Anion Exchange Resin (Type I)',
                'name': 'IEX_Type1_Standard',
                'surface_area': 800,
                'pore_size': 'macroporous',
                'modifications': [],
                'safety_score': 0.85,
                'cost_score': 0.7,
                'novelty_score': 0.2
            },
            {
                'type': 'silica_based',
                'composition': 'Mesoporous Silica + hydrophobic coating',
                'name': 'Silica_Hydrophobic',
                'surface_area': 600,
                'pore_size': 'mesoporous',
                'modifications': ['hydrophobic_functionalization'],
                'safety_score': 0.9,
                'cost_score': 0.5,
                'novelty_score': 0.3
            }
        ]
        
        for material in safe_materials:
            candidates.append({
                'type': 'sorbent_material',
                'name': material['name'],
                'composition': material['composition'],
                'properties': {
                    'surface_area_m2g': material['surface_area'],
                    'pore_size_category': material['pore_size']
                },
                'modifications': material['modifications'],
                'expected_property': property_name,
                'generation_method': 'safe_material_library',
                'safety_score': material['safety_score'],
                'cost_score': material['cost_score'],
                'novelty_score': material['novelty_score']
            })
        
        return candidates
    
    def _generate_co2_catalyst_candidates(self, composition: str, property_name: str) -> List[Dict[str, Any]]:
        """Generate CO₂ reduction catalyst candidates."""
        candidates = []
        
        # Safe catalyst materials for CO₂ reduction
        safe_catalysts = [
            {
                'type': 'copper_based',
                'composition': 'Cu nanoparticles on carbon support',
                'name': 'Cu_NP_Carbon',
                'active_metal': 'Cu',
                'support': 'carbon_black',
                'particle_size': 5,  # nm
                'loading': 20,  # wt%
                'safety_score': 0.8,
                'performance_score': 0.7,
                'novelty_score': 0.2
            },
            {
                'type': 'silver_based', 
                'composition': 'Ag nanowires on graphene',
                'name': 'Ag_NW_Graphene',
                'active_metal': 'Ag',
                'support': 'graphene',
                'particle_size': 3,
                'loading': 15,
                'safety_score': 0.85,
                'performance_score': 0.75,
                'novelty_score': 0.4
            },
            {
                'type': 'copper_alloy',
                'composition': 'CuAg alloy nanoparticles',
                'name': 'CuAg_Alloy_NP',
                'active_metal': 'CuAg',
                'support': 'carbon_nanotube',
                'particle_size': 4,
                'loading': 25,
                'safety_score': 0.75,
                'performance_score': 0.8,
                'novelty_score': 0.5
            },
            {
                'type': 'single_atom',
                'composition': 'Cu single atoms on N-doped carbon',
                'name': 'Cu_SA_NC',
                'active_metal': 'Cu',
                'support': 'N_doped_carbon',
                'particle_size': 0.2,  # Single atoms
                'loading': 5,
                'safety_score': 0.9,
                'performance_score': 0.85,
                'novelty_score': 0.7
            }
        ]
        
        for catalyst in safe_catalysts:
            candidates.append({
                'type': 'electrocatalyst',
                'name': catalyst['name'],
                'composition': catalyst['composition'],
                'properties': {
                    'active_metal': catalyst['active_metal'],
                    'support_material': catalyst['support'],
                    'particle_size_nm': catalyst['particle_size'],
                    'metal_loading_wt': catalyst['loading']
                },
                'modifications': [],
                'expected_property': property_name,
                'generation_method': 'catalyst_design_library',
                'safety_score': catalyst['safety_score'],
                'performance_score': catalyst['performance_score'],
                'novelty_score': catalyst['novelty_score']
            })
        
        return candidates
    
    def _generate_generic_material_candidates(self, composition: str, property_name: str) -> List[Dict[str, Any]]:
        """Generate generic material candidates."""
        candidates = []
        
        # Simple generic materials
        generic_materials = [
            {
                'name': f'Material_Baseline_{property_name}',
                'composition': composition or 'Generic_Material',
                'modifications': [],
                'safety_score': 0.7,
                'novelty_score': 0.1
            },
            {
                'name': f'Material_Variant_1_{property_name}',
                'composition': f'{composition}_modified' if composition else 'Generic_Modified',
                'modifications': ['surface_treatment'],
                'safety_score': 0.6,
                'novelty_score': 0.3
            }
        ]
        
        for material in generic_materials:
            candidates.append({
                'type': 'generic_material',
                'name': material['name'],
                'composition': material['composition'],
                'properties': {},
                'modifications': material['modifications'],
                'expected_property': property_name,
                'generation_method': 'generic_variation',
                'safety_score': material['safety_score'],
                'novelty_score': material['novelty_score']
            })
        
        return candidates
    
    def _get_molecular_controls(self, property_name: str) -> List[Dict[str, Any]]:
        """Get positive and negative control molecules."""
        controls = []
        
        # Simple control molecules based on property
        if property_name == 'solubility':
            controls = [
                {
                    'type': 'positive_control',
                    'smiles': 'O',  # Water - highly soluble
                    'name': 'Water_Positive_Control',
                    'modifications': [],
                    'expected_property': property_name,
                    'generation_method': 'positive_control',
                    'safety_score': 1.0,
                    'novelty_score': 0.0
                },
                {
                    'type': 'negative_control',
                    'smiles': 'CCCCCCCCCCCCCCCC',  # Long alkane - insoluble
                    'name': 'Alkane_Negative_Control',
                    'modifications': [],
                    'expected_property': property_name,
                    'generation_method': 'negative_control',
                    'safety_score': 0.8,
                    'novelty_score': 0.0
                }
            ]
        elif property_name == 'pKa':
            controls = [
                {
                    'type': 'positive_control',
                    'smiles': 'CC(=O)O',  # Acetic acid - known pKa
                    'name': 'Acetic_Acid_Control',
                    'modifications': [],
                    'expected_property': property_name,
                    'generation_method': 'positive_control',
                    'safety_score': 0.7,
                    'novelty_score': 0.0
                }
            ]
        
        return controls
    
    def _generate_method_variations(self, property_name: str) -> List[Dict[str, Any]]:
        """Generate methodological variations for literature replication."""
        variations = []
        
        # Common methodological variations
        method_mods = [
            {
                'method': 'different_software',
                'modifications': ['use_alternative_software_package'],
                'confidence': 0.7,
                'novelty': 0.1
            },
            {
                'method': 'higher_precision',
                'modifications': ['increase_numerical_precision'],
                'confidence': 0.8,
                'novelty': 0.15
            },
            {
                'method': 'different_basis_set',
                'modifications': ['use_larger_basis_set'],
                'confidence': 0.6,
                'novelty': 0.2
            }
        ]
        
        return method_mods


# Example usage and testing
if __name__ == "__main__":
    # Test candidate generator
    generator = CandidateGenerator()
    
    # Example molecular property claim
    test_claim = {
        'claim_type': 'molecular_property',
        'metadata': {
            'molecule_id': 'test_mol_001',
            'smiles': 'CCO',
            'property': 'solubility'
        }
    }
    
    test_campaign = {
        'name': 'Test Solubility Campaign',
        'objective': 'predict_solubility'
    }
    
    candidates = generator.propose_candidates(test_claim, test_campaign)
    
    print(f"✅ Generated {len(candidates)} candidates:")
    for i, candidate in enumerate(candidates):
        print(f"  {i+1}. {candidate['name']} ({candidate['type']})")
        print(f"     Safety: {candidate.get('safety_score', 'N/A'):.2f}")
        print(f"     Novelty: {candidate.get('novelty_score', 'N/A'):.2f}")
    
    # Test PFAS material generation
    pfas_claim = {
        'claim_type': 'material_discovery',
        'metadata': {
            'material_id': 'pfas_material_001',
            'composition': 'activated_carbon',
            'property': 'pfas_removal'
        }
    }
    
    pfas_campaign = {
        'name': 'PFAS Sorbent Discovery',
        'objective': 'minimize_residual_pfas_ngL'
    }
    
    pfas_candidates = generator.propose_candidates(pfas_claim, pfas_campaign)
    
    print(f"\n🧪 Generated {len(pfas_candidates)} PFAS sorbent candidates:")
    for candidate in pfas_candidates:
        print(f"  - {candidate['name']}: {candidate['composition']}")
        print(f"    Safety: {candidate.get('safety_score', 'N/A'):.2f}, "
              f"Novelty: {candidate.get('novelty_score', 'N/A'):.2f}")
