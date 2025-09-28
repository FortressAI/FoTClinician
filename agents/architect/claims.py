#!/usr/bin/env python3
"""
FoTChemistry Architect Agent - Claim Generation

The Architect agent converts discovery signals into structured, testable claims
with uncertainty quantification and virtue weighting according to FoT principles.
"""

import logging
import uuid
from typing import Dict, List, Any, Optional
from datetime import datetime
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class Claim:
    """A structured scientific claim with FoT metadata."""
    id: str
    signal_source: str
    claim_type: str
    hypothesis: str
    testable_prediction: str
    uncertainty: float  # 0.0 to 1.0
    virtue_weight: float  # 0.0 to 1.0
    estimated_cost: float  # Computational/experimental cost
    priority: str  # "high", "medium", "low"
    collapse_rules: Dict[str, Any]
    metadata: Dict[str, Any]
    created: datetime


class ClaimFactory:
    """Factory for generating testable claims from discovery signals."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize claim factory with configuration."""
        self.config = config or {}
        self.claim_templates = {
            'molecular_property': self._generate_property_claim,
            'reaction_performance': self._generate_reaction_claim,
            'material_discovery': self._generate_material_claim,
            'literature_validation': self._generate_literature_claim
        }
    
    def make_claim(self, signal: Dict[str, Any], campaign: Dict[str, Any]) -> Optional[Claim]:
        """
        Generate a testable claim from a discovery signal.
        
        Args:
            signal: Discovery signal from Scout agent
            campaign: Campaign configuration
            
        Returns:
            Structured claim ready for testing
        """
        try:
            signal_type = signal.get('signal_type', 'unknown')
            claim_type = self._infer_claim_type(signal, campaign)
            
            logger.debug(f"🏗️ Generating {claim_type} claim from {signal_type} signal")
            
            # Route to appropriate claim generator
            if claim_type in self.claim_templates:
                generator = self.claim_templates[claim_type]
                claim = generator(signal, campaign)
                
                if claim:
                    logger.info(f"✅ Generated claim: {claim.id}")
                    return claim
                else:
                    logger.warning(f"⚠️ Failed to generate {claim_type} claim")
                    return None
            else:
                logger.warning(f"⚠️ Unknown claim type: {claim_type}")
                return None
                
        except Exception as e:
            logger.error(f"❌ Claim generation failed: {e}")
            return None
    
    def _infer_claim_type(self, signal: Dict[str, Any], campaign: Dict[str, Any]) -> str:
        """Infer the type of claim to generate based on signal and campaign."""
        signal_type = signal.get('signal_type', '')
        campaign_objective = campaign.get('objective', '')
        
        # Map signal types to claim types
        if 'molecule' in signal_type or 'molecular' in campaign_objective:
            return 'molecular_property'
        elif 'reaction' in signal_type or 'reaction' in campaign_objective:
            return 'reaction_performance'
        elif 'material' in signal_type or 'catalyst' in campaign_objective:
            return 'material_discovery'
        elif 'literature' in signal_type or 'paper' in signal_type:
            return 'literature_validation'
        else:
            # Default to molecular property for chemistry
            return 'molecular_property'
    
    def _generate_property_claim(self, signal: Dict[str, Any], campaign: Dict[str, Any]) -> Optional[Claim]:
        """Generate a molecular property prediction claim."""
        content = signal.get('content', {})
        
        # Extract molecular information
        molecule_id = content.get('id', 'unknown')
        smiles = content.get('smiles', '')
        
        if not smiles:
            logger.warning("No SMILES found in signal, cannot generate property claim")
            return None
        
        # Determine property to predict based on campaign
        objective = campaign.get('objective', '')
        if 'solubility' in objective:
            property_name = 'solubility'
            property_unit = 'mg/mL'
            testable_prediction = f"Molecule {molecule_id} has aqueous solubility >1 mg/mL"
        elif 'permeability' in objective:
            property_name = 'permeability'
            property_unit = 'cm/s'
            testable_prediction = f"Molecule {molecule_id} has membrane permeability >1e-6 cm/s"
        elif 'pka' in objective.lower():
            property_name = 'pKa'
            property_unit = 'pH units'
            testable_prediction = f"Molecule {molecule_id} has pKa between 6-8"
        else:
            property_name = 'activity'
            property_unit = 'relative'
            testable_prediction = f"Molecule {molecule_id} shows target activity >50%"
        
        # Calculate uncertainty based on molecular complexity
        molecular_weight = content.get('molecular_weight', 300)  # Default estimate
        complexity_factor = min(molecular_weight / 500, 1.0)  # Higher MW = higher uncertainty
        base_uncertainty = 0.3
        uncertainty = base_uncertainty + (0.4 * complexity_factor)
        
        # Calculate virtue weight
        virtue_weight = self._calculate_virtue_weight(content, campaign)
        
        # Estimate computational cost
        estimated_cost = self._estimate_property_cost(property_name, complexity_factor)
        
        # Generate collapse rules
        collapse_rules = self._generate_property_collapse_rules(property_name, campaign)
        
        claim = Claim(
            id=str(uuid.uuid4()),
            signal_source=signal.get('source_id', 'unknown'),
            claim_type='molecular_property',
            hypothesis=f"Molecule {molecule_id} ({smiles}) has predictable {property_name}",
            testable_prediction=testable_prediction,
            uncertainty=min(uncertainty, 1.0),
            virtue_weight=virtue_weight,
            estimated_cost=estimated_cost,
            priority=self._calculate_priority(virtue_weight, uncertainty, estimated_cost),
            collapse_rules=collapse_rules,
            metadata={
                'molecule_id': molecule_id,
                'smiles': smiles,
                'property': property_name,
                'unit': property_unit,
                'signal_type': signal.get('signal_type'),
                'campaign': campaign.get('name')
            },
            created=datetime.now()
        )
        
        return claim
    
    def _generate_reaction_claim(self, signal: Dict[str, Any], campaign: Dict[str, Any]) -> Optional[Claim]:
        """Generate a reaction performance claim."""
        content = signal.get('content', {})
        
        reaction_id = content.get('id', 'unknown')
        reaction_smiles = content.get('reaction_smiles', '')
        
        if not reaction_smiles:
            logger.warning("No reaction SMILES found, cannot generate reaction claim")
            return None
        
        # Generate testable prediction based on campaign objective
        objective = campaign.get('objective', '')
        if 'yield' in objective:
            testable_prediction = f"Reaction {reaction_id} achieves >80% yield under optimized conditions"
            property_name = 'yield'
        elif 'selectivity' in objective:
            testable_prediction = f"Reaction {reaction_id} shows >90% selectivity for target product"
            property_name = 'selectivity'
        else:
            testable_prediction = f"Reaction {reaction_id} proceeds successfully under standard conditions"
            property_name = 'feasibility'
        
        # Calculate uncertainty based on reaction complexity
        reactant_count = reaction_smiles.count('.') + 1
        complexity_factor = min(reactant_count / 5, 1.0)
        uncertainty = 0.25 + (0.3 * complexity_factor)
        
        virtue_weight = self._calculate_virtue_weight(content, campaign)
        estimated_cost = self._estimate_reaction_cost(complexity_factor)
        collapse_rules = self._generate_reaction_collapse_rules(property_name, campaign)
        
        claim = Claim(
            id=str(uuid.uuid4()),
            signal_source=signal.get('source_id', 'unknown'),
            claim_type='reaction_performance',
            hypothesis=f"Reaction {reaction_id} shows predictable {property_name}",
            testable_prediction=testable_prediction,
            uncertainty=min(uncertainty, 1.0),
            virtue_weight=virtue_weight,
            estimated_cost=estimated_cost,
            priority=self._calculate_priority(virtue_weight, uncertainty, estimated_cost),
            collapse_rules=collapse_rules,
            metadata={
                'reaction_id': reaction_id,
                'reaction_smiles': reaction_smiles,
                'property': property_name,
                'signal_type': signal.get('signal_type'),
                'campaign': campaign.get('name')
            },
            created=datetime.now()
        )
        
        return claim
    
    def _generate_material_claim(self, signal: Dict[str, Any], campaign: Dict[str, Any]) -> Optional[Claim]:
        """Generate a material discovery claim."""
        content = signal.get('content', {})
        
        material_id = content.get('id', 'unknown')
        composition = content.get('composition', 'unknown')
        
        # Generate claim based on campaign objective
        objective = campaign.get('objective', '')
        if 'pfas' in objective.lower():
            testable_prediction = f"Material {material_id} removes PFAS to <10 ng/L"
            property_name = 'pfas_removal'
        elif 'co2' in objective.lower():
            testable_prediction = f"Material {material_id} achieves >85% CO₂ reduction efficiency"
            property_name = 'co2_reduction'
        elif 'catalyst' in objective:
            testable_prediction = f"Material {material_id} shows >50% improvement in catalytic activity"
            property_name = 'catalytic_activity'
        else:
            testable_prediction = f"Material {material_id} meets performance targets"
            property_name = 'performance'
        
        # Calculate uncertainty based on material novelty
        novelty_score = content.get('novelty_score', 0.5)
        uncertainty = 0.3 + (0.4 * novelty_score)
        
        virtue_weight = self._calculate_virtue_weight(content, campaign)
        estimated_cost = self._estimate_material_cost(composition)
        collapse_rules = self._generate_material_collapse_rules(property_name, campaign)
        
        claim = Claim(
            id=str(uuid.uuid4()),
            signal_source=signal.get('source_id', 'unknown'),
            claim_type='material_discovery',
            hypothesis=f"Material {material_id} ({composition}) shows superior {property_name}",
            testable_prediction=testable_prediction,
            uncertainty=min(uncertainty, 1.0),
            virtue_weight=virtue_weight,
            estimated_cost=estimated_cost,
            priority=self._calculate_priority(virtue_weight, uncertainty, estimated_cost),
            collapse_rules=collapse_rules,
            metadata={
                'material_id': material_id,
                'composition': composition,
                'property': property_name,
                'signal_type': signal.get('signal_type'),
                'campaign': campaign.get('name')
            },
            created=datetime.now()
        )
        
        return claim
    
    def _generate_literature_claim(self, signal: Dict[str, Any], campaign: Dict[str, Any]) -> Optional[Claim]:
        """Generate a literature validation claim."""
        content = signal.get('content', {})
        
        paper_id = content.get('id', 'unknown')
        title = content.get('title', 'Unknown paper')
        
        # Extract key claim from paper (simplified)
        if 'pka' in title.lower():
            testable_prediction = "Literature pKa values can be reproduced within 0.3 units"
            property_name = 'pka_reproducibility'
        elif 'logp' in title.lower():
            testable_prediction = "Literature logP values can be reproduced within 0.5 units"
            property_name = 'logp_reproducibility'
        elif 'yield' in title.lower():
            testable_prediction = "Literature reaction yields can be reproduced within 10%"
            property_name = 'yield_reproducibility'
        else:
            testable_prediction = "Literature claims can be independently validated"
            property_name = 'general_reproducibility'
        
        # Literature validation typically has high uncertainty initially
        uncertainty = 0.6
        virtue_weight = 0.8  # High virtue for reproducibility efforts
        estimated_cost = 50  # Moderate cost for replication
        
        collapse_rules = {
            'success_criteria': {
                'replication_success': {'>=': 0.8},
                'statistical_agreement': {'>=': 0.7}
            },
            'uncertainty_max': 0.3,
            'min_replications': 3
        }
        
        claim = Claim(
            id=str(uuid.uuid4()),
            signal_source=signal.get('source_id', 'unknown'),
            claim_type='literature_validation',
            hypothesis=f"Claims in paper {paper_id} are reproducible",
            testable_prediction=testable_prediction,
            uncertainty=uncertainty,
            virtue_weight=virtue_weight,
            estimated_cost=estimated_cost,
            priority='high',  # Reproducibility is always high priority
            collapse_rules=collapse_rules,
            metadata={
                'paper_id': paper_id,
                'title': title,
                'property': property_name,
                'signal_type': signal.get('signal_type'),
                'campaign': campaign.get('name')
            },
            created=datetime.now()
        )
        
        return claim
    
    def _calculate_virtue_weight(self, content: Dict[str, Any], campaign: Dict[str, Any]) -> float:
        """Calculate FoT virtue weight for a claim."""
        virtue_weights = campaign.get('virtue_weighting', {
            'Beneficence': 0.4,
            'Prudence': 0.3,
            'Honesty': 0.2,
            'Temperance': 0.1
        })
        
        # Estimate virtue components based on content
        beneficence = 0.5  # Default neutral
        prudence = 0.5
        honesty = 0.5
        temperance = 0.5
        
        # Adjust based on safety indicators
        if content.get('safety_score', 0) > 0.7:
            prudence += 0.3
        if content.get('toxicity_score', 0) > 0.7:
            prudence -= 0.3
        
        # Adjust based on environmental impact
        if content.get('environmental_score', 0) > 0.7:
            beneficence += 0.2
            temperance += 0.2
        
        # Adjust based on data quality
        if content.get('data_quality', 0) > 0.8:
            honesty += 0.3
        
        # Normalize to 0-1 range
        beneficence = max(0, min(1, beneficence))
        prudence = max(0, min(1, prudence))
        honesty = max(0, min(1, honesty))
        temperance = max(0, min(1, temperance))
        
        # Calculate weighted virtue score
        virtue_score = (
            virtue_weights.get('Beneficence', 0.4) * beneficence +
            virtue_weights.get('Prudence', 0.3) * prudence +
            virtue_weights.get('Honesty', 0.2) * honesty +
            virtue_weights.get('Temperance', 0.1) * temperance
        )
        
        return virtue_score
    
    def _calculate_priority(self, virtue_weight: float, uncertainty: float, cost: float) -> str:
        """Calculate claim priority based on virtue, uncertainty, and cost."""
        # High virtue, low uncertainty, low cost = high priority
        priority_score = virtue_weight * (1 - uncertainty) / max(cost / 100, 0.1)
        
        if priority_score > 2.0:
            return 'high'
        elif priority_score > 1.0:
            return 'medium'
        else:
            return 'low'
    
    def _estimate_property_cost(self, property_name: str, complexity: float) -> float:
        """Estimate computational cost for property prediction."""
        base_costs = {
            'solubility': 10,
            'permeability': 15,
            'pKa': 20,
            'activity': 25
        }
        
        base_cost = base_costs.get(property_name, 15)
        complexity_multiplier = 1 + complexity
        
        return base_cost * complexity_multiplier
    
    def _estimate_reaction_cost(self, complexity: float) -> float:
        """Estimate cost for reaction prediction."""
        base_cost = 30
        complexity_multiplier = 1 + (2 * complexity)  # Reactions more sensitive to complexity
        
        return base_cost * complexity_multiplier
    
    def _estimate_material_cost(self, composition: str) -> float:
        """Estimate cost for material characterization."""
        # Simple heuristic based on composition complexity
        if 'MOF' in composition or 'framework' in composition:
            return 100  # Complex materials
        elif 'carbon' in composition or 'oxide' in composition:
            return 50   # Moderate complexity
        else:
            return 75   # Default
    
    def _generate_property_collapse_rules(self, property_name: str, campaign: Dict[str, Any]) -> Dict[str, Any]:
        """Generate collapse rules for property claims."""
        base_rules = {
            'uncertainty_max': 0.2,
            'min_replications': 2,
            'statistical_significance': 0.05
        }
        
        # Property-specific success criteria
        if property_name == 'solubility':
            success_criteria = {'prediction_error': {'<=': 0.5}}  # log units
        elif property_name == 'pKa':
            success_criteria = {'prediction_error': {'<=': 0.3}}  # pH units
        elif property_name == 'permeability':
            success_criteria = {'prediction_error': {'<=': 0.5}}  # log units
        else:
            success_criteria = {'accuracy': {'>=': 0.7}}
        
        rules = {**base_rules, 'success_criteria': success_criteria}
        
        # Override with campaign-specific rules
        campaign_rules = campaign.get('collapse_rules', {})
        rules.update(campaign_rules)
        
        return rules
    
    def _generate_reaction_collapse_rules(self, property_name: str, campaign: Dict[str, Any]) -> Dict[str, Any]:
        """Generate collapse rules for reaction claims."""
        if property_name == 'yield':
            success_criteria = {'yield_prediction_error': {'<=': 10}}  # %
        elif property_name == 'selectivity':
            success_criteria = {'selectivity_prediction_error': {'<=': 15}}  # %
        else:
            success_criteria = {'success_rate': {'>=': 0.8}}
        
        rules = {
            'success_criteria': success_criteria,
            'uncertainty_max': 0.25,
            'min_replications': 2
        }
        
        # Override with campaign rules
        campaign_rules = campaign.get('collapse_rules', {})
        rules.update(campaign_rules)
        
        return rules
    
    def _generate_material_collapse_rules(self, property_name: str, campaign: Dict[str, Any]) -> Dict[str, Any]:
        """Generate collapse rules for material claims."""
        # Use campaign-specific rules for materials
        campaign_rules = campaign.get('collapse_rules', {})
        
        if not campaign_rules:
            # Default rules
            rules = {
                'success_criteria': {'performance_target': {'>=': 0.8}},
                'uncertainty_max': 0.3,
                'min_replications': 2
            }
        else:
            rules = campaign_rules.copy()
        
        return rules


# Example usage and testing
if __name__ == "__main__":
    # Test claim factory
    factory = ClaimFactory()
    
    # Example signal
    test_signal = {
        'source_id': 'test_source_001',
        'signal_type': 'new_molecule',
        'content': {
            'id': 'mol_001',
            'smiles': 'CCO',
            'molecular_weight': 46.07,
            'safety_score': 0.8,
            'data_quality': 0.9
        },
        'relevance_score': 0.8
    }
    
    # Example campaign
    test_campaign = {
        'name': 'Test Campaign',
        'objective': 'predict_solubility',
        'virtue_weighting': {
            'Beneficence': 0.5,
            'Prudence': 0.3,
            'Honesty': 0.2
        },
        'collapse_rules': {
            'success_criteria': {'accuracy': {'>=': 0.8}},
            'uncertainty_max': 0.2
        }
    }
    
    claim = factory.make_claim(test_signal, test_campaign)
    if claim:
        print(f"✅ Generated test claim: {claim.hypothesis}")
        print(f"🎯 Prediction: {claim.testable_prediction}")
        print(f"📊 Uncertainty: {claim.uncertainty:.2f}")
        print(f"⚖️ Virtue weight: {claim.virtue_weight:.2f}")
        print(f"💰 Estimated cost: {claim.estimated_cost:.1f}")
        print(f"🏷️ Priority: {claim.priority}")
    else:
        print("❌ Failed to generate claim")
