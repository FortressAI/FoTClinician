#!/usr/bin/env python3
"""
Mass Balance Validation Agent

Validates that chemical reactions satisfy mass balance conservation laws
according to Field of Truth principles.
"""

from typing import Dict, List, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class Molecule:
    """Simple molecule representation."""
    smiles: str
    molecular_formula: str
    molecular_weight: float


@dataclass
class Reaction:
    """Chemical reaction representation."""
    reactants: List[Molecule]
    products: List[Molecule]
    stoichiometry: Dict[str, float]  # SMILES -> stoichiometric coefficient


@dataclass
class ValidationResult:
    """Result of mass balance validation."""
    is_valid: bool
    mass_difference: float
    confidence: float
    errors: List[str]


class MassBalanceValidator:
    """Agent for validating chemical reaction mass balance."""
    
    def __init__(self, tolerance: float = 0.01):
        """
        Initialize mass balance validator.
        
        Args:
            tolerance: Acceptable mass difference (fraction, e.g., 0.01 = 1%)
        """
        self.tolerance = tolerance
    
    def validate_reaction(self, reaction: Reaction) -> ValidationResult:
        """
        Validate that reaction satisfies mass balance.
        
        Args:
            reaction: Chemical reaction to validate
            
        Returns:
            Validation result with pass/fail and diagnostics
        """
        logger.info("Validating mass balance for reaction")
        
        try:
            # Calculate total mass of reactants
            reactant_mass = self._calculate_total_mass(
                reaction.reactants, 
                reaction.stoichiometry,
                is_reactant=True
            )
            
            # Calculate total mass of products  
            product_mass = self._calculate_total_mass(
                reaction.products,
                reaction.stoichiometry, 
                is_reactant=False
            )
            
            # Calculate mass difference
            mass_diff = abs(reactant_mass - product_mass)
            relative_diff = mass_diff / max(reactant_mass, product_mass)
            
            # Determine if within tolerance
            is_valid = relative_diff <= self.tolerance
            
            # Calculate confidence based on precision
            confidence = max(0.0, 1.0 - (relative_diff / self.tolerance))
            
            errors = []
            if not is_valid:
                errors.append(f"Mass imbalance: {relative_diff:.3%} > tolerance {self.tolerance:.1%}")
            
            result = ValidationResult(
                is_valid=is_valid,
                mass_difference=relative_diff,
                confidence=confidence,
                errors=errors
            )
            
            logger.info(f"Mass balance validation: {'PASS' if is_valid else 'FAIL'}")
            return result
            
        except Exception as e:
            logger.error(f"Mass balance validation failed: {e}")
            return ValidationResult(
                is_valid=False,
                mass_difference=float('inf'),
                confidence=0.0,
                errors=[str(e)]
            )
    
    def _calculate_total_mass(self, molecules: List[Molecule], 
                             stoichiometry: Dict[str, float], 
                             is_reactant: bool) -> float:
        """
        Calculate total mass for reactants or products.
        
        Args:
            molecules: List of molecules
            stoichiometry: Stoichiometric coefficients
            is_reactant: True for reactants, False for products
            
        Returns:
            Total mass in atomic mass units
        """
        total_mass = 0.0
        
        for molecule in molecules:
            # Get stoichiometric coefficient
            coeff = stoichiometry.get(molecule.smiles, 1.0)
            
            # Reactants are consumed (negative), products formed (positive)
            if is_reactant:
                coeff = abs(coeff)  # Take absolute value for mass calculation
            
            mass_contribution = molecule.molecular_weight * coeff
            total_mass += mass_contribution
            
            logger.debug(f"Molecule {molecule.smiles}: MW={molecule.molecular_weight:.3f}, "
                        f"coeff={coeff}, contribution={mass_contribution:.3f}")
        
        return total_mass
    
    def validate_elemental_balance(self, reaction: Reaction) -> ValidationResult:
        """
        Validate elemental balance (C, H, N, O, etc.) for reaction.
        
        Args:
            reaction: Chemical reaction to validate
            
        Returns:
            Validation result for elemental balance
        """
        # Placeholder for elemental balance validation
        # Would require molecular formula parsing
        logger.warning("Elemental balance validation not yet implemented")
        
        return ValidationResult(
            is_valid=True,
            mass_difference=0.0,
            confidence=0.5,  # Lower confidence since not fully implemented
            errors=["Elemental balance validation is placeholder"]
        )


if __name__ == "__main__":
    # Example usage
    validator = MassBalanceValidator(tolerance=0.01)
    
    # Create example reaction: CH4 + 2O2 -> CO2 + 2H2O
    methane = Molecule("C", "CH4", 16.043)
    oxygen = Molecule("O=O", "O2", 31.998)
    co2 = Molecule("O=C=O", "CO2", 44.009)
    water = Molecule("O", "H2O", 18.015)
    
    reaction = Reaction(
        reactants=[methane, oxygen],
        products=[co2, water],
        stoichiometry={
            "C": 1.0,      # 1 CH4
            "O=O": 2.0,    # 2 O2
            "O=C=O": 1.0,  # 1 CO2
            "O": 2.0       # 2 H2O
        }
    )
    
    result = validator.validate_reaction(reaction)
    print(f"Validation result: {result}")
