"""
Simple molecular visualization without RDKit dependencies
Creates basic 2D molecular representations from SMILES
"""

import re
import math
import random
from typing import List, Tuple, Dict

def smiles_to_simple_svg(smiles: str, width: int = 400, height: int = 400) -> str:
    """
    Convert SMILES to a simple SVG representation
    This is a basic implementation for when RDKit is not available
    """
    if not smiles:
        return create_error_svg("No SMILES provided", width, height)
    
    try:
        # Parse basic SMILES components
        atoms = extract_atoms(smiles)
        bonds = extract_bonds(smiles)
        
        if not atoms:
            return create_error_svg("Could not parse SMILES", width, height)
        
        # Generate simple 2D coordinates
        coordinates = generate_2d_coordinates(atoms, bonds, width, height)
        
        # Create SVG
        svg = create_molecule_svg(atoms, bonds, coordinates, width, height, smiles)
        return svg
        
    except Exception as e:
        return create_error_svg(f"Error: {str(e)}", width, height)

def extract_atoms(smiles: str) -> List[str]:
    """Extract atoms from SMILES string"""
    # Remove bonds and brackets for simple parsing
    cleaned = re.sub(r'[=\-#\(\)\[\]0-9@+\-]', '', smiles)
    
    atoms = []
    i = 0
    while i < len(cleaned):
        if cleaned[i].isupper():
            atom = cleaned[i]
            # Check for two-letter atoms
            if i + 1 < len(cleaned) and cleaned[i + 1].islower():
                atom += cleaned[i + 1]
                i += 1
            atoms.append(atom)
        i += 1
    
    # If no atoms found, try a different approach
    if not atoms:
        # Count carbons implicitly
        carbon_count = len([c for c in smiles if c.isalnum() and c.upper() not in 'NOSPFCL'])
        atoms = ['C'] * max(carbon_count, 1)
    
    return atoms

def extract_bonds(smiles: str) -> List[Tuple[int, int, str]]:
    """Extract bond information (simplified)"""
    bonds = []
    # This is a very simplified bond extraction
    # In reality, SMILES parsing is much more complex
    for i in range(len(smiles) - 1):
        if smiles[i].isalnum() and smiles[i + 1].isalnum():
            bonds.append((i, i + 1, 'single'))
        elif smiles[i] == '=' and i > 0 and i < len(smiles) - 1:
            bonds.append((i - 1, i + 1, 'double'))
    
    return bonds

def generate_2d_coordinates(atoms: List[str], bonds: List[Tuple[int, int, str]], width: int, height: int) -> Dict[int, Tuple[float, float]]:
    """Generate simple 2D coordinates for atoms"""
    coordinates = {}
    
    if len(atoms) == 1:
        coordinates[0] = (width // 2, height // 2)
        return coordinates
    
    # Create a simple circular or linear layout
    center_x, center_y = width // 2, height // 2
    radius = min(width, height) // 3
    
    if len(atoms) <= 6:
        # Circular layout for small molecules
        for i, atom in enumerate(atoms):
            angle = 2 * math.pi * i / len(atoms)
            x = center_x + radius * math.cos(angle)
            y = center_y + radius * math.sin(angle)
            coordinates[i] = (x, y)
    else:
        # Linear layout for larger molecules
        spacing = min(width, height) / (len(atoms) + 1)
        for i, atom in enumerate(atoms):
            x = spacing * (i + 1)
            y = center_y + (random.random() - 0.5) * height * 0.3
            coordinates[i] = (x, y)
    
    return coordinates

def create_molecule_svg(atoms: List[str], bonds: List[Tuple[int, int, str]], coordinates: Dict[int, Tuple[float, float]], width: int, height: int, smiles: str) -> str:
    """Create SVG representation of the molecule"""
    
    svg_lines = [
        f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">',
        f'<rect width="{width}" height="{height}" fill="white" stroke="#ddd" stroke-width="1"/>',
        f'<title>Molecular Structure: {smiles}</title>'
    ]
    
    # Add bonds
    for bond in bonds:
        if len(bond) >= 3 and bond[0] in coordinates and bond[1] in coordinates:
            x1, y1 = coordinates[bond[0]]
            x2, y2 = coordinates[bond[1]]
            
            stroke_width = 2 if bond[2] == 'single' else 3
            svg_lines.append(
                f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                f'stroke="#333" stroke-width="{stroke_width}"/>'
            )
    
    # Add atoms
    for i, atom in enumerate(atoms):
        if i in coordinates:
            x, y = coordinates[i]
            
            # Atom colors
            color = get_atom_color(atom)
            
            # Don't show carbon atoms unless they're special
            if atom.upper() != 'C' or len(atoms) < 3:
                # Circle for atom
                svg_lines.append(
                    f'<circle cx="{x:.1f}" cy="{y:.1f}" r="12" fill="{color}" stroke="#333" stroke-width="1"/>'
                )
                
                # Atom label
                svg_lines.append(
                    f'<text x="{x:.1f}" y="{y + 4:.1f}" text-anchor="middle" '
                    f'font-family="Arial, sans-serif" font-size="12" font-weight="bold" fill="white">{atom}</text>'
                )
    
    # Add SMILES string at bottom
    svg_lines.append(
        f'<text x="{width // 2}" y="{height - 10}" text-anchor="middle" '
        f'font-family="monospace" font-size="14" fill="#666">{smiles}</text>'
    )
    
    svg_lines.append('</svg>')
    
    return '\n'.join(svg_lines)

def get_atom_color(atom: str) -> str:
    """Get color for atom type"""
    colors = {
        'C': '#333333',
        'N': '#3050F8',
        'O': '#FF0D0D',
        'S': '#FFFF30',
        'P': '#FF8000',
        'F': '#90E050',
        'Cl': '#1FF01F',
        'Br': '#A62929',
        'I': '#940094'
    }
    return colors.get(atom.upper(), '#666666')

def create_error_svg(message: str, width: int, height: int) -> str:
    """Create error SVG when parsing fails"""
    return f'''
    <svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
        <rect width="{width}" height="{height}" fill="#f8f9fa" stroke="#ddd" stroke-width="1"/>
        <text x="{width // 2}" y="{height // 2 - 10}" text-anchor="middle" 
              font-family="Arial, sans-serif" font-size="14" fill="#666">
            Unable to render molecular structure
        </text>
        <text x="{width // 2}" y="{height // 2 + 10}" text-anchor="middle" 
              font-family="monospace" font-size="12" fill="#999">
            {message}
        </text>
    </svg>
    '''

def create_3d_fallback_html(smiles: str) -> str:
    """Create 3D visualization fallback"""
    return f'''
    <div style="width: 400px; height: 400px; border: 1px solid #ddd; border-radius: 8px; 
                background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%); 
                display: flex; flex-direction: column; align-items: center; justify-content: center;
                font-family: Arial, sans-serif;">
        <div style="font-size: 48px; margin-bottom: 20px;">🧬</div>
        <div style="font-size: 18px; font-weight: bold; margin-bottom: 10px; color: #333;">
            3D Molecular Structure
        </div>
        <div style="font-size: 14px; color: #666; margin-bottom: 15px; text-align: center; padding: 0 20px;">
            Interactive 3D visualization requires additional chemistry packages
        </div>
        <div style="font-family: monospace; background: white; padding: 8px 12px; 
                    border-radius: 4px; border: 1px solid #ddd; font-size: 12px; color: #333;">
            {smiles}
        </div>
    </div>
    '''
