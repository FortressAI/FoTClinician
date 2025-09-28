#!/usr/bin/env python3
"""
Test script to reproduce the visualization issues
"""

import json
import sys
sys.path.append('.')

# Test loading data and rendering like Streamlit does
def test_2d_rendering():
    print("🧪 Testing 2D molecular visualization...")
    
    # Load the data
    with open('cloud_data_snapshot_with_viz.json', 'r') as f:
        data = json.load(f)
    
    # Get first molecule (CCCO - ethanol)
    molecule = data['discoveries'][0]
    print(f"Testing molecule: {molecule['smiles']}")
    
    # Test the 2D rendering logic
    if 'visualization_2d' in molecule:
        viz_data = molecule['visualization_2d']
        if viz_data.get('svg'):
            svg_content = viz_data['svg']
            print(f"SVG content length: {len(svg_content)}")
            
            # Check what our validation logic does
            if 'bond-' in svg_content or 'rdkit' in svg_content.lower():
                print("✅ SVG passes validation checks")
                
                # Show first and last 200 chars to see structure
                print("\nFirst 200 chars of SVG:")
                print(svg_content[:200])
                print("\nLast 200 chars of SVG:")
                print(svg_content[-200:])
                
                # Check for common HTML rendering issues
                if '<rect' in svg_content and '<path' in svg_content:
                    print("✅ SVG has rect and path elements (bonds)")
                else:
                    print("❌ SVG missing essential elements")
                    
                return svg_content
            else:
                print("❌ SVG fails validation")
        else:
            print("❌ No SVG data found")
    else:
        print("❌ No visualization_2d found")
    
    return None

def test_3d_rendering():
    print("\n🧬 Testing 3D molecular visualization...")
    
    # Load the data
    with open('cloud_data_snapshot_with_viz.json', 'r') as f:
        data = json.load(f)
    
    # Get first molecule
    molecule = data['discoveries'][0]
    print(f"Testing molecule: {molecule['smiles']}")
    
    if 'visualization_3d' in molecule:
        viz_data = molecule['visualization_3d']
        if viz_data.get('molblock'):
            molblock = viz_data['molblock']
            print(f"MOL block length: {len(molblock)}")
            print("First 5 lines of MOL block:")
            print('\n'.join(molblock.split('\n')[:5]))
            return molblock
        else:
            print("❌ No molblock data found")
    else:
        print("❌ No visualization_3d found")
    
    return None

if __name__ == "__main__":
    svg_result = test_2d_rendering()
    molblock_result = test_3d_rendering()
    
    print(f"\n📊 Results:")
    print(f"2D SVG available: {svg_result is not None}")
    print(f"3D MOL block available: {molblock_result is not None}")
