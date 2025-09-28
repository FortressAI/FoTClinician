#!/usr/bin/env python3
"""
Test the actual Streamlit rendering functions
"""

import json
import sys
sys.path.append('.')

# Mock streamlit functions for testing
class MockStreamlit:
    def warning(self, msg):
        print(f"WARNING: {msg}")
    
    def info(self, msg):
        print(f"INFO: {msg}")
    
    def error(self, msg):
        print(f"ERROR: {msg}")

# Import the actual render functions from streamlit_app
try:
    # Set up mock streamlit before importing
    sys.modules['streamlit'] = MockStreamlit()
    
    # Now import our render functions
    from streamlit_app import render_molecule_2d, render_molecule_3d
    print("✅ Successfully imported render functions")
except Exception as e:
    print(f"❌ Failed to import render functions: {e}")
    exit(1)

def test_actual_render_functions():
    print("\n🧪 Testing actual Streamlit render functions...")
    
    # Load the data
    with open('cloud_data_snapshot_with_viz.json', 'r') as f:
        data = json.load(f)
    
    # Get first molecule
    molecule = data['discoveries'][0]
    print(f"Testing molecule: {molecule['smiles']}")
    
    # Test 2D rendering
    print("\n📊 Testing render_molecule_2d...")
    try:
        svg_result = render_molecule_2d(molecule)
        if svg_result:
            print(f"✅ 2D render returned SVG, length: {len(svg_result)}")
            print(f"   Contains bonds: {'bond-' in svg_result}")
            print(f"   Properly closed: {'</svg>' in svg_result}")
            
            # Show structure
            print("First 150 chars:")
            print(svg_result[:150])
        else:
            print("❌ 2D render returned None")
    except Exception as e:
        print(f"❌ 2D render failed: {e}")
    
    # Test 3D rendering
    print("\n📊 Testing render_molecule_3d...")
    try:
        result_3d = render_molecule_3d(molecule)
        print(f"3D render result: {result_3d}")
    except Exception as e:
        print(f"❌ 3D render failed: {e}")

if __name__ == "__main__":
    test_actual_render_functions()
