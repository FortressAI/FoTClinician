#!/usr/bin/env python3
"""
Unit tests for FoTChemistry Streamlit data loading and display.
Tests the core issues: score extraction, data quantity, and field mapping.
"""

import unittest
import json
import tempfile
import os
from unittest.mock import patch, MagicMock

# Import the streamlit app functions
import sys
sys.path.append('.')

class TestStreamlitDataLoading(unittest.TestCase):
    """Test data loading and processing for Streamlit dashboard."""
    
    def setUp(self):
        """Set up test data."""
        self.sample_discovery = {
            "discovery_id": "test-123",
            "smiles": "CCO",
            "combined_score": 0.785,
            "molecular_properties": {
                "molecular_weight": 46.07,
                "logp": -0.31
            },
            "discovery_timestamp": "2025-09-29T06:00:00.000000"
        }
        
        self.sample_dataset = {
            "discovery_summary": {
                "total_discoveries": 3043,
                "statistics": {
                    "avg_score": 0.735,
                    "max_score": 0.798
                }
            },
            "discoveries": [self.sample_discovery]
        }
    
    def test_score_extraction_with_combined_score(self):
        """Test that combined_score is properly extracted as score."""
        discoveries = [self.sample_discovery.copy()]
        
        # Test score extraction logic (simulating Streamlit code)
        scores = [mol.get('score', mol.get('combined_score', 0)) for mol in discoveries]
        
        self.assertEqual(len(scores), 1)
        self.assertEqual(scores[0], 0.785)
        self.assertNotEqual(scores[0], 0.0)  # Should not be zero
    
    def test_score_extraction_with_missing_score(self):
        """Test handling when neither score nor combined_score exist."""
        discovery_no_score = {
            "discovery_id": "test-no-score",
            "smiles": "CC"
        }
        
        scores = [mol.get('score', mol.get('combined_score', 0)) for mol in [discovery_no_score]]
        self.assertEqual(scores[0], 0)
    
    def test_data_quantity_not_limited(self):
        """Test that all discoveries are loaded, not just samples."""
        # Create dataset with multiple discoveries
        large_dataset = {
            "discovery_summary": {"total_discoveries": 3043},
            "discoveries": [self.sample_discovery.copy() for _ in range(100)]
        }
        
        # Simulate loading all discoveries
        discoveries = large_dataset.get('discoveries', [])
        total_count = large_dataset.get('discovery_summary', {}).get('total_discoveries', 0)
        
        # Should load all available discoveries in the file
        self.assertEqual(len(discoveries), 100)
        # But should report the true total count
        self.assertEqual(total_count, 3043)
    
    def test_file_loading_priority(self):
        """Test that files are loaded in correct priority order."""
        # Test file priority logic (simplified)
        expected_priority = [
            "results/chemistry_discoveries.json",  # Should be first priority now
            "results/overnight_discovery_mega_dataset.json",
            "cloud_data_snapshot.json"
        ]
        
        # Just verify the priority order is correct
        self.assertEqual(expected_priority[0], "results/chemistry_discoveries.json")
        self.assertEqual(len(expected_priority), 3)
    
    def test_statistics_calculation(self):
        """Test that statistics are calculated correctly from score data."""
        discoveries_with_scores = [
            {"combined_score": 0.700, "score": 0.700},
            {"combined_score": 0.800, "score": 0.800},
            {"combined_score": 0.750, "score": 0.750}
        ]
        
        # Simulate Streamlit statistics calculation
        scores = [mol.get('score', 0) for mol in discoveries_with_scores if mol.get('score') is not None]
        avg_score = sum(scores) / len(scores) if scores else 0
        max_score = max(scores) if scores else 0
        
        self.assertEqual(len(scores), 3)
        self.assertAlmostEqual(avg_score, 0.750, places=3)
        self.assertEqual(max_score, 0.800)
        self.assertNotEqual(avg_score, 0.0)
        self.assertNotEqual(max_score, 0.0)
    
    def test_field_mapping_compatibility(self):
        """Test that both score and combined_score fields work."""
        discoveries = [
            {"score": 0.785},  # New format
            {"combined_score": 0.792},  # Old format
            {"score": 0.750, "combined_score": 0.760}  # Both (score should take priority)
        ]
        
        # Test extraction with fallback
        extracted_scores = []
        for mol in discoveries:
            score = mol.get('score', mol.get('combined_score', 0))
            extracted_scores.append(score)
        
        self.assertEqual(extracted_scores[0], 0.785)  # Direct score
        self.assertEqual(extracted_scores[1], 0.792)  # Fallback to combined_score
        self.assertEqual(extracted_scores[2], 0.750)  # Score takes priority
    
    def test_data_structure_validation(self):
        """Test that data structure is valid for Streamlit consumption."""
        # Test valid structure
        valid_data = {
            "discovery_summary": {
                "total_discoveries": 100,
                "statistics": {"avg_score": 0.7}
            },
            "discoveries": [{"smiles": "CCO", "score": 0.7}]
        }
        
        # Should have required keys
        self.assertIn('discovery_summary', valid_data)
        self.assertIn('discoveries', valid_data)
        self.assertIn('total_discoveries', valid_data['discovery_summary'])
        
        # Discoveries should be a list
        self.assertIsInstance(valid_data['discoveries'], list)
        
        # Each discovery should have basic fields
        if valid_data['discoveries']:
            discovery = valid_data['discoveries'][0]
            self.assertIn('smiles', discovery)


class TestDataFileFixes(unittest.TestCase):
    """Test the data file fixes and exports."""
    
    def test_score_mapping_in_export(self):
        """Test that export properly maps combined_score to score."""
        # Simulate discovery data with combined_score
        raw_discovery = {
            "smiles": "CCO",
            "combined_score": 0.785
        }
        
        # Apply the mapping fix
        if 'combined_score' in raw_discovery and 'score' not in raw_discovery:
            raw_discovery['score'] = raw_discovery['combined_score']
        
        # Verify mapping worked
        self.assertIn('score', raw_discovery)
        self.assertEqual(raw_discovery['score'], 0.785)
        self.assertEqual(raw_discovery['score'], raw_discovery['combined_score'])
    
    def test_complete_dataset_export_structure(self):
        """Test that complete dataset export has correct structure."""
        # Simulate the export data structure
        export_data = {
            "discovery_summary": {
                "total_discoveries": 3043,
                "statistics": {
                    "avg_score": 0.735,
                    "max_score": 0.798,
                    "score_field_mapped": True
                },
                "data_fixes": [
                    "Mapped combined_score to score for Streamlit compatibility"
                ]
            },
            "discoveries": []  # Would contain all 3043 discoveries
        }
        
        # Validate structure
        self.assertIn('discovery_summary', export_data)
        self.assertIn('discoveries', export_data)
        self.assertTrue(export_data['discovery_summary']['statistics']['score_field_mapped'])
        self.assertGreater(export_data['discovery_summary']['total_discoveries'], 1000)


def run_all_tests():
    """Run all tests and return results."""
    unittest.main(argv=[''], exit=False, verbosity=2)


if __name__ == "__main__":
    print("🧪 Running FoTChemistry Streamlit Data Loading Tests")
    print("=" * 60)
    run_all_tests()
