#!/usr/bin/env python3
"""
Ethics Agent - Safety Guard Module

Handles safety screening and ethical validation.
"""

from typing import Dict, List, Any, Optional

class EthicsGate:
    """Ethics and safety validation gate"""
    
    def __init__(self):
        # Default safety lists
        self.default_blocklist = [
            'explosives', 'nerve_agents', 'carcinogens', 
            'hazard_class_1', 'toxic_heavy_metals'
        ]
        
        self.default_allowlist = [
            'activated_carbons', 'silicas', 'aluminas', 
            'safe_ionic_liquids', 'food_grade_materials'
        ]
    
    def screen_candidate(self, candidate: Dict[str, Any], ethics_config: Dict[str, Any]) -> Dict[str, Any]:
        """Screen candidate for ethical/safety concerns"""
        
        blocklist = ethics_config.get('blocklists', self.default_blocklist)
        allowlist = ethics_config.get('allowlists', self.default_allowlist)
        
        # Simulate safety screening
        candidate_type = candidate.get('type', 'unknown')
        candidate_desc = candidate.get('description', '').lower()
        
        # Check blocklist
        blocked = False
        block_reason = None
        
        for blocked_item in blocklist:
            if blocked_item.lower() in candidate_desc:
                blocked = True
                block_reason = f"Contains blocked material: {blocked_item}"
                break
        
        # Check allowlist (if present, must match)
        if allowlist and not blocked:
            allowed = False
            for allowed_item in allowlist:
                if allowed_item.lower() in candidate_desc:
                    allowed = True
                    break
            
            if not allowed:
                blocked = True
                block_reason = "Not in allowlist"
        
        # Ethics score (higher = more ethical)
        ethics_score = 0.9 if not blocked else 0.1
        
        return {
            'candidate': candidate,
            'passed': not blocked,
            'ethics_score': ethics_score,
            'block_reason': block_reason,
            'screening_agent': 'ethics_guard'
        }
    
    def batch_screen(self, candidates: List[Dict[str, Any]], ethics_config: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Screen multiple candidates"""
        results = []
        
        for candidate in candidates:
            result = self.screen_candidate(candidate, ethics_config)
            results.append(result)
            
        return results
