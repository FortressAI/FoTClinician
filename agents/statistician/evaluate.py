#!/usr/bin/env python3
"""
FoTChemistry Statistician Agent - Truth Collapse Evaluation

The Statistician agent applies Field of Truth collapse rules to determine
when claims should be collapsed to "Truth", "Refuted", or "Needs Evidence".

Uses Pareto virtue optimization and Bayesian evidence evaluation.
"""

import numpy as np
import logging
from typing import Dict, List, Any, Tuple, Optional
from dataclasses import dataclass
from scipy import stats
import json

logger = logging.getLogger(__name__)


@dataclass
class CollapseVerdict:
    """Result of claim evaluation and potential truth collapse."""
    status: str  # "truth", "refuted", "needs_evidence", "insufficient_data"
    confidence: float  # 0.0 to 1.0
    evidence_strength: float
    virtue_score: float
    replication_count: int
    statistical_summary: Dict[str, Any]
    reasoning: str
    recommendation: str


def pareto_wave_collapse(front: List[Tuple[float, ...]]) -> np.ndarray:
    """
    Collapse Pareto front to stable center of mass.
    
    For FoT virtue vectors (Beneficence, Prudence, Honesty, Temperance),
    find the non-dominated solutions and return their centroid.
    
    Args:
        front: List of virtue tuples
        
    Returns:
        Collapsed virtue vector as numpy array
    """
    if not front:
        return np.array([0.5, 0.5, 0.5, 0.5])  # Default neutral virtues
    
    arr = np.array(front)
    
    # Find non-dominated solutions
    non_dominated = []
    for i, a in enumerate(arr):
        # Check if any other solution dominates this one
        dominated = False
        for j, b in enumerate(arr):
            if i != j and np.all(b >= a) and np.any(b > a):
                dominated = True
                break
        
        if not dominated:
            non_dominated.append(tuple(a))
    
    # Return center of mass of non-dominated set
    if non_dominated:
        return np.mean(np.array(non_dominated), axis=0)
    else:
        return np.mean(arr, axis=0)


class CollapseRules:
    """
    FoT truth collapse rules engine.
    
    Determines when scientific claims should collapse to accepted truth
    based on statistical evidence, replication, and virtue alignment.
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize collapse rules with configuration."""
        self.config = config or {}
        self.default_thresholds = {
            'min_confidence': 0.8,
            'min_replications': 2,
            'max_uncertainty': 0.2,
            'statistical_significance': 0.05,
            'effect_size_threshold': 0.3,
            'virtue_threshold': 0.6
        }
        self.thresholds = {**self.default_thresholds, **self.config.get('thresholds', {})}
    
    def judge_claim(self, claim: Dict[str, Any], candidate: Dict[str, Any], 
                   measurement_result: Dict[str, Any]) -> CollapseVerdict:
        """
        Judge a claim and determine if it should collapse to truth.
        
        Args:
            claim: Original claim being tested
            candidate: Specific candidate being evaluated
            measurement_result: Results from measurement agent
            
        Returns:
            CollapseVerdict with recommendation
        """
        logger.info(f"⚖️ Evaluating claim: {claim.get('id', 'unknown')}")
        
        try:
            # Extract collapse rules from claim
            collapse_rules = claim.get('collapse_rules', {})
            success_criteria = collapse_rules.get('success_criteria', {})
            
            # Get measurement metrics
            metrics = measurement_result.get('metrics', {})
            uncertainty = measurement_result.get('uncertainty', 1.0)
            replications = measurement_result.get('replications', [])
            
            # Evaluate success criteria
            criteria_met = self._evaluate_success_criteria(metrics, success_criteria)
            
            # Statistical evaluation
            statistical_summary = self._statistical_evaluation(
                metrics, replications, collapse_rules
            )
            
            # Virtue evaluation
            virtue_history = measurement_result.get('virtue_vector_history', [])
            virtue_score = self._evaluate_virtues(virtue_history, candidate)
            
            # Evidence strength assessment
            evidence_strength = self._calculate_evidence_strength(
                metrics, replications, statistical_summary
            )
            
            # Replication assessment
            replication_count = len(replications)
            replication_agreement = self._assess_replication_agreement(replications)
            
            # Overall confidence calculation
            confidence = self._calculate_overall_confidence(
                criteria_met, statistical_summary, virtue_score, 
                evidence_strength, replication_agreement
            )
            
            # Determine collapse decision
            status, reasoning, recommendation = self._make_collapse_decision(
                criteria_met, confidence, uncertainty, replication_count,
                statistical_summary, virtue_score
            )
            
            verdict = CollapseVerdict(
                status=status,
                confidence=confidence,
                evidence_strength=evidence_strength,
                virtue_score=virtue_score,
                replication_count=replication_count,
                statistical_summary=statistical_summary,
                reasoning=reasoning,
                recommendation=recommendation
            )
            
            logger.info(f"⚖️ Verdict: {status} (confidence: {confidence:.2f})")
            return verdict
            
        except Exception as e:
            logger.error(f"❌ Claim evaluation failed: {e}")
            return CollapseVerdict(
                status="insufficient_data",
                confidence=0.0,
                evidence_strength=0.0,
                virtue_score=0.0,
                replication_count=0,
                statistical_summary={},
                reasoning=f"Evaluation failed: {str(e)}",
                recommendation="Fix evaluation pipeline and retry"
            )
    
    def _evaluate_success_criteria(self, metrics: Dict[str, Any], 
                                 criteria: Dict[str, Any]) -> Dict[str, bool]:
        """Evaluate whether metrics meet success criteria."""
        results = {}
        
        for criterion, condition in criteria.items():
            if criterion not in metrics:
                results[criterion] = False
                continue
            
            value = metrics[criterion]
            
            if isinstance(condition, dict):
                # Handle operators like {">=": 0.85}
                for operator, threshold in condition.items():
                    if operator == ">=":
                        results[criterion] = value >= threshold
                    elif operator == "<=":
                        results[criterion] = value <= threshold
                    elif operator == ">":
                        results[criterion] = value > threshold
                    elif operator == "<":
                        results[criterion] = value < threshold
                    elif operator == "==":
                        results[criterion] = abs(value - threshold) < 1e-6
                    else:
                        logger.warning(f"Unknown operator: {operator}")
                        results[criterion] = False
                    break  # Only process first operator
            else:
                # Direct comparison
                results[criterion] = value >= condition
        
        met_count = sum(results.values())
        total_count = len(results)
        
        logger.info(f"📊 Success criteria: {met_count}/{total_count} met")
        return results
    
    def _statistical_evaluation(self, metrics: Dict[str, Any], 
                              replications: List[Dict[str, Any]],
                              collapse_rules: Dict[str, Any]) -> Dict[str, Any]:
        """Perform statistical evaluation of results."""
        summary = {
            'sample_size': len(replications),
            'mean_values': {},
            'std_values': {},
            'confidence_intervals': {},
            'p_values': {},
            'effect_sizes': {}
        }
        
        if len(replications) < 2:
            summary['sufficient_data'] = False
            return summary
        
        summary['sufficient_data'] = True
        
        # Extract values for each metric across replications
        for metric_name in metrics.keys():
            values = []
            for rep in replications:
                if metric_name in rep.get('metrics', {}):
                    values.append(rep['metrics'][metric_name])
            
            if len(values) >= 2:
                # Basic statistics
                mean_val = np.mean(values)
                std_val = np.std(values, ddof=1)
                
                summary['mean_values'][metric_name] = mean_val
                summary['std_values'][metric_name] = std_val
                
                # Confidence interval
                if len(values) >= 3:
                    ci = stats.t.interval(0.95, len(values)-1, 
                                        loc=mean_val, 
                                        scale=stats.sem(values))
                    summary['confidence_intervals'][metric_name] = ci
                
                # One-sample t-test against threshold
                thresholds = collapse_rules.get('success_criteria', {})
                if metric_name in thresholds:
                    threshold_dict = thresholds[metric_name]
                    if isinstance(threshold_dict, dict):
                        threshold = list(threshold_dict.values())[0]
                        t_stat, p_val = stats.ttest_1samp(values, threshold)
                        summary['p_values'][metric_name] = p_val
                        
                        # Effect size (Cohen's d)
                        effect_size = (mean_val - threshold) / std_val if std_val > 0 else 0
                        summary['effect_sizes'][metric_name] = effect_size
        
        return summary
    
    def _evaluate_virtues(self, virtue_history: List[List[float]], 
                         candidate: Dict[str, Any]) -> float:
        """Evaluate virtue alignment using Pareto collapse."""
        if not virtue_history:
            # Default virtue assessment based on candidate properties
            return self._default_virtue_assessment(candidate)
        
        # Apply Pareto wave collapse to virtue history
        collapsed_virtues = pareto_wave_collapse(virtue_history)
        
        # Overall virtue score as weighted average
        virtue_weights = [0.4, 0.3, 0.2, 0.1]  # Beneficence, Prudence, Honesty, Temperance
        virtue_score = np.dot(collapsed_virtues, virtue_weights)
        
        logger.debug(f"🎯 Virtue assessment: {virtue_score:.3f}")
        return float(virtue_score)
    
    def _default_virtue_assessment(self, candidate: Dict[str, Any]) -> float:
        """Default virtue assessment when no history available."""
        # Simple heuristic based on candidate properties
        score = 0.5  # Neutral baseline
        
        # Boost for safety indicators
        if candidate.get('safety_score', 0) > 0.7:
            score += 0.1
        
        # Boost for environmental friendliness
        if candidate.get('environmental_score', 0) > 0.7:
            score += 0.1
        
        # Penalty for high toxicity
        if candidate.get('toxicity_score', 0) > 0.7:
            score -= 0.2
        
        return max(0.0, min(1.0, score))
    
    def _calculate_evidence_strength(self, metrics: Dict[str, Any],
                                   replications: List[Dict[str, Any]],
                                   statistical_summary: Dict[str, Any]) -> float:
        """Calculate overall evidence strength."""
        strength = 0.0
        
        # Sample size contribution
        n = len(replications)
        if n >= 5:
            strength += 0.3
        elif n >= 3:
            strength += 0.2
        elif n >= 2:
            strength += 0.1
        
        # Statistical significance contribution
        p_values = statistical_summary.get('p_values', {})
        significant_count = sum(1 for p in p_values.values() if p < 0.05)
        if p_values:
            strength += 0.3 * (significant_count / len(p_values))
        
        # Effect size contribution
        effect_sizes = statistical_summary.get('effect_sizes', {})
        large_effects = sum(1 for es in effect_sizes.values() if abs(es) > 0.8)
        if effect_sizes:
            strength += 0.2 * (large_effects / len(effect_sizes))
        
        # Consistency contribution (low standard deviation)
        std_values = statistical_summary.get('std_values', {})
        if std_values:
            mean_cv = np.mean([std/abs(mean) for mean, std in 
                             zip(statistical_summary.get('mean_values', {}).values(),
                                 std_values.values()) if mean != 0])
            if mean_cv < 0.1:  # CV < 10%
                strength += 0.2
            elif mean_cv < 0.2:  # CV < 20%
                strength += 0.1
        
        return min(1.0, strength)
    
    def _assess_replication_agreement(self, replications: List[Dict[str, Any]]) -> float:
        """Assess agreement between independent replications."""
        if len(replications) < 2:
            return 0.0
        
        # Calculate coefficient of variation for key metrics
        agreements = []
        
        for replication in replications:
            metrics = replication.get('metrics', {})
            for metric_name, value in metrics.items():
                # Collect values for this metric across replications
                values = []
                for rep in replications:
                    if metric_name in rep.get('metrics', {}):
                        values.append(rep['metrics'][metric_name])
                
                if len(values) >= 2:
                    cv = np.std(values) / np.mean(values) if np.mean(values) != 0 else 1.0
                    agreement = max(0.0, 1.0 - cv)  # Lower CV = higher agreement
                    agreements.append(agreement)
        
        return np.mean(agreements) if agreements else 0.0
    
    def _calculate_overall_confidence(self, criteria_met: Dict[str, bool],
                                    statistical_summary: Dict[str, Any],
                                    virtue_score: float,
                                    evidence_strength: float,
                                    replication_agreement: float) -> float:
        """Calculate overall confidence in the claim."""
        
        # Success criteria contribution (40%)
        criteria_score = sum(criteria_met.values()) / max(len(criteria_met), 1)
        
        # Statistical evidence contribution (30%) 
        statistical_score = 0.0
        if statistical_summary.get('sufficient_data', False):
            # Factor in p-values and effect sizes
            p_values = statistical_summary.get('p_values', {})
            effect_sizes = statistical_summary.get('effect_sizes', {})
            
            if p_values:
                significant_fraction = sum(1 for p in p_values.values() if p < 0.05) / len(p_values)
                statistical_score += 0.5 * significant_fraction
            
            if effect_sizes:
                large_effect_fraction = sum(1 for es in effect_sizes.values() if abs(es) > 0.5) / len(effect_sizes)
                statistical_score += 0.5 * large_effect_fraction
        
        # Replication agreement contribution (20%)
        replication_score = replication_agreement
        
        # Virtue alignment contribution (10%)
        virtue_contribution = virtue_score
        
        # Weighted combination
        confidence = (0.4 * criteria_score + 
                     0.3 * statistical_score +
                     0.2 * replication_score + 
                     0.1 * virtue_contribution)
        
        # Evidence strength modifier
        confidence *= (0.5 + 0.5 * evidence_strength)
        
        return min(1.0, confidence)
    
    def _make_collapse_decision(self, criteria_met: Dict[str, bool],
                              confidence: float, uncertainty: float,
                              replication_count: int,
                              statistical_summary: Dict[str, Any],
                              virtue_score: float) -> Tuple[str, str, str]:
        """Make the final collapse decision."""
        
        # Check minimum requirements
        min_confidence = self.thresholds['min_confidence']
        min_replications = self.thresholds['min_replications']
        max_uncertainty = self.thresholds['max_uncertainty']
        virtue_threshold = self.thresholds['virtue_threshold']
        
        # All success criteria must be met
        all_criteria_met = all(criteria_met.values()) if criteria_met else False
        
        # Determine status
        if (all_criteria_met and 
            confidence >= min_confidence and 
            uncertainty <= max_uncertainty and
            replication_count >= min_replications and
            virtue_score >= virtue_threshold):
            
            status = "truth"
            reasoning = (f"All criteria met with high confidence ({confidence:.2f}), "
                        f"low uncertainty ({uncertainty:.2f}), and sufficient replications ({replication_count})")
            recommendation = "Collapse to truth and publish as validated discovery"
            
        elif not any(criteria_met.values()) and confidence < 0.3:
            status = "refuted"
            reasoning = f"No success criteria met and low confidence ({confidence:.2f})"
            recommendation = "Mark as refuted and investigate alternative approaches"
            
        elif replication_count < min_replications:
            status = "needs_evidence"
            reasoning = f"Insufficient replications ({replication_count} < {min_replications})"
            recommendation = "Seek additional independent replications"
            
        elif uncertainty > max_uncertainty:
            status = "needs_evidence"
            reasoning = f"High uncertainty ({uncertainty:.2f} > {max_uncertainty})"
            recommendation = "Improve measurement precision and reduce uncertainty"
            
        elif confidence < min_confidence:
            status = "needs_evidence"
            reasoning = f"Low confidence ({confidence:.2f} < {min_confidence})"
            recommendation = "Gather more evidence or refine experimental design"
            
        elif virtue_score < virtue_threshold:
            status = "needs_evidence"
            reasoning = f"Low virtue score ({virtue_score:.2f} < {virtue_threshold})"
            recommendation = "Address ethical/safety concerns before proceeding"
            
        else:
            status = "needs_evidence"
            reasoning = "Some criteria met but insufficient for truth collapse"
            recommendation = "Continue investigation with focused improvements"
        
        return status, reasoning, recommendation


# Example usage and testing
if __name__ == "__main__":
    # Test the collapse rules engine
    rules = CollapseRules()
    
    # Example claim
    test_claim = {
        "id": "test_claim_001",
        "collapse_rules": {
            "success_criteria": {
                "removal_efficiency": {">=": 95},
                "stability_hours": {">=": 24}
            }
        }
    }
    
    # Example measurement result
    test_result = {
        "metrics": {
            "removal_efficiency": 97.5,
            "stability_hours": 36
        },
        "uncertainty": 0.15,
        "replications": [
            {"metrics": {"removal_efficiency": 96.8, "stability_hours": 34}},
            {"metrics": {"removal_efficiency": 98.1, "stability_hours": 38}},
            {"metrics": {"removal_efficiency": 97.6, "stability_hours": 36}}
        ],
        "virtue_vector_history": [[0.8, 0.7, 0.9, 0.6], [0.85, 0.75, 0.88, 0.65]]
    }
    
    test_candidate = {"safety_score": 0.8, "environmental_score": 0.9}
    
    verdict = rules.judge_claim(test_claim, test_candidate, test_result)
    print(f"🎯 Test verdict: {verdict.status} (confidence: {verdict.confidence:.2f})")
    print(f"📝 Reasoning: {verdict.reasoning}")
    print(f"💡 Recommendation: {verdict.recommendation}")
