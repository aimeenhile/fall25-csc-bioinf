# bio_codon/motifs/thresholds.py

"""Approximate calculation of appropriate thresholds for motif finding."""

import math
import numpy as np
from typing import Optional, List, Dict, Tuple, Union
from .matrix import PositionSpecificScoringMatrix
from collections import defaultdict
from .matrix import Background


class ScoreDistribution:
    """
    Class representing approximate score distribution for a given motif.
    
    Provides a simplified method for threshold calculation (actual dynamic 
    programming is omitted for porting simplicity).
    """
    min_score: float
    interval: float
    n_points: int
    ic: float # Information content / mean
    step: float
    mo_density: List[float] # Motif occurrence density
    bg_density: List[float] # Background density

    # Store the actual distribution as a dictionary of score -> probability
    # FIX: Use Dict[float, float] for distribution
    distribution: Dict[float, float] 

    def __init__(self, pssm: PositionSpecificScoringMatrix, background: Background, precision: int = 10000):
        """
        Initialize the distribution calculator.
        
        Uses dynamic programming (convolution) to compute the score distribution.
        """
        self.pssm = pssm
        self.background = background
        self.precision = precision
        
        # Calculate scores and probabilities for each column
        # list of {score_int: prob} dictionaries for each column
        self.column_distributions: List[Dict[int, float]] = [] 
        
        self.min_score: float = 0.0
        self.max_score: float = 0.0
        
        # Find the total background probability of the alphabet used
        sum_bg = sum(self.background.get(base, 0.0) for base in pssm.alphabet)
        if sum_bg == 0.0:
            raise ValueError("Total background probability is zero for the motif's alphabet.")
        
        for j in range(pssm.length):
            col_dist: Dict[int, float] = defaultdict(float)
            
            for base in pssm.alphabet:
                score = pssm[base, j]
                # Scale score to an integer point for DP binning
                score_int = int(round(score * precision)) 
                
                # Probability of this score is proportional to the background probability of the base
                prob = self.background.get(base, 0.0) / sum_bg
                
                col_dist[score_int] = prob
                
                # Update min/max scores
                if j == 0 and base == pssm.alphabet[0]:
                    self.min_score = self.max_score = score
                else:
                    self.min_score = min(self.min_score, score)
                    self.max_score = max(self.max_score, score)
            
            self.column_distributions.append(col_dist)

        # Total min/max score across all columns
        self.min_score = int(round(self.min_score * precision)) * pssm.length / precision
        self.max_score = int(round(self.max_score * precision)) * pssm.length / precision
        
        # --- Dynamic Programming / Convolution ---
        # Initialize the main distribution with the first column's distribution
        # The keys are now scaled integer scores (score * precision)
        self.distribution = defaultdict(float)
        
        if self.column_distributions:
            self.distribution = self.column_distributions[0]
            
            for j in range(1, len(self.column_distributions)):
                next_col_dist = self.column_distributions[j]
                new_dist: Dict[int, float] = defaultdict(float)
                
                # Convolve current distribution with the next column's distribution
                for score1, prob1 in self.distribution.items():
                    for score2, prob2 in next_col_dist.items():
                        new_score = score1 + score2
                        new_prob = prob1 * prob2
                        new_dist[new_score] += new_prob
                
                self.distribution = new_dist
        
        # Convert integer scores back to float scores
        final_distribution: Dict[float, float] = {}
        total_prob = 0.0
        for score_int, prob in self.distribution.items():
            score_float = score_int / precision
            final_distribution[score_float] = prob
            total_prob += prob

        # Re-normalize the final distribution (should be close to 1.0)
        if total_prob > 0.0:
             if abs(total_prob - 1.0) > 1e-6:
                 for score in final_distribution:
                     final_distribution[score] /= total_prob
        
        self.distribution = final_distribution
        
        # Calculate Information Content (IC) / Mean Score
        self.ic = sum(score * prob for score, prob in self.distribution.items())
        
        # Set interval for threshold finding
        self.interval = 1.0 / precision


    def _get_target_score(self, target_value: float, is_fpr: bool) -> float:
        """Helper to find the score threshold corresponding to a target FPR or FNR."""
        if not (0.0 < target_value <= 1.0):
            raise ValueError("Target value (FPR/FNR) must be between 0.0 and 1.0.")
            
        # Get scores and probabilities as sorted lists
        sorted_scores = sorted(self.distribution.keys())
        probabilities = [self.distribution[s] for s in sorted_scores]
        
        if is_fpr:
            # FPR (False Positive Rate): Cumulative probability of scores >= threshold
            # Iterate high to low
            sorted_scores.reverse()
            probabilities.reverse()
            
        # FNR (False Negative Rate): Cumulative probability of scores < threshold
        # Iterate low to high (already sorted if is_fpr is False)

        cumulative_prob = 0.0
        
        # The threshold is the score *at or above which* the cumulative probability first exceeds the target.
        for i, score in enumerate(sorted_scores):
            prob = probabilities[i]
            cumulative_prob += prob
            
            if cumulative_prob >= target_value:
                # This score (and everything higher, if FPR) or lower (if FNR)
                # is the threshold.
                return score
        
        # Should not happen if total probability is 1.0
        if is_fpr:
            return self.min_score # If target is 1.0, return min score
        else:
            return self.max_score # If target is 1.0, return max score

    def threshold_fpr(self, target_fpr: float) -> float:
        """Approximate the log-odds threshold which makes the type I error (false positive rate)."""
        return self._get_target_score(target_fpr, is_fpr=True)
    
    def threshold_fnr(self, target_fnr: float) -> float:
        """Approximate the log-odds threshold which makes the type II error (false negative rate)."""
        return self._get_target_score(target_fnr, is_fpr=False)

    def threshold_balanced(self, balance_factor: float = 1.0) -> float:
        """Approximate log-odds threshold making FNR equal to FPR times balance_factor."""
        
        sorted_scores = sorted(self.distribution.keys())
        
        # Pre-calculate cumulative probabilities (low-to-high) for FNR
        fnr_prob = np.cumsum([self.distribution.get(s, 0.0) for s in sorted_scores])
        # Pre-calculate cumulative probabilities (high-to-low) for FPR
        fpr_prob_raw = [self.distribution.get(s, 0.0) for s in reversed(sorted_scores)]
        # Need to reverse the cumsum to map back to the original sorted_scores order
        fpr_prob = np.cumsum(fpr_prob_raw)[::-1] 
        
        min_diff = float('inf')
        best_score_index = -1
        
        # Iterate through scores in ascending order
        for i in range(len(sorted_scores)):
            # FPR at threshold s: P(score >= s) -> fpr_prob[i]
            current_fpr = fpr_prob[i] 
            
            # FNR at threshold s: P(score < s) -> sum of all probabilities up to sorted_scores[i-1]
            current_fnr = fnr_prob[i-1] if i > 0 else 0.0
            
            # Find the score 's' that minimizes |FPR - balance_factor * FNR|
            diff = abs(current_fpr - balance_factor * current_fnr)
            if diff < min_diff:
                min_diff = diff
                best_score_index = i
            
        if best_score_index != -1:
            return sorted_scores[best_score_index]
        else:
            return self.ic # Fallback to mean/IC
