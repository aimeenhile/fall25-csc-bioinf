# bio_codon/motifs/thresholds.py

"""Approximate calculation of appropriate thresholds for motif finding."""

import math
import numpy as np
from typing import Optional, List, Dict, Union, Tuple
from .matrix import PositionSpecificScoringMatrix

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
    distribution: Dict[int, float] 

    def __init__(self, pssm: PositionSpecificScoringMatrix, background: Dict[str, float], precision: int = 10000):
        """
        Initialize the distribution calculator.
        
        Uses dynamic programming to compute the score distribution.
        """
        self.pssm = pssm
        self.background = background
        self.precision = precision
        
        # Calculate scores and probabilities for each column
        # list of {score_int: prob} dictionaries for each column
        self.column_distributions: List[Dict[int, float]] = [] 
        
        self.min_score: float = 0.0
        self.max_score: float = 0.0
        
        for j in range(pssm.length):
            col_scores: Dict[int, float] = {}
            for base in pssm.alphabet:
                pssm_score = pssm[base, j]
                bg_prob = background.get(base, 0.25) # Default to 0.25 if not specified
                
                # Rescale the score to an integer for dynamic programming, based on precision
                # E.g., Score of 1.25 with precision 100 becomes 125
                score_int = int(round(pssm_score * precision))
                
                col_scores[score_int] = col_scores.get(score_int, 0.0) + bg_prob

            self.column_distributions.append(col_scores)
            
            # Keep track of min/max score for the whole PSSM
            column_min = min(col_scores.keys())
            column_max = max(col_scores.keys())
            self.min_score += column_min
            self.max_score += column_max

        # Dynamic Programming to compute the full distribution
        # Start with the distribution of the first column
        self.distribution: Dict[int, float] = self.column_distributions[0]
        
        for i in range(1, pssm.length):
            next_distribution: Dict[int, float] = defaultdict(float)
            current_col = self.column_distributions[i]
            
            # Convolve the current distribution with the next column's distribution
            for score1, prob1 in self.distribution.items():
                for score2, prob2 in current_col.items():
                    total_score = score1 + score2
                    total_prob = prob1 * prob2
                    next_distribution[total_score] += total_prob
            
            self.distribution = next_distribution

        # Normalize the scores back to float space 
        self.min_score /= precision
        self.max_score /= precision
        self.step = 1.0 / precision # Score unit is 1/precision
        self.n_points = len(self.distribution) # Number of discrete points

    def _get_target_score(self, target_rate: float, is_fpr: bool) -> float:
        """Internal helper to find the score corresponding to a target cumulative rate (FPR or FNR)."""
        
        sorted_scores = sorted(self.distribution.keys())
        
        # Pre-calculate cumulative probabilities
        if is_fpr:
            # FPR (False Positive Rate): P(S >= s) -> Sum of probabilities from high score down
            scores_to_use = list(reversed(sorted_scores))
        else:
            # FNR (False Negative Rate): P(S <= s) -> Sum of probabilities from low score up
            scores_to_use = sorted_scores
            
        cumulative_prob = 0.0
        target_score_int = sorted_scores[0] # Fallback
        
        for score in scores_to_use:
            prob = self.distribution.get(score, 0.0)
            cumulative_prob += prob
            
            if cumulative_prob >= target_rate:
                target_score_int = score
                break
        
        return target_score_int / self.precision

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
        fpr_prob = np.cumsum([self.distribution.get(s, 0.0) for s in reversed(sorted_scores)])[::-1]
        
        min_diff = float('inf')
        best_score_index = -1
        
        for i in range(len(sorted_scores)):
            current_fpr = fpr_prob[i]
            current_fnr = fnr_prob[i]
            
            # Find the score 's' that minimizes |FPR - balance_factor * FNR|
            diff = abs(current_fpr - balance_factor * current_fnr)
            if diff < min_diff:
                min_diff = diff
                best_score_index = i
            
        if best_score_index != -1:
            return sorted_scores[best_score_index] / self.precision
        else:
            return self.min_score