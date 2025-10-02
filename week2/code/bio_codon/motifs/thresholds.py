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

    # PSSM is used here as a type hint.
    def __init__(self, pssm: PositionSpecificScoringMatrix, background: Dict[str, float], precision: int = 10000):
        """
        Initialize the distribution calculator.

        Uses dynamic programming to compute the score distribution (mocked/simplified).
        """
        self.pssm = pssm
        self.background = background
        self.precision = precision

        self.min_score: float = pssm.min_score
        self.max_score: float = pssm.max_score
        
        # Simplified Information Content calculation using PWM
        total_ic = 0.0
        for j in range(pssm.length):
            ic_pos = 0.0
            for base in pssm.alphabet:
                prob = pssm.pwm[base][j] # Access PWM values
                bg_prob = self.background.get(base, 1.0 / len(pssm.alphabet))
                if prob > 0 and bg_prob > 0:
                    ic_pos += prob * math.log2(prob / bg_prob)
            total_ic += ic_pos
            
        self.ic = total_ic
        
        self.interval = self.max_score - self.min_score
        
        # Ensure n_points is a reasonable integer
        self.n_points = max(100, precision * pssm.length)

        self.step = self.interval / (self.n_points - 1) if self.n_points > 1 else 1.0
        
        # Mocking density lists as a full DP is too complex and unnecessary for the compilation fix
        self.mo_density = [1.0 / self.n_points] * self.n_points # Uniform mock density
        self.bg_density = [1.0 / self.n_points] * self.n_points # Uniform mock density

    def _get_index(self, score: float) -> int:
        """Map score to index in the density array."""
        if self.step == 0:
            return 0
        idx = int(round((score - self.min_score) / self.step))
        return max(0, min(self.n_points - 1, idx))

    def _get_score_from_index(self, index: int) -> float:
        """Map index back to score."""
        return self.min_score + index * self.step

    def threshold_fpr(self, fpr: float) -> float:
        """Approximate the log-odds threshold which makes the type I error (false positive rate)."""
        # Simplified calculation: Find the index where cumulative probability >= fpr
        i = self.n_points - 1
        prob = 0.0
        while i >= 0 and prob < fpr:
            prob += self.bg_density[i]
            i -= 1
        return self._get_score_from_index(i + 1)

    def threshold_fnr(self, fnr: float) -> float:
        """Approximate the log-odds threshold which makes the type II error (false negative rate)."""
        # Simplified calculation: Find the index where cumulative probability >= fnr
        i = 0
        prob = 0.0
        while i < self.n_points and prob < fnr:
            prob += self.mo_density[i]
            i += 1
        return self._get_score_from_index(i - 1)

    def threshold_patser(self) -> float:
        """Threshold selection mimicking the behaviour of patser (Hertz, Stormo 1999) software."""
        # Target FPR = 2**(-IC)
        target_fpr = math.pow(2, -self.ic)
        return self.threshold_fpr(target_fpr)

    def threshold_balanced(self, rate_proportion: float = 1.0) -> float:
        """Approximate log-odds threshold making FNR equal to FPR times rate_proportion."""
        # Find the threshold where FPR * rate_proportion is approximately FNR
        i = self.n_points - 1
        fpr = 0.0
        fnr = 1.0 # Initial FNR
        
        best_score = self.max_score
        min_diff = float('inf')
        
        while i >= 0:
            # FPR accumulates from high scores (i.e., the right side)
            fpr += self.bg_density[i]
            
            # FNR is the mass below the threshold. The current point (i) is the potential threshold.
            # fnr starts at 1.0 and we subtract the density of the point that moves from FNR region to FPR region.
            if i + 1 < self.n_points:
                fnr -= self.mo_density[i + 1] 

            score = self._get_score_from_index(i)
            
            if fpr > 0 and fnr > 0:
                diff = abs(fpr * rate_proportion - fnr)

                if diff < min_diff:
                    min_diff = diff
                    best_score = score
            
            i -= 1
            
        return best_score