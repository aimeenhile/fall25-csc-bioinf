# bio_codon/motifs/thresholds.py

"""Approximate calculation of appropriate thresholds for motif finding."""

import math
import numpy as np
from typing import list, dict
from .matrix import PositionWeightMatrix

class ScoreDistribution:
    """
    Class representing approximate score distribution for a given motif.
    
    Provides a simplified method for threshold calculation (actual dynamic 
    programming is omitted for porting simplicity).
    """
    pwm: PositionWeightMatrix
    min_score: float
    max_score: float
    interval: float
    n_points: int
    step: float

    def __init__(self, pwm: PositionWeightMatrix, precision: int = 1000):
        """Initializes the ScoreDistribution based on the PWM."""
        self.pwm = pwm
        
        self.min_score = pwm.min_score
        self.max_score = pwm.max_score
        self.interval = self.max_score - self.min_score
        # Precision fields are kept for structure, but not used in the simplified calculation
        self.n_points = precision * pwm.length
        self.step = self.interval / float(self.n_points - 1) if self.n_points > 1 else 0.0

    def __str__(self) -> str:
        """String representation."""
        return f"<ScoreDistribution for PWM of length {self.pwm.length}, range [{self.min_score:.2f}, {self.max_score:.2f}]>"

    def threshold_for_p_value(self, p_value: float) -> float:
        """
        Simplified approximation: Finds the score threshold corresponding to a given P-value.
        
        Maps the negative logarithm of the P-value linearly onto the score range.
        """
        if p_value <= 0.0:
            return self.max_score
        if p_value >= 1.0:
            return self.min_score
            
        # Cap log_p calculation for stability (e.g., P=1e-8 gives 8.0)
        log_p_capped: float = -math.log10(max(p_value, 1e-8))
        
        # Normalize the log scale
        norm_log_p: float = min(log_p_capped / 8.0, 1.0) 
        
        # Linearly map the normalized value to the score interval
        return self.min_score + norm_log_p * self.interval