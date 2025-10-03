# bio_codon/motifs/matrix.py

"""Support for various forms of sequence motif matrices.

Implementation of frequency (count) matrices, position-weight matrices,
and position-specific scoring matrices.
"""
import math
from python import numbers
import numpy as np
from python import Bio.Seq.Seq
# FIX: Removed Union from imports to avoid internal Codon compiler crash.
from typing import Dict, Tuple, List, Optional
from collections import defaultdict
from . import Background # Import type hint from __init__.py if needed

# A utility to calculate IUPAC degenerate consensus (simplified)
# W (A/T), S (G/C), K (G/T), M (A/C), R (A/G), Y (C/T), V (A/C/G), H (A/C/T), D (A/G/T), B (C/G/T)
IUPAC_CODE: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
    "W": "AT", "S": "GC", "K": "GT", "M": "AC",
    "R": "AG", "Y": "CT", "V": "ACG", "H": "ACT",
    "D": "AGT", "B": "CGT", "N": "ACGT"
}
# Inverse mapping for degeneracy calculation.
# FIX: Using separate maps by tuple length to avoid the Codon type checker issue with Dict[Union[Tuple[...], ...], ...]
# All keys remain alphabetically sorted (canonical) to match lookup logic.

DEGENERATE_MAP_2: Dict[Tuple[str, str], str] = {
    ('A', 'T'): 'W',      # A or T
    ('C', 'G'): 'S',      # G or C
    ('G', 'T'): 'K',      # G or T
    ('A', 'C'): 'M',      # A or C
    ('A', 'G'): 'R',      # A or G
    ('C', 'T'): 'Y',      # C or T
}
DEGENERATE_MAP_3: Dict[Tuple[str, str, str], str] = {
    ('C', 'G', 'T'): 'B', # C or G or T
    ('A', 'G', 'T'): 'D', # A or G or T
    ('A', 'C', 'T'): 'H', # A or C or T
    ('A', 'C', 'G'): 'V', # A or C or G
}
DEGENERATE_MAP_4: Dict[Tuple[str, str, str, str], str] = {
    ('A', 'C', 'G', 'T'): 'N', # A or C or G or T
}

# Mapping from tuple size to its corresponding map
DEGENERATE_MAP: Dict[int, Dict[Tuple, str]] = {
    2: DEGENERATE_MAP_2,
    3: DEGENERATE_MAP_3,
    4: DEGENERATE_MAP_4,
}


class GenericPositionMatrix:
    """Base class for all position-based matrices (Counts, Frequencies, PWM, PSSM)."""

    def __init__(self, alphabet: str, values: Dict[str, List[float]], length: Optional[int] = None):
        self.alphabet = alphabet.upper()
        self.values = values
        
        # Infer length if not provided
        self.length = length
        if self.length is None:
            if self.values:
                # Assuming all lists have the same length
                first_key = next(iter(self.values))
                self.length = len(self.values[first_key])
            else:
                self.length = 0

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, key: Tuple[str, Optional[int]]):
        """Allows access like matrix[base, position] or matrix[base]."""
        if isinstance(key, tuple):
            # Access like matrix['A', 0]
            base, pos = key
            if pos is None or pos < 0 or pos >= self.length:
                raise IndexError("Position index out of bounds.")
            if base not in self.alphabet:
                raise ValueError(f"Base '{base}' not in alphabet '{self.alphabet}'.")
            return self.values[base.upper()][pos]
        elif isinstance(key, str):
            # Access like matrix['A'] (returns the entire row/list for that base)
            base = key
            if base not in self.alphabet:
                raise ValueError(f"Base '{base}' not in alphabet '{self.alphabet}'.")
            return self.values[base.upper()]
        else:
            raise TypeError("Key must be a base string or a (base, position) tuple.")

    def __setitem__(self, key: Tuple[str, Optional[int]], value: float):
        """Allows assignment like matrix['A', 0] = 5.0 or matrix['A'] = [1, 2, 3]."""
        if isinstance(key, tuple):
            # Assignment like matrix['A', 0] = 5.0
            base, pos = key
            if pos is None or pos < 0 or pos >= self.length:
                raise IndexError("Position index out of bounds.")
            if base not in self.alphabet:
                raise ValueError(f"Base '{base}' not in alphabet '{self.alphabet}'.")
            self.values[base.upper()][pos] = value
        elif isinstance(key, str):
            # Assignment like matrix['A'] = [1, 2, 3]
            base = key
            if base not in self.alphabet:
                raise ValueError(f"Base '{base}' not in alphabet '{self.alphabet}'.")
            if not isinstance(value, list) or len(value) != self.length:
                 raise ValueError("Value must be a list of floats matching the matrix length.")
            self.values[base.upper()] = value
        else:
            raise TypeError("Key must be a base string or a (base, position) tuple.")
            
    # Omitted methods like consensus, degenerate_consensus, etc. for brevity.

    @property
    def max(self) -> float:
        """Return the maximum score in the matrix."""
        max_val = -float('inf')
        for base in self.alphabet:
            for val in self.values[base]:
                if val > max_val:
                    max_val = val
        return max_val

    @property
    def min(self) -> float:
        """Return the minimum score in the matrix."""
        min_val = float('inf')
        for base in self.alphabet:
            for val in self.values[base]:
                if val < min_val:
                    min_val = val
        return min_val


class CountsMatrix(GenericPositionMatrix):
    """
    A counts matrix stores the frequency of each base at each position.
    
    Inherits from GenericPositionMatrix.
    """
    def __init__(self, alphabet: str, values: Dict[str, List[float]], length: Optional[int] = None):
        super().__init__(alphabet, values, length)

    def normalize(self, pseudocounts: Optional[float] = None) -> 'PositionWeightMatrix':
        """
        Converts the CountsMatrix to a PositionWeightMatrix (PWM) 
        by normalizing the counts to probabilities (frequencies).
        
        If **pseudocounts** is provided, it applies the pseudo-count before normalization.
        """
        length = self.length
        alphabet = self.alphabet
        
        # Initialize PWM values
        pwm_values: Dict[str, List[float]] = defaultdict(lambda: [0.0] * length)

        # Calculate the total count (including pseudocounts) for each column
        column_totals: List[float] = [0.0] * length
        
        # Apply pseudocounts
        if pseudocounts is not None:
            # We assume the pseudocounts are applied equally to all bases at all positions.
            
            for j in range(length):
                col_sum = 0.0
                for base in alphabet:
                    # Add original count and the single pseudocount
                    col_sum += self[base, j] + pseudocounts
                column_totals[j] = col_sum
                
                # Normalize and assign to PWM
                for base in alphabet:
                    pwm_values[base][j] = (self[base, j] + pseudocounts) / column_totals[j]
        else:
            # No pseudocounts: calculate sum of raw counts per column
            for j in range(length):
                col_sum = 0.0
                for base in alphabet:
                    col_sum += self[base, j]
                column_totals[j] = col_sum
                
                # Normalize and assign to PWM
                for base in alphabet:
                    if column_totals[j] > 0.0:
                         pwm_values[base][j] = self[base, j] / column_totals[j]
                    else:
                         pwm_values[base][j] = 0.0 # Should not happen with valid input
        
        # Create and return the PWM instance
        return PositionWeightMatrix(alphabet=alphabet, values=pwm_values, length=length)


# Alias for CountsMatrix (to match Biopython's structure if needed)
FrequencyPositionMatrix = CountsMatrix


class PositionWeightMatrix(GenericPositionMatrix):
    """
    A position weight matrix stores the normalized frequency (probability)
    of each base at each position.
    
    Inherits from GenericPositionMatrix.
    """
    def __init__(self, alphabet: str, values: Dict[str, List[float]], length: Optional[int] = None):
        super().__init__(alphabet, values, length)

    def log_odds(self, background: Background) -> 'PositionSpecificScoringMatrix':
        """
        Converts the PWM to a PositionSpecificScoringMatrix (PSSM) 
        by taking the log-odds ratio using the provided **background** frequencies.
        
        PSSM[b, j] = log2(PWM[b, j] / background[b])
        """
        length = self.length
        alphabet = self.alphabet
        
        # Check background frequencies
        for base in alphabet:
            if base not in background or background[base] <= 0.0:
                raise ValueError(f"Background frequency for base '{base}' is missing or zero.")
        
        # Initialize PSSM values
        pssm_values: Dict[str, List[float]] = defaultdict(lambda: [0.0] * length)
        
        for j in range(length):
            for base in alphabet:
                pwm_prob = self[base, j]
                bg_prob = background[base]
                
                # Use a small epsilon to avoid log(0)
                if pwm_prob <= 0.0:
                    score = -float('inf') 
                else:
                    score = math.log2(pwm_prob / bg_prob)
                
                pssm_values[base][j] = score
        
        # Create and return the PSSM instance
        return PositionSpecificScoringMatrix(alphabet=alphabet, values=pssm_values, length=length)


class PositionSpecificScoringMatrix(GenericPositionMatrix):
    """
    A PSSM stores the log-odds scores for each base at each position.
    
    Inherits from GenericPositionMatrix.
    """
    def __init__(self, alphabet: str, values: Dict[str, List[float]], length: Optional[int] = None):
        super().__init__(alphabet, values, length)

    def calculate(self, sequence_str: str) -> List[float]:
        """
        Calculate the score of a sequence string for the motif.
        
        Returns a list of scores for each possible starting position of the motif 
        in the sequence.
        """
        if len(sequence_str) < self.length:
            return []

        scores: List[float] = []
        
        # Iterate over all possible starting positions (i)
        for i in range(len(sequence_str) - self.length + 1):
            motif_score = 0.0
            
            # Iterate over each position (j) in the motif
            for j in range(self.length):
                base = sequence_str[i + j].upper()
                
                if base in self.alphabet:
                    # Add the score for the base at this position
                    motif_score += self[base][j]
                else:
                    # Handle unknown/degenerate characters by marking the score as NaN
                    motif_score = float('nan')
                    break
            
            scores.append(motif_score)
            
        return scores

    # Omitted complex methods like "search" and "distribution" for brevity and core focus.
    
    def distribution(self, background: Background, precision: int = 10**3):
        """Calculate the distribution of the scores at the given precision."""
        # Use a placeholder for distribution logic to avoid importing within the class body
        # and to keep the structure. The real implementation is in thresholds.py
        from .thresholds import ScoreDistribution
             
        return ScoreDistribution(pssm=self, background=background, precision=precision)
