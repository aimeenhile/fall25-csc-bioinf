# bio_codon/motifs/matrix.py

"""Support for various forms of sequence motif matrices.

Implementation of frequency (count) matrices, position-weight matrices,
and position-specific scoring matrices.
"""
import math
from python import numbers
import numpy as np
from python.Bio.Seq import Seq
from typing import Dict, Tuple, List, Optional, Union, Sequence, Callable

# A utility to calculate IUPAC degenerate consensus (simplified)
# W (A/T), S (G/C), K (G/T), M (A/C), R (A/G), Y (C/T), V (A/C/G), H (A/C/T), D (A/G/T), B (C/G/T)
IUPAC_CODE: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
    "W": "AT", "S": "GC", "K": "GT", "M": "AC",
    "R": "AG", "Y": "CT", "V": "ACG", "H": "ACT",
    "D": "AGT", "B": "CGT", "N": "ACGT"
}
# Inverse mapping for degeneracy calculation (simplified for common cases)
DEGENERATE_MAP: Dict[Tuple[str, ...], str] = {
    ('A', 'T'): 'W', ('G', 'C'): 'S', ('A', 'G'): 'R', ('C', 'T'): 'Y',
    ('A', 'C', 'G'): 'V',
}

def _get_score(matrix, base, position) -> float:
    """Safely get score, returning 0.0 or a very low score for missing bases (PRIVATE)."""
    try:
        return matrix[base][position]
    except (KeyError, IndexError):
        # A common practice in PSSMs is to return a very low score for unknown bases (like 'N')
        # Here we'll return the lowest possible score in the matrix or 0.0 if not defined
        min_score = min(min(col) for col in matrix.values()) if matrix.values() else -100.0
        return min_score


class GenericPositionMatrix(dict):
    """Base class for the support of position matrix operations."""
    
    length: int
    alphabet: str

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the class."""
        self.length = 0
        self.alphabet = alphabet
        
        for letter in alphabet:
            if letter not in values:
                raise ValueError(f"Missing count data for alphabet letter: {letter}")
                
            col_data = values[letter]
            if self.length == 0:
                self.length = len(col_data)
            elif self.length != len(col_data):
                raise ValueError("Data has inconsistent lengths")
            
            # Cast any numpy floats into Python floats:
            self[letter] = [float(_) for _ in col_data]

    def __str__(self) -> str:
        """Return a string containing bases and counts/weights of the alphabet in the Matrix."""
        # Max width for position index
        index_width = len(str(self.length - 1)) + 2
        
        # Header row: '    0 1 2 3...'
        words = [f"{i:>{index_width}}" for i in range(self.length)]
        line = " " * len(self.alphabet) + " " + " ".join(words)
        lines = [line]
        
        # Data rows: 'A 1.00 0.00 1.00...'
        for letter in self.alphabet:
            words = [f"{value:{index_width}.2f}" for value in self[letter]]
            line = f"{letter} " + " ".join(words)
            lines.append(line)
            
        return "\n".join(lines)
        
    def __getitem__(self, key: Union[str, Tuple[str, int], Tuple[int, object]]):
        """Allow access like m.counts['A', 3] or m.counts[:, 3].
        
        The type hint Tuple[int, object] is used to represent a slice like m[:, 3] 
        where the slice(None) object is interpreted as 'object' by Codon.
        """
        if isinstance(key, str):
            # Access like m['A'] -> returns list of floats for base A
            return super().__getitem__(key)
        
        if isinstance(key, tuple):
            if len(key) == 2:
                base_or_index, pos = key
                
                # m.counts[base, position]
                if isinstance(base_or_index, str):
                    if base_or_index not in self.alphabet:
                        raise KeyError(f"Base '{base_or_index}' not in alphabet '{self.alphabet}'")
                    if not (0 <= pos < self.length):
                        raise IndexError(f"Position index {pos} out of range [0, {self.length-1}]")
                    return self[base_or_index][pos]
                
                # m.counts[base_index, position] - simplified access
                elif isinstance(base_or_index, int):
                    if not (0 <= base_or_index < len(self.alphabet)):
                        raise IndexError(f"Base index {base_or_index} out of range [0, {len(self.alphabet)-1}]")
                    base = self.alphabet[base_or_index]
                    return self[base][pos]
                
                # m.counts[:, position] - access to a column (slice(None) is of type object)
                elif base_or_index == slice(None) and isinstance(pos, int):
                    if not (0 <= pos < self.length):
                        raise IndexError(f"Position index {pos} out of range [0, {self.length-1}]")
                    # Returns column as a dictionary
                    return {base: self[base][pos] for base in self.alphabet}
                    
        raise KeyError("Invalid key for PositionMatrix access.")


class CountsMatrix(GenericPositionMatrix):
    """Counts matrix (frequency matrix)."""

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the class."""
        super().__init__(alphabet, values)

    def normalize(self, pseudocounts: float = 0.0) -> 'PositionWeightMatrix':
        """Normalize the counts to obtain a Position Weight Matrix (PWM)."""
        pwm_values: Dict[str, List[float]] = defaultdict(list)
        
        # Total number of sequences/instances
        num_instances: float = sum(self[base][0] for base in self.alphabet)
        
        # Effective sequence count after pseudocounts
        N = num_instances
        if N == 0:
            N = 1.0 
             
        denominator = N + pseudocounts * len(self.alphabet)

        for i in range(self.length):
            for base in self.alphabet:
                count_with_pc = self[base][i] + pseudocounts
                freq = count_with_pc / denominator
                pwm_values[base].append(freq)
                
        return PositionWeightMatrix(self.alphabet, pwm_values)


class PositionWeightMatrix(GenericPositionMatrix):
    """Position Weight Matrix (PWM)."""

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the class."""
        super().__init__(alphabet, values)

    def log_odds(self, background: Optional[Dict[str, float]] = None) -> 'PositionSpecificScoringMatrix':
        """Calculate the PSSM (Position Specific Scoring Matrix) from the PWM."""
        if background is None:
            # Uniform background
            background = dict.fromkeys(self.alphabet, 1.0 / len(self.alphabet))
            
        pssm_values: Dict[str, List[float]] = defaultdict(list)
        
        for base in self.alphabet:
            bg_prob = background.get(base, 0.0)
            if bg_prob == 0.0:
                raise ValueError(f"Background probability for base '{base}' is zero.")
            
            for i in range(self.length):
                pwm_val = self[base][i]
                if pwm_val == 0.0:
                    # Replace 0.0 with a very small number to avoid log(0)
                    log_odds_score = math.log2(1e-10 / bg_prob)
                else:
                    log_odds_score = math.log2(pwm_val / bg_prob)
                    
                pssm_values[base].append(log_odds_score)
                
        return PositionSpecificScoringMatrix(self.alphabet, pssm_values, background)


class PositionSpecificScoringMatrix(GenericPositionMatrix):
    """Position Specific Scoring Matrix (PSSM)."""
    
    background: Dict[str, float]
    
    def __init__(self, alphabet: str, values: Dict[str, List[float]], background: Dict[str, float]):
        """Initialize the class."""
        super().__init__(alphabet, values)
        self.background = background
        
    @property
    def min_score(self) -> float:
        """Return the minimum possible score of a sequence with this PSSM."""
        return sum(min(scores) for scores in self.values())

    @property
    def max_score(self) -> float:
        """Return the maximum possible score of a sequence with this PSSM."""
        return sum(max(scores) for scores in self.values())

    def calculate(self, sequence: Seq) -> List[float]:
        """Calculate the score for all positions in the given sequence."""
        if len(sequence) < self.length:
            return [] # Sequence is too short
            
        scores: List[float] = []
        for start in range(len(sequence) - self.length + 1):
            subsequence = sequence[start:start + self.length]
            score = 0.0
            
            for i in range(self.length):
                base = str(subsequence)[i]
                # If the base is not in the alphabet, use a placeholder score
                score += _get_score(self, base, i)
            
            scores.append(score)
            
        return scores

    def distribution(self, background: Optional[Dict[str, float]] = None, precision: int = 10**3):
        """Calculate the score distribution for the PSSM."""
        from .thresholds import ScoreDistribution
        
        if background is None:
            background = self.background
            
        # The ScoreDistribution init will calculate the actual distribution
        return ScoreDistribution(pssm=self, background=background, precision=precision)
        
    def __sub__(self, other: 'PositionSpecificScoringMatrix') -> float:
        """Returns the correlation between two PSSMs. The two PSSMs are aligned in an optimal way."""
        
        total_max_correlation = -float('inf')
        
        # Range includes negative offsets (other starts before self)
        # Offset is defined as: other starts at self[offset]
        for offset in range(-(other.length - 1), self.length):
            
            if offset >= 0:
                # self.cc(other, offset): other starts at self[offset]
                correlation = self.cc(other, offset)
            else:
                # offset < 0 means other starts -offset positions *before* self.
                # This is equivalent to self starting at other[-offset]
                # other.cc(self, -offset)
                correlation = other.cc(self, -offset)

            total_max_correlation = max(total_max_correlation, correlation)

        return total_max_correlation if total_max_correlation != -float('inf') else 0.0

    def cc(self, other: 'PositionSpecificScoringMatrix', offset: int) -> float:
        """Return the similarity score based on pearson correlation at the given offset."""
        
        # Overlapping range in 'self' coordinates
        overlap_start_self = max(0, offset)
        overlap_end_self = min(self.length, other.length + offset)
        
        overlap_length = overlap_end_self - overlap_start_self
        
        if overlap_length <= 0:
            return 0.0

        letters = self.alphabet
        
        # Number of terms in the sum: (overlap_length) * (len(letters))
        norm_factor = overlap_length * len(letters)
        
        sx = 0.0
        sy = 0.0
        sxx = 0.0
        sxy = 0.0
        syy = 0.0

        for i in range(overlap_length):
            pos_self = overlap_start_self + i
            pos_other = pos_self - offset
            
            for letter in letters:
                x = self[letter, pos_self]
                y = other[letter, pos_other]
                
                sx += x
                sy += y
                sxx += x * x
                sxy += x * y
                syy += y * y
                
        # Calculate expected values (normalized by the number of terms)
        Ex = sx / norm_factor
        Ey = sy / norm_factor
        Exx = sxx / norm_factor
        Exy = sxy / norm_factor
        Eyy = syy / norm_factor
        
        # Pearson correlation formula
        numerator = Exy - Ex * Ey
        denominator_sq_x = Exx - Ex * Ex
        denominator_sq_y = Eyy - Ey * Ey
        
        # Handle cases where variance is zero
        if denominator_sq_x < 1e-9 or denominator_sq_y < 1e-9:
            return 0.0 if abs(numerator) < 1e-9 else 1.0 if numerator > 0 else -1.0
            
        denominator = math.sqrt(denominator_sq_x * denominator_sq_y)
        
        if denominator == 0.0:
            return 0.0
            
        correlation = numerator / denominator
        
        # Clip to [-1, 1] due to potential floating point errors
        return max(-1.0, min(1.0, correlation))