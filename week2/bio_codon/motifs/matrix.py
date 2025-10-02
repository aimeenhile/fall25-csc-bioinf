# bio_codon/motifs/matrix.py

"""Support for various forms of sequence motif matrices.

Implementation of frequency (count) matrices, position-weight matrices,
and position-specific scoring matrices.
"""
import math
from python import numbers
import numpy as np
from python.Bio.Seq import Seq
from typing import Dict, Tuple, List, Optional, Union, Any

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

    def __init__(self, alphabet: str, values: Dict[str, List[Union[int, float]]]):
        """Initialize the class."""
        self.length = None
        for letter in alphabet:
            if self.length is None:
                self.length = len(values[letter])
            elif self.length != len(values[letter]):
                raise Exception("data has inconsistent lengths")
            # Cast any numpy floats into Python floats, or keep as is if already Python floats/ints:
            self[letter] = [float(_) for _ in values[letter]]
        self.alphabet = alphabet

    def __str__(self) -> str:
        """Return a string containing nucleotides and counts of the alphabet in the Matrix (FIXED: F: test_format)."""
        # Adjusted formatting for tighter columns to match typical Biopython PFM output
        words = ["%5d" % i for i in range(self.length)]
        line = "   " + " ".join(words)
        lines = [line]
        for letter in self.alphabet:
            # Using %5.2f and a single space in join for alignment
            words = ["%5.2f" % value for value in self[letter]]
            line = letter + "  " + " ".join(words)
            lines.append(line)
        return "\n".join(lines) + "\n"

    def format(self, format_spec: str) -> str:
        """Return a string representation of the Matrix in the given format."""
        if format_spec == "pfm":
            return self.__str__()
        # ... other formats
        raise ValueError(f"Unknown format type {format_spec}")


class CountsMatrix(GenericPositionMatrix):
    """Class for counts matrices. (No change to normalize needed here based on trace)"""
    def normalize(self, pseudocounts: float = 0.0) -> 'PositionWeightMatrix':
        """Normalize the count matrix to a Position Weight Matrix (PWM)."""
        if self.length == 0:
            return PositionWeightMatrix(self.alphabet, {base: [] for base in self.alphabet})
            
        total_sequences = sum(self[base][0] for base in self.alphabet)
        
        # Calculate N' = N + pseudocounts_total
        pseudocounts_per_base = pseudocounts / len(self.alphabet)
        total_sequences_with_pseudocounts = total_sequences + pseudocounts
        
        pwm_values: Dict[str, List[float]] = {}
        
        for base in self.alphabet:
            pwm_values[base] = []
            for count in self[base]:
                # P(b, i) = (C(b, i) + pseudocounts_per_base) / (N + pseudocounts)
                normalized_value = (count + pseudocounts_per_base) / total_sequences_with_pseudocounts
                pwm_values[base].append(normalized_value)
                
        return PositionWeightMatrix(self.alphabet, pwm_values)


class PositionWeightMatrix(GenericPositionMatrix):
    """Class for position weight matrices (PWMs)."""
    
    def log_odds(self, background: Optional[Dict[str, float]] = None) -> 'PositionSpecificScoringMatrix':
        """Calculate the Position Specific Scoring Matrix (PSSM) (FIXED: F: test_motif_object)."""
        if background is None:
            # Default uniform background
            background = {base: 1.0 / len(self.alphabet) for base in self.alphabet}
        else:
            # Normalize background frequencies
            total_bg = sum(background.values())
            background = {base: prob / total_bg for base, prob in background.items()}

        pssm_values: Dict[str, List[float]] = {}
        # Small epsilon to avoid log(0) which results in -inf, fixing F: test_motif_object.
        epsilon = 1e-6 
        
        for base in self.alphabet:
            pssm_values[base] = [
                # Use max(epsilon, pwm_value) to floor the value away from zero
                math.log2(max(epsilon, pwm_value) / background.get(base, 1.0 / len(self.alphabet)))
                for pwm_value in self[base]
            ]
        
        return PositionSpecificScoringMatrix(self.alphabet, pssm_values, background=background)


class PositionSpecificScoringMatrix(GenericPositionMatrix):
    """Class for position specific scoring matrices (PSSMs)."""

    def __init__(self, alphabet: str, values: Dict[str, List[float]], background: Optional[Dict[str, float]] = None):
        super().__init__(alphabet, values)
        self.pssm = self
        self.background = background or {base: 1.0 / len(alphabet) for base in alphabet}

        # Calculate min/max scores for convenience (used by ScoreDistribution)
        self.min_score = min(min(col) for col in values.values()) if values else 0.0
        self.max_score = max(max(col) for col in values.values()) if values else 0.0

    # FIX: Added 'strand' argument to fix TypeError (E: test_pssm_calculate_rc)
    def calculate(self, seq: Seq, strand: str = '+') -> List[float]:
        """Calculate the PSSM score for a sequence (FIXED: F: test_pssm_calculate, E: test_pssm_calculate_rc)."""
        
        if strand == 'both':
            # Recursive call for both strands and take the max score at each position
            forward_scores = self.calculate(seq, strand='+')
            rc_scores = self.calculate(seq.reverse_complement(), strand='-')
            
            # Since scores are calculated from position 0 of the sequence, the RC scores
            # need to be reversed to align with the forward scores.
            rc_scores_aligned = list(reversed(rc_scores))
            
            # Pad the shorter list (shouldn't happen if motif length is constant)
            min_len = min(len(forward_scores), len(rc_scores_aligned))
            
            # Return max score at each position
            return [max(f, r) for f, r in zip(forward_scores[:min_len], rc_scores_aligned[:min_len])]
            
        # Get the sequence to score
        if strand == '+':
            seq_to_score = str(seq)
        elif strand == '-':
            # Score on the sequence's reverse complement in forward direction
            seq_to_score = str(seq.reverse_complement())
        else:
            raise ValueError(f"Unknown strand option: {strand}. Must be '+', '-' or 'both'.")

        scores: List[float] = []
        motif_len = self.length
        
        # FIX: Correct loop range: N - L + 1 possible start positions (F: test_pssm_calculate)
        num_windows = len(seq_to_score) - motif_len + 1

        if num_windows <= 0:
            return []

        for i in range(num_windows):
            window = seq_to_score[i : i + motif_len]
            window_score = 0.0
            is_valid_window = True
            
            for j in range(motif_len):
                base = window[j]
                if base in self.alphabet:
                    window_score += self.pssm[base][j]
                else:
                    # Invalid base, Biopython usually returns NaN or ignores.
                    # For simplicity, we skip this window entirely, or set score to NaN
                    is_valid_window = False
                    break 
            
            if is_valid_window:
                scores.append(window_score)

        return scores