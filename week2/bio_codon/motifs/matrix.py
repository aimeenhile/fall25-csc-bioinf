# bio_codon/motifs/matrix.py

"""Support for various forms of sequence motif matrices.

Implementation of frequency (count) matrices, position-weight matrices,
and position-specific scoring matrices.
"""
import math
from python import numbers
import numpy as np
from python import Bio.Seq.Seq
from typing import Dict, Tuple, List, Optional, Union

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
    ('C', 'G'): 'S',      # C or G
    ('A', 'G'): 'R',      # A or G
    ('C', 'T'): 'Y',      # C or T
    ('G', 'T'): 'K',      # G or T
    ('A', 'C'): 'M',      # A or C
}

DEGENERATE_MAP_3: Dict[Tuple[str, str, str], str] = {
    ('A', 'C', 'G'): 'V', # A, C, or G
    ('A', 'C', 'T'): 'H', # A, C, or T
    ('A', 'G', 'T'): 'D', # A, G, or T
    ('C', 'G', 'T'): 'B', # C, G, or T
}

DEGENERATE_MAP_4: Dict[Tuple[str, str, str, str], str] = {
    ('A', 'C', 'G', 'T'): 'N', # Any base
}


def _get_background_prob(background: Optional[Dict[str, float]], base: str) -> float:
    """Safely get background probability for a base."""
    if background is None:
        return 0.25 # Default uniform background
    return background.get(base, 0.0)

class GenericPositionMatrix(dict):
    """Base class for the support of position matrix operations."""

    length: int
    alphabet: str

    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the class."""
        self.length = None
        for letter in alphabet:
            if self.length is None:
                if isinstance(values[letter], np.ndarray):
                    self.length = values[letter].size
                else:
                    self.length = len(values[letter])
            elif self.length != (values[letter].size if isinstance(values[letter], np.ndarray) else len(values[letter])):
                raise Exception("data has inconsistent lengths")
            # Cast any numpy floats into Python floats/standard list:
            self[letter] = [float(_) for _ in values[letter]]
        self.alphabet = alphabet

    def __str__(self) -> str:
        """Return a string containing nucleotides and counts of the alphabet in the Matrix."""
        words = [f"{i:6d}" for i in range(self.length)]
        line = "   " + " ".join(words)
        lines = [line]
        for letter in self.alphabet:
            words = [f"{value:6.2f}" for value in self[letter]]
            lines.append(f"{letter}: " + " ".join(words))
        return "\n".join(lines)

    # FIX: Simplify the Union hint for key to avoid the 'union already sealed' compiler crash.
    def __getitem__(self, key: Union[str, Tuple[object, int]]):
        """Allow access like m.counts['A', 3] or m.counts[:, 3].
        
        The Union hint on the key argument handles str for row access (m['A'])
        and Tuple for cell/column access (m['A', 3] or m[:, 3]).
        """
        if isinstance(key, str):
            # Access like m['A'] - returns the list of values for base 'A'
            return super().__getitem__(key)
        
        if isinstance(key, tuple):
            if len(key) == 2:
                base_or_index, pos = key
                
                # m.counts[base, position]
                if isinstance(base_or_index, str):
                    return self[base_or_index][pos]
                
                # m.counts[base_index, position] - simplified access
                elif isinstance(base_or_index, int):
                    base = self.alphabet[base_or_index]
                    return self[base][pos]
                
                # m.counts[:, position] - access to a column (requires pos to be an int)
                elif base_or_index == slice(None):
                    if isinstance(pos, int):
                        # Returns column as a dictionary
                        return {base: self[base][pos] for base in self.alphabet}
                    
        raise KeyError("Invalid key for PositionMatrix access.")


    def _calculate_consensus(self) -> str:
        """Calculate the consensus sequence (PRIVATE)."""
        # NOTE: self._consensus is not defined in __init__, so we'll treat this as
        # a standard property calculation. For proper state, it should be initialized.
        # Assuming the caller/Motif handles caching if needed.
        
        consensus_seq: List[str] = []
        for i in range(self.length or 0):
            column = self[:, i]
            # Find the base with the maximum value
            max_value: float = -1.0
            max_base: str = ''
            for base, value in column.items():
                if value > max_value:
                    max_value = value
                    max_base = base
                elif value == max_value:
                    # Tie-breaking: Biopython usually just takes the first one found
                    pass
            consensus_seq.append(max_base)
        
        # NOTE: Skipping internal attribute caching (self._consensus = ...) for simplicity
        # and to match the limited snippet context.
        return "".join(consensus_seq)

    @property
    def consensus(self) -> str:
        """Return the consensus sequence."""
        return self._calculate_consensus()

    def _calculate_degenerate_consensus(self) -> str:
        """Calculate the degenerate consensus sequence (PRIVATE)."""
        # NOTE: self._degenerate_consensus is not defined in __init__.
        
        degenerate_consensus_seq: List[str] = []
        # Calculate the threshold for degenerate consensus based on max value
        # Simplified: bases with value > (max_value / 2) are included in degeneracy
        for i in range(self.length or 0):
            column = self[:, i]
            max_value: float = max(column.values())
            
            # Bases that contribute significantly to the position (e.g., > 50% of max value)
            # The list of bases MUST be sorted to match the canonical order of keys in the maps
            significant_bases: List[str] = sorted([
                base for base, value in column.items()
                if value >= max_value * 0.5
            ])
            
            # Convert the list of bases to an IUPAC degenerate code
            base_tuple = tuple(significant_bases)
            
            if len(base_tuple) == 1:
                degenerate_consensus_seq.append(base_tuple[0])
            else:
                code: Optional[str] = None
                
                # Search the correct map based on the tuple length
                # We use specific map types (DEGENERATE_MAP_2, _3, _4) to avoid the previous compiler error.
                if len(base_tuple) == 2:
                    code = DEGENERATE_MAP_2.get(base_tuple) # type: ignore
                elif len(base_tuple) == 3:
                    code = DEGENERATE_MAP_3.get(base_tuple) # type: ignore
                elif len(base_tuple) == 4:
                    code = DEGENERATE_MAP_4.get(base_tuple) # type: ignore
                
                # Append the degenerate code or fall back to 'N'
                if code:
                    degenerate_consensus_seq.append(code)
                else:
                    degenerate_consensus_seq.append('N')
        
        # NOTE: Skipping internal attribute caching (self._degenerate_consensus = ...)
        return "".join(degenerate_consensus_seq)

    @property
    def degenerate_consensus(self) -> str:
        """Return the degenerate consensus sequence."""
        return self._calculate_degenerate_consensus()

    def reverse_complement(self, alphabet: str = "ACGT") -> 'GenericPositionMatrix':
        """Return the reverse complement of the motif."""
        if self.alphabet != alphabet:
            raise ValueError("Alphabet mismatch for reverse complement calculation.")
        
        if alphabet == "ACGT":
            # DNA alphabet complement mapping
            complement_map: Dict[str, str] = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        elif alphabet == "ACGU":
            # RNA alphabet complement mapping
            complement_map = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
        else:
            raise ValueError("Unsupported alphabet for reverse complement.")

        rev_comp_values: Dict[str, List[float]] = {}
        length = self.length or 0
        
        for base in alphabet:
            complement_base = complement_map.get(base)
            if complement_base is None:
                raise ValueError(f"Base {base} not supported in complement mapping.")
            
            # Values are reversed and then mapped to the complement base
            rev_comp_values[base] = [
                self[complement_base][length - 1 - i]
                for i in range(length)
            ]
            
        return GenericPositionMatrix(alphabet, rev_comp_values)


class CountsMatrix(GenericPositionMatrix):
    """Counts Matrix."""

    def normalize(self, pseudocounts: float = 0.0) -> 'PositionWeightMatrix':
        """Normalize the matrix to get a PositionWeightMatrix."""
        pwm_values: Dict[str, List[float]] = {}
        
        # Calculate column sums (N) for normalization
        column_sums = [sum(self[base][j] for base in self.alphabet) for j in range(self.length)]
        
        for base in self.alphabet:
            counts = self[base]
            freqs = []
            for j in range(self.length):
                N = column_sums[j]
                # Formula: (count + pseudocount) / (N + len(alphabet) * pseudocount)
                denominator = N + len(self.alphabet) * pseudocount
                if denominator == 0:
                    # Avoid division by zero if N is 0 and pseudocount is 0
                    freq = 0.0
                else:
                    freq = (counts[j] + pseudocount) / denominator
                freqs.append(freq)
            pwm_values[base] = freqs
            
        return PositionWeightMatrix(self.alphabet, pwm_values)

    def pssm(self, background: Optional[Dict[str, float]] = None) -> 'PositionSpecificScoringMatrix':
        """Calculate the PSSM (Requires normalization first)."""
        # This implementation simply wraps normalize then log_odds.
        return self.normalize(pseudocounts=0.0).log_odds(background)

class PositionWeightMatrix(GenericPositionMatrix):
    """Position Weight Matrix (PWM)."""

    def log_odds(self, background: Optional[Dict[str, float]] = None) -> 'PositionSpecificScoringMatrix':
        """Calculate the PSSM (Position-Specific Scoring Matrix) from the PWM."""
        pssm_values: Dict[str, List[float]] = {}
        
        for base in self.alphabet:
            probs = self[base]
            bg_prob = _get_background_prob(background, base)
            
            if bg_prob == 0.0:
                 # Default to log(P_ij / 0.25) if no background is provided,
                 # or if background is provided but the base is not in it.
                 bg_prob = 0.25 

            scores = []
            for prob in probs:
                if prob == 0.0:
                    # Log(0) is negative infinity. Use a small epsilon to avoid math domain error.
                    score = math.log2(1e-10 / bg_prob) # Using a tiny probability
                else:
                    score = math.log2(prob / bg_prob)
                scores.append(score)
            pssm_values[base] = scores
            
        return PositionSpecificScoringMatrix(self.alphabet, pssm_values)


class PositionSpecificScoringMatrix(GenericPositionMatrix):
    """Position Specific Scoring Matrix (PSSM)."""
    
    # Pre-calculate min/max score for convenience
    min: float
    max: float
    
    def __init__(self, alphabet: str, values: Dict[str, List[float]]):
        """Initialize the class."""
        super().__init__(alphabet, values)
        
        self.min = 0.0
        self.max = 0.0
        
        if self.length > 0:
            all_scores = [self[base][j] for base in self.alphabet for j in range(self.length)]
            self.min = min(all_scores)
            self.max = max(all_scores)

    def calculate(self, sequence: Seq) -> List[float]:
        """Calculate the PSSM score for all positions in the given sequence."""
        sequence_str = str(sequence).upper()
        if len(sequence_str) < self.length:
            return []

        scores: List[float] = []
        
        # Iterate over all possible starting positions (i)
        for i in range(len(sequence_str) - self.length + 1):
            motif_score = 0.0
            is_valid = True
            
            # Iterate over each position (j) in the motif
            for j in range(self.length):
                base = sequence_str[i + j]
                
                if base in self.alphabet:
                    # Add the score for the base at this position
                    motif_score += self[base][j]
                else:
                    # Handle unknown/degenerate characters by marking the score as NaN
                    motif_score = float('nan')
                    is_valid = False
                    break
            
            scores.append(motif_score)
            
        return scores

    # Omitted complex methods like "search" and "distribution" for brevity and core focus.
    
    def distribution(self, background: Optional[Dict[str, float]] = None, precision: int = 10**3):
        """Calculate the distribution of the scores at the given precision."""
        # Use a placeholder for distribution logic to avoid importing within the class body
        # and to keep the structure. The real implementation is in thresholds.py
        from .thresholds import ScoreDistribution
        return ScoreDistribution(pssm=self, background=background, precision=precision)