# bio_codon/motifs/__init__.py

"""Tools for sequence motif analysis.

from typing import Any, Optional, List, Dict, Tuple, Union

from python import warnings
from python.urllib.parse import urlencode
from python.urllib.request import Request, urlopen
import numpy as np

from python.Bio.Align import Alignment
from python.Bio.Seq import Seq

from .matrix import CountsMatrix, PositionWeightMatrix, PositionSpecificScoringMatrix
from . import minimal
from . import thresholds
"""

from typing import list, dict, Optional, Any
import math
from .matrix import (
    CountMatrix,
    FrequencyPositionMatrix,
    PositionWeightMatrix,
    PositionSpecificScoringMatrix,
    DNA_ALPHABET
)
from .thresholds import ScoreDistribution
from . import minimal 


class Instances:
    """Class containing a list of sequences that made the motifs."""
    sequences: list[str]
    length: int

    def __init__(self, sequences: list[str]):
        if not sequences:
            self.sequences = []
            self.length = 0
        else:
            self.sequences = sequences
            self.length = len(sequences[0])

    def reverse_complement(self) -> 'Instances':
        """Returns the reverse complement of all contained instances."""
        # Simple mapping for DNA bases and unknown 'N'
        rev_comp_map: dict[str, str] = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        rev_comp_sequences: list[str] = []
        
        for seq in self.sequences:
            rev_seq: list[str] = []
            # Manually reverse complement the sequence
            for i in range(len(seq) - 1, -1, -1):
                base_char: str = seq[i]
                rev_seq.append(rev_comp_map.get(str(base_char).upper(), 'N'))
            rev_comp_sequences.append("".join(rev_seq))
            
        return Instances(rev_comp_sequences)
        
    def __len__(self) -> int:
        return len(self.sequences)
        
    def __str__(self) -> str:
        return f"<Instances: {len(self.sequences)} sequences of length {self.length}>"


class Motif:
    """A class representing sequence motifs."""

    name: str
    instances: Instances
    counts: CountMatrix
    length: int
    
    # Configuration properties with internal caches/defaults
    _pseudocounts: float = 0.5
    _background: dict[str, float] = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    _freq_matrix: Optional[FrequencyPositionMatrix] = None
    _pwm: Optional[PositionWeightMatrix] = None

    def __init__(self, sequences: list[str], name: str = "Unnamed Motif"):
        self.name = name
        self.instances = Instances(sequences)
        self.counts = CountMatrix(sequences)
        self.length = self.counts.length
    
    # Getter/Setter for pseudocounts, clearing cache on change
    @property
    def pseudocounts(self) -> float:
        return self._pseudocounts

    @pseudocounts.setter
    def pseudocounts(self, value: float):
        self._pseudocounts = value
        # Clear cache when settings change
        self._freq_matrix = None
        self._pwm = None

    # Getter/Setter for background, clearing cache on change
    @property
    def background(self) -> dict[str, float]:
        return self._background
        
    @background.setter
    def background(self, value: dict[str, float]):
        self._background = value
        # Clear cache when settings change
        self._pwm = None

    def _get_freq_matrix(self) -> FrequencyPositionMatrix:
        """Computes and caches the frequency matrix using current pseudocounts."""
        if self._freq_matrix is None:
            self._freq_matrix = FrequencyPositionMatrix(
                self.counts, 
                pseudocount=self._pseudocounts
            )
        return self._freq_matrix

    @property
    def consensus(self) -> str:
        """Returns the consensus sequence based on frequency."""
        return self._get_freq_matrix().consensus()

    @property
    def pwm(self) -> PositionWeightMatrix:
        """Returns the PositionWeightMatrix (log-odds scores)."""
        if self._pwm is None:
            freq_mat = self._get_freq_matrix()
            self._pwm = PositionWeightMatrix(freq_mat, self._background)
        return self._pwm

    @property
    def pssm(self) -> PositionSpecificScoringMatrix:
        """Returns the PSSM (an alias for PWM)."""
        return self.pwm

    def reverse_complement(self) -> 'Motif':
        """Returns a new Motif instance that is the reverse complement."""
        rc_instances = self.instances.reverse_complement()
        rc_motif = Motif(rc_instances.sequences, name=f"RC of {self.name}")
        # Copy settings to the new motif
        rc_motif.pseudocounts = self._pseudocounts
        rc_motif.background = self._background
        return rc_motif
        
    def __len__(self) -> int:
        return self.length

    def __str__(self) -> str:
        return f"{self.name} (L={self.length}, N={len(self.instances)}) Consensus: {self.consensus}"
        
    def __format__(self, format_spec: str) -> str:
        """Allows formatting using f-string syntax (e.g., f'{motif:minimal}')."""
        if format_spec == 'minimal':
            return write(self, format='minimal')
        return str(self)

    
# --- Public I/O Functions ---

def create(instances: list[str], name: str = "Motif from Instances") -> Motif:
    """Factory function to create a Motif from instances."""
    return Motif(instances, name=name)

def read(data: list[str], format: str = 'minimal') -> list[Motif]:
    """
    Reads motifs from a list of strings (lines) in the specified format.
    Only 'minimal' format is supported in this Codon port.
    """
    if format == 'minimal':
        # Use the minimal I/O module
        record = minimal.read(data)
        return record.motifs
    else:
        raise ValueError(f"Unknown or unsupported format '{format}'. Only 'minimal' is implemented.")

def parse(data: str, format: str = 'minimal') -> list[Motif]:
    """Parses a string containing motif data by splitting it into lines."""
    lines: list[str] = data.split('\n')
    return read(lines, format)

def write(motif: Motif, format: str = 'minimal') -> str:
    """Converts a Motif object into a string representation in the given format."""
    if format == 'minimal':
        # Simple minimal format: Name + Consensus
        return f">{motif.name}\n{motif.consensus}"
    else:
        raise ValueError(f"Unknown or unsupported format '{format}'. Only 'minimal' is implemented.")

# --- Example Usage ---

def main():
    print("--- Codon Bio.motifs Port Demonstration ---")
    
    # Sequences that define a TATA-like box (consensus TACAAT)
    sequences: list[str] = [
        "TACAAT",
        "TACAAG",
        "CACAAT",
        "GAGAAG",
        "AACAAG",
        "TACAAC",
        "TACAAA",
    ]

    # 1. Create and Display Motif
    my_motif = create(sequences, name="TATA-like Box")
    print("\n--- 1. Full Motif Data ---")
    print(str(my_motif.counts)) 
    print(str(my_motif.pwm)) 
    print(f"Consensus: {my_motif.consensus}")
    
    # 2. Score and Threshold
    test_pwm = my_motif.pwm
    good_sequence: str = "TACAAT"
    poor_sequence: str = "GGGGGG"
    
    good_score: float = test_pwm.calculate_score(good_sequence)
    poor_score: float = test_pwm.calculate_score(poor_sequence)
    
    dist = ScoreDistribution(test_pwm)
    p_value: float = 0.001
    threshold: float = dist.threshold_for_p_value(p_value)
    
    print("\n--- 2. Scoring and Threshold ---")
    print(f"Score for consensus '{good_sequence}': {good_score:.4f}")
    print(f"Score for non-site '{poor_sequence}': {poor_score:.4f}")
    print(f"Approximate score threshold for P={p_value}: {threshold:.4f}")
    
    # 3. I/O Test
    minimal_output = write(my_motif, format='minimal')
    print("\n--- 3. Minimal Format Output ---")
    print(minimal_output)
    
    # Read back (simulating a file read)
    lines_to_read: list[str] = [
        ">MY_NEW_MOTIF_1",
        "ATGC",
        "ATTT",
        ">MY_NEW_MOTIF_2",
        "GGCA",
        "GGGG"
    ]
    read_motifs: list[Motif] = read(lines_to_read, format='minimal')
    print("\n--- 4. Reading Minimal Format (2 motifs found) ---")
    for m in read_motifs:
        print(f"Read Motif: {m.name}, Consensus: {m.consensus}")

# The entry point of the Codon program
if __name__ == "__main__":
    main()