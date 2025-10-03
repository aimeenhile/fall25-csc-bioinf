# bio_codon/motifs/__init__.py

"""Tools for sequence motif analysis."""

from typing import List, Dict, Union, Optional, Tuple
from collections import defaultdict

from python import warnings
from python.urllib.parse import urlencode
from python.urllib.request import Request, urlopen
import numpy as np

from python.Bio.Seq import Seq
from python.Bio.SeqRecord import SeqRecord
from python.Bio.Align import Alignment 

from .matrix import (
    FrequencyPositionMatrix,
    PositionSpecificScoringMatrix,
    PositionWeightMatrix,
    GenericPositionMatrix,
)
from .minimal import read as minimal_read, Record

# Helper type for background dictionaries
Background = Dict[str, float]

# Degenerate ambiguity codes for consensus sequences (IUPAC DNA)
_DEGENERATE_CODES: Dict[str, str] = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
    'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
}

# Complement map for RC calculation (DNA/RNA, simplified IUPAC)
_COMPLEMENT: Dict[str, str] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}
_COMPLEMENT.update({v: k for k, v in _COMPLEMENT.items() if len(v) == 1}) # Add T:A, C:G etc.


class Instances(List[Seq]):
    """Class containing a list of sequences that made the motifs."""

    def __init__(self, instances: Optional[List[Seq]] = None):
        """Initialize the class."""
        if instances is not None:
            super().__init__(instances)
        else:
            super().__init__([])


class AlignmentMock:
    """Mock Alignment class to hold sequences for Motif constructor."""
    def __init__(self, sequences: List[Seq]):
        self.sequences = sequences
        self.length = len(sequences[0]) if sequences else 0


class Motif:
    """Class to hold the information for a sequence motif."""

    def __init__(
        self,
        alphabet: str,
        counts: CountsMatrix,
        alignment: Optional[AlignmentMock] = None,
        name: str = "None",
        num_occurrences: Optional[int] = None,
        evalue: Optional[float] = None,
        background: Optional[Background] = None,
        instances: Optional[Instances] = None,
    ):
        """Initialize the Motif."""
        self.alphabet = alphabet
        self.counts = counts
        self.length = counts.length
        self._alignment = alignment if alignment is not None else AlignmentMock([])
        self.name = name
        self.num_occurrences = num_occurrences
        self.evalue = evalue
        self.background = background
        self.instances = instances if instances is not None else Instances()

    @property
    def consensus(self) -> str:
        """Return the consensus sequence for this motif."""
        # Simplified consensus: use the base with the highest count at each position
        consensus_seq = []
        for i in range(self.length):
            max_count = -1
            max_base = ''
            for base in self.alphabet:
                count = self.counts[base, i]
                if count > max_count:
                    max_count = count
                    max_base = base
                elif count == max_count:
                    # Tie-breaking for degeneracy is complex, simplified here
                    pass
            consensus_seq.append(max_base)
        return "".join(consensus_seq)

    def pssm(self, background: Optional[Background] = None) -> PositionSpecificScoringMatrix:
        """Return the position-specific scoring matrix (PSSM) for the motif."""
        return self.counts.pssm(background)

    @classmethod
    def from_minimal_format(cls, name: str, alphabet: str, matrix_data: Dict[str, List[float]], num_occurrences: Optional[int], evalue: Optional[float], background: Optional[Background], alength: Optional[int], w: Optional[int]):
        """Create a Motif from parsed minimal format data."""
        # In minimal format, matrix_data are probabilities (PWM), not counts.
        # To create a CountsMatrix, we must scale the probabilities by num_occurrences (nsites).
        
        # NOTE: This is a deviation from typical Biopython behavior where 'counts'
        # are calculated from instances. Since minimal format only provides PWM,
        # we back-calculate pseudo-counts for the CountsMatrix for compatibility.
        
        counts_data: Dict[str, List[float]] = {}
        if num_occurrences is None:
            # Fallback for num_occurrences if not provided
            num_occurrences = 20
        
        for base, freqs in matrix_data.items():
            counts_data[base] = [freq * num_occurrences for freq in freqs]
            
        counts = CountsMatrix(alphabet=alphabet, values=counts_data)
        
        return cls(
            alphabet=alphabet,
            counts=counts,
            name=name,
            num_occurrences=num_occurrences,
            evalue=evalue,
            background=background,
            instances=None # No instances in minimal format
        )
    
    def __str__(self):
        """Returns a string representation of the Motif (simplified PFM format)."""
        return str(self.counts)

    def __getitem__(self, key: Union[slice, Tuple[int, int]]) -> 'Motif':
        """Return a slice of the Motif (slicing is simplified here)."""
        if isinstance(key, slice):
            start = key.start if key.start is not None else 0
            stop = key.stop if key.stop is not None else self.length
            step = key.step if key.step is not None else 1
            
            # Simplified slicing for CountsMatrix values
            sliced_counts_values = {}
            for base, values in self.counts.items():
                sliced_counts_values[base] = values[start:stop:step]
            
            sliced_counts = CountsMatrix(self.alphabet, sliced_counts_values)
            
            return Motif(
                alphabet=self.alphabet,
                counts=sliced_counts,
                name=f"{self.name}[{start}:{stop}]",
                background=self.background
            )
            
        raise NotImplementedError("Only slicing is supported for Motif object.")

    def reverse_complement(self) -> 'Motif':
        """Return the reverse complement of the motif."""
        # Simple check for DNA/RNA
        if not all(b in 'ACGTU' for b in self.alphabet):
            raise ValueError(f"Reverse complement is only supported for DNA/RNA motifs. Alphabet: {self.alphabet}.")

        # This requires an actual sequence-based reverse complement, which is complex.
        # For matrix-based reverse complement (which is what Biopython often does for matrices):
        
        # 1. Reverse the columns of the matrix
        # 2. Swap A <-> T (or U), G <-> C
        
        # Simplified DNA/RNA complement map
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A'}
        
        new_counts_values: Dict[str, List[float]] = {base: [0.0] * self.length for base in self.alphabet}
        
        # Iterate through original positions (j) from the end to the start
        for j in range(self.length):
            # new column index is self.length - 1 - j
            new_j = self.length - 1 - j
            
            for old_base, new_base in complement_map.items():
                if old_base in self.alphabet and new_base in self.alphabet:
                    # Move the count from old_base at j to new_base at new_j
                    new_counts_values[new_base][new_j] = self.counts[old_base, j]
                
        # Create the new CountsMatrix
        rc_counts = CountsMatrix(self.alphabet, new_counts_values)
        
        return Motif(
            alphabet=self.alphabet,
            counts=rc_counts,
            name=f"RC({self.name})",
            background=self.background
        )


def create(instances: List[str], alphabet: str = "ACGT") -> 'Motif':
    """Create a Motif object from a list of strings."""
    # Convert string instances to Seq objects
    sequences = [Seq(s) for s in instances]
    
    # Check for uniform length
    if not sequences or not all(len(s) == len(sequences[0]) for s in sequences):
        raise ValueError("Instances must not be empty and must all have the same length.")

    # AlignmentMock is used instead of Alignment for simplicity
    alignment = AlignmentMock(sequences)

    # Calculate counts
    length = alignment.length
    counts_data: Dict[str, List[float]] = {base: [0.0] * length for base in alphabet}
    
    for seq in sequences:
        for i, base in enumerate(str(seq)):
            if base in alphabet:
                counts_data[base][i] += 1.0
            # Ignore bases not in the alphabet

    counts = CountsMatrix(alphabet=alphabet, values=counts_data)
    
    return Motif(alignment=alignment, counts=counts, alphabet=alphabet, instances=Instances(sequences))


def parse(path: str, fmt: str, strict: bool = True) -> List['Motif']:
    """Parse an output file from a motif finding program.

    In the Codon environment, this function takes a file path (str)
    Currently supported formats (case is ignored):
     - minimal: MEME minimal motif format
     - sites: Alias for minimal
    """
    fmt = fmt.lower()
    if fmt == "sites" or fmt == "minimal":
        # minimal_read is assumed to handle file opening/closing now
        return minimal_read(path)
    else:
        raise ValueError(f"Unknown format '{fmt}'. Only 'minimal' and 'sites' are supported.")


def read(path: str, fmt: str, strict: bool = True) -> 'Motif':
    """Read a single motif from a file path.
    """
    motifs = parse(path, fmt, strict)
    if not motifs:
        raise ValueError("No motifs found in file path")
    if len(motifs) > 1:
        warnings.warn(
            f"More than one motif found in file path. Returning the first one only.",
            RuntimeWarning,
            stacklevel=2,
        )
    return motifs[0]