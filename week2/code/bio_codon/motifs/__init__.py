# bio_codon/motifs/__init__.py

"""Tools for sequence motif analysis.
"""

from typing import Optional, List, Dict, Tuple, Union, TextIO
from collections import defaultdict
import warnings
from python.urllib.parse import urlencode
from python.urllib.request import Request, urlopen
import numpy as np

from python.Bio.Align import Alignment
from python.Bio.Seq import Seq

from .matrix import CountsMatrix, PositionWeightMatrix, PositionSpecificScoringMatrix
from . import minimal
from . import thresholds


from .matrix import (
    FrequencyPositionMatrix,
    PositionSpecificScoringMatrix,
    PositionWeightMatrix,
    GenericPositionMatrix,
)
from .minimal import read as minimal_read, Record
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Helper type for background dictionaries
Background = Dict[str, float]

# Degenerate ambiguity codes for consensus sequences (IUPAC DNA)
_DEGENERATE_CODES: Dict[str, str] = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
    'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
}

class Instances(List[Seq]):
    """Class containing a list of sequences that made the motifs."""

    def __init__(self, instances: Optional[List[Seq]] = None):
        """Initialize the class."""
        super().__init__(instances or [])
        
    def __len__(self) -> int:
        return len(self.sequences)
        
    def __str__(self) -> str:
        return f"<Instances: {len(self.sequences)} sequences of length {self.length}>"


class AlignmentMock:
    """Mock class to replicate necessary attributes of Bio.Align.Alignment for Motif creation."""
    def __init__(self, sequences: List[Seq]):
        self.sequences = sequences
        self.length = len(sequences[0]) if sequences else 0
        self.column_annotations = {} # Simplification

    def get_alignment_length(self) -> int:
        return self.length


class Motif:
    """Class for sequence motif analysis."""

    alphabet: str
    sequence: list[str]
    alignment: AlignmentMock

    def __init__(self, alignment: AlignmentMock, alphabet: str):
        """Initialize the Motif object."""
        self._alignment = alignment
        self._alphabet = alphabet

        # 1. CountsMatrix
        self.counts = self._calculate_counts()

        # 2. PositionWeightMatrix (PWM) - Normalized Counts
        # Default to no pseudocounts for basic PWM
        self.pwm: PositionWeightMatrix = self.counts.normalize(pseudocounts=0.0)

        # 3. PositionSpecificScoringMatrix (PSSM) - Log odds
        # Default to uniform background
        background = dict.fromkeys(self._alphabet, 1.0 / len(self._alphabet))
        self.pssm: PositionSpecificScoringMatrix = self.pwm.log_odds(background)

        # Other properties
        self.name: Optional[str] = None
        self.altname: Optional[str] = None
        self.evalue: Optional[float] = None
        self.num_occurrences: Optional[int] = None
        self.instances = Instances(alignment.sequences)

    def _calculate_counts(self) -> CountsMatrix:
        """Calculate the CountsMatrix from the alignment."""
        sequences = self._alignment.sequences
        length = self._alignment.get_alignment_length()
        counts: Dict[str, List[float]] = {base: [0.0] * length for base in self._alphabet}

        for i in range(length):
            column_bases = [str(seq)[i] for seq in sequences]
            for base in column_bases:
                if base in self._alphabet:
                    counts[base][i] += 1.0

        return CountsMatrix(self._alphabet, counts)

    @property
    def consensus(self) -> str:
        """Returns the consensus sequence for this motif."""
        consensus_seq: List[str] = []
        for i in range(self.counts.length):
            column = {base: self.counts[base][i] for base in self._alphabet}
            max_count = -1.0
            consensus_base: Optional[str] = None
            for base, count in column.items():
                if count > max_count:
                    max_count = count
                    consensus_base = base
                elif count == max_count and consensus_base:
                    # Tie-breaking rule: use the base that comes first in the alphabet
                    if self._alphabet.index(base) < self._alphabet.index(consensus_base):
                         consensus_base = base

            if consensus_base is not None:
                consensus_seq.append(consensus_base)
            else:
                consensus_seq.append('N') # Should not happen with valid alphabet

        return "".join(consensus_seq)

    def __len__(self) -> int:
        """Return the length of the Motif."""
        return self.counts.length

    def __getitem__(self, key: Union[int, slice]) -> 'Motif':
        """Return a slice of the Motif."""
        if isinstance(key, slice):
            start = key.start if key.start is not None else 0
            stop = key.stop if key.stop is not None else len(self)
            step = key.step if key.step is not None else 1

            new_sequences: List[Seq] = []
            for seq in self._alignment.sequences:
                new_sequences.append(seq[start:stop:step])

            rc_alignment = AlignmentMock(new_sequences)
            return Motif(rc_alignment, self._alphabet)

        raise TypeError("Motif slicing only supports slice objects.")

    def reverse_complement(self) -> 'Motif':
        """Return the reverse complement of the Motif."""
        if 'T' in self._alphabet and 'A' in self._alphabet and 'G' in self._alphabet and 'C' in self._alphabet:
            # Assuming standard DNA/RNA alphabet for RC
            pass
        else:
            warnings.warn(f"Reverse complement not defined for non-ACGT alphabet ({self._alphabet}).")

        rc_sequences = [s.reverse_complement() for s in self._alignment.sequences]
        rc_alignment = AlignmentMock(rc_sequences)
        return Motif(rc_alignment, self._alphabet)

    
# --- Public I/O Functions ---

def create(instances: List[str], alphabet: str = "ACGT") -> 'Motif':
    """Create a Motif object from a list of strings."""
    # Convert string instances to Seq objects
    sequences = [Seq(s) for s in instances]

    # Check for uniform length
    if not sequences or not all(len(s) == len(sequences[0]) for s in sequences):
        raise ValueError("Instances must not be empty and must all have the same length.")

    alignment = AlignmentMock(sequences) # Use AlignmentMock
    return Motif(alignment=alignment, alphabet=alphabet)


def parse(handle: TextIO, fmt: str, strict: bool = True) -> List['Motif']:
    """Parse an output file from a motif finding program."""
    fmt = fmt.lower()
    if fmt == "sites":
        # Simplified 'sites' format for testing, assumes minimal style
        return minimal.read(handle)
    if fmt == "minimal":
        return minimal.read(handle)
    else:
        raise ValueError(f"Unknown format '{fmt}'. Only 'minimal' and 'sites' are supported.")

def read(handle: TextIO, fmt: str, strict: bool = True) -> 'Motif':
    """Read a single motif from a handle."""
    motifs = parse(handle, fmt, strict)
    if not motifs:
        raise ValueError("No motifs found in handle")
    if len(motifs) > 1:
        warnings.warn("Found more than one motif, returning the first.")
    return motifs[0]