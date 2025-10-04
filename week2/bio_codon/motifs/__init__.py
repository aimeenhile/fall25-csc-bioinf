"""Tools for sequence motif analysis."""

from typing import List, Dict, Optional, Tuple, Union
from collections import defaultdict

from python import warnings 
import numpy as np

from python import Bio.Seq.Seq as Seq
from python import Bio.SeqRecord.SeqRecord as SeqRec
from python import Bio.Align.Alignment as Alignment

from .matrix import CountsMatrix
# from .minimal import read as minimal_read, Record

# Degenerate ambiguity codes for consensus sequences (IUPAC DNA)
_DEGENERATE_CODES: Dict[str, str] = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
    'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
}

# Complement map for RC calculation (DNA/RNA, simplified IUPAC)
_COMPLEMENT: Dict[str, str] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}

class Motif:
    pass

class Record:
    """
    A container for motifs read from a single file, replacing Bio.motifs.Record.
    It holds global file information (version, alphabet, background) and a list of motifs.
    """
    version: str
    alphabet: str
    background: Dict[str, float]
    motif: List['Motif']

    def __init__(self):
        self.version = ""
        self.alphabet = ""
        self.background: Dict[str, float] = {}
        self.motifs: List['Motif'] = []

    def __len__(self) -> int:
        return len(self.motifs)

    def __iter__(self):
        return iter(self.motifs)

    def __getitem__(self, index: int) -> 'Motif':
        return self.motifs[index]


class Instances:
    """
    A container for the sequences used to create the motif.
    """

    sequence: List['Seq']
    length: int
    alphabet: str

    def __init__(self, sequences: List['Seq']):
        self.sequences = sequences
        self.length = len(sequences[0]) if sequences else 0
        self.alphabet = sequences[0].alphabet if sequences and hasattr(sequences[0], 'alphabet') else 'ACGT' # Default
        
    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, index: int) -> Seq:
        return self.sequences[index]

    def __iter__(self):
        return iter(self.sequences)

    def __str__(self) -> str:
        return f"Instances of {self.length} sequences"


class Motif:
    """
    Represents a sequence motif.
    
    This class is the core representation, holding the counts matrix,
    alphabet information, and providing properties for consensus,
    reverse complement, etc.
    """

    # alignment: Optional[Alignment]
    # counts: Optinal[CountsMatrix]
    alphabet: str
    # instance: Optional[Instances]
    name: str 
    length: int

    def __init__(self, alignment: Optional['Alignment'], counts: optinal[CountsMatrix], alphabet: str, instances: Optional['Instances'], name: str, length: int, pwm_data: Optional[Dict[str, List[float]]] = None):
        self.alignment = alignment
        self.counts = counts
        self.alphabet = alphabet
        self.instances = instances
        self.name = name
        self.length = length
        self.pwm_data = pwm_data

        # If created from minimal format (pwm_data provided but counts is None)
        if self.counts is None and self.pwm_data is not None:
            self.counts = self._pwm_to_counts(self.pwm_data)
        
    def _pwm_to_counts(self, pwm_data: Dict[str, List[float]]) -> CountsMatrix:
        """Converts raw PWM data (floats/probabilities) into a CountsMatrix.
        
        placeholder to allow the CountsMatrix to hold the data 
        for compatibility with other Motif methods, treating probabilities as counts.
        """
        # Ensure the values are floats, not just ints
        float_values = {base: [float(v) for v in vals] for base, vals in pwm_data.items()}
        return CountsMatrix(alphabet=self.alphabet, values=float_values, length=self.length)

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, index: slice) -> 'Motif':
        """Allows slicing the motif to create a sub-motif."""
        if self.counts is None:
             raise ValueError("Motif has no counts matrix to slice.")

        sliced_counts_data = {
            base: self.counts[base][index] for base in self.counts.alphabet
        }
        
        # Determine the new length
        start, stop, step = index.indices(self.length)
        new_length = len(range(start, stop, step))

        sliced_counts = CountsMatrix(alphabet=self.alphabet, values=sliced_counts_data, length=new_length)
        
        return Motif(alignment=None, counts=sliced_counts, alphabet=self.alphabet, 
                     instances=None, name=f"{self.name}_slice", length=new_length)
        
        return Motif(alignment=None, counts=sliced_counts, alphabet=self.alphabet, instances=Instances([]))

    @property
    def consensus(self) -> str:
        """Returns the consensus sequence based on the majority base at each position."""
        consensus_seq: List[str] = []

        matrix = self.counts

        for j in range(self.length):
            max_value = -1.0
            consensus_base = 'N'
            for base in self.alphabet:
                value = matrix[base][j]
                if value > max_value:
                    max_value = value
                    consensus_base = base
                elif value == max_value and max_value > 0:
                    pass 
            consensus_seq.append(consensus_base)
        return "".join(consensus_seq)

    def reverse_complement(self) -> 'Motif':
        """Returns the reverse complement of the motif."""
        if self.counts is None:
             raise ValueError("Motif has no counts matrix to reverse complement.")

        new_counts_values: Dict[str, List[float]] = defaultdict(lambda: [0.0] * self.length)

        for new_pos in range(self.length):
            old_pos = self.length - 1 - new_pos
            for old_base in self.alphabet:
                new_base = _COMPLEMENT.get(old_base, old_base)
                
                new_counts_values[new_base][new_pos] = self.counts[old_base][old_pos]
                
        rc_counts = CountsMatrix(alphabet=self.alphabet, values=new_counts_values, length=self.length)

        return Motif(alignment=None, counts=rc_counts, alphabet=self.alphabet, 
                instances=None, name=self.name + "_RC", length=self.length)


def create(alignment: List[Seq], alphabet: Optional[str] = None) -> 'Motif':
    """
    Create a Motif object from a list of sequences.
    """

    if not alignment:
        raise ValueError("Cannot create motif from empty sequence list.")

    # Determine alphabet if not provided
    if alphabet is None:
        alphabet = alignment[0].alphabet.letters if hasattr(alignment[0].alphabet, 'letters') else 'ACGT'
        # Handle RNA U/T distinction
        if 'U' in str(alignment[0]).upper():
            alphabet = alphabet.replace('T', 'U')

    length = len(alignment[0])
    
    # Initialize counts matrix data
    counts_data: Dict[str, List[float]] = {base: [0.0] * length for base in alphabet}

    for seq in alignment:
        if len(seq) != length:
            raise ValueError("All sequences must have the same length.")
        
        for i, base in enumerate(str(seq).upper()):
            if base in alphabet:
                counts_data[base][i] += 1.0

    counts = CountsMatrix(alphabet=alphabet, values=counts_data, length=length)
    
    return Motif(alignment=None, counts=counts, alphabet=alphabet, instances=Instances(alignment))


# Import after classes are defined
from .minimal import read as minimal_read


def parse(path: str, fmt: str, strict: bool = True) -> List['Motif']:
    """Parse an output file from a motif finding program.

    Currently supported formats (case is ignored):
      - minimal: MEME minimal motif format
      - sites: Alias for minimal
    """
    fmt = fmt.lower()
    if fmt == "sites" or fmt == "minimal":
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
            stacklevel=2 
        )
        return motifs[0] 
    return motifs[0]
