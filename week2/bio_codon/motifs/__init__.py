"""Tools for sequence motif analysis."""

# FIX: Removed Union from imports and simplified type hints in functions.
from typing import List, Dict, Optional, Tuple
from collections import defaultdict

from python import warnings 
import numpy as np

from python import Bio.Seq.Seq
from python import Bio.SeqRecord.SeqRecord
from python import Bio.Align.Alignment 

from .matrix import (
    FrequencyPositionMatrix,
    PositionSpecificScoringMatrix,
    PositionWeightMatrix,
    CountsMatrix,
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


class Instances:
    """
    A container for the sequences used to create the motif.
    
    This replaces Bio.motifs.Instances and is used internally by Motif.
    """
    def __init__(self, sequences: List[Bio.Seq.Seq]):
        self.sequences = sequences
        self.length = len(sequences[0]) if sequences else 0
        self.alphabet = sequences[0].alphabet if sequences and hasattr(sequences[0], 'alphabet') else 'ACGT' # Default
        
    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, index: int) -> Bio.Seq.Seq:
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
    def __init__(self, alignment: Optional[Bio.Align.Alignment], counts: CountsMatrix, alphabet: str, instances: Instances):
        self.alignment = alignment
        self.counts = counts
        self.alphabet = alphabet
        self.instances = instances
        self.length = counts.length
        
    def __len__(self) -> int:
        return self.length

    def __getitem__(self, index: slice) -> 'Motif':
        """Allows slicing the motif to create a sub-motif."""
        # Create a new CountsMatrix from the sliced columns
        sliced_counts_data = {
            base: col[index] for base, col in self.counts.values.items()
        }
        
        # Determine the new length
        start, stop, step = index.indices(self.length)
        new_length = len(range(start, stop, step))
        
        # FIX: Ensure length is passed to the constructor
        sliced_counts = CountsMatrix(alphabet=self.alphabet, values=sliced_counts_data, length=new_length)
        
        # Slicing is difficult for alignment/instances, so we return a new motif 
        # based only on the sliced counts, setting alignment/instances to None.
        return Motif(alignment=None, counts=sliced_counts, alphabet=self.alphabet, instances=Instances([]))

    def __str__(self) -> str:
        return f"Motif of length {self.length} and alphabet {self.alphabet}\nConsensus: {self.consensus}"

    @property
    def consensus(self) -> str:
        """Returns the consensus sequence based on the majority base at each position."""
        consensus_seq: List[str] = []
        for j in range(self.length):
            # Find the base with the maximum count at position j
            max_count = -1.0
            consensus_base = 'N'
            for base in self.alphabet:
                count = self.counts[base][j]
                if count > max_count:
                    max_count = count
                    consensus_base = base
                elif count == max_count and max_count > 0:
                    # Tie: use the one that comes first in the alphabet or choose the IUPAC code
                    # For simplicity, we currently choose the alphabetically first base in a tie.
                    pass 
            consensus_seq.append(consensus_base)
        return "".join(consensus_seq)

    def reverse_complement(self) -> 'Motif':
        """Returns the reverse complement of the motif."""
        
        # A simpler way: create a new, empty CountsMatrix structure.
        new_counts_values: Dict[str, List[float]] = defaultdict(lambda: [0.0] * self.length)

        for new_pos in range(self.length):
            old_pos = self.length - 1 - new_pos
            for old_base in self.alphabet:
                new_base = _COMPLEMENT.get(old_base, old_base)
                
                # The count of new_base at new_pos is the count of old_base at old_pos
                new_counts_values[new_base][new_pos] = self.counts[old_base][old_pos]
                
        rc_counts = CountsMatrix(alphabet=self.alphabet, values=new_counts_values, length=self.length)

        # Create the new Motif instance
        return Motif(alignment=None, counts=rc_counts, alphabet=self.alphabet, instances=Instances([]))


# FIX: Simplified the type hint to List[Bio.Seq.Seq] to avoid the Codon "union already sealed" crash.
def create(alignment: List[Bio.Seq.Seq], alphabet: Optional[str] = None) -> 'Motif':
    """
    Create a Motif object from a list of sequences.
    
    This function replaces Bio.motifs.create and simplifies its logic
    to handle only a List of Bio.Seq.Seq objects.
    """
    if not alignment:
        raise ValueError("Cannot create motif from empty sequence list.")

    # Determine alphabet if not provided
    if alphabet is None:
        alphabet = alignment[0].alphabet.letters if hasattr(alignment[0].alphabet, 'letters') else 'ACGT'
        # Handle RNA U/T distinction if possible
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
            # Ignore bases not in the alphabet

    # FIX: Ensure length is passed to the constructor
    counts = CountsMatrix(alphabet=alphabet, values=counts_data, length=length)
    
    return Motif(alignment=None, counts=counts, alphabet=alphabet, instances=Instances(alignment))


def parse(path: str, fmt: str, strict: bool = True) -> List['Motif']:
    """Parse an output file from a motif finding program.

    Accepts a file **path** (str) instead of a file handle for Codon compatibility.

    Currently supported formats (case is ignored):
      - minimal: MEME minimal motif format
      - sites: Alias for minimal
    """
    fmt = fmt.lower()
    if fmt == "sites" or fmt == "minimal":
        # minimal_read must be adapted to accept a path string instead of a file handle
        return minimal_read(path)
    else:
        raise ValueError(f"Unknown format '{fmt}'. Only 'minimal' and 'sites' are supported.")


def read(path: str, fmt: str, strict: bool = True) -> 'Motif':
    """Read a single motif from a file path.

    Accepts a file **path** (str) instead of a file handle for Codon compatibility.
    """
    motifs = parse(path, fmt, strict)
    if not motifs:
        raise ValueError("No motifs found in file path")
    if len(motifs) > 1:
        # Using standard warnings mechanism
        warnings.warn(
            f"More than one motif found in file path. Returning the first one only.",
            RuntimeWarning,
            stacklevel=2 # stacklevel is 2 because this is called from 'read'
        )
        return motifs[0] # Return the first one
    return motifs[0]
