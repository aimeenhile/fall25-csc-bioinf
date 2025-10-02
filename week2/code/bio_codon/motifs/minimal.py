# bio_codon/motifs/minimal.py

"""Module for the support of MEME minimal motif format."""

from typing import list, dict, Any

class Record:
    """Represents a record containing a list of motifs parsed from input."""
    motifs: list[Any] 

    def __init__(self, motifs: list[Any]):
        self.motifs = motifs
        
    def __len__(self) -> int:
        return len(self.motifs)
        
    def __getitem__(self, index: int) -> Any:
        return self.motifs[index]

def read(data: list[str]) -> Record:
    """
    Reads motifs from a list of strings (lines) in a minimal format.
    
    It expects sequences to be aligned and separated by optional '> Name' lines.
    """
    # Import create function from the parent package's __init__.py
    from . import create 

    if not data:
        return Record([])
        
    name: str = "Unnamed Motif"
    sequences: list[str] = []
    first_seq_len: int = -1
    motifs: list[Any] = []

    def process_current_motif():
        """Creates a Motif object from current sequences and resets state."""
        nonlocal sequences, name, first_seq_len
        if sequences:
            # Create the Motif object using the gathered sequences
            motif = create(sequences, name=name)
            motifs.append(motif)
        
        # Reset state for the next motif
        sequences = []
        first_seq_len = -1
        name = "Unnamed Motif" 

    for line in data:
        stripped_line = line.strip()
        if not stripped_line:
            continue
            
        if stripped_line.startswith('>'):
            # Found a new motif header: finalize the previous one and start a new one
            process_current_motif()
            name = stripped_line[1:].strip()
        else:
            # Sequence data
            seq = stripped_line.upper().replace(" ", "") 
            if first_seq_len == -1:
                first_seq_len = len(seq)
                sequences.append(seq)
            elif len(seq) == first_seq_len:
                sequences.append(seq)
            else:
                # Skip unaligned sequences
                continue
    
    # Process the last motif after the loop finishes
    process_current_motif()

    return Record(motifs)