# bio_codon/motifs/minimal.py

"""Module for the support of MEME minimal motif format."""

from typing import TextIO, List, Optional, Dict, Tuple, Union
from . import Motif, Record # Import Record here
from python.Bio.Seq import Seq
from collections import defaultdict
import re

class Record(List[Motif]):
    """Container for multiple motifs read from a single file."""
    def __init__(self, motifs: Optional[List[Motif]] = None):
        """Initialize the class."""
        super().__init__(motifs or [])
        self.version: Optional[str] = None
        self.alphabet: Optional[str] = None
        self.strands: Optional[str] = None
        self.background: Dict[str, float] = {}
        # Allows access by motif name
        self._motif_name_map: Dict[str, Motif] = {}

    def __setitem__(self, key: int, value: Motif):
        """Set item and update name map."""
        if value.name:
            self._motif_name_map[value.name] = value
        super().__setitem__(key, value)

    def append(self, motif: Motif):
        """Append a motif and update name map."""
        if motif.name:
            self._motif_name_map[motif.name] = motif
        super().append(motif)

    def __getitem__(self, key: Union[int, str]) -> Motif:
        """Access motif by index or name."""
        if isinstance(key, int):
            return super().__getitem__(key)
        elif isinstance(key, str):
            if key in self._motif_name_map:
                return self._motif_name_map[key]
            else:
                raise KeyError(f"Motif with name '{key}' not found.")
        else:
            raise TypeError("Index must be an integer or motif name (string).")


def _read_version(record: Record, handle: TextIO):
    """Read MEME version (PRIVATE)."""
    for line in handle:
        if line.startswith("MEME version "):
            record.version = line.split()[-1]
            return
    raise ValueError("MEME version not found")


def _read_alphabet(record: Record, handle: TextIO):
    """Read alphabet (PRIVATE)."""
    for line in handle:
        if line.startswith("ALPHABET="):
            record.alphabet = line.split("=")[1].strip()
            # Also read strands line if present
            for next_line in handle:
                if next_line.startswith("strands:"):
                    record.strands = next_line.split(":")[-1].strip()
                    return
                elif next_line.strip() != "":
                    # If the next line is not strands and not empty, we are done
                    return
            return # If EOF reached


def _read_background(record: Record, handle: TextIO):
    """Read background frequencies (PRIVATE)."""
    for line in handle:
        if line.startswith("Background letter frequencies"):
            parts = line.split()
            # Frequencies start after the fourth word (f.e. 'frequencies')
            for i in range(4, len(parts), 2):
                base = parts[i]
                freq = float(parts[i + 1])
                record.background[base] = freq
            return
    raise ValueError("Background letter frequencies not found")


def _read_motif_info(line: str) -> Tuple[Optional[int], int, float]:
    """Read motif attributes (PRIVATE)."""
    # Example line: letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009
    
    # The "nsites= source sites" will default to 20 if it is not provided.
    num_occurrences = (
        int(line.split("nsites=")[1].split()[0]) if line.find("nsites=") != -1 else 20
    )
    # Length can be infered later if it is not provided.
    length = int(line.split("w=")[1].split()[0]) if line.find("w=") != -1 else None
    # E-value will default to zero if it is not provided.
    evalue = float(line.split("E=")[1].split()[0]) if line.find("E=") != -1 else 0.0
    return length, num_occurrences, evalue


def read(handle: TextIO) -> Record:
    """Parse the text output of the MEME program into a meme.Record object."""
    motif_number = 0
    record = Record()
    _read_version(record, handle)
    _read_alphabet(record, handle)
    _read_background(record, handle)

    while True:
        # 1. Find the MOTIF line
        motif_name = None
        for line in handle:
            if line.startswith("MOTIF"):
                motif_name = line.split()[-1]
                break
        else:
            # End of file
            return record

        # 2. Find the matrix header line
        matrix_header_line = ""
        for line in handle:
            if line.startswith("letter-probability matrix:"):
                matrix_header_line = line
                break
        else:
            raise ValueError(f"Matrix header not found after MOTIF {motif_name}")

        length, num_occurrences, evalue = _read_motif_info(matrix_header_line)

        # 3. Read the matrix itself
        if not record.alphabet:
            raise ValueError("Alphabet not defined in MEME file.")

        alphabet = record.alphabet
        counts: Dict[str, List[float]] = defaultdict(list)
        
        base_order = list(alphabet)
        
        matrix_lines: List[str] = []
        for line in handle:
            line = line.strip()
            if not line:
                break
            matrix_lines.append(line)
            
        if not matrix_lines:
            raise ValueError(f"No matrix found for MOTIF {motif_name}")

        # Process matrix lines
        for line in matrix_lines:
            values = line.split()
            if len(values) != len(base_order):
                raise ValueError(f"Expected {len(base_order)} values, got {len(values)} in matrix line: {line}")

            for base, prob_str in zip(base_order, values):
                prob = float(prob_str)
                count = prob * num_occurrences
                counts[base].append(count)

        # Infer length if not provided
        if length is None:
            length = len(counts[base_order[0]])
            
        # Create a mock alignment with a single sequence of length `length`
        from . import AlignmentMock
        dummy_sequence_str = "".join(base_order[0] for _ in range(length))
        dummy_sequence = Seq(dummy_sequence_str)
        dummy_alignment = AlignmentMock([dummy_sequence])

        # Create a new Motif object with the mock alignment
        motif = Motif(dummy_alignment, alphabet)
        
        # Manually overwrite the calculated counts with the parsed counts
        from .matrix import CountsMatrix
        final_counts: Dict[str, List[float]] = {
            base: [c for c in counts[base]] 
            for base in base_order
        }
        motif.counts = CountsMatrix(alphabet, final_counts)
        
        # Recalculate PWM and PSSM based on the new counts
        motif.pwm = motif.counts.normalize(pseudocounts=0.0)
        
        background_for_pssm = record.background or dict.fromkeys(alphabet, 1.0 / len(alphabet))
        motif.pssm = motif.pwm.log_odds(background_for_pssm)
        
        # Set metadata
        motif.name = motif_name
        motif.altname = f"MEME {motif_number+1}"
        motif.evalue = evalue
        motif.num_occurrences = num_occurrences
        from . import Instances
        motif.instances = Instances([]) # Minimal MEME has no instances

        record.append(motif)
        motif_number += 1