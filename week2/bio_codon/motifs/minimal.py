# bio_codon/motifs/minimal.py

"""Module for the support of MEME minimal motif format."""

from typing import List, Optional, Dict, Tuple, Union
from .__init__ import Motif, Record 
from python.Bio.Seq import Seq
from collections import defaultdict
import re


def read(path: str) -> List[Motif]:
    """Parse the text output of the MEME program into a meme.Record object (takes path in Codon).

    In the Codon environment, this function takes a file path (str)
    and handles file opening/closing internally.
    """
    motif_number = 0
    record = Record()
    
    # Use 'with open(path)' to handle the file
    try:
        with open(path, 'r') as handle:
            _read_version(record, handle)
            _read_alphabet(record, handle)
            _read_background(record, handle)

            while True:
                # Find the next MOTIF line
                line = ""
                for line in handle:
                    if line.startswith("MOTIF"):
                        break
                else:
                    return record
                
                # Parse the MOTIF line for info and the next 'letter-probability matrix:' line
                name, alength, w, nsites, E, matrix_start_line = _read_motif_metadata(handle, line)
                
                # Read the matrix data from the handle
                matrix_data = _read_matrix(handle, name, alength, w, record.alphabet)
                
                # Create the Motif object
                motif = Motif.from_minimal_format(
                    name=name,
                    alphabet=record.alphabet,
                    matrix_data=matrix_data,
                    num_occurrences=nsites,
                    evalue=E,
                    background=record.background,
                    alength=alength,
                    w=w
                )
                record.append(motif)
    except FileNotFoundError:
        print(f"Error: Minimal format file not found at path: {path}")
        return []
    except Exception as e:
        print(f"An error occurred during minimal file parsing: {e}")
        return []

# Helper functions that operate on the open file handle (internal to minimal.py)

def _read_version(record, handle):
    """Read MEME version (PRIVATE)."""
    # Assuming the first line is always 'MEME version X'
    try:
        line = next(handle).strip()
        if not line.startswith("MEME version"):
            raise ValueError("File does not start with MEME version tag.")
        record.version = line.split()[-1]
    except StopIteration:
        raise ValueError("File is empty or truncated.")

def _read_alphabet(record, handle):
    """Read ALPHABET (PRIVATE)."""
    for line in handle:
        if line.startswith("ALPHABET="):
            record.alphabet = line.strip().split('=')[1].strip()
            break
    else:
        raise ValueError("ALPHABET not found.")

def _read_background(record, handle):
    """Read Background letter frequencies (PRIVATE)."""
    record.background = defaultdict(float)
    for line in handle:
        if line.startswith("Background letter frequencies"):
            break
    else:
        return # Background is optional

    # Read the actual frequencies on the next line
    try:
        line = next(handle).strip()
        parts = line.split()
        for i in range(0, len(parts), 2):
            base = parts[i]
            freq = float(parts[i+1])
            record.background[base] = freq
    except StopIteration:
        # File ended abruptly after the header
        pass


def _read_motif_metadata(handle, motif_line: str) -> Tuple[str, Optional[int], Optional[int], Optional[int], Optional[float], str]:
    """Read motif's name and metadata (PRIVATE)."""
    # motif_line example: MOTIF KRP
    name = motif_line.split()[1]

    # Find the 'letter-probability matrix:' line
    matrix_start_line = ""
    for line in handle:
        if line.startswith("letter-probability matrix:"):
            matrix_start_line = line.strip()
            break
    
    if not matrix_start_line:
        raise ValueError(f"Could not find matrix for motif {name}")

    # The "nsites= source sites" will default to 20 if it is not provided.
    num_occurrences = None
    if "nsites=" in matrix_start_line:
        num_occurrences = int(matrix_start_line.split("nsites=")[1].split()[0])
    
    # Length (w)
    length = None
    if "w=" in matrix_start_line:
        length = int(matrix_start_line.split("w=")[1].split()[0])
        
    # E-value will default to zero if it is not provided.
    evalue = None
    if "E=" in matrix_start_line:
        evalue = float(matrix_start_line.split("E=")[1].split()[0])
        
    return name, length, num_occurrences, evalue, matrix_start_line


def _read_matrix(handle, name: str, alength: Optional[int], w: Optional[int], alphabet: str) -> Dict[str, List[float]]:
    """Read the probability matrix (PWM) data (PRIVATE)."""
    matrix_data: Dict[str, List[float]] = defaultdict(list)
    
    # The matrix starts on the line after 'letter-probability matrix:'
    
    expected_matrix_lines = w if w is not None else 1000 # Heuristic max if w is unknown
    lines_read = 0
    
    for line in handle:
        line = line.strip()
        if not line: # Empty line marks the end of the matrix
            break
            
        parts = line.split()
        if not parts:
            break
            
        if len(parts) != len(alphabet):
            raise ValueError(f"Motif {name} matrix line has {len(parts)} columns but alphabet has {len(alphabet)} bases.")
            
        for i, base in enumerate(alphabet):
            try:
                matrix_data[base].append(float(parts[i]))
            except ValueError as e:
                raise ValueError(f"Error parsing matrix value '{parts[i]}' for motif {name}: {e}")
                
        lines_read += 1
        if w is not None and lines_read >= w:
            break
            
    # Final check on length
    if w is not None and lines_read != w:
        raise ValueError(f"Motif {name}: Expected length W={w} but read {lines_read} rows.")

    return matrix_data
