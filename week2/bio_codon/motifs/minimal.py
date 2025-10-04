# bio_codon/motifs/minimal.py

"""Module for the support of MEME minimal motif format."""

from typing import List, Optional, Dict, Tuple, Union
from . import Motif, Record
from python import Bio.Seq.Seq
from collections import defaultdict
import re


def read(path: str) -> List[Motif]:
    """Parse the text output of the MEME program into a list of Motifs """

    record = Record()
    
    # Read all lines from the file path
    try:
        with open(path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File not found at path: {path}")
        return []
    except Exception as e:
        print(f"Error reading file: {e}")
        return []
        
    line_num = 0
    total_lines = len(lines)

    # Parse global headers
    try:
        line_num = _read_version(record, lines, line_num)
        line_num = _read_alphabet(record, lines, line_num)
        line_num = _read_background(record, lines, line_num)
    except ValueError as e:
        print(f"Error during header parsing: {e}")
        return []
    
    # Parse Motifs
    while line_num < total_lines:
        
        # Find the next MOTIF line
        motif_line = ""
        while line_num < total_lines:
            line = lines[line_num].strip()
            line_num += 1
            if line.startswith("MOTIF"):
                motif_line = line
                break
        else:
            break 
            
        try:
            # Parse the MOTIF line and the subsequent metadata line
            name, alength, w, nsites, E, line_num = _read_motif_metadata(lines, line_num, motif_line)
            
            # Read the matrix data
            matrix_data, line_num = _read_matrix(lines, line_num, name, w, record.alphabet)

            # Determine length, using 'w' if present, otherwise the actual number of rows read
            motif_length = w if w is not None else len(list(matrix_data.values())[0])
            
            # Create the Motif object 
            motif = Motif(
                alignment=None,
                counts=None, 
                alphabet=record.alphabet,
                instances=None, 
                name=name,
                length=motif_length,
                pwm_data=matrix_data
            )
            
            record.motifs.append(motif)

        except ValueError as e:
            print(f"Error parsing motif starting near line {line_num}: {e}")
            break 
            
    return record.motifs

# --- PRIVATE functions ---

def _read_version(record: Record, lines: List[str], line_num: int) -> int:
    # first line contains the version
    if line_num < len(lines):
        header = lines[line_num].strip()
        if header.startswith("MEME version"):
            record.version = header
        return line_num + 1
    return line_num


def _read_alphabet(record: Record, lines: List[str], line_num: int) -> int:
    """Find and set the alphabet from the file (PRIVATE)."""
    while line_num < len(lines):
        line = lines[line_num].strip()
        line_num += 1
        if line.startswith("ALPHABET="):
            # Extract bases from ALPHABET=A C G T...
            record.alphabet = "".join(line.split("=")[1].strip().split())
            return line_num
            
    # If not found by EOF
    raise ValueError("MEME minimal file is missing the ALPHABET definition.")


def _read_background(record: Record, lines: List[str], line_num: int) -> int:
    """Find and set the background probabilities (PRIVATE)."""
    while line_num < len(lines):
        line = lines[line_num].strip()
        line_num += 1
            
        if line.startswith("Background letter frequencies"):
            if line_num >= len(lines):
                return line_num
                 
            data_line = lines[line_num].strip()
            line_num += 1 
            
            # Data line looks like: A 0.303 C 0.183 G 0.209 T 0.306
            parts = data_line.split()
            
            background: Dict[str, float] = {}
            for i in range(0, len(parts), 2):
                base = parts[i]
                try:
                    freq = float(parts[i+1])
                    background[base] = freq
                except (IndexError, ValueError):
                    continue 
            record.background = background
            return line_num

    return line_num


def _read_motif_metadata(lines: List[str], line_num: int, motif_line: str) -> Tuple[str, Optional[int], Optional[int], Optional[int], Optional[float], int]:
    """Read MOTIF line and find the 'letter-probability matrix' line"""
    
    # Example motif_line: "MOTIF KRP"
    parts = motif_line.split()
    if len(parts) < 2 or parts[0] != "MOTIF":
        raise ValueError(f"Expected 'MOTIF' line, got: {motif_line.strip()}")
        
    name = parts[1]

    # Read the next line contains the matrix metadata
    if line_num >= len(lines):
        raise ValueError(f"File ended unexpectedly after MOTIF {name}")

    matrix_header_line = lines[line_num].strip()
    line_num += 1

    # Ensure it's the expected matrix header line
    if "letter-probability matrix:" not in matrix_header_line:
        raise ValueError(f"Expected 'letter-probability matrix:' line after MOTIF {name}, but got: {matrix_header_line}")
    
    # Initialization
    alength: Optional[int] = None
    w: Optional[int] = None
    nsites: Optional[int] = None
    E: Optional[float] = None

    # Example metadata line: "letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009"

    # Split by the colon and then by spaces, ignoring the initial text
    parts = matrix_header_line.split(":")[-1].strip().split()
    
    # Parse key=value pairs (e.g., alength= 4, w= 19)
    for part in parts:
        if part.startswith("alength="):
            alength = int(part.split("=")[1].strip())
        elif part.startswith("w="):
            w = int(part.split("=")[1].strip())
        elif part.startswith("nsites="):
            nsites = int(part.split("=")[1].strip())
        elif part.startswith("E="):
            E = float(part.split("=")[1].strip())
            
    return name, alength, w, nsites, E, line_num


def _read_matrix(lines: List[str], line_num: int, name: str, w: Optional[int], alphabet: str) -> Tuple[Dict[str, List[float]], int]:
    """Read the probability matrix (PWM) data.
    
    Assumes the file line_num is currently positioned right before the first line of the matrix data.
    """
    matrix_data: Dict[str, List[float]] = defaultdict(list)
    
    lines_read = 0
    
    # The loop immediately starts reading matrix rows
    while line_num < len(lines):
        line = lines[line_num].strip()
        line_num += 1
        
        # Matrix ends with an empty line or the start of the next block
        if not line or line.startswith("URL") or line.startswith("MOTIF"): 
            line_num -= 1
            break
            
        parts = line.split()
        if not parts:
            continue
            
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
        raise ValueError(f"Motif {name} expected width {w} but found {lines_read} rows.")
        
    return matrix_data, line_num