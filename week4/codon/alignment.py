import numpy as np

MATCH: int = 3
MISMATCH: int = -3
GAP: int = -2
GAP_OPEN: int = -5
GAP_EXTENSION: int = -1
      

def global_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Find an alignment with maximum alignment score """
    n: int = len(s1)
    m: int = len(s2)

    return 0

def local_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Substrings of s1 and s2 whose best global alignment score is maximized """
    n: int = len(s1)
    m: int = len(s2)
    
    return 0

def fitting_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    n: int = len(s1)
    m: int = len(s2)

    return 0

def affine_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap_open: int = GAP_OPEN, gap_extend: int = GAP_EXTENSION):
    n: int = len(s1)
    m: int = len(s2)

    return 0