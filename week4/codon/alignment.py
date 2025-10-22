import numpy as np

MATCH: int = 3
MISMATCH: int = -3
GAP: int = -2
GAP_OPEN: int = -5
gAP_EXTENSION: int = -1

# Helper functions
def score(a: str, b: str) -> int:
    """ Calculate alignment score """
    if a == '-' or b == '-':
        return GAP
    elif a == b:
        return MATCH
    else:
        return MISMATCH
      

def global_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Find an alignment with maximum alignment score """
    n: int = len(s1)
    m: int = len(s2)
    D = np.full((n+1, m+1), -1000, dtype=int)
    P = np.full((n+1, m+1), -1, dtype=int)

    def dp(i, j):
        if D[i, j] != -np.inf:
            return D[i, j]
        if i == 0 and j == 0:
            # base case: empty strings
            D[i, j] = 0
            P[i, j] = -1
            return 0
        if i == 0:
            # gap in s1
            D[i, j] = j * gap
            P[i, j] = 2  # left
            return D[i, j]
        if j == 0:
            # gap in s2
            D[i, j] = i * gap
            P[i, j] = 1  # up
            return D[i, j]

        diag: int = dp(i-1, j-1) + score(s1[i-1], s2[j-1], match, mismatch, gap)
        up: int = dp(i-1, j) + gap
        left: int = dp(i, j-1) + gap

        D[i, j] = max(diag, up, left)
        if D[i, j] == diag:
            P[i, j] = 0
        elif D[i, j] == up:
            P[i, j] = 1
        else:
            P[i, j] = 2

        return D[i, j]

    score_max: int = dp(n, m)

    # Traceback
    align1: str = ""
    align2: str = ""
    i: int = n
    j: int = m

    while i > 0 or j > 0:
        if P[i, j] == 0:
            align1 = s1[i - 1] + align1
            align2 = s2[j - 1] + align2
            i -= 1
            j -= 1
        elif P[i, j] == 1:
            align1 = s1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        elif P[i, j] == 2:
            align1 = '-' + align1
            align2 = s2[j - 1] + align2
            j -= 1
        else:
            break

    return D[n, m], align1, align2


def local_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Substrings of s1 and s2 whose best global alignment score is maximized """
    n: int = len(s1)
    m: int = len(s2)
    D = np.full((n+1, m+1), -1000, dtype=int)
    P = np.full((n+1, m+1), -1, dtype=int)
    
    return 0

def fitting_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    n: int = len(s1)
    m: int = len(s2)
    D = np.full((n+1, m+1), -1000, dtype=int)
    P = np.full((n+1, m+1), -1, dtype=int)

    return 0

def affine_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap_open: int = GAP_OPEN, gap_extend: int = GAP_EXTENSION):
    n: int = len(s1)
    m: int = len(s2)
    D = np.full((n+1, m+1), -1000, dtype=int)
    P = np.full((n+1, m+1), -1, dtype=int)

    return 0