import numpy as np

# Scoring parameters 
MATCH = 3
MISMATCH = -3
GAP = -2
GAP_OPEN = -5
GAP_EXTENSION = -1


# Helper functions
def score(a: str, b: str) -> int:
    """ Calculate alignment score """
    if a == '-' or b == '-':
        return GAP
    elif a == b:
        return MATCH
    else:
        return MISMATCH
      
def bwt(text):
    """Return BWT string and suffix array of text"""
    text = text + "$"
    rotations = sorted([(text[i:] + text[:i], i) for i in range(len(text))])
    bwt_str = "".join(rot[0][-1] for rot in rotations)
    sa = [rot[1] for rot in rotations]
    return bwt_str, sa

def build_first_occurrence(bwt_str):
    """Return first occurrence dictionary for backward search"""
    first_col = sorted(bwt_str)
    first_occ = {}
    for i, c in enumerate(first_col):
        if c not in first_occ:
            first_occ[c] = i
    return first_occ

def build_count_matrix(bwt_str):
    """Build count array: count[c][i] = # of c in bwt_str[0:i]"""
    import collections
    counts = {c: [0]*(len(bwt_str)+1) for c in set(bwt_str)}
    for i, ch in enumerate(bwt_str):
        for c in counts:
            counts[c][i+1] = counts[c][i] + (1 if ch == c else 0)
    return counts

def backward_search(pattern, bwt_str, first_occ, counts):
    """
    Return suffix array range where pattern occurs in text.
    """
    l = 0
    r = len(bwt_str)-1
    for char in reversed(pattern):
        if char not in counts:
            return -1, -1
        l = first_occ[char] + counts[char][l]
        r = first_occ[char] + counts[char][r+1] - 1
        if r < l:
            return -1, -1
    return l, r 

def global_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Find an alignment with maximum alignment score """
    n = len(s1)
    m = len(s2)
    D = np.zeros((n+1, m+1), dtype=int)
    P = np.zeros((n+1, m+1), dtype=int) # diag=0, up=1, left=2

    # Initialization
    for i in range(1, n+1):
        D[i,0] = i * gap
        P[i,0] = 1  # up
    for j in range(1, m+1):
        D[0,j] = j * gap
        P[0,j] = 2  # left

    # Fill DP
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = D[i-1,j-1] + score(s1[i-1], s2[j-1])
            up = D[i-1,j] + gap
            left = D[i,j-1] + gap
            D[i,j] = max(diag, up, left)
            if D[i,j] == diag:
                P[i,j] = 0
            elif D[i,j] == up:
                P[i,j] = 1
            else:
                P[i,j] = 2

    # Backtrace
    align1 = ""
    align2 = ""
    i = n
    j = m

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
        else:
            align1 = '-' + align1
            align2 = s2[j - 1] + align2
            j -= 1

    return D[n, m], align1, align2


def local_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Substrings of s1 and s2 whose best global alignment score is maximized """
    n = len(s1)
    m = len(s2)
    V = np.zeros((n+1, m+1), dtype=int)
    P = np.zeros((n+1, m+1), dtype=int) # diag=0, up=1, left=2

    max_i = 0
    max_j = 0

    # Fill DP
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = V[i-1,j-1] + score(s1[i-1], s2[j-1])
            up = V[i-1,j] + gap
            left = V[i,j-1] + gap
            V[i,j] = max(diag, up, left, 0)
            if V[i,j] == diag:
                P[i,j] = 0
            elif V[i,j] == up:
                P[i,j] = 1
            elif V[i,j] == left:
                P[i,j] = 2
            else:
                P[i,j] = -1 
            if V[i,j] > V[max_i,max_j]:
                max_i, max_j = i, j

    # Backtrace: start from max score
    i = max_i
    j = max_j
    align1 = ""
    align2 = ""

    while i > 0 and j > 0 and V[i,j] != 0:
        if P[i, j] == 0:
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i -= 1
            j -= 1
        elif P[i, j] == 1:
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            j -= 1

    return V[max_i,max_j], align1, align2
    

def fitting_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Fitting alignment of string s1 to string s2 using BWT-inspired candidate positions """
    n = len(s1)
    m = len(s2)
    D = np.full((n+1, m+1), -np.inf, dtype=int)
    P = np.full((n+1, m+1), -1, dtype=int)
    
    # Initialization
    D[0, :] = 0

    # Fill DP table
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = D[i-1, j-1] + score(v[i-1], w[j-1])
            up = D[i-1, j] + GAP
            left = D[i, j-1] + GAP
            D[i,j] = max(diag, up, left)

            if D[i,j] == diag:
                P[i,j] = 0
            elif D[i,j] == up:
                P[i,j] = 1
            else:
                P[i,j] = 2

    # maximum score is in last column of any row
    i = np.argmax(D[:, m])
    max_score = D[i, m]

    # Backtrace
    j = m
    align1 = ""
    align2 = ""
    while j > 0:
        if i > 0 and P[i, j] == 0:
            aligned_v = v[i-1] + aligned_v
            aligned_w = w[j-1] + aligned_w
            i -= 1
            j -= 1
        elif i > 0 and P[i, j] == 1:
            aligned_v = v[i-1] + aligned_v
            aligned_w = '-' + aligned_w
            i -= 1
        else:  # P[i, j] == 2
            aligned_v = '-' + aligned_v
            aligned_w = w[j-1] + aligned_w
            j -= 1

    return max_score, align1, align2

def affine_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap_open: int = GAP_OPEN, gap_extend: int = GAP_EXTENSION):
    n = len(s1)
    m = len(s2)

    # DP matrices
    lower = np.full((n+1, m+1), -np.inf) # insertion (gap in s2)
    middle = np.full((n+1, m+1), -np.inf) # matches/mismatches
    upper = np.full((n+1, m+1), -np.inf) # deletion (gap in s1)

    # Pointer matrices (0=middle, 1=lower, 2=upper)
    ptr_middle = np.zeros((n+1, m+1), dtype=int)
    ptr_lower = np.zeros((n+1, m+1), dtype=int)
    ptr_upper = np.zeros((n+1, m+1), dtype=int)

    middle[0,0] = 0
    for i in range(1, n+1):
        lower[i,0] = gap_open + (i-1)*gap_extend
        ptr_lower[i,0] = 1  # up
    for j in range(1, m+1):
        upper[0,j] = gap_open + (j-1)*gap_extend
        ptr_upper[0,j] = 2  # left

    # Fill DP
    for i in range(1, n+1):
        for j in range(1, m+1):
            # Lower: gap in s2 (vertical)
            val1 = lower[i-1,j] + gap_extend
            val2 = middle[i-1,j] + gap_open
            lower[i,j] = max(val1, val2)
            ptr_lower[i,j] = 1 if lower[i,j]==val1 else 0

            # Upper: gap in s1 (horizontal)
            val1 = upper[i,j-1] + gap_extend
            val2 = middle[i,j-1] + gap_open
            upper[i,j] = max(val1, val2)
            ptr_upper[i,j] = 2 if upper[i,j]==val1 else 0

            # Middle: match/mismatch
            diag = middle[i-1,j-1] + score(s1[i-1], s2[j-1])
            low = lower[i,j]
            up = upper[i,j]
            middle[i,j] = max(diag, low, up)
            if middle[i,j] == diag:
                ptr_middle[i,j] = 0
            elif middle[i,j] == low:
                ptr_middle[i,j] = 1
            else:
                ptr_middle[i,j] = 2

    # Traceback
    i = n
    j = m
    align1 = ""
    align2 = ""
    current = 0  # start from middle

    while i > 0 or j > 0:
        if current == 0:  
            # middle
            if ptr_middle[i,j] == 0:
                align1 = s1[i-1] + align1
                align2 = s2[j-1] + align2
                i -= 1
                j -= 1
                current = 0
            elif ptr_middle[i,j] == 1: 
                current = 1
            else:
                current = 2
        elif current == 1:  
            # lower
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            i -= 1
            current = ptr_lower[i+1,j]
        elif current == 2:  
            # upper
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            j -= 1
            current = ptr_upper[i,j+1]

    return middle[n,m], align, align