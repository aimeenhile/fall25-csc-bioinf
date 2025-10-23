import numpy as np

# Scoring parameters 
MATCH = 3
MISMATCH = -3
GAP = -2
GAP_OPEN = -5
GAP_EXTENSION = -1


# Helper functions
def score(a: str, b: str, match = MATCH, mismatch = MISMATCH) -> int:
    """ Calculate alignment score """
    if a == b:
        return MATCH
    else:
        return MISMATCH


def global_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Find an alignment with maximum alignment score """
    n = len(s1)
    m = len(s2)
    D = np.zeros((n+1, m+1), dtype=float)
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
            diag = D[i-1,j-1] + score(s1[i-1], s2[j-1], match, mismatch)
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
    i = n
    j = m
    align1_list = []
    align2_list = []

    while i > 0 or j > 0:
        if P[i, j] == 0:
            align1_list.append(s1[i-1])
            align2_list.append(s2[j-1])
            i -= 1
            j -= 1
        elif P[i, j] == 1:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            i -= 1
        else:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            j -= 1

        if i == 0 and j > 0:
            while j > 0:
                align1_list.append('-')
                align2_list.append(s2[j-1])
                j -= 1
            break
        if j == 0 and i > 0:
            while i > 0:
                align1_list.append(s1[i-1])
                align2_list.append('-')
                i -= 1
            break

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    return D[n, m], align1, align2


def local_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Substrings of s1 and s2 whose best global alignment score is maximized """
    n = len(s1)
    m = len(s2)
    V = np.zeros((n+1, m+1), dtype=float)
    P = np.zeros((n+1, m+1), dtype=int) # diag=0, up=1, left=2

    max_score = 0.0
    max_i = 0
    max_j = 0

    # Fill DP
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = V[i-1,j-1] + score(s1[i-1], s2[j-1], match, mismatch)
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

            if V[i,j] > max_score:
                max_score = V[i,j]
                max_i, max_j = i, j

    # Backtrace: start from max score
    i = max_i
    j = max_j
    align1_list = []
    align2_list = []

    while i > 0 and j > 0 and V[i,j] != 0:
        if P[i, j] == 0:
            align1_list.append(s1[i-1])
            align2_list.append(s2[j-1])
            i -= 1
            j -= 1
        elif P[i, j] == 1:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            i -= 1
        elif P[i, j] == 2:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            j -= 1
        else: # P[i,j] == -1 (hit a 0)
            break

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    return V[max_i,max_j], align1, align2
    

def fitting_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Fitting alignment of string s1 to string s2 """
    n = len(s1)
    m = len(s2)
    D = np.zeros((n+1, m+1), dtype=float)
    P = np.full((n+1, m+1), -1, dtype=int)
    
    # Initialization
    for i in range(1, n+1):
        D[i, 0] = i * gap
        P[i, 0] = 1 # up

    # Fill DP table
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = D[i-1, j-1] + score(s1[i-1], s2[j-1], match, mismatch)
            up = D[i-1, j] + gap
            left = D[i, j-1] + gap
            D[i,j] = max(diag, up, left)

            if D[i,j] == diag:
                P[i,j] = 0
            elif D[i,j] == up:
                P[i,j] = 1
            else:
                P[i,j] = 2

    # maximum score is in last row
    j = np.argmax(D[n, :])
    max_score = D[n, j]
    i = n 

    # Backtrace
    align1_list = []
    align2_list = []

    while i > 0:
        if P[i, j] == 0:
            align1_list.append(s1[i-1])
            align2_list.append(s2[j-1])
            i -= 1
            j -= 1
        elif P[i, j] == 1:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            i -= 1
        elif P[i, j] == 2:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            j -= 1
        else:
            break

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    return max_score, align1, align2

def affine_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap_open: int = GAP_OPEN, gap_extend: int = GAP_EXTENSION):
    n = len(s1)
    m = len(s2)

    # DP matrices
    lower = np.full((n+1, m+1), -np.inf, dtype=float) # insertion (gap in s2)
    middle = np.full((n+1, m+1), -np.inf, dtype=float) # matches/mismatches
    upper = np.full((n+1, m+1), -np.inf, dtype=float) # deletion (gap in s1)

    # Pointer matrices (0=middle, 1=lower, 2=upper)
    ptr_middle = np.zeros((n+1, m+1), dtype=int)
    ptr_lower = np.zeros((n+1, m+1), dtype=int)
    ptr_upper = np.zeros((n+1, m+1), dtype=int)

    # Initialization
    middle[0,0] = 0.0

    for i in range(1, n+1):
        lower[i,0] = gap_open + (i-1) * gap_extend
        ptr_lower[i,0] = 1  # From lower[i-1, 0]
    for j in range(1, m+1):
        upper[0,j] = gap_open + (j-1) * gap_extend
        ptr_upper[0,j] = 2  # From upper[0, j-1]

    # Fill DP
    for i in range(1, n+1):
        for j in range(1, m+1):
            # Lower: gap in s2 (vertical)
            val1 = lower[i-1,j] + gap_extend
            val2 = middle[i-1,j] + gap_open
            lower[i,j] = max(val1, val2)
            ptr_lower[i,j] = 1 if lower[i,j] == val1 else 0

            # Upper: gap in s1 (horizontal)
            val1 = upper[i,j-1] + gap_extend
            val2 = middle[i,j-1] + gap_open
            upper[i,j] = max(val1, val2)
            ptr_upper[i,j] = 2 if upper[i,j] == val1 else 0

            # Middle: match/mismatch
            diag = middle[i-1,j-1] + score(s1[i-1], s2[j-1], match, mismatch)
            low = lower[i,j]
            up = upper[i,j]

            middle[i,j] = max(diag, low, up)
            if middle[i,j] == diag:
                ptr_middle[i,j] = 0
            elif middle[i,j] == low:
                ptr_middle[i,j] = 1
            else:
                ptr_middle[i,j] = 2

    final_score = max(middle[n,m], lower[n,m], upper[n,m])

    # Starting state for traceback 
    if final_score == middle[n,m]:
        # Start in middle
        current = 0 
    elif final_score == lower[n,m]:
        # Start in lower
        current = 1 
    else: 
        # Start in upper
        current = 2 

    # Traceback
    i = n
    j = m
    align1_list = []
    align2_list = []
    current = 0  # start from middle

    while i > 0 or j > 0:
        if current == 0:  
            # middle
            if ptr_middle[i,j] == 0:
                align1_list.append(s1[i-1])
                align2_list.append(s2[j-1])
                i -= 1
                j -= 1
                current = 0
            elif ptr_middle[i,j] == 1: 
                current = 1
            else:
                current = 2
        elif current == 1:  
            # lower
            align1_list.append(s1[i-1])
            align2_list.append('-')
            current = ptr_lower[i,j]
            i -= 1
        elif current == 2:  
            # upper
            align1_list.append('-')
            align2_list.append(s2[j-1])
            current = ptr_upper[i,j]
            j -= 1

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    return middle[n,m], align1, align2