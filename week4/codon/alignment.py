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

    D_prev = np.ndarray[int, 1]((m + 1,))
    D_curr = np.ndarray[int, 1]((m+1,))

    for j in range(m+1):
        D_prev[j] = 0
        D_curr[j] = 0

    P = np.ndarray[int, 2]((n+1, m+1)) # diag=0, up=1, left=2

    for i in range(n+1):
        for j in range(m+1):
            P[i,j] = -1

    # Initialization
    for j in range(1, m+1):
        D_prev[j] = j * gap
        P[0,j] = 2  # left

    # Fill DP
    for i in range(1, n+1):
        D_curr[0] = i * gap
        P[i,0] = 1  # up
        s1_i = s1[i-1]
        for j in range(1, m+1):
            s2_j = s2[j-1]
            diag = D_prev[j-1] + (match if s1_i == s2_j else mismatch)
            up = D_prev[j] + gap
            left = D_curr[j-1] + gap

            val = max(diag, up, left)
            D_curr[j] = val

            if val == diag:
                P[i,j] = 0
            elif val == up:
                P[i,j] = 1
            else:
                P[i,j] = 2

        D_prev, D_curr = D_curr, D_prev

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

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    max_score = D_prev[m]

    return max_score, align1, align2


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