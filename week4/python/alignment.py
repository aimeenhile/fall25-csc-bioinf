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

    # Only two rows at a time
    prev = np.zeros(m+1, dtype=int)
    curr = np.zeros(m+1, dtype=int)
    # Pointers for traceback (store indices for backtrace)
    trace = [ [0]*(m+1) for _ in range(n+1) ]

    # Initialization
    for j in range(1, m+1):
        prev[j] = j * gap
        trace[0][j] = 2  # left
    trace[0][0] = -1

    for i in range(1, n+1):
        curr[0] = i * gap
        trace[i][0] = 1  # up
        for j in range(1, m+1):
            diag = prev[j-1] + score(s1[i-1], s2[j-1], match, mismatch)
            up = prev[j] + gap
            left = curr[j-1] + gap
            best = max(diag, up, left)
            curr[j] = best
            if best == diag:
                trace[i][j] = 0
            elif best == up:
                trace[i][j] = 1
            else:
                trace[i][j] = 2
        prev, curr = curr, prev

    # Traceback
    i, j = n, m
    align1_list, align2_list = [], []
    while i > 0 or j > 0:
        if trace[i][j] == 0:
            align1_list.append(s1[i-1])
            align2_list.append(s2[j-1])
            i -= 1
            j -= 1
        elif trace[i][j] == 1:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            i -= 1
        else:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            j -= 1

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))
    return prev[m], align1, align2


def local_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Substrings of s1 and s2 whose best global alignment score is maximized """
    n = len(s1)
    m = len(s2)

    # Only two rows needed for DP
    prev = np.zeros(m+1, dtype=int)
    curr = np.zeros(m+1, dtype=int)
    trace = [ [0]*(m+1) for _ in range(n+1) ]  # pointers: diag=0, up=1, left=2, start=None=-1

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n+1):
        curr[0] = 0
        for j in range(1, m+1):
            diag = prev[j-1] + score(s1[i-1], s2[j-1], match, mismatch)
            up = prev[j] + gap
            left = curr[j-1] + gap
            best = max(diag, up, left, 0)
            curr[j] = best

            # Trace pointer
            if best == 0:
                trace[i][j] = -1
            elif best == diag:
                trace[i][j] = 0
            elif best == up:
                trace[i][j] = 1
            else:
                trace[i][j] = 2

            # Track maximum score
            if best > max_score:
                max_score = best
                max_pos = (i, j)
        prev, curr = curr, prev

    # Traceback from max_score
    i, j = max_pos
    align1_list, align2_list = [], []

    while i > 0 and j > 0 and trace[i][j] != -1:
        if trace[i][j] == 0:
            align1_list.append(s1[i-1])
            align2_list.append(s2[j-1])
            i -= 1
            j -= 1
        elif trace[i][j] == 1:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            i -= 1
        elif trace[i][j] == 2:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            j -= 1

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    return max_score, align1, align2
    

def fitting_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Fitting alignment of string s1 to string s2 using BWT-inspired candidate positions """
    n = len(s1)
    m = len(s2)

    # Only store two rows at a time
    prev = np.zeros(m+1, dtype=int)
    curr = np.zeros(m+1, dtype=int)
    # To remember which j in last column gives max score
    end_i = 0

    # Initialization: first row is all zeros (fits anywhere in s2)
    prev[:] = 0

    # Fill DP table row by row
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = prev[j-1] + score(s1[i-1], s2[j-1], match, mismatch)
            up = prev[j] + gap
            left = curr[j-1] + gap
            curr[j] = max(diag, up, left)
        prev, curr = curr, prev

    # Find maximum in last column (fitting alignment ends at s2[-1])
    j = m
    i = np.argmax(prev)  # row index of max score in last column
    max_score = prev[i]

    # Recompute DP for traceback submatrix (small enough)
    sub_s1 = s1[:i]
    sub_s2 = s2
    D = np.zeros((len(sub_s1)+1, len(sub_s2)+1), dtype=int)
    P = np.zeros((len(sub_s1)+1, len(sub_s2)+1), dtype=int)  # diag=0, up=1, left=2

    for ii in range(1, len(sub_s1)+1):
        for jj in range(1, len(sub_s2)+1):
            diag = D[ii-1,jj-1] + score(sub_s1[ii-1], sub_s2[jj-1], match, mismatch)
            up = D[ii-1,jj] + gap
            left = D[ii,jj-1] + gap
            D[ii,jj] = max(diag, up, left)
            if D[ii,jj] == diag:
                P[ii,jj] = 0
            elif D[ii,jj] == up:
                P[ii,jj] = 1
            else:
                P[ii,jj] = 2

    # Traceback from (i, m)
    align1_list = []
    align2_list = []
    ii = len(sub_s1)
    jj = len(sub_s2)
    while ii > 0 and jj > 0:
        if P[ii,jj] == 0:
            align1_list.append(sub_s1[ii-1])
            align2_list.append(sub_s2[jj-1])
            ii -= 1
            jj -= 1
        elif P[ii,jj] == 1:
            align1_list.append(sub_s1[ii-1])
            align2_list.append('-')
            ii -= 1
        else:  # left
            align1_list.append('-')
            align2_list.append(sub_s2[jj-1])
            jj -= 1

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    return max_score, align1, align2


def affine_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap_open: int = GAP_OPEN, gap_extend: int = GAP_EXTENSION):
    n = len(s1)
    m = len(s2)

    # Initialize matrices: only two rows for each
    lower_prev = np.full(m+1, -np.inf, dtype=float)
    lower_curr = np.full(m+1, -np.inf, dtype=float)
    upper_prev = np.full(m+1, -np.inf, dtype=float)
    upper_curr = np.full(m+1, -np.inf, dtype=float)
    middle_prev = np.full(m+1, -np.inf, dtype=float)
    middle_curr = np.full(m+1, -np.inf, dtype=float)

    # Trace pointers for reconstruction
    trace = [ [0]*(m+1) for _ in range(n+1) ]

    middle_prev[0] = 0
    lower_prev[0] = -np.inf
    upper_prev[0] = -np.inf
    trace[0][0] = -1

    # Fill first row
    for j in range(1, m+1):
        upper_prev[j] = gap_open + (j-1)*gap_extend
        middle_prev[j] = upper_prev[j]
        trace[0][j] = 2  # left

    for i in range(1, n+1):
        lower_curr[0] = gap_open + (i-1)*gap_extend
        middle_curr[0] = lower_curr[0]
        upper_curr[0] = -np.inf
        trace[i][0] = 1  # up

        for j in range(1, m+1):
            # Lower (gap in s2)
            lower_curr[j] = max(lower_prev[j] + gap_extend, middle_prev[j] + gap_open)
            # Upper (gap in s1)
            upper_curr[j] = max(upper_curr[j-1] + gap_extend, middle_curr[j-1] + gap_open)
            # Middle (match/mismatch)
            diag = middle_prev[j-1] + score(s1[i-1], s2[j-1], match, mismatch)
            middle_curr[j] = max(diag, lower_curr[j], upper_curr[j])

            # Trace pointer for reconstruction
            if middle_curr[j] == diag:
                trace[i][j] = 0
            elif middle_curr[j] == lower_curr[j]:
                trace[i][j] = 1
            else:
                trace[i][j] = 2

        # Swap rows
        lower_prev, lower_curr = lower_curr, lower_prev
        upper_prev, upper_curr = upper_curr, upper_prev
        middle_prev, middle_curr = middle_curr, middle_prev

    # Traceback from bottom-right
    i, j = n, m
    align1_list, align2_list = [], []
    while i > 0 or j > 0:
        if trace[i][j] == 0:
            align1_list.append(s1[i-1])
            align2_list.append(s2[j-1])
            i -= 1
            j -= 1
        elif trace[i][j] == 1:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            i -= 1
        else:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            j -= 1

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    return middle_prev[m], align1, align2