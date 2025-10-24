import numpy as np

# Scoring parameters 
MATCH = 3
MISMATCH = -3
GAP = -2
GAP_OPEN = -5
GAP_EXTENSION = -1


def global_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Find an alignment with maximum alignment score """

    def nw_score(X, Y):
        """ Needleman-Wunsch score vector using only two rows (for Hirschberg) """
        m = len(Y)
        prev = np.zeros(m + 1, dtype=float)
        for j in range(1, m + 1):
            prev[j] = j * gap
        for i in range(1, len(X) + 1):
            curr = np.zeros(m + 1, dtype=float)
            curr[0] = i * gap
            for j in range(1, m + 1):
                diag = prev[j-1] + (match if X[i-1] == Y[j-1] else mismatch)
                up = prev[j] + gap
                left = curr[j-1] + gap
                curr[j] = max(diag, up, left)
            prev = curr
        return prev

    def hirschberg(X, Y):
        """ Hirschberg recursive alignment """
        n = len(X)
        m = len(Y)
        if n == 0:
            return '-'*m, Y
        elif m == 0:
            return X, '-'*n
        elif n == 1 or m == 1:
            # Small problem: compute full DP and traceback
            D = np.zeros((n+1, m+1), dtype=float)
            for i in range(1, n+1):
                D[i,0] = i * gap
            for j in range(1, m+1):
                D[0,j] = j * gap
            for i in range(1, n+1):
                for j in range(1, m+1):
                    diag = D[i-1,j-1] + (match if X[i-1] == Y[j-1] else mismatch)
                    up = D[i-1,j] + gap
                    left = D[i,j-1] + gap
                    D[i,j] = max(diag, up, left)
            # Traceback
            i, j = n, m
            align_X, align_Y = [], []
            while i > 0 or j > 0:
                if i > 0 and j > 0 and D[i,j] == D[i-1,j-1] + (match if X[i-1] == Y[j-1] else mismatch):
                    align_X.append(X[i-1])
                    align_Y.append(Y[j-1])
                    i -= 1
                    j -= 1
                elif i > 0 and D[i,j] == D[i-1,j] + gap:
                    align_X.append(X[i-1])
                    align_Y.append('-')
                    i -= 1
                else:
                    align_X.append('-')
                    align_Y.append(Y[j-1])
                    j -= 1
            return ''.join(reversed(align_X)), ''.join(reversed(align_Y))
        else:
            # Divide X in half
            i_mid = n // 2
            score_L = nw_score(X[:i_mid], Y)
            score_R = nw_score(X[i_mid:][::-1], Y[::-1])
            j_split = np.argmax(score_L + score_R[::-1])
            # Recurse
            left_X, left_Y = hirschberg(X[:i_mid], Y[:j_split])
            right_X, right_Y = hirschberg(X[i_mid:], Y[j_split:])
            return left_X + right_X, left_Y + right_Y

    # Compute alignment
    align1, align2 = hirschberg(s1, s2)

    # Compute final score
    score = 0.0
    for a, b in zip(align1, align2):
        if a == '-' or b == '-':
            score += gap
        else:
            score += match if a == b else mismatch

    return score, align1, align2

"""
    n = len(s1)
    m = len(s2)
    D_prev = np.zeros(m+1, dtype=float)
    D_curr = np.zeros(m+1, dtype=float)
    P = np.zeros((n+1, m+1), dtype=np.int8) # diag=0, up=1, left=2

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
            D_curr[j] = max(diag, up, left)

            if D_curr[j] == diag:
                P[i,j] = 0
            elif D_curr[j] == up:
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

    score = D_prev[m]

    return score, align1, align2
"""


def local_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Substrings of s1 and s2 whose best global alignment score is maximized """
    n = len(s1)
    m = len(s2)

    prev = np.zeros(m+1, dtype=int)
    curr = np.zeros(m+1, dtype=int)

    max_score = 0
    max_pos = (0, 0)

    # Fill DP with rolling rows
    for i in range(1, n+1):
        curr[0] = 0
        for j in range(1, m+1):
            s1_i = s1[i-1]
            s2_j = s2[j-1]
            diag = prev[j-1] + (match if s1_i == s2_j else mismatch)
            up = prev[j] + gap
            left = curr[j-1] + gap
            curr[j] = max(0, diag, up, left)
            
            if curr[j] > max_score:
                max_score = curr[j]
                max_pos = (i, j)
        
        prev, curr = curr, prev  

    # Traceback for small window
    i, j = max_pos
    align1 = []
    align2 = []
    score_sw = max_score

    while score_sw > 0 and i > 0 and j > 0:
        diag = (match if s1[i-1] == s2[j-1] else mismatch)
        if prev[j-1] + diag == score_sw:
            align1.append(s1[i-1])
            align2.append(s2[j-1])
            i -= 1
            j -= 1
            score_sw -= diag
        elif prev[j] + gap == score_sw:
            align1.append(s1[i-1])
            align2.append('-')
            i -= 1
            score_sw -= gap
        else:
            align1.append('-')
            align2.append(s2[j-1])
            j -= 1
            score_sw -= gap

    align1 = ''.join(reversed(align1))
    align2 = ''.join(reversed(align2))

    return max_score, align1, align2    

    """
    V_prev = np.zeros(m + 1, dtype=float)
    V_curr = np.zeros(m + 1, dtype=float)
    P = np.full((n+1, m+1), -1, dtype=np.int8) # diag=0, up=1, left=2

    max_score = 0.0
    max_i = 0
    max_j = 0

    # Fill DP
    for i in range(1, n+1):
        V_curr[0] = 0.0
        s1_i = s1[i-1]
        for j in range(1, m+1):
            s2_j = s2[j-1]
            diag = V_prev[j-1] + (match if s1_i == s2_j else mismatch)
            up = V_prev[j] + gap
            left = V_curr[j-1] + gap

            val = max(diag, up, left, 0)
            V_curr[j] = val
            if val == 0.0:
                P[i, j] = -1
            else:
                if val == diag:
                    P[i,j] = 0
                elif val == up:
                    P[i,j] = 1
                else:
                    P[i,j] = 2

            if val > max_score:
                max_score = val
                max_i, max_j = i, j

        V_prev, V_curr = V_curr, V_prev

    # Backtrace: start from max score
    i = max_i
    j = max_j
    align1_list = []
    align2_list = []

    while i > 0 and j > 0 and P[i,j] != -1:
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

    return max_score, align1, align2          
    """
        
    

def fitting_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Fitting alignment of string s1 to string s2 """
    n = len(s1)
    m = len(s2)
    D_prev = np.zeros(m + 1, dtype=np.float32)
    D_curr = np.zeros(m + 1, dtype=np.float32)
    P = np.full((n+1, m+1), -1, dtype=np.int8)

    # Fill DP table
    for i in range(1, n+1):
        D[i, 0] = i * gap
        s1_i = s1[i-1]
        P[i, 0] = 1 # up
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

    # maximum score is in last row
    j = int(np.argmax(D_prev))
    max_score = float(D_prev[j])
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
            if val1 >= val2:
                lower[i,j] = val1
                ptr_lower[i,j] = 1
            else:
                lower[i,j] = val2
                ptr_lower[i,j] = 0

            # Upper: gap in s1 (horizontal)
            val1 = upper[i,j-1] + gap_extend
            val2 = middle[i,j-1] + gap_open
            if val1 >= val2:
                upper[i,j] = val1
                ptr_upper[i,j] = 2
            else:
                upper[i,j] = val2
                ptr_upper[i,j] = 0

            # Middle: match/mismatch
            s1_i = s1[i-1]
            s2_j = s2[j-1]
            diag = middle[i-1, j-1] + (match if s1_i == s2_j else mismatch)
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

    return final_score, align1, align2