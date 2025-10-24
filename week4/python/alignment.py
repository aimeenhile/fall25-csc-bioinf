import numpy as np

# Scoring parameters 
MATCH = 3
MISMATCH = -3
GAP = -2
GAP_OPEN = -5
GAP_EXTENSION = -1


def global_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Find an alignment with maximum alignment score """
    n = len(s1)
    m = len(s2)
    D_prev = np.zeros(m+1, dtype=np.float32)
    D_curr = np.zeros(m+1, dtype=np.float32)
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
    n = len(s1)
    m = len (s2)

    V_prev = np.zeros(m + 1, dtype=np.float32)
    V_curr = np.zeros(m + 1, dtype=np.float32)
    P = np.full((n+1, m+1), -1, dtype=np.int8) # diag=0, up=1, left=2

    max_score = 0
    max_i = 0
    max_j = 0

    # Fill DP
    for i in range(1, n+1):
        V_curr[0] = 0
        s1_i = s1[i-1]
        for j in range(1, m+1):
            s2_j = s2[j-1]
            diag = V_prev[j-1] + (match if s1_i == s2_j else mismatch)
            up = V_prev[j] + gap
            left = V_curr[j-1] + gap

            val = max(diag, up, left, 0)
            V_curr[j] = val
            if val == 0:
                P[i, j] = -1
            elif val == diag:
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
        

def fitting_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Fitting alignment of string s1 to string s2 """
    n = len(s1) 
    m = len(s2)

    if n <= m:
        free_start_s1, free_end_s1 = True, True
        free_start_s2, free_end_s2 = False, False
    else:
        free_start_s1, free_end_s1 = False, False
        free_start_s2, free_end_s2 = True, True
        
    prev = np.zeros(m + 1, dtype=np.float32)
    curr = np.zeros(m + 1, dtype=np.float32)
    P = np.full((n+1, m+1), -1, dtype=np.int8)

    # Initialization
    # Row 0
    for j in range(m + 1):
        if free_start_s1:
            prev[j] = 0
        else:
            prev[j] = j * gap
    # Col 0 
    if not free_start_s2:
        for i in range(1, n + 1):
            P[i, 0] = 1
    else:
        for i in range(1, n + 1):
            P[i, 0] = -1  

    # Fill DP
    for i in range(1, n + 1):
        s1_i = s1[i-1]

        if free_start_s2:
            curr[0] = 0
        else:
            curr[0] = prev[0] + gap
            P[i, 0] = 1
        
        for j in range(1, m + 1):
            s2_j = s2[j-1]
            diag = prev[j-1] + (match if s1_i == s2_j else mismatch)
            up = prev[j] + gap
            left = curr[j-1] + gap

            val = max(diag, up, left)
            curr[j] = val

            if val == diag:
                P[i,j] = 0
            elif val == up:
                P[i,j] = 1
            else:
                P[i,j] = 2

        prev, curr = curr, prev

    # Traceback start pos
    if free_end_s2:
        j = int(np.argmax(prev))
        max_score = prev[j]
        i = n
    elif free_end_s1:
        # find max over last column
        col_vals = np.array([prev[-1] if k == n else P[k, m] for k in range(n + 1)])
        i = int(np.argmax(col_vals))
        j = m
        max_score = col_vals[i]
    else:
        i, j = n, m
        max_score = prev[m]

    # Traceback
    align1_list = []
    align2_list = []

    while i > 0:
        pos = P[i,j]
        if pos == 0:
            align1_list.append(s1[i-1])
            align2_list.append(s2[j-1])
            i -= 1
            j -= 1
        elif pos == 1:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            i -= 1
        elif pos == 2:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            j -= 1
        else:
            break
        
        if (free_start_s1 and i == 0) or (free_start_s2 and j == 0):
            break
        
    while j > 0 and not free_start_s2:
        align1_list.append('-')
        align2_list.append(s2[j-1])
        j -= 1
    while i > 0 and not free_start_s1:
        align1_list.append(s1[i - 1])
        align2_list.append('-')
        i -= 1

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    return max_score, align1, align2


def affine_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap_open: int = GAP_OPEN, gap_extend: int = GAP_EXTENSION):

    def affine_forward(a, b):
        """Compute last row of affine DP (forward)."""
        n, m = len(a), len(b)
        M_prev = np.full(m + 1, -np.inf)
        X_prev = np.full(m + 1, -np.inf)
        Y_prev = np.full(m + 1, -np.inf)
        M_prev[0] = 0
        for j in range(1, m + 1):
            Y_prev[j] = gap_open + (j - 1) * gap_extend
        for i in range(1, n + 1):
            M_curr = np.full(m + 1, -np.inf)
            X_curr = np.full(m + 1, -np.inf)
            Y_curr = np.full(m + 1, -np.inf)
            X_curr[0] = gap_open + (i - 1) * gap_extend
            for j in range(1, m + 1):
                sc = match if a[i - 1] == b[j - 1] else mismatch
                X_curr[j] = max(X_prev[j] + gap_extend, M_prev[j] + gap_open)
                Y_curr[j] = max(Y_curr[j - 1] + gap_extend, M_curr[j - 1] + gap_open)
                M_curr[j] = max(M_prev[j - 1] + sc, X_curr[j], Y_curr[j])
            M_prev, X_prev, Y_prev = M_curr, X_curr, Y_curr
        return np.maximum.reduce([M_prev, X_prev, Y_prev])

    def affine_reverse(a, b):
        """Compute first row of affine DP backwards (for suffix alignment)."""
        n, m = len(a), len(b)
        M_prev = np.full(m + 1, -np.inf)
        X_prev = np.full(m + 1, -np.inf)
        Y_prev = np.full(m + 1, -np.inf)
        M_prev[m] = 0
        for j in range(m - 1, -1, -1):
            Y_prev[j] = gap_open + (m - j - 1) * gap_extend
        for i in range(n - 1, -1, -1):
            M_curr = np.full(m + 1, -np.inf)
            X_curr = np.full(m + 1, -np.inf)
            Y_curr = np.full(m + 1, -np.inf)
            X_curr[m] = gap_open + (n - i - 1) * gap_extend
            for j in range(m - 1, -1, -1):
                sc = match if a[i] == b[j] else mismatch
                X_curr[j] = max(X_prev[j] + gap_extend, M_prev[j] + gap_open)
                Y_curr[j] = max(Y_curr[j + 1] + gap_extend, M_curr[j + 1] + gap_open)
                M_curr[j] = max(M_prev[j + 1] + sc, X_curr[j], Y_curr[j])
            M_prev, X_prev, Y_prev = M_curr, X_curr, Y_curr
        return np.maximum.reduce([M_prev, X_prev, Y_prev])

    def divide_and_conquer(a, b):
        """Recursive alignment reconstruction."""
        n, m = len(a), len(b)
        if n == 0:
            return "-" * m, b
        if m == 0:
            return a, "-" * n
        if n == 1 or m == 1:
            # Fallback to simple DP for small subproblems
            return needleman_wunsch_small(a, b)

        mid = n // 2
        score_left = affine_forward(a[:mid], b)
        score_right = affine_reverse(a[mid:], b)
        split = int(np.argmax(score_left + score_right))

        left_a, left_b = divide_and_conquer(a[:mid], b[:split])
        right_a, right_b = divide_and_conquer(a[mid:], b[split:])
        return left_a + right_a, left_b + right_b

    def needleman_wunsch_small(a, b):
        """Simple O(n*m) DP for small subproblems."""
        n, m = len(a), len(b)
        dp = np.zeros((n + 1, m + 1))
        for i in range(1, n + 1):
            dp[i, 0] = gap_open + (i - 1) * gap_extend
        for j in range(1, m + 1):
            dp[0, j] = gap_open + (j - 1) * gap_extend
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                sc = match if a[i - 1] == b[j - 1] else mismatch
                dp[i, j] = max(
                    dp[i - 1, j - 1] + sc,
                    dp[i - 1, j] + (gap_extend if dp[i - 2, j] != 0 else gap_open),
                    dp[i, j - 1] + (gap_extend if dp[i, j - 2] != 0 else gap_open)
                )
        # Traceback for small case
        i, j = n, m
        a_aln, b_aln = [], []
        while i > 0 or j > 0:
            if i > 0 and j > 0 and dp[i, j] == dp[i - 1, j - 1] + (match if a[i - 1] == b[j - 1] else mismatch):
                a_aln.append(a[i - 1])
                b_aln.append(b[j - 1])
                i -= 1
                j -= 1
            elif i > 0 and dp[i, j] == dp[i - 1, j] + gap_extend or dp[i, j] == dp[i - 1, j] + gap_open:
                a_aln.append(a[i - 1])
                b_aln.append('-')
                i -= 1
            else:
                a_aln.append('-')
                b_aln.append(b[j - 1])
                j -= 1
        return ''.join(reversed(a_aln)), ''.join(reversed(b_aln))

    aligned_s1, aligned_s2 = divide_and_conquer(s1, s2)
    score = affine_forward(s1, s2)[-1]

    return score, aligned_s1, aligned_s2

"""
    n = len(s1)
    m = len(s2)

    # DP matrices
    lower = np.full((n+1, m+1), -np.inf, dtype=float) # insertion (gap in s2)
    middle = np.full((n+1, m+1), -np.inf, dtype=float) # matches/mismatches
    upper = np.full((n+1, m+1), -np.inf, dtype=float) # deletion (gap in s1)

    # Pointer matrices (0=middle, 1=lower, 2=upper)
    ptr_middle = np.zeros((n+1, m+1), dtype=np.int8)
    ptr_lower = np.zeros((n+1, m+1), dtype=np.int8)
    ptr_upper = np.zeros((n+1, m+1), dtype=np.int8)

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
"""