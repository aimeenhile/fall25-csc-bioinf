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

    n = len(s1)
    m = len(s2)

    lower_prev = np.full(m+1, -np.inf, dtype=float)
    lower_curr = np.full(m+1, -np.inf, dtype=float)
    upper_prev = np.full(m+1, -np.inf, dtype=float)

    upper_curr = np.full(m+1, -np.inf, dtype=float)
    middle_prev = np.full(m+1, -np.inf, dtype=float)
    middle_curr = np.full(m+1, -np.inf, dtype=float)

    # Initialization
    middle_prev[0] = 0.0
    for j in range(1, m+1):
        upper_prev[j] = gap_open + (j-1)*gap_extend
        middle_prev[j] = -np.inf
        lower_prev[j] = -np.inf

    # Keep traceback pointers for small window
    trace_middle = np.zeros((n+1, m+1), dtype=np.int8)
    trace_lower = np.zeros((n+1, m+1), dtype=np.int8)
    trace_upper = np.zeros((n+1, m+1), dtype=np.int8)

    for i in range(1, n+1):
        s1_i = s1[i-1]
        # first column
        lower_curr[0] = gap_open + (i-1)*gap_extend
        middle_curr[0] = -np.inf
        upper_curr[0] = -np.inf

        for j in range(1, m+1):
            s2_j = s2[j-1]

            # lower (vertical gap)
            val1 = lower_prev[j] + gap_extend
            val2 = middle_prev[j] + gap_open
            if val1 >= val2:
                lower_curr[j] = val1
                trace_lower[i,j] = 1  # came from lower
            else:
                lower_curr[j] = val2
                trace_lower[i,j] = 0  # came from middle

            # upper (horizontal gap)
            val1 = upper_curr[j-1] + gap_extend
            val2 = middle_curr[j-1] + gap_open
            if val1 >= val2:
                upper_curr[j] = val1
                trace_upper[i,j] = 2  # came from upper
            else:
                upper_curr[j] = val2
                trace_upper[i,j] = 0  # came from middle

            # middle (match/mismatch)
            diag = middle_prev[j-1] + (match if s1_i == s2_j else mismatch)
            low = lower_curr[j]
            up = upper_curr[j]
            middle_curr[j] = max(diag, low, up)

            if middle_curr[j] == diag:
                trace_middle[i,j] = 0
            elif middle_curr[j] == low:
                trace_middle[i,j] = 1
            else:
                trace_middle[i,j] = 2

        # swap rows
        lower_prev, lower_curr = lower_curr, lower_prev
        upper_prev, upper_curr = upper_curr, upper_prev
        middle_prev, middle_curr = middle_curr, middle_prev

    # final score
    final_score = max(middle_prev[m], lower_prev[m], upper_prev[m])

    # Traceback
    i = n
    j = m
    # Starting state for traceback 
    if final_score == middle_prev[m]:
        # middle
        current = 0
    elif final_score == lower_prev[m]:
        # lower
        current = 1
    else:
        # upper
        current = 2

    align1_list = []
    align2_list = []

    while i > 0 or j > 0:
        if current == 0:
            # middle
            if trace_middle[i,j] == 0:
                align1_list.append(s1[i-1])
                align2_list.append(s2[j-1])
                i -= 1
                j -= 1
                current = 0
            elif trace_middle[i,j] == 1:
                current = 1
            else:
                current = 2
        elif current == 1:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            current = trace_lower[i,j]
            i -= 1
        elif current == 2:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            current = trace_upper[i,j]
            j -= 1

        if i == 0 and j > 0:
            align1_list.extend(['-'] * j)
            align2_list.extend(list(s2[:j][::-1]))
            break
        if j == 0 and i > 0:
            align1_list.extend(list(s1[:i][::-1]))
            align2_list.extend(['-']*i)
            break

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    return final_score, align1, align2