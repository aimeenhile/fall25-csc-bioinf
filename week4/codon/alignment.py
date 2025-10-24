import numpy as np
from numpy import ndarray, int8, float64

MATCH: int = 3
MISMATCH: int = -3
GAP: int = -2
GAP_OPEN: int = -5
GAP_EXTENSION: int = -1
      

def global_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Find an alignment with maximum alignment score """
    n: int = len(s1)
    m: int = len(s2)

    D_prev: ndarray[int, 1] = np.ndarray((m+1,), dtype=int8, ndim=1)
    D_curr: ndarray[int, 1] = np.ndarray((m+1,), dtype=int8, ndim=1)
    for j in range(m+1):
        D_prev[j] = 0
        D_curr[j] = 0

    P: ndarray[int8, 2] = np.ndarray((n+1, m+1), dtype=int8, ndim=2)
    for i in range(n+1):
        for j in range(m+1):
            P[i, j] = -1
            
    # Initialization
    for j in range(1, m+1):
        D_prev[j] = j * gap
        P[0][j] = 2  # left

    # Fill DP
    for i in range(1, n+1):
        D_curr[0] = i * gap
        P[i][0] = 1  # up
        s1_i: str = s1[i-1]
        for j in range(1, m+1):
            s2_j: str = s2[j-1]
            diag = D_prev[j-1] + (match if s1_i == s2_j else mismatch)
            up = D_prev[j] + gap
            left = D_curr[j-1] + gap

            val = diag
            pos: int = 0
            if up > val:
                val = up
                pos = 1
            if left > val:
                val = left
                pos = 2

            D_curr[j] = val
            P[i][j] = pos

        tmp = D_prev
        D_prev = D_curr
        D_curr = tmp

    # Backtrace
    i: int = n
    j: int = m
    align1_list: list[str] = []
    align2_list: list[str] = []

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            pos: int = P[i][j]
            if pos == 0:
                align1_list.append(s1[i-1])
                align2_list.append(s2[j-1])
                i -= 1
                j -= 1
            elif pos == 1:
                align1_list.append(s1[i-1])
                align2_list.append('-')
                i -= 1
            else:
                align1_list.append('-')
                align2_list.append(s2[j-1])
                j -= 1
        elif i > 0:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            i -= 1
        elif j > 0:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            j -= 1

    align1_list.reverse()
    align2_list.reverse()
    align1: str = ''.join(align1_list)
    align2: str = ''.join(align2_list)

    max_score: int = D_prev[m]

    return max_score, align1, align2


def local_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Substrings of s1 and s2 whose best global alignment score is maximized """
    n: int = len(s1)
    m: int = len(s2)
    
    V_prev: list[int] = [0] * (m + 1)
    V_curr: list[int] = [0] * (m + 1)

    P: list[list[int]] = [[-1 for _ in range(m+1)] for _ in range(n+1)]

    max_score: int = 0
    max_i: int = 0
    max_j: int = 0

    for i in range(1, n+1):
        V_curr[0] = 0
        s1_i: str = s1[i-1]
        for j in range(1, m+1):
            s2_j: str = s2[j-1]
            diag = V_prev[j-1] + (match if s1_i == s2_j else mismatch)
            up = V_prev[j] + gap
            left = V_curr[j-1] + gap

            val = diag
            pos: int = 0
            if up > val:
                val = up
                pos = 1
            if left > val:
                val = left
                pos = 2
            if val < 0:
                val = 0
                pos = -1

            V_curr[j] = val
            P[i][j] = pos

            if val > max_score:
                max_score = val
                max_i = i
                max_j = j

        tmp = V_prev
        V_prev = V_curr
        V_curr = tmp

    # Backtrace
    i: int = max_i
    j: int = max_j
    align1_list: list[str] = []
    align2_list: list[str] = []

    while i > 0 and j > 0 and P[i][j] != -1:
        pos: int = P[i][j]
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

    align1_list.reverse()
    align2_list.reverse()
    align1: str = ''.join(align1_list)
    align2: str = ''.join(align2_list)

    return max_score, align1, align2

def fitting_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    n: int = len(s1)
    m: int = len(s2)

    prev: list[int] = [0] * (m + 1)
    curr: list[int] = [0] * (m + 1)

    P: list[list[int]] = [[-1 for _ in range(m+1)] for _ in range(n+1)]

    # Fill DP
    for i in range(1, n + 1):
        s1_i: str = s1[i-1]
        curr[0] = prev[0] + gap
        P[i][0] = 1 # up pointer

        for j in range(1, m + 1):
            s2_j: str = s2[j-1]
            diag = prev[j-1] + (match if s1_i == s2_j else mismatch)
            up = prev[j] + gap
            left = curr[j-1] + gap

            val = diag
            pos: int = 0
            if up > val:
                val = up
                pos = 1
            if left > val:
                val = left
                pos = 2

            curr[j] = val
            P[i][j] = pos

        tmp = prev
        prev = curr
        curr = tmp

    # Traceback start
    j: int = 0
    max_score: int = prev[0]
    for jj in range(1, m + 1):
        if prev[jj] > max_score:
            max_score = prev[jj]
            j = jj
    i: int = n

    # Traceback
    align1_list: list[str] = []
    align2_list: list[str] = []

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            pos: int = P[i][j]
            if pos == 0:
                align1_list.append(s1[i-1])
                align2_list.append(s2[j-1])
                i -= 1
                j -= 1
            elif pos == 1:
                align1_list.append(s1[i-1])
                align2_list.append('-')
                i -= 1
            else:
                align1_list.append('-')
                align2_list.append(s2[j-1])
                j -= 1
        elif i > 0:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            i -= 1
        elif j > 0:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            j -= 1

    align1_list.reverse()
    align2_list.reverse()
    align1: str = ''.join(align1_list)
    align2: str = ''.join(align2_list)

    return max_score, align1, align2

def affine_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap_open: int = GAP_OPEN, gap_extend: int = GAP_EXTENSION):
    n: int = len(s1)
    m: int = len(s2)

    NEG_INF: float = -1e9

    lower_prev: list[float] = [NEG_INF] * (m + 1)
    lower_curr: list[float] = [NEG_INF] * (m + 1)
    upper_prev: list[float] = [NEG_INF] * (m + 1)
    upper_curr: list[float] = [NEG_INF] * (m + 1)
    middle_prev: list[float] = [NEG_INF] * (m + 1)
    middle_curr: list[float] = [NEG_INF] * (m + 1)

    middle_prev[0] = 0.0
    for j in range(1, m + 1):
        upper_prev[j] = gap_open + (j - 1) * gap_extend

    trace_lower: list[list[int]] = [[0 for _ in range(m+1)] for _ in range(n+1)]
    trace_upper: list[list[int]] = [[0 for _ in range(m+1)] for _ in range(n+1)]
    trace_middle: list[list[int]] = [[0 for _ in range(m+1)] for _ in range(n+1)]

    # fill DP
    for i in range(1, n+1):
        s1_i: str = s1[i-1]
        lower_curr[0] = gap_open + (i-1)*gap_extend
        middle_curr[0] = NEG_INF
        upper_curr[0] = NEG_INF

        for j in range(1, m+1):
            s2_j: str = s2[j-1]

            # lower
            val1 = lower_prev[j] + gap_extend
            val2 = middle_prev[j] + gap_open
            if val1 >= val2:
                lower_curr[j] = val1
                trace_lower[i][j] = 1
            else:
                lower_curr[j] = val2
                trace_lower[i][j] = 0

            # upper
            val1 = upper_curr[j-1] + gap_extend
            val2 = middle_curr[j-1] + gap_open
            if val1 >= val2:
                upper_curr[j] = val1
                trace_upper[i][j] = 2
            else:
                upper_curr[j] = val2
                trace_upper[i][j] = 0

            # middle
            diag = middle_prev[j-1] + (match if s1_i == s2_j else mismatch)
            low = lower_curr[j]
            up = upper_curr[j]
            best = diag
            pos: int = 0
            if low > best:
                best = low
                pos = 1
            if up > best:
                best = up
                pos = 2
            middle_curr[j] = best
            trace_middle[i][j] = pos

        tmp = lower_prev
        lower_prev = lower_curr
        lower_curr = tmp

        tmp = upper_prev
        upper_prev = upper_curr
        upper_curr = tmp

        tmp = middle_prev
        middle_prev = middle_curr
        middle_curr = tmp

    # final score
    final_score: float = middle_prev[m]
    if lower_prev[m] > final_score:
        final_score = lower_prev[m]
        current = 1
    if upper_prev[m] > final_score:
        final_score = upper_prev[m]
        current = 2
    else:
        current = 0

    i: int = n
    j: int = m
    align1_list: list[str] = []
    align2_list: list[str] = []

    while i > 0 or j > 0:
        if current == 0:
            tm = int(trace_middle[i][j])
            if tm == 0:
                align1_list.append(s1[i-1])
                align2_list.append(s2[j-1])
                i -= 1
                j -= 1
                current = 0
            elif tm == 1:
                current = 1
            else:
                current = 2
        elif current == 1:
            align1_list.append(s1[i-1])
            align2_list.append('-')
            current = trace_lower[i][j]
            i -= 1
        else:
            align1_list.append('-')
            align2_list.append(s2[j-1])
            current = trace_upper[i][j]
            j -= 1

        if i == 0 and j > 0:
            for jj in range(j - 1, -1, -1):
                align1_list.append('-')
                align2_list.append(s2[jj])
            break
        if j == 0 and i > 0:
            for ii in range(i - 1, -1, -1):
                align1_list.append(s1[ii])
                align2_list.append('-')
            break

    align1_list.reverse()
    align2_list.reverse()
    align1: str = ''.join(align1_list)
    align2: str = ''.join(align2_list)

    return final_score, align1, align2