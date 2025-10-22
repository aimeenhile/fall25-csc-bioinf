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

    def nw_score(a, b):
        """Compute last row of Needleman-Wunsch DP in linear space."""
        n = len(a)
        m = len(b)
        prev = np.arange(m+1) * gap
        for i in range(1, n+1):
            curr = np.zeros(m+1, dtype=int)
            curr[0] = i * gap
            for j in range(1, m+1):
                diag = prev[j-1] + score(a[i-1], b[j-1], match, mismatch)
                up = prev[j] + gap
                left = curr[j-1] + gap
                curr[j] = max(diag, up, left)
            prev = curr
        return prev

    def hirschberg(a, b):
        """Recursively compute global alignment."""
        if len(a) == 0:
            return '-'*len(b), b
        if len(b) == 0:
            return a, '-'*len(a)
        if len(a) == 1 or len(b) == 1:
            # small enough to compute full DP
            n = len(a)
            m = len(b)
            D = np.zeros((n+1, m+1), dtype=int)
            P = np.zeros((n+1, m+1), dtype=int)
            for i in range(1, n+1):
                D[i,0] = i*gap
            for j in range(1, m+1):
                D[0,j] = j*gap
            for i in range(1, n+1):
                for j in range(1, m+1):
                    diag = D[i-1,j-1] + score(a[i-1], b[j-1], match, mismatch)
                    up = D[i-1,j] + gap
                    left = D[i,j-1] + gap
                    D[i,j] = max(diag, up, left)
            # Backtrace
            i,j = n,m
            align_a, align_b = [], []
            while i>0 or j>0:
                if i>0 and j>0 and D[i,j] == D[i-1,j-1] + score(a[i-1],b[j-1],match,mismatch):
                    align_a.append(a[i-1]); align_b.append(b[j-1]); i-=1; j-=1
                elif i>0 and D[i,j] == D[i-1,j] + gap:
                    align_a.append(a[i-1]); align_b.append('-'); i-=1
                else:
                    align_a.append('-'); align_b.append(b[j-1]); j-=1
            return ''.join(reversed(align_a)), ''.join(reversed(align_b))
        else:
            # Divide
            xlen = len(a)//2
            scoreL = nw_score(a[:xlen], b)
            scoreR = nw_score(a[xlen:][::-1], b[::-1])
            sums = scoreL + scoreR[::-1]
            k = np.argmax(sums)
            left_a, left_b = hirschberg(a[:xlen], b[:k])
            right_a, right_b = hirschberg(a[xlen:], b[k:])
            return left_a+right_a, left_b+right_b

    align1, align2 = hirschberg(s1, s2)
    
    # Compute final score
    final_score = 0
    for c1, c2 in zip(align1, align2):
        if c1 == '-' or c2 == '-':
            final_score += gap
        else:
            final_score += score(c1, c2, match, mismatch)
    
    return final_score, align1, align2


def local_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Substrings of s1 and s2 whose best global alignment score is maximized """
    n = len(s1)
    m = len(s2)

    # Compute last row and find max score position
    prev = np.zeros(m+1, dtype=int)
    curr = np.zeros(m+1, dtype=int)
    max_score = 0
    max_pos = (0,0)

    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = prev[j-1] + score(s1[i-1], s2[j-1], match, mismatch)
            up = prev[j] + gap
            left = curr[j-1] + gap
            curr[j] = max(0, diag, up, left)
            if curr[j] > max_score:
                max_score = curr[j]
                max_pos = (i, j)
        prev, curr = curr, prev

    # Recursive Hirschberg-style traceback for local alignment
    def traceback(si, sj):
        """Reconstruct alignment from position si, sj backwards to first 0 score."""
        if si == 0 or sj == 0:
            return "", ""
        # Recompute DP for submatrix ending at si,sj
        sub_s1 = s1[:si]
        sub_s2 = s2[:sj]
        dp = np.zeros((len(sub_s1)+1, len(sub_s2)+1), dtype=int)
        for i in range(1, len(sub_s1)+1):
            for j in range(1, len(sub_s2)+1):
                diag = dp[i-1,j-1] + score(sub_s1[i-1], sub_s2[j-1], match, mismatch)
                up = dp[i-1,j] + gap
                left = dp[i,j-1] + gap
                dp[i,j] = max(0, diag, up, left)
        # Start traceback
        i,j = len(sub_s1), len(sub_s2)
        while i>0 and j>0 and dp[i,j]>0:
            if dp[i,j] == dp[i-1,j-1] + score(sub_s1[i-1], sub_s2[j-1], match, mismatch):
                i-=1; j-=1
            elif dp[i,j] == dp[i-1,j] + gap:
                i-=1
            else:
                j-=1
        # Recurse if needed
        return traceback(i,j) + (sub_s1[i:], sub_s2[j:])

    # Get final alignment
    align1, align2 = traceback(*max_pos)
    align1 = ''.join(align1)
    align2 = ''.join(align2)

    return max_score, align1, align2
    

def fitting_alignment(s1: str, s2: str, match: int = MATCH, mismatch: int = MISMATCH, gap: int = GAP):
    """ Fitting alignment of string s1 to string s2 using BWT-inspired candidate positions """
    n = len(s1)
    m = len(s2)
    prev = np.zeros(m+1, dtype=int)
    curr = np.zeros(m+1, dtype=int)

    # Fill DP
    for i in range(1,n+1):
        for j in range(1,m+1):
            diag = prev[j-1]+score(s1[i-1],s2[j-1],match,mismatch)
            up = prev[j]+gap
            left = curr[j-1]+gap
            curr[j] = max(diag, up, left)
        prev, curr = curr, prev

    # Maximum in last column
    last_col = prev
    j = m
    i = np.argmax(last_col)
    max_score = last_col[i]

    # traceback (approximate, for memory)
    align1, align2 = [], []
    ci, cj = i, j
    while cj>0:
        if ci>0:
            align1.append(s1[ci-1])
            align2.append(s2[cj-1])
            ci -= 1
            cj -= 1
        else:
            align1.append('-')
            align2.append(s2[cj-1])
            cj -= 1
    return max_score, ''.join(reversed(align1)), ''.join(reversed(align2))


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

    middle[0,0] = 0.0
    lower[0,0] = 0.0
    upper[0,0] = 0.0

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
            i -= 1
            current = ptr_lower[i,j]
        elif current == 2:  
            # upper
            align1_list.append('-')
            align2_list.append(s2[j-1])
            j -= 1
            current = ptr_upper[i,j]

    align1 = ''.join(reversed(align1_list))
    align2 = ''.join(reversed(align2_list))

    return middle[n,m], align1, align2