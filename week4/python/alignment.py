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