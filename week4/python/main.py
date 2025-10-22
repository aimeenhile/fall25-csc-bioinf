from utils import read_data
from alignment import global_alignment, local_alignment, fitting_alignment, affine_alignment
import sys
import os
import time

ALIGNMENTS = {
    "global": global_alignment,
    "local": local_alignment,
    "semi-global": fitting_alignment,
    "affine-global": affine_alignment
}

if __name__ == "__main__":
    argv = sys.argv
    mt_human, mt_orang, q1, t1 = read_data(os.path.join('./', argv[1]))

    query = q1 
    target = t1 

    # Datasets
    datasets = [("mt_human", mt_human, mt_orang)]
    for i in range(len(queries)):
        datasets.append((f"q{i+1}", query[i], target[i]))

    # Alignment methods
    methods = [
        ("global", global_alignment),
        ("local", local_alignment),
        ("semi_global", fitting_alignment),
        ("affine", affine_alignment)
    ]

    # Scoring parameters
    MATCH = 3
    MISMATCH = -3
    GAP = -2
    GAP_OPEN = -5
    gAP_EXTENSION = -1

    # Run all methods on all datasets and measure runtime
    for name, s1, s2 in datasets:
        for method, align in methods:
            start = time.time()
            if method == "affine":
                alignment = align(s1, s2, match=MATCH, mismatch=MISMATCH, gap_open=GAP_OPEN, gap_extend=gAP_EXTENSION)
            else:
                alignment = align(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP)
            end = time.time()
            runtime_ms = int((end - start) * 1000)

            print(f"{method}-{align}\tpython\t{runtime_ms}ms")
    