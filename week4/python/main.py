from utils import read_data
from alignment import global_alignment, local_alignment, fitting_alignment, affine_alignment
import sys
import os
import time

sys.setrecursionlimit(2000000)

# Scoring parameters
MATCH = 3
MISMATCH = -3
GAP = -2
GAP_OPEN = -5
GAP_EXTENSION = -1

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python main.py <data_directory_path>")
        sys.exit(1)
        
    argv = sys.argv
    data_path = os.path.join(os.path.dirname(__file__), '..', argv[1])

    print(f"Loading data from: {data_path}")
    mt_human, mt_orang, q1, t1 = read_data(data_path)

    query = q1 
    target = t1 

    # Datasets
    datasets = [("mt_human", mt_human[0], mt_orang[0])]
    for i in range(len(query)):
        datasets.append((f"q{i+1}", query[i], target[i]))

    # Alignment methods
    methods = [
        ("global", global_alignment),
        ("local", local_alignment),
        ("semi_global", fitting_alignment),
        ("affine", affine_alignment)
    ]

    # Run all methods on all datasets and measure runtime
    for name, s1, s2 in datasets:
        for method, align_function in methods:
            start = time.perf_counter()
            score_val, align1, align2 = (0, "", "")

            try:
                if method == "affine":
                    score_val, align1, align2 = align_function(s1, s2, match=MATCH, mismatch=MISMATCH, gap_open=GAP_OPEN, gap_extend=GAP_EXTENSION)
                else:
                    score_val, align1, align2 = align_function(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP)
            except Exception as e:
                print(f"Error running {method} on {name}: {e}")
                # Print stack trace for debugging
                import traceback
                traceback.print_exc()
                continue

            end = time.perf_counter()
            runtime_ms = int((end - start) * 1000)  

            print(f"{method}-{name}\tpython\t{runtime_ms}ms")

            # Print score and alignment
            print(f"Score: {score_val}")
            print(f"Alignment 1: {align1}")
            print(f"Alignment 2: {align2}\n")
    