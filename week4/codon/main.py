from utils import read_data
from alignment import global_alignment, local_alignment, fitting_alignment, affine_alignment
from python import sys
from python import os
import time 
from typing import List, Tuple

# Scoring parameters
MATCH = 3
MISMATCH = -3
GAP = -2
GAP_OPEN = -5
GAP_EXTENSION = -1

def run_alignment(method: str, s1: str, s2: str) -> Tuple[int, str, str]:
    """Dispatch alignment based on method name."""
    if method == "global":
        return global_alignment(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP)
    elif method == "local":
        return local_alignment(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP)
    elif method == "semi_global":
        return fitting_alignment(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP)
    elif method == "affine":
        return affine_alignment(
            s1, s2, match=MATCH, mismatch=MISMATCH, gap_open=GAP_OPEN, gap_extend=GAP_EXTENSION
        )
    else:
        raise Exception(f"Unknown alignment method: {method}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: codon main.py <data_directory_path>")
        sys.exit(1)
        
    argv = sys.argv
    data_path = str(os.path.join('..', argv[1]))

    print(f"Loading data from: {data_path}")
    mt_human, mt_orang, q1, t1 = read_data(data_path)

    mt_human_cod: list[str] = list(mt_human)
    mt_orang_cod: list[str] = list(mt_orang)
    q1_cod: list[str] = list(q1)
    t1_cod: list[str] = list(t1)

    query = q1_cod 
    target = t1_cod

    # Datasets
    datasets: list[tuple[str, str, str]] = [("mt_human", mt_orang_cod[0], mt_human_cod[0])]
    for i in range(len(query)):
        datasets.append((f"q{i+1}", query[i], target[i]))

    # Alignment methods
    methods: list[str] = ["global", "local", "semi_global", "affine"]

    # Run all methods on all datasets and measure runtime
    for name, s1, s2 in datasets:
        for method in methods:
            start = time.perf_counter()
            score_val: int = 0
            align1: str = ""
            align2: str = ""

            try:
                score_val, align1, align2 = run_alignment(method, s1, s2)
            except Exception as e:
                print(f"Error running {method} on {name}: {e}")
                raise e
            
            end = time.perf_counter()
            runtime_ms: int = int((end - start) * 1000)  

            print(f"{method}-{name}\tcodon\t{runtime_ms}ms")

            # Print score and alignment
            print(f"Score: {score_val}  Length s1:{len(align1)} s2:{len(align2)}")
            print(f"Alignment 1: {align1}")
            print(f"Alignment 2: {align2}\n")