from utils import read_data
from alignment import global_alignment, local_alignment, fitting_alignment, affine_alignment
from python import sys
from python import os
import time 
from typing import Callable, List, Tuple

# Scoring parameters
MATCH = 3
MISMATCH = -3
GAP = -2
GAP_OPEN = -5
GAP_EXTENSION = -1


def global_wrap(s1: str, s2: str):
    return global_alignment(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP)

def local_wrap(s1: str, s2: str):
    return local_alignment(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP)

def fitting_wrap(s1: str, s2: str):
    return fitting_alignment(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP)

def affine_wrap(s1: str, s2: str):
    return affine_alignment(s1, s2, match=MATCH, mismatch=MISMATCH, gap_open=GAP_OPEN, gap_extend=GAP_EXTENSION)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: codon main.py <data_directory_path>")
        sys.exit(1)
        
    argv = sys.argv
    data_path = str(os.path.join('..', argv[1]))

    print(f"Loading data from: {data_path}")
    mt_human, mt_orang, q1, t1 = read_data(data_path)

    query = q1 
    target = t1 

    # Datasets
    datasets = [("mt_human", mt_orang[0], mt_human[0])]
    for i in range(len(query)):
        datasets.append((f"q{i+1}", query[i], target[i]))

    # Define a unified function type
    AlignFunction = Callable[[str, str], Tuple[int, str, str]]

    # Alignment methods
    methods: list = [
        ("global", global_wrap),
        ("local", local_wrap),
        ("semi_global", fitting_wrap),
        ("affine", affine_wrap)
    ]

    # Run all methods on all datasets and measure runtime
    for name, s1, s2 in datasets:
        for method, align_function in methods:
            start = time.perf_counter()
            score_val: int = 0
            align1: str = ""
            align2: str = ""

            try:
                if method == "affine":
                    score_val, align1, align2 = align_function(s1, s2, match=MATCH, mismatch=MISMATCH, gap_open=GAP_OPEN, gap_extend=GAP_EXTENSION)
                else:
                    score_val, align1, align2 = align_function(s1, s2, match=MATCH, mismatch=MISMATCH, gap=GAP)
            except Exception as e:
                print(f"Error running {method} on {name}: {e}")
                continue

            end = time.perf_counter()
            runtime_ms: int = int((end - start) * 1000)  

            print(f"{method}-{name}\tcodon\t{runtime_ms}ms")

            # Print score and alignment
            print(f"Score: {score_val}  Length s1:{len(align1)} s2:{len(align2)}")
            print(f"Alignment 1: {align1}")
            print(f"Alignment 2: {align2}\n")