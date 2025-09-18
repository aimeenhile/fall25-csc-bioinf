from dbg import DBG
from utils import read_data
import sys
import os
import traceback
from typing import List, Dict
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
from datetime import datetime

sys.setrecursionlimit(1000000)

def compute_N50_from_lengths(lengths: List[int]) -> str:
    """Compute N50 from a list of contig lengths (in memory)."""
    if not lengths:
        return "NA"
    lengths.sort(reverse=True)
    total_len = sum(lengths)
    cum_len = 0
    for l in lengths:
        cum_len += l
        if cum_len >= total_len / 2:
            return str(l)
    return "NA"

def process_dataset(dataset_path: str, dataset_name: str) -> Dict[str, str]:
    """Process one dataset: build DBG, generate contigs in memory, compute N50, and measure runtime."""
    start_time = time.time()
    try:
        print(f"Processing {dataset_name}...")
        # Read reads
        short1, short2, long1 = read_data(dataset_path)
        dbg = DBG(k=25, data_list=[short1, short2, long1])

        # Generate contigs in memory (up to 20 longest contigs)
        contigs: List[str] = []
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            contigs.append(c)


        # Compute N50 in memory
        contig_lengths = [len(c) for c in contigs]
        N50 = compute_N50_from_lengths(contig_lengths)

        end_time = time.time()
        runtime_sec = int(end_time - start_time)
        runtime_str = f"{runtime_sec // 3600}:{(runtime_sec % 3600) // 60:02d}:{runtime_sec % 60:02d}"

        # Print runtime immediately
        print(f"Processed {dataset_name} in {runtime_str}, N50={N50}")

        metrics: Dict[str, str] = {
            "Dataset": dataset_name,
            "Rank": "NA",  
            "Submission_Time": datetime.now().strftime("%Y/%m/%d %I:%M:%S%p"),
            "Submission_Count": "NA",
            "Genome_Fraction(%)": "NA",
            "Duplication ratio": "NA",
            "N50": N50,
            "Misassemblies": "NA",
            "Mismatches per 100kbp": "NA",
            "Runtime_sec": str(runtime_sec)
        }

    except Exception:
        print(f"Error processing {dataset_name}")
        traceback.print_exc()
        metrics = {
            "Dataset": dataset_name,
            "Rank": "NA",
            "Submission_Time": datetime.now().strftime("%Y/%m/%d %I:%M:%S%p"),
            "Submission_Count": "NA",
            "Genome_Fraction(%)": "NA",
            "Duplication ratio": "NA",
            "N50": "NA",
            "Misassemblies": "NA",
            "Mismatches per 100kbp": "NA",
            "Runtime_sec": "0"
        }

    return metrics

def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python main.py <data_root>")
        return

    data_root: str = sys.argv[1]
    results: List[Dict[str, str]] = []

    # List all datasets
    datasets: List[str] = sorted(
        [d for d in os.listdir(data_root) if os.path.isdir(os.path.join(data_root, d))]
    )

    # Parallel processing
    max_threads = 1
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        future_to_dataset = {
            executor.submit(process_dataset, os.path.join(data_root, d), d): d for d in datasets
        }

        for future in as_completed(future_to_dataset):
            try:
                metrics = future.result()
                results.append(metrics)
            except Exception:
                dataset_name = future_to_dataset[future]
                print(f"Unhandled error processing {dataset_name}")
                traceback.print_exc()

    # Rank by N50 descending (NA treated as 0)
    results.sort(key=lambda x: int(x["N50"]) if x["N50"] != "NA" else 0, reverse=True)
    for rank, res in enumerate(results, 1):
        res["Rank"] = str(rank)

    # Print Markdown table (compatible with evaluate.sh)
    header = [
        "Rank",
        "Dataset",
        "Submission_Time",
        "Submission_Count",
        "Genome_Fraction(%)",
        "Duplication ratio",
        "N50",
        "Misassemblies",
        "Mismatches per 100kbp",
        "Runtime_sec", 
    ]
    print("| " + " | ".join(header) + " |")
    print("|" + "|".join(["---"] * len(header)) + "|")
    for res in results:
        print(
            f"| {res['Rank']} | {res['Dataset']} | {res['Submission_Time']} | "
            f"{res['Submission_Count']} | {res['Genome_Fraction(%)']} | "
            f"{res['Duplication ratio']} | {res['N50']} | {res['Misassemblies']} | "
            f"{res['Mismatches per 100kbp']} | {res['Runtime_sec']} |"
        )

if __name__ == "__main__":
    main()