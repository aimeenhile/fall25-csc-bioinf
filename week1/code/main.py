from dbg_kmer_as_key import DBG
from utils import read_data
import sys
from typing import List, Dict
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

def format_hms(seconds: float) -> str:
    """Convert seconds to h:mm:ss format."""
    seconds = int(seconds)
    hh = seconds // 3600
    mm = (seconds % 3600) // 60
    ss = seconds % 60
    return f"{hh}:{mm:02}:{ss:02}"

def process_dataset(dataset_path: str, dataset_name: str) -> Dict[str, str]:
    """Process a dataset entirely in memory and return metrics."""
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
        }

    except Exception as e:
        print(f"Error processing {dataset_name}: {e}")
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
        }
    end_time = time.time()
    runtime_hms = format_hms(end_time - start_time)
    metrics["Runtime"] = runtime_hms

    # Print per-dataset runtime immediately
    print(f"{dataset_name} completed in {runtime_hms}")

    return metrics

def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python main.py <data_root>")
        return

    data_root: str = sys.argv[1]
    datasets: List[str] = sorted(["data1", "data2", "data3", "data4"])
    results: List[Dict[str, str]] = []

    for dataset in datasets:
        dataset_path: str = f"{data_root}/{dataset}"
        metrics = process_dataset(dataset_path, dataset)
        results.append(metrics)

    ## Rank by N50 descending (NA treated as 0)
    results.sort(key=lambda x: int(x["N50"]) if x["N50"] != "NA" else 0, reverse=True)
    for rank, res in enumerate(results, 1):
        res["Rank"] = str(rank)

    # Print table 
    header: List[str] = [
        "Rank", "Dataset", "Submission_Time", "Submission_Count",
        "Genome_Fraction(%)", "Duplication ratio", "N50",
        "Misassemblies", "Mismatches per 100kbp", "Runtime"
    ]
    
    print("| " + " | ".join(header) + " |")
    print("|" + "|".join(["---"]*len(header)) + "|")
    for res in results:
        print(
            f"| {res['Rank']} | {res['Dataset']} | {res['Submission_Time']} | {res['Submission_Count']} | "
            f"{res['Genome_Fraction(%)']} | {res['Duplication ratio']} | {res['N50']} | "
            f"{res['Misassemblies']} | {res['Mismatches per 100kbp']} | {res['Runtime']} |"
        )

if __name__ == "__main__":
    main()