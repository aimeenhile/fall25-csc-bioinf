from dbg import DBG
from utils import read_data
from typing import List, Optional
import sys
import fs

def compute_N50(contig_file: str) -> str:
    lengths: List[int] = []
    seq_len: int = 0

    with open(contig_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_len > 0:
                    lengths.append(seq_len)
                seq_len = 0
            else:
                seq_len += len(line)
        if seq_len > 0:
            lengths.append(seq_len)

    if not lengths:
        return "NA"

    lengths.sort(reverse=True)
    total_len: int = sum(lengths)
    cum_len: int = 0
    for l in lengths:
        cum_len += l
        if cum_len >= total_len // 2:
            return str(l)
    return "NA"


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: codon run -release code/main.codon.py <data_path>")
        return

    data_root: str = sys.argv[1]
    results: List[Dict[str, str]] = []

    datasets: List[str] = fs.listdir(data_root)  
    datasets.sort()

    for dataset in datasets:
        dataset_path: str = fs.joinpath(data_root, dataset)
        if not fs.isdir(dataset_path):
            continue

        print(f"Processing {dataset}...")
        short1, short2, long1 = read_data(dataset_path)
        dbg = DBG(k=25, data_list=[short1, short2, long1])

        contig_file: str = fs.joinpath(dataset_path, "contig.fasta")
        with open(contig_file, "w") as f:
            for i in range(20):
                c: Optional[str] = dbg.get_longest_contig()
                if c is None:
                    break
                f.write(f">contig_{i}\n{c}\n")

        N50: str = compute_N50(contig_file)

        results.append({
            "Dataset": dataset,
            "Genome_Fraction(%)": "NA",
            "Duplication ratio": "NA",
            "N50": N50,
            "Misassemblies": "NA",
            "Mismatches per 100kbp": "NA"
        })

    # Rank by N50 descending (NA treated as 0)
    results.sort(key=lambda x: int(x["N50"]) if x["N50"] != "NA" else 0, reverse=True)
    for rank, res in enumerate(results, 1):
        res["Rank"] = rank

    # Print Markdown table
    header: List[str] = ["Rank", "Dataset", "Genome_Fraction(%)", "Duplication ratio", "N50", "Misassemblies", "Mismatches per 100kbp"]
    print("| " + " | ".join(header) + " |")
    print("|" + "|".join(["---"]*len(header)) + "|")
    for res in results:
        print(f"| {res['Rank']} | {res['Dataset']} | {res['Genome_Fraction(%)']} | "
              f"{res['Duplication ratio']} | {res['N50']} | {res['Misassemblies']} | "
              f"{res['Mismatches per 100kbp']} |")
if __name__ == "__main__":
    main()

