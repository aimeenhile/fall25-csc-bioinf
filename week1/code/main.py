from dbg import DBG
from utils import read_data
import sys
import os
import traceback

sys.setrecursionlimit(1000000)

def compute_N50(contig_file):
    """
    Compute N50 from contig lengths 
    """
    lengths = []
    with open(contig_file) as f:
        seq_len = 0
        for line in f:
            if line.startswith(">"):
                if seq_len > 0:
                    lengths.append(seq_len)
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if seq_len > 0:
            lengths.append(seq_len)

    if not lengths:
        return "NA"

    lengths.sort(reverse=True)
    total_len = sum(lengths)
    cum_len = 0
    for l in lengths:
        cum_len += l
        if cum_len >= total_len / 2:
            return l
    return "NA"

if __name__ == "__main__":
    try:
        data_root = os.path.join("./", "data")
        results = []

        for dataset in sorted(os.listdir(data_root)):
            dataset_path = os.path.join(data_root, dataset)
            if not os.path.isdir(dataset_path):
                continue

            print(f"Processing {dataset}...")
            short1, short2, long1 = read_data(dataset_path)
            dbg = DBG(k=25, data_list=[short1, short2, long1])

            contig_file = os.path.join(dataset_path, "contig.fasta")
            with open(contig_file, "w") as f:
                for i in range(20):
                    c = dbg.get_longest_contig()
                    if c is None:
                        break
                    f.write(f">contig_{i}\n{c}\n")

            # Compute N50 (reproducible)
            N50 = compute_N50(contig_file)

            # Other metrics require a reference genome, mark as NA
            metrics = {
                "Genome_Fraction(%)": "NA",
                "Duplication ratio": "NA",
                "N50": N50,
                "Misassemblies": "NA",
                "Mismatches per 100kbp": "NA"
            }

            results.append({
                "Dataset": dataset,
                **metrics
            })

        # Rank by N50 descending (NA treated as 0)
        results.sort(key=lambda x: float(x["N50"]) if x["N50"] != "NA" else 0, reverse=True)
        for rank, res in enumerate(results, 1):
            res["Rank"] = rank

        # Print Markdown table
        header = ["Rank", "Dataset", "Genome_Fraction(%)", "Duplication ratio", "N50", "Misassemblies", "Mismatches per 100kbp"]
        print("| " + " | ".join(header) + " |")
        print("|" + "|".join(["---"]*len(header)) + "|")
        for res in results:
            print(f"| {res['Rank']} | {res['Dataset']} | {res['Genome_Fraction(%)']} | {res['Duplication ratio']} | {res['N50']} | {res['Misassemblies']} | {res['Mismatches per 100kbp']} |")

    except Exception:
        traceback.print_exc()
        sys.exit(1)
