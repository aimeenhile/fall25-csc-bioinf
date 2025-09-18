from dbg import DBG
from utils import read_data
import sys
import os
import traceback
import subprocess

sys.setrecursionlimit(1000000)

# === CONFIG ===
DATA_ROOT = '../data'       # path to data folders (data1, data2, ...)
K = 25
TOP_CONTIGS = 20
QUAST_PATH = 'quast'        # ensure quast is in PATH

def run_assembler(dataset_path):
    """Build DBG and output contigs."""
    short1, short2, long1 = read_data(dataset_path)
    dbg = DBG(k=K, data_list=[short1, short2, long1])

    contig_file = os.path.join(dataset_path, 'contig.fasta')
    with open(contig_file, 'w') as f:
        for i in range(TOP_CONTIGS):
            c = dbg.get_longest_contig()
            if c is None:
                break
            f.write(f'>contig_{i}\n')
            f.write(c + '\n')
    return contig_file

def run_quast(contig_file, dataset_name):
    """Run QUAST on contigs and return parsed stats."""
    quast_out = os.path.join(os.path.dirname(contig_file), 'quast_output')
    cmd = [
        QUAST_PATH,
        '-o', quast_out,
        '--threads', '4',
        contig_file
    ]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # parse report.tsv
    report_file = os.path.join(quast_out, 'report.tsv')
    stats = {}
    if os.path.exists(report_file):
        with open(report_file, 'r') as f:
            headers = f.readline().strip().split('\t')
            values = f.readline().strip().split('\t')
            stats = dict(zip(headers, values))
    stats['Dataset'] = dataset_name
    return stats

def compute_rank(stats_list):
    """Rank contigs by N50 descending (higher N50 -> better)."""
    sorted_list = sorted(stats_list, key=lambda x: float(x.get('N50', 0)), reverse=True)
    for i, s in enumerate(sorted_list, 1):
        s['Rank'] = i
    return sorted_list

def print_markdown_table(stats_list):
    """Print the final Markdown table."""
    headers = ['Rank','Dataset','Submission Count','Genome Fraction (%)','Duplication ratio','N50','Misassemblies','Mismatches per 100 kbp']
    print('| ' + ' | '.join(headers) + ' |')
    print('|' + '---|'*len(headers))
    for s in stats_list:
        print('| {Rank} | {Dataset} | {# contigs} | {Genome Fraction} | {Duplication ratio} | {N50} | {# misassemblies} | {# mismatches per 100 kbp} |'.format(
            Rank=s.get('Rank', ''),
            Dataset=s.get('Dataset',''),
            **s
        ))

if __name__ == "__main__":
    try:
        if not os.path.exists(DATA_ROOT):
            raise FileNotFoundError(f"Data root '{DATA_ROOT}' not found. Please adjust DATA_ROOT path.")
        
        stats_all = []
        for dataset in sorted(os.listdir(DATA_ROOT)):
            dataset_path = os.path.join(DATA_ROOT, dataset)
            if not os.path.isdir(dataset_path):
                continue
            print(f'Processing {dataset} ...')
            contig_file = run_assembler(dataset_path)
            stats = run_quast(contig_file, dataset)
            # add submission count as number of contigs
            stats['# contigs'] = sum(1 for line in open(contig_file) if line.startswith('>'))
            stats_all.append(stats)

        ranked_stats = compute_rank(stats_all)
        print('\n## Genome Assembly Results\n')
        print_markdown_table(ranked_stats)

    except Exception:
        traceback.print_exc()
        sys.exit(1)
