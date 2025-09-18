from dbg import DBG
from utils import read_data
from typing import Optional
import sys

def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: codon run -release code/main.codon.py <data_path>")
        return

    data_path = sys.argv[1]
    data_list = read_data(data_path)

    k = 25
    dbg = DBG(k, data_list)

    print("Longest contigs from DBG:")
    while True:
        contig: Optional[str] = dbg.get_longest_contig()
        if contig is None:
            break
        print(contig)

if __name__ == "__main__":
    main()

