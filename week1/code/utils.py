from typing import List

def read_fasta(path: str, name: str) -> List[str]:
    data: List[str] = []
    seq: str = ""
    full_path: str = path + "/" + name  
    with open(full_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    data.append(seq)
                    seq = ""
            else:
                seq += line
        if seq:
            data.append(seq)
    print(name, len(data), len(data[0]) if data else 0)
    return data

def read_data(path: str) -> List[List[str]]:
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1 = read_fasta(path, "long.fasta")
    return [short1, short2, long1]
