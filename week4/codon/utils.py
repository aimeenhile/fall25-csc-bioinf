from python import os

def os_path_join(*args):
    return '/' .join(args)
with open(str(os.path.join(path, name)), 'r') as f: ...

def read_fasta(path: str, name: str) -> list[str]:
    """
    Reads a FASTA file and returns a list of sequences (strings).
    Each sequence corresponds to one FASTA entry.
    """
    data: list[str] = []
    sequence: str = ""

    with open(str(os.path.join(path, name)), 'r') as f:
        for line in f:
            line = str(line.strip())
            if line[0] == ">":
                if sequence:
                    data.append(sequence)
                    sequence = ""
            else:
                sequence += line
        if sequence:  # save last sequence
            data.append(sequence)

    print(name, len(data), [len(s) for s in data])

    return data


def read_data(path):

    mt_human = read_fasta(path, "MT-human.fa")
    mt_orang = read_fasta(path, "MT-orang.fa")
    q1 = read_fasta(path, "q1.fa")
    t1 = read_fasta(path, "t1.fa")
    
    return mt_human, mt_orang, q1, t1