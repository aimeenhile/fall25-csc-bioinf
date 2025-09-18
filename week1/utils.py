import os


def read_fasta(path, name):
    data = []
    with open(os.path.join(path, name), 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line[0] != '>':
                data.append(line)
    print(name, len(data), len(data[0]))
    # print('Sample:', data[0])
    return data


def read_data(path):
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1 = read_fasta(path, "long.fasta")
    return short1, short2, long1
