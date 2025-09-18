import copy
#from matplotlib import pyplot as plt


def reverse_complement(key):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    key = list(key[::-1])
    for i in range(len(key)):
        key[i] = complement[key[i]]
    return ''.join(key)


class Node:
    def __init__(self, kmer):
        self._children = set()
        self._count = 0
        self.kmer = kmer
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def add_child(self, kmer):
        self._children.add(kmer)

    def increase(self):
        self._count += 1

    def reset(self):
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def get_count(self):
        return self._count

    def get_children(self):
        return list(self._children)

    def remove_children(self, target):
        self._children = self._children - target


class DBG:
    def __init__(self, k, data_list):
        self.k = k
        self.nodes = {}
        # private
        self.kmer2idx = {}
        self.kmer_count = 0
        # build
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list):
        # check data list
        assert len(data_list) > 0
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list):
        for data in data_list:
            for original in data:
                rc = reverse_complement(original)
                for i in range(len(original) - self.k - 1):
                    self._add_arc(original[i: i + self.k], original[i + 1: i + 1 + self.k])
                    self._add_arc(rc[i: i + self.k], rc[i + 1: i + 1 + self.k])

    def show_count_distribution(self):
        count = [0] * 30
        for idx in self.nodes:
            count[self.nodes[idx].get_count()] += 1
        print(count[0:10])
        # plt.plot(count)
        # plt.show()

    def _add_node(self, kmer):
        if kmer not in self.kmer2idx:
            self.kmer2idx[kmer] = self.kmer_count
            self.nodes[self.kmer_count] = Node(kmer)
            self.kmer_count += 1
        idx = self.kmer2idx[kmer]
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1, kmer2):
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, child):
        return self.nodes[child].get_count()

    def _get_sorted_children(self, idx):
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, idx):
        if not self.nodes[idx].visited:
            self.nodes[idx].visited = True
            children = self._get_sorted_children(idx)
            max_depth, max_child = 0, None
            for child in children:
                depth = self._get_depth(child)
                if depth > max_depth:
                    max_depth, max_child = depth, child
            self.nodes[idx].depth, self.nodes[idx].max_depth_child = max_depth + 1, max_child
        return self.nodes[idx].depth

    def _reset(self):
        for idx in self.nodes.keys():
            self.nodes[idx].reset()

    def _get_longest_path(self):
        max_depth, max_idx = 0, None
        for idx in self.nodes.keys():
            depth = self._get_depth(idx)
            if depth > max_depth:
                max_depth, max_idx = depth, idx

        path = []
        while max_idx is not None:
            path.append(max_idx)
            max_idx = self.nodes[max_idx].max_depth_child
        return path

    def _delete_path(self, path):
        for idx in path:
            del self.nodes[idx]
        path_set = set(path)
        for idx in self.nodes.keys():
            self.nodes[idx].remove_children(path_set)

    def _concat_path(self, path):
        if len(path) < 1:
            return None
        concat = copy.copy(self.nodes[path[0]].kmer)
        for i in range(1, len(path)):
            concat += self.nodes[path[i]].kmer[-1]
        return concat

    def get_longest_contig(self):
        # reset params in nodes for getting longest path
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig
