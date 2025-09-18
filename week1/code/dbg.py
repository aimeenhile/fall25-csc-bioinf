from typing import List, Dict, Optional, Set

def reverse_complement(seq: str) -> str:
    complement: Dict[str, str] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[c] for c in reversed(seq))

class Node:
    kmer: str
    _children: Set[int]
    _count: int
    visited: bool
    depth: int
    max_depth_child: Optional[int]

    def __init__(self, kmer: str) -> None:
        self.kmer = kmer
        self._children = set()
        self._count = 0
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def add_child(self, idx: int) -> None:
        self._children.add(idx)

    def increase(self) -> None:
        self._count += 1

    def reset(self) -> None:
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def get_count(self) -> int:
        return self._count

    def get_children(self) -> List[int]:
        return list(self._children)

    def remove_children(self, target: Set[int]) -> None:
        self._children -= target

class DBG:
    k: int
    nodes: Dict[int, Node]
    kmer2idx: Dict[str, int]
    kmer_count: int

    def __init__(self, k: int, data_list: List[List[str]]) -> None:
        self.k = k
        self.nodes = {}
        self.kmer2idx = {}
        self.kmer_count = 0
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list: List[List[str]]) -> None:
        if len(data_list) == 0 or len(data_list[0]) == 0:
            raise ValueError("Data list cannot be empty")
        if self.k > min(len(seq) for data in data_list for seq in data):
            raise ValueError("k-mer length larger than sequence")

    def _build(self, data_list: List[List[str]]) -> None:
        for data in data_list:
            for seq in data:
                if len(seq) < self.k + 1:
                    continue
                rc = reverse_complement(seq)
                for i in range(len(seq) - self.k):
                    kmer1 = seq[i:i+self.k]
                    kmer2 = seq[i+1:i+1+self.k]
                    self._add_arc(kmer1, kmer2)
                    rc_kmer1 = rc[i:i+self.k]
                    rc_kmer2 = rc[i+1:i+1+self.k]
                    self._add_arc(rc_kmer1, rc_kmer2)

    def _add_node(self, kmer: str) -> int:
        if kmer not in self.kmer2idx:
            idx = self.kmer_count
            self.kmer2idx[kmer] = idx
            self.nodes[idx] = Node(kmer)
            self.kmer_count += 1
        idx = self.kmer2idx[kmer]
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1: str, kmer2: str) -> None:
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, idx: int) -> int:
        return self.nodes[idx].get_count()

    def _get_sorted_children(self, idx: int) -> List[int]:
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, idx: int) -> int:
        node = self.nodes[idx]
        if not node.visited:
            node.visited = True
            max_depth = 0
            max_child = None
            for child in self._get_sorted_children(idx):
                d = self._get_depth(child)
                if d > max_depth:
                    max_depth = d
                    max_child = child
            node.depth = max_depth + 1
            node.max_depth_child = max_child
        return node.depth

    def _reset(self) -> None:
        for node in self.nodes.values():
            node.reset()

    def _get_longest_path(self) -> List[int]:
        self._reset()
        max_depth = 0
        max_idx: Optional[int] = None
        for idx in self.nodes:
            d = self._get_depth(idx)
            if d > max_depth:
                max_depth = d
                max_idx = idx
        path: List[int] = []
        while max_idx is not None:
            path.append(max_idx)
            max_idx = self.nodes[max_idx].max_depth_child
        return path

    def _delete_path(self, path: List[int]) -> None:
        path_set = set(path)
        for idx in path:
            del self.nodes[idx]
        for node in self.nodes.values():
            node.remove_children(path_set)

    def _concat_path(self, path: List[int]) -> Optional[str]:
        if not path:
            return None
        concat = self.nodes[path[0]].kmer
        for idx in path[1:]:
            concat += self.nodes[idx].kmer[-1]
        return concat

    def get_longest_contig(self) -> Optional[str]:
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig
