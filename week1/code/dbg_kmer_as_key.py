import copy
from typing import List, Optional, Dict, Set

def reverse_complement(key: str) -> str:
    complement: Dict[str, str] = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    key_list: List[str] = list(key[::-1])
    for i in range(len(key_list)):
        key_list[i] = complement[key_list[i]]
    return ''.join(key_list)


class Node:
    _children: Set[str]
    _count: int
    visited: bool
    depth: int
    max_depth_child: Optional[str]

    def __init__(self) -> None:
        self._children = set()
        self._count = 0
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def add_child(self, kmer: str) -> None:
        if kmer is not None:
            self._children.add(kmer)

    def increase(self) -> None:
        self._count += 1

    def reset(self) -> None:
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def get_count(self) -> int:
        return self._count

    def get_children(self) -> List[str]:
        return list(self._children)

    def remove_children(self, target: Set[str]) -> None:
        self._children -= target


class DBG:
    k: int
    nodes: Dict[str, Node]

    def __init__(self, k: int, data_list: List[List[str]]) -> None:
        self.k = k
        self.nodes = {}
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list: List[List[str]]) -> None:
        assert len(data_list) > 0, "data_list must not be empty"
        assert self.k <= len(data_list[0][0]), "k-mer size larger than read length"

    def _build(self, data_list: List[List[str]]) -> None:
        for data in data_list:
            for original in data:
                rc: str = reverse_complement(original)
                for i in range(len(original) - self.k):
                    kmer1: str = original[i:i+self.k]
                    kmer2: str = original[i+1:i+1+self.k]
                    rc1: str = rc[i:i+self.k]
                    rc2: str = rc[i+1:i+1+self.k]
                    self._add_arc(kmer1, kmer2)
                    self._add_arc(rc1, rc2)

    def _add_node(self, kmer: Optional[str]) -> None:
        if kmer is None or len(kmer) != self.k:
            return
        if kmer not in self.nodes:
            self.nodes[kmer] = Node()
        self.nodes[kmer].increase()

    def _add_arc(self, kmer1: Optional[str], kmer2: Optional[str]) -> None:
        self._add_node(kmer1)
        self._add_node(kmer2)
        if kmer1 is not None and kmer2 is not None:
            self.nodes[kmer1].add_child(kmer2)

    def _get_count(self, child: str) -> int:
        return self.nodes[child].get_count()

    def _get_sorted_children(self, kmer: str) -> List[str]:
        children: List[str] = self.nodes[kmer].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, kmer: str) -> int:
        node: Node = self.nodes[kmer]
        if not node.visited:
            node.visited = True
            children: List[str] = self._get_sorted_children(kmer)
            max_depth: int = 0
            max_child: Optional[str] = None
            for child in children:
                depth: int = self._get_depth(child)
                if depth > max_depth:
                    max_depth, max_child = depth, child
            node.depth = max_depth + 1
            node.max_depth_child = max_child
        return node.depth

    def _reset(self) -> None:
        for node in self.nodes.values():
            node.reset()

    def _get_longest_path(self) -> List[str]:
        max_depth: int = 0
        max_kmer: Optional[str] = None
        for kmer in self.nodes.keys():
            depth: int = self._get_depth(kmer)
            if depth > max_depth:
                max_depth, max_kmer = depth, kmer
        path: List[str] = []
        while max_kmer is not None:
            path.append(max_kmer)
            max_kmer = self.nodes[max_kmer].max_depth_child
        return path

    def _delete_path(self, path: List[str]) -> None:
        path_set: Set[str] = set(path)
        for kmer in path:
            if kmer in self.nodes:
                del self.nodes[kmer]
        for kmer in self.nodes.keys():
            self.nodes[kmer].remove_children(path_set)

    def _concat_path(self, path: List[str]) -> Optional[str]:
        if not path:
            return None
        concat: str = copy.copy(path[0])
        for i in range(1, len(path)):
            concat += path[i][-1]
        return concat

    def get_longest_contig(self) -> Optional[str]:
        self._reset()
        path: List[str] = self._get_longest_path()
        contig: Optional[str] = self._concat_path(path)
        self._delete_path(path)
        return contig
