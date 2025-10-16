#import numpy as np
from typing import List, Optional
import math
import copy
from python import numpy as pnp

MAX_FLOAT = pnp.finfo(pnp.float64).max

# --- TREE ---

@extend
class set:
    def __hash__(self):
        MAX = int.MAX
        MASK = 2 * MAX + 1
        n = len(self)
        h = 1927868237 * (n + 1)
        h &= MASK
        for x in self:
            hx = hash(x)
            h ^= (hx ^ (hx << 16) ^ 89869747) * 3644798167
            h &= MASK
        h = h * 69069 + 907133923
        h &= MASK
        if h > MAX:
            h -= MASK + 1
        if h == -1:
            h = 590923713
        return h


class TreeError(Static[Exception]):
    """Exception used for tree topology related errors."""
    pass


class TreeNode:

    def __init__(self, children: List[TreeNode] = None, distances: List[float] = None, index: int = None):
        """
        If index is provided -> leaf node.
        Otherwise -> intermediate node, children and distances must be provided.
        """
        self._is_root: bool = False
        self._distance: float = 0.0
        self._parent: Optional[TreeNode] = None
        self._children: List[TreeNode] = []
        self._index: int = -1

        if index is not None:
            # Leaf node
            self._index = int(index)
            if children is not None or distances is not None:
                raise TypeError("Leaf node cannot have children or distances")
        else:
            # Internal node
            if children is None or distances is None:
                raise TypeError("Internal node requires children and distances")
            if len(children) == 0:
                raise ValueError("Internal node must have at least one child")
            if len(children) != len(distances):
                raise ValueError("Number of children must match number of distances")
            self._children = [c for c in children]
            for child, d in zip(children, distances):
                child._set_parent(self, float(d))

    def _set_parent(self, parent: Optional[TreeNode], distance: float) -> None:
        if self._parent is not None or self._is_root:
            raise TreeError("Node already has a parent")
        self._parent = parent
        self._distance = float(distance)

    def is_leaf(self) -> bool:
        return self._index != -1

    def is_root(self) -> bool:
        return bool(self._is_root)

    def as_root(self) -> None:
        if self._parent is not None:
            raise TreeError("Node has parent, cannot be a root node")
        self._is_root = True

    @property
    def index(self):
        return None if self._index == -1 else self._index

    @property
    def children(self):
        return self._children

    @property
    def parent(self):
        return self._parent

    @property
    def distance(self):
        return None if self._parent is None else self._distance

    def copy(self):
        if self.is_leaf():
            return TreeNode(index=self._index)
        else:
            children_clones = [child.copy() for child in self._children if child is not None]
            distances_f = [float(child.distance) if child is not None else 0.0 for child in self._children]
            return TreeNode(children_clones, distances_f)

    def get_leaves(self):
        """
        Return List of leaf nodes (direct or indirect).
        """
        leaf_list: List[TreeNode] = []
        self._collect_leaves(leaf_list)
        return leaf_list

    def _collect_leaves(self, leaves: List[TreeNode]):
        if self.is_leaf():
            leaves.append(self)
        else:
            for c in self._children:
                c._collect_leaves(leaves)

    def get_indices(self):
        leaves = self.get_leaves()
        return pnp.array([leaf._index for leaf in leaves], dtype=pnp.int64)

    def get_leaf_count(self):
        return _get_leaf_count(self)

    def to_newick(self, labels: Optional[List[str]] = None,
                  include_distance: bool = True,
                  round_distance: int = None) -> str:
        if self.is_leaf():
            if labels is not None:
                #lbls = list(labels)
                label = labels[self._index]
                illegal_chars: List[str] = [",", ":", ";", "(", ")"]
                for ch in illegal_chars:
                    if ch in label:
                        raise ValueError(f"Label '{label}' contains illegal character '{ch}'")
            else:
                label = str(self._index)
            if include_distance:
                if round_distance is None:
                    return f"{label}:{self._distance}"
                else:
                    return f"{label}:{self._distance:.{round_distance}f}"
            else:
                return f"{label}"
        else:
            child_strings: List[str] = [child.to_newick(labels, include_distance, round_distance)
                             for child in self._children]
            if include_distance:
                if round_distance is None:
                    return f"({','.join(child_strings)}):{self._distance}"
                else:
                    return f"({','.join(child_strings)}):{self._distance:.{round_distance}f}"
            else:
                return f"({','.join(child_strings)})"

    @staticmethod
    def from_newick(newick: str, labels: Optional[List[str]] = None):
        # remove whitespace
        s = "".join(newick.split())
        if len(newick) == 0:
            raise ValueError("Newick string is empty")

        if s[0] != "(":
            # Leaf node
            if ":" in s:
                parts = s.split(":")
                label = parts[0]
                distance = float(parts[1])
            else:
                label = s
                distance = 0.0
            if labels is None:
                try:
                    idx = int(label)
                except:
                    idx = int(float(label))
            else:
                idx = labels.index(label)
            return TreeNode(index=idx), distance

        # Internal node
        # Find matching parentheses for top-level split
        level = 0
        split_indices: List[int] = []
        for i, ch in enumerate(s):
            if ch == "(":
                level += 1
            elif ch == ")":
                level -= 1
            elif ch == "," and level == 1:
                split_indices.append(i)
            if level < 0:
                raise ValueError("Mismatched parentheses")

        children: List[TreeNode] = []
        distances: List[float] = []
        start = 1
        for idx in split_indices + [len(s)-1]:
            sub = s[start:idx]
            child, d = TreeNode.from_newick(sub, labels)
            children.append(child)
            distances.append(d)
            start = idx + 1

        # Parse distance after closing parenthesis
        remaining = s[len(s)-1:]
        if remaining.startswith(":"):
            distance = float(remaining[1:])
        else:
            distance = 0.0

        return TreeNode(children=children, distances=distances), distance

        """
        substart = -1
        substop = -1
        # find first '(' and last ')'
        for i, ch in enumerate(newick):
            if ch == "(":
                substart = i
                break
            if ch == ")":
                raise ValueError("Bracket closed before it was opened")
        for i in range(len(newick)-1, -1, -1):
            ch = newick[i]
            if ch == ")":
                substop = i + 1
                break
            if ch == "(":
                raise ValueError("Bracket was opened but not closed")

        if substart == -1 and substop == -1:
            # leaf
            label_and_distance = newick
            try:
                label, dist_s = label_and_distance.split(":")
                distance = float(dist_s)
            except ValueError:
                distance = 0.0
                label = label_and_distance
            #idx = int(float(label)) if labels is None else labels.index(label)
            #return TreeNode(index=idx), distance
            if labels is None:
                try:
                    idx = int(label)  # normal integer string
                except:
                    idx = int(float(label))  # handles "29.0"
            else:
                # lookup in labels List
                idx = labels.index(label)

            return TreeNode(index=idx), distance
        
        else:
            # intermediate node
            if substop == len(newick):
                distance = 0.0
            else:
                label_and_distance = newick[substop:]
                try:
                    _, dist_s = label_and_distance.split(":")
                    distance = float(dist_s)
                except ValueError:
                    distance = 0.0
            sub = newick[substart+1:substop-1]
            if len(sub) == 0:
                raise ValueError("Intermediate node must at least have one child")
            # split top-level commas
            comma_pos: List[int] = []
            level = 0
            for i, ch in enumerate(sub):
                if ch == "(":
                    level += 1
                elif ch == ")":
                    level -= 1
                elif ch == ",":
                    if level == 0:
                        comma_pos.append(i)
                if level < 0:
                    raise ValueError("Bracket closed before it was opened")
            children: List[TreeNode] = []
            distances: List[float] = []
            if len(comma_pos) != 0:
                start = 0
                for pos in comma_pos:
                    child_s = sub[start:pos]
                    child, d = TreeNode.from_newick(child_s, labels=labels)
                    children.append(child)
                    distances.append(d)
                    start = pos + 1
                # last
                child_s = sub[start:]
                child, d = TreeNode.from_newick(child_s, labels=labels)
                children.append(child)
                distances.append(d)
            else:
                child, d = TreeNode.from_newick(sub, labels=labels)
                children.append(child)
                distances.append(d)
            return TreeNode(children, distances), distance
        """



    def lowest_common_ancestor(self, node: TreeNode):
        self_path = _create_path_to_root(self)
        other_path = _create_path_to_root(node)
        lca: TreeNode = None
        min_len = min(len(self_path), len(other_path))
        for i in range(1, min_len + 1):
            if self_path[-i] is other_path[-i]:
                lca = self_path[-i]
            else:
                break
        return lca

    def distance_to(self, node: TreeNode, topological: bool = False) -> float:
        lca = self.lowest_common_ancestor(node)
        if lca is None:
            raise TreeError("The nodes do not have a common ancestor")
        distance: float = 0.0
        current = self
        while current is not lca:
            if topological:
                distance += 1.0
            else:
                distance += float(current._distance)
            current = current._parent  # type: ignore[assignment]
        current = node
        while current is not lca:
            if topological:
                distance += 1.0
            else:
                distance += float(current._distance)
            current = current._parent  # type: ignore[assignment]
        return distance

    def __eq__(self, other: object):
        if not isinstance(other, TreeNode):
            return False
        node: TreeNode = other
        if self._distance != node._distance:
            return False
        if self._index != -1:
            return self._index == node._index
        else:
            # order of children not important
            return set(self._children) == set(node._children)

    def __hash__(self):
        children_set = set(self._children) if len(self._children) > 0 else None
        return hash((self._index, children_set, self._distance))


# --- Helper functions ---

def _get_leaves(node: TreeNode, leaf_list: List[TreeNode]):
    if node._index == -1:
        for child in node._children:
            _get_leaves(child, leaf_list)
    else:
        leaf_list.append(node)



def _get_leaf_count(node: TreeNode):
    if node._index == -1:
        count = 0
        for child in node._children:
            count += _get_leaf_count(child)
        return count
    else:
        return 1


def _create_path_to_root(node: TreeNode):
    path: List[TreeNode] = []
    current: TreeNode = node
    while current is not None:
        path.append(current)
        current = current._parent
    return path


class Tree:

    def __init__(self, root: TreeNode):
        root.as_root()
        self._root: TreeNode = root
        leaves_unsorted = self._root.get_leaves()
        leaf_count = len(leaves_unsorted)
        indices = pnp.array([leaf._index for leaf in leaves_unsorted], dtype=pnp.int64)
        self._leaves: Optional[List[TreeNode]] = [None for _ in range(leaf_count)]  # type: List[TreeNode]
        for i in range(len(indices)):
            idx = int(indices[i])
            if idx >= leaf_count or idx < 0:
                raise TreeError("The tree's indices are out of range")
            self._leaves[idx] = leaves_unsorted[i]

    def __copy_create__(self):
        return Tree(self._root.copy())

    @property
    def root(self):
        return self._root

    @property
    def leaves(self):
        return copy.copy(self._leaves)

    def __len__(self):
        return len(self._leaves)

    def get_distance(self, index1: int, index2: int, topological: bool = False):
        return self._leaves[index1].distance_to(self._leaves[index2], topological)

    def to_newick(self, labels: List[str] = None, include_distance: bool = True,
                  round_distance: int = None):
        return self._root.to_newick(labels, include_distance, round_distance) + ";"

    @staticmethod
    def from_newick(newick: str, labels: List[str] = None):
        s = newick.strip()
        if len(s) == 0:
            raise ValueError("Newick string is empty")
        if s.endswith(";"):
            s = s[:-1]
        root, _ = TreeNode.from_newick(s, labels)

        return Tree(root)

    def __str__(self):
        return self.to_newick()

    def __eq__(self, other: object):
        if not isinstance(other, Tree):
            return False
        return self._root == other._root

    def __hash__(self):
        return hash(self._root)


def _find_min_pair_triangular(mat: pnp.ndarray[float, 2], mask: List[bool]):
    """
    Finds indices (i,j) with i>j of minimum mat[i,j] among entries where mask[i] and mask[j] are True.
    """
    dist_min = float(MAX_FLOAT)
    i_min = -1
    j_min = -1
    n = mat.shape[0]
    for i in range(n):
        if not mask[i]:
            continue
        for j in range(i):
            if not mask[j]:
                continue
            d = float(mat[i, j])
            if d < dist_min:
                dist_min = d
                i_min = i
                j_min = j
    return i_min, j_min


# --- UPGMA ---

def upgma(distances: pnp.ndarray):
    """
    distances: square numpy array (any dtype) -> converted to float64 inside
    Returns a Tree constructed with the UPGMA algorithm.
    """
    if distances.shape[0] != distances.shape[1] or not pnp.allclose(distances.T, distances):
        raise ValueError("Distance matrix must be symmetric")
    if pnp.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")

    n0 = distances.shape[0]
    if n0 < 2:
        raise ValueError("At least 2 nodes are required")

    D = distances.astype(pnp.float64, copy=True)

    nodes: List[TreeNode] = [TreeNode(index=i) for i in range(n0)]
    cluster_size: List[int] = [1 for _ in range(n0)]
    node_heights: List[float] = [0.0 for _ in range(n0)]
    active: List[bool] = [True for _ in range(n0)]
    remaining = n0

    while remaining > 1:
        # find min pair
        i_min, j_min = _find_min_pair_triangular(D, active)
        if i_min == -1 or j_min == -1:
            break

        # sanity check
        # assert active[i_min] and active[j_min], "Merging inactive nodes!"

        dist_min = float(D[i_min, j_min])
        height = dist_min / 2.0

        child_i = nodes[i_min].copy()
        child_j = nodes[j_min].copy()

        nodes[i_min] = TreeNode(
        children=[child_i, child_j],
        distances=[float(height - node_heights[i_min]),
                   float(height - node_heights[j_min])]
        )
        node_heights[i_min] = height
        active[j_min] = False
        #nodes[j_min] = None

        remaining -= 1

        # update distances: arithmetic mean weighted by cluster sizes
        for k in range(n0):
            if k == i_min or not active[k]:
                continue
            mean = ((D[i_min, k] * cluster_size[i_min] +
                 D[j_min, k] * cluster_size[j_min]) /
                (cluster_size[i_min] + cluster_size[j_min]))
            D[i_min, k] = D[k, i_min] = mean

        # update cluster size
        cluster_size[i_min] += cluster_size[j_min]

    # find the last active node to be root
    for idx in range(n0-1, -1, -1):
        if active[idx]:
            return Tree(nodes[idx])

    # fallback
    raise RuntimeError("UPGMA failed")


# --- NEIGHBOUR JOINING ---

def neighbor_joining(distances: pnp.ndarray):
    """
    distances: square numpy array -> converted to float64 inside
    Returns a Tree constructed with the Neighbor-Joining algorithm.
    """
    if distances.shape[0] != distances.shape[1] or not pnp.allclose(distances.T, distances):
        raise ValueError("Distance matrix must be symmetric")
    if pnp.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")

    n0 = distances.shape[0]
    if n0 < 3:
        raise ValueError("At least 3 nodes are required")

    D = distances.astype(pnp.float64, copy=True)
    nodes: List[TreeNode] = [TreeNode(index=i) for i in range(n0)]
    active: List[bool] = [True for _ in range(n0)]
    remaining = n0

    while remaining > 3:
        # total distances per row 
        total = pnp.zeros(D.shape[0], dtype=pnp.float64)
        for i in range(D.shape[0]):
            if active[i]:
                total[i] = sum(D[i, j] for j in range(D.shape[0]) if active[j])

        # compute Q-matrix 
        Q = pnp.full(D.shape, float(MAX_FLOAT), dtype=pnp.float64)
        for i in range(D.shape[0]):
            if not active[i]:
                continue
            for j in range(i):
                if not active[j]:
                    continue
                Q[i, j] = (remaining - 2) * D[i, j] - total[i] - total[j]
                Q[j, i] = Q[i, j]

        # find min pair in Q
        i_min, j_min = _find_min_pair_triangular(Q, active)
        if i_min == -1 or j_min == -1:
            break

        dist_ij = float(D[i_min, j_min])
        delta = 0.0
        if remaining > 2:
            delta = (total[i_min] - total[j_min]) / (remaining - 2)
        limb_i = 0.5 * (dist_ij + delta)
        limb_j = 0.5 * (dist_ij - delta)

        child_i = nodes[i_min].copy()
        child_j = nodes[j_min].copy()

        # create new node 
        new_node = TreeNode(children=[child_i, child_j], distances=[float(limb_i), float(limb_j)])

        # update nodes and active array
        nodes[i_min] = new_node
        active[j_min] = False
        #nodes[j_min] = None  

        # update distances for new cluster i_min 
        mask = [k for k in range(D.shape[0]) if active[k] and (k != i_min or not active[k])]
        for k in mask:
            D[i_min, k] = D[k, i_min] = 0.5 * (D[i_min, k] + D[j_min, k] - dist_ij)

        remaining -= 1

    # final combine remaining nodes into a root
    # collect active nodes
    final_nodes = [nodes[i] for i in range(D.shape[0]) if active[i]]
    final_idx = [i for i in range(D.shape[0]) if active[i]]

    if len(final_nodes) == 2:
        a, b = final_nodes
        ia, ib = final_idx
        d = float(D[ia, ib])
        root = TreeNode(children=[a, b], distances=[d/2.0, d/2.0])
    elif len(final_nodes) == 3:
        a, b, c = final_nodes
        ia, ib, ic = final_idx
        da = 0.5 * (D[ia, ic] + D[ia, ib] - D[ib, ic])
        db = 0.5 * (D[ib, ia] + D[ib, ic] - D[ia, ic])
        dc = 0.5 * (D[ic, ia] + D[ic, ib] - D[ia, ib])
        root = TreeNode(children=[a, b, c], distances=[float(da), float(db), float(dc)])
    else:
        raise RuntimeError("Neighbor-Joining failed")

    return Tree(root)