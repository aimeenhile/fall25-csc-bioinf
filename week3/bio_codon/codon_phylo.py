import numpy as np
from typing import List, Tuple, Optional
import math
import copy
from python import numpy as pnp


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

    def __init__(self, children: Optional[List[TreeNode]] = None, distances: Optional[List[float]] = None, index: Optional[int] = None):
        """
        If index is provided -> leaf node.
        Otherwise -> intermediate node, children and distances must be provided.
        """
        self._is_root: bool = False
        self._distance: float = 0.0
        self._parent: Optional[TreeNode] = None
        self._children: Optional[List[TreeNode]] = [] 
        self._index: int = -1

        if index is None:
            # intermediate node -> need children and distances
            if children is None or distances is None:
                raise TypeError(
                    "Either reference index (for terminal node) or "
                    "child nodes including the distance (for intermediate node) must be set"
                )
            child_list = [c for c in children]
            dist_list = [float(d) for d in distances]
            if len(child_list) == 0:
                raise TreeError("Intermediate nodes must at least contain one child node")
            if len(child_list) != len(dist_list):
                raise ValueError("The number of children must equal the number of distances")
            # ensure distinct child objects
            for i in range(len(child_list)):
                for j in range(len(child_list)):
                    if i != j and child_list[i] is child_list[j]:
                        raise TreeError("Two child nodes cannot be the same object")
            # set fields
            self._index = -1
            self._children = [c for c in child_list]
            # assign parent & distances to children
            for child, d in zip(child_list, dist_list):
                child._set_parent(self, float(d))
        else:
            # leaf node
            if children is not None or distances is not None:
                raise TypeError("Reference index and child nodes are mutually exclusive")
            if index < 0:
                raise ValueError("Index cannot be negative")
            self._index = int(index)
            self._children = []

    # internal
    def _set_parent(self, parent: Optional[TreeNode], distance: float) -> None:
        if self._parent is not None or self._is_root:
            raise TreeError("Node already has a parent")
        self._parent = parent
        self._distance = float(distance)
        return None

    # public API used in tests
    def is_leaf(self) -> bool:
        return False if self._index == -1 else True

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
            distances = [child.distance for child in self._children]
            children_clones = [child.copy() for child in self._children]
            # distances might contain None if parent's distance missing; ensure floats
            distances_f = [float(d) if d is not None else 0.0 for d in distances]
            return TreeNode(children_clones, distances_f)

    def get_leaves(self):
        """
        Return List of leaf nodes (direct or indirect).
        """
        leaf_list: list["TreeNode"] = []
        _get_leaves(self, leaf_list)
        return leaf_list

    def get_indices(self):
        leaves = self.get_leaves()
        return np.array([leaf._index for leaf in leaves], dtype=np.int32)

    def get_leaf_count(self):
        return _get_leaf_count(self)

    def to_newick(self, labels: Optional[list[str]] = None,
                  include_distance: bool = True,
                  round_distance: Optional[int] = None):
        if self.is_leaf():
            if labels is not None:
                lbls = list(labels)
                label = lbls[self._index]
                illegal_chars = [",", ":", ";", "(", ")"]
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
            child_strings = [child.to_newick(labels, include_distance, round_distance)
                             for child in self._children]
            if include_distance:
                if round_distance is None:
                    return f"({','.join(child_strings)}):{self._distance}"
                else:
                    return f"({','.join(child_strings)}):{self._distance:.{round_distance}f}"
            else:
                return f"({','.join(child_strings)})"

    @staticmethod
    def from_newick(newick: str, labels: Optional[list[str]] = None):
        # remove whitespace
        newick = "".join(newick.split())
        if len(newick) == 0:
            raise ValueError("Newick string is empty")

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
            children: list["TreeNode"] = []
            distances: list[float] = []
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

    def lowest_common_ancestor(self, node: "TreeNode"):
        self_path = _create_path_to_root(self)
        other_path = _create_path_to_root(node)
        lca: Optional["TreeNode"] = None
        min_len = min(len(self_path), len(other_path))
        for i in range(1, min_len + 1):
            if self_path[-i] is other_path[-i]:
                lca = self_path[-i]
            else:
                break
        return lca

    def distance_to(self, node: "TreeNode", topological: bool = False) -> float:
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
        node: "TreeNode" = other
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

def _get_leaves(node: "TreeNode", leaf_list: List[TreeNode]):
    if node._index == -1:
        for child in node._children:
            _get_leaves(child, leaf_list)
    else:
        leaf_list.append(node)


def _get_leaf_count(node: "TreeNode"):
    if node._index == -1:
        count = 0
        for child in node._children:
            count += _get_leaf_count(child)
        return count
    else:
        return 1


def _create_path_to_root(node: "TreeNode"):
    path: list["TreeNode"] = []
    current: Optional["TreeNode"] = node
    while current is not None:
        path.append(current)
        current = current._parent
    return path


class Tree:

    def __init__(self, root: TreeNode):
        root.as_root()
        self._root: TreeNode = root
        # gather leaves and place them at positions equal to leaf.index
        leaves_unsorted = self._root.get_leaves()
        leaf_count = len(leaves_unsorted)
        indices = np.array([leaf._index for leaf in leaves_unsorted], dtype=np.int32)
        self._leaves: list[TreeNode] = [None] * leaf_count  # type: ignore[assignment]
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

    def to_newick(self, labels: Optional[list[str]] = None, include_distance: bool = True,
                  round_distance: Optional[int] = None):
        return self._root.to_newick(labels, include_distance, round_distance) + ";"

    @staticmethod
    def from_newick(newick: str, labels: Optional[list[str]] = None):
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


MAX_FLOAT = np.finfo(np.float64).max


def _find_min_pair_triangular(mat: np.ndarray, mask: list[bool]):
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

def upgma(distances: np.ndarray):
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

    D = distances.astype(np.float64, copy=True)

    # nodes: current nodes (TreeNode)
    nodes: list[TreeNode] = [TreeNode(index=i) for i in range(n0)]
    # cluster sizes
    cluster_size: list[int] = [1 for _ in range(n0)]
    # node heights
    node_heights: list[float] = [0.0 for _ in range(n0)]
    # track whether position is active
    active: list[bool] = [True for _ in range(n0)]
    remaining = n0

    while remaining > 1:
        # find min pair
        i_min, j_min = _find_min_pair_triangular(D, active)
        if i_min == -1 or j_min == -1:
            break

        dist_min = float(D[i_min, j_min])
        height = dist_min / 2.0

        # create new node at i_min, mark j_min inactive
        child_i = nodes[i_min]
        child_j = nodes[j_min]
        child_dist_i = float(height - node_heights[i_min])
        child_dist_j = float(height - node_heights[j_min])
        nodes[i_min] = TreeNode(children=[child_i, child_j], distances=[child_dist_i, child_dist_j])
        node_heights[i_min] = height
        nodes[j_min] = None  
        active[j_min] = False

        # update distances: arithmetic mean weighted by cluster sizes
        for k in range(n0):
            if not active[k] or k == i_min:
                continue
            left = float(D[i_min, k]) * float(cluster_size[i_min])
            right = float(D[j_min, k]) * float(cluster_size[j_min])
            denom = float(cluster_size[i_min] + cluster_size[j_min])
            mean = 0.0
            if denom != 0.0:
                mean = (left + right) / denom
            D[i_min, k] = mean
            D[k, i_min] = mean

        # update cluster size
        cluster_size[i_min] = cluster_size[i_min] + cluster_size[j_min]
        remaining -= 1

    # find the last active node to be root
    for idx in range(n0-1, -1, -1):
        if active[idx]:
            root_node = nodes[idx]
            return Tree(root_node)

    # fallback
    raise RuntimeError("UPGMA failed to construct a tree")


# --- NEIGHBOUR JOINING ---

def neighbor_joining(distances: np.ndarray):
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

    D = distances.astype(np.float64, copy=True)
    nodes: list[TreeNode] = [TreeNode(index=i) for i in range(n0)]
    active: list[bool] = [True for _ in range(n0)]
    remaining = n0

    while remaining > 3:
        # total distances per row 
        total = np.zeros(D.shape[0], dtype=np.float64)
        for i in range(D.shape[0]):
            if not active[i]:
                total[i] = 0.0
                continue
            s = 0.0
            for j in range(D.shape[0]):
                if not active[j]:
                    continue
                s += float(D[i, j])
            total[i] = s

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
        denom = float(remaining - 2)
        if denom == 0.0:
            delta = 0.0
        else:
            delta = (total[i_min] - total[j_min]) / denom
        limb_i = 0.5 * (dist_ij + delta)
        limb_j = 0.5 * (dist_ij - delta)

        # create new node at i_min, mark j_min inactive
        new_node = TreeNode(children=[nodes[i_min], nodes[j_min]], distances=[float(limb_i), float(limb_j)])
        nodes[i_min] = new_node
        nodes[j_min] = None  
        active[j_min] = False

        # compute new distances to replace row/col i_min
        # mask of remaining indices excluding i_min and j_min
        mask = [k for k in range(D.shape[0]) if active[k] and k != i_min]
        new_row = []
        for k in mask:
            d = 0.5 * (D[i_min, k] + D[j_min, k] - dist_ij)
            new_row.append(d)

        # build new D' matrix in-place: set distances for i_min to masked entries
        # and set symmetric counterparts
        for idx_k, k in enumerate(mask):
            D[i_min, k] = new_row[idx_k]
            D[k, i_min] = new_row[idx_k]

        remaining -= 1

    # final combine remaining nodes into a root
    # collect active nodes
    final_nodes: list[TreeNode] = []
    final_idx: list[int] = []

    for idx in range(D.shape[0]):
        if active[idx]:
            final_nodes.append(nodes[idx])
            final_idx.append(idx)

    if len(final_nodes) == 2:
        d = float(D[final_idx[0], final_idx[1]])
        root = TreeNode(children=[final_nodes[0], final_nodes[1]], distances=[d/2.0, d/2.0])
    else:
        # len == 3
        a, b, c = final_nodes
        ia, ib, ic = final_idx
        da = 0.5 * (D[ia, ic] + D[ia, ib] - D[ib, ic])
        db = 0.5 * (D[ib, ia] + D[ib, ic] - D[ia, ic])
        dc = 0.5 * (D[ic, ia] + D[ic, ib] - D[ia, ib])
        root = TreeNode(children=[a, b, c], distances=[float(da), float(db), float(dc)])

    return Tree(root)