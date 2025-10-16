from typing import List, Optional, Tuple
import math
import copy
from python import numpy as pnp

MAX_FLOAT: float = 1.7976931348623157e+308 


# --- TREE ---

@extend
class set:
    def __hash__(self):
        #MAX = int.MAX
        MAX = 2**63 - 1
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
    _is_root: bool
    _distance: float
    _parent: Optional[TreeNode]
    _children: List[TreeNode]
    _index: int

    def __init__(self, children, distances: List[float] = None, index: int = -1):
        """
        If index is provided -> leaf node.
        Otherwise -> intermediate node, children and distances must be provided.
        """
        self._is_root: bool = False
        self._distance: float = 0.0
        self._parent = None
        self._children = [] 
        self._index = index

        if children is not None and distances is not None:
            if len(children) != len(distances):
                raise ValueError("Children and distances must have same length")
            for c, d in zip(children, distances):
                c._set_parent(self, float(d))
                self._children.append(c)
        
    @staticmethod
    def leaf(index: int) -> TreeNode:
        node = TreeNode(index=index)
        return node

    @staticmethod
    def internal(children: List[TreeNode], distances: List[float]) -> TreeNode:
        node = TreeNode(children=children, distances=distances)
        node._index = -1
        node._children = [c for c in children]  
        for c, d in zip(node._children, distances):
            c._set_parent(node, float(d))
        return node

    def _set_parent(self, parent: TreeNode, distance: float) -> None:
        if self._parent is not None or self._is_root:
            raise TreeError("Node already has a parent")
        self._parent = parent
        self._distance = distance

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
    def children(self) -> List[TreeNode]:
        return [c for c in self._children]

    @property
    def parent(self):
        return self._parent

    @property
    def distance(self):
        return None if self._parent is None else self._distance

    visited = set()
    def copy(self, visited=None):
        if visited is None:
            visited = set()
        if id(self) in visited:
            raise TreeError("Cycle detected in tree")
        visited.add(id(self))
        if self.is_leaf():
            return TreeNode(index=self._index)
        children_clones: list[TreeNode] = [child.copy(visited) for child in self._children]
        distances_f: list[float] = [float(child.distance) for child in self._children]
        return TreeNode(children_clones, distances_f)

    def get_leaves(self):
        """
        Return List of leaf nodes (direct or indirect).
        """
        leaf_list: List[Optional[TreeNode]] = []
        self._collect_leaves(leaf_list)
        return leaf_list

    def _collect_leaves(self, leaves: List[Optional[TreeNode]]):
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
                    distance_str = str(round(self._distance, round_distance))
                    return f"{label}:{distance_str}"
            else:
                return f"{label}"
        else:
            child_strings: List[str] = [child.to_newick(labels, include_distance, round_distance)
                             for child in self._children]
            joined_children = ",".join(child_strings)
            if include_distance:
                if round_distance is None:
                    return f"({joined_children}):{self._distance}"
                else:
                    distance_str = str(round(self._distance, round_distance))
                    return f"({joined_children}):{distance_str}"
            else:
                return f"({joined_children})"

    @staticmethod
    def from_newick(newick: str, labels: Optional[List[str]] = None) -> TreeNode:
        # remove whitespace
        s = "".join(newick.split())
        if len(newick) == 0:
            raise ValueError("Newick string is empty")
        if s.endswith(";"):
            s = s[:-1]

        stack: List[Tuple[List[Optional[TreeNode]], List[float]]] = []
        children: List[Optional[TreeNode]] = []
        distances: List[float] = []
        i = 0
        start = 0
        while i < len(s):
            c = s[i]
            if c == "(":
                # start new subtree
                stack.append((children, distances))
                children = []
                distances = []
                start = i + 1
            elif c == "," or c == ")":
                # parse node between start:i
                if start < i:
                    node_str = s[start:i]
                    if node_str != "":
                        node, dist = TreeNode._parse_leaf(node_str, labels)
                        children.append(node)
                        distances.append(dist)
                start = i + 1
                if c == ")":
                    # finish current subtree
                    parent_children = [i for i in children]
                    parent_distances = [i for i in distances]
                    if stack:
                        children, distances = stack.pop()
                    else:
                        children = []
                        distances = []

                    # parse distance after ')' if any
                    j = i + 1
                    dist_val = 0.0
                    if j < len(s) and s[j] == ":":
                        j += 1
                        end = j
                        while end < len(s) and (s[end].isdigit() or s[end] in ".eE-"):
                            end += 1
                        dist_val = float(s[j:end])
                        i = end - 1

                    # create internal node, children/distances must be lists
                    node = TreeNode(children=parent_children, distances=parent_distances)
                    children.append(node)
                    distances.append(dist_val)
            i += 1

        # handle final leaf if any
        if start < len(s):
            node_str = s[start:]
            if node_str != "":
                node, dist = TreeNode._parse_leaf(node_str, labels)
                children.append(node)
                distances.append(dist)

        if len(children) != 1:
            raise ValueError("Malformed Newick string")

        return children[0]

    @staticmethod
    def _parse_leaf(s: str, labels: Optional[List[str]]) -> Tuple[TreeNode, float]:
        """
        Parse a leaf string from Newick format.
        Returns a TreeNode (leaf) and its distance.
        """
        s = s.strip()
        if len(s) == 0:
            raise ValueError("Empty leaf string")

        if ":" in s:
            parts = s.split(":")
            if len(parts) != 2:
                raise ValueError(f"Malformed leaf: {s}")
            label_str, dist_str = parts
            try:
                dist = float(dist_str)
            except:
                raise ValueError(f"Invalid distance '{dist_str}' in leaf '{s}'")
        else:
            label_str = s
            dist = 0.0

        label_str = label_str.strip()

        # determine index
        if labels is None:
            try:
                idx = int(label_str)  # normal integer
            except:
                try:
                    idx = int(float(label_str))  # handles "12.0"
                except:
                    raise ValueError(f"Invalid numeric leaf label: '{label_str}'")
        else:
            try:
                idx = labels.index(label_str)
            except ValueError:
                raise ValueError(f"Label '{label_str}' not found in labels list")

        # create leaf node
        leaf = TreeNode(index=idx)
        return leaf, dist


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
    _root: TreeNode
    _leaves: List[Optional[TreeNode]]

    def __init__(self, root: TreeNode):
        root.as_root()
        self._root: TreeNode = root
        leaves_unsorted = self._root.get_leaves()
        leaf_count = len(leaves_unsorted)
        indices = pnp.array([leaf._index for leaf in leaves_unsorted], dtype=pnp.int64)
        self._leaves: List[Optional[TreeNode]] = [None for _ in range(leaf_count)]  
        for i in range(len(indices)):
            idx: int = int(indices[i])
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

        root = TreeNode.from_newick(s, labels) 
        return Tree(root)

    def __str__(self):
        return self.to_newick()


def _find_min_pair_triangular(mat, mask: list[bool]) -> Tuple[int, int]: 
    """
    Finds indices (i,j) with i>j of minimum mat[i,j] among entries where mask[i] and mask[j] are True.
    """
    dist_min: float = MAX_FLOAT
    i_min: int = -1
    j_min: int = -1
    n: int = mat.shape[0]
    for i in range(n):
        if not mask[i]:
            continue
        for j in range(i):
            if not mask[j]:
                continue
            d = mat[i, j]
            if d < dist_min:
                dist_min = d
                i_min = i
                j_min = j
    return int(i_min), int(j_min)


# --- UPGMA ---

def upgma(distances): 
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

    n0: int = distances.shape[0]
    if n0 < 2:
        raise ValueError("At least 2 nodes are required")

    D = pnp.array(distances, dtype=pnp.float64)

    #nodes: List[TreeNode] = [TreeNode(index=i) for i in range(n0)]
    nodes: List[Optional[TreeNode]] = []
    for i in range(n0):
        nodes.append(TreeNode.leaf(i))
    cluster_size: List[int] = [1 for _ in range(n0)]
    node_heights: List[float] = [0.0 for _ in range(n0)]
    active: List[bool] = [True for _ in range(n0)]
    remaining: int = n0

    while remaining > 1:
        # find min pair
        i_min, j_min = _find_min_pair_triangular(D, active)
        if i_min == -1 or j_min == -1:
            raise RuntimeError("UPGMA failed: no valid pair found")

        dist_min: float = D[i_min, j_min]
        height: float = dist_min / 2.0

        child_i = nodes[i_min]
        child_j = nodes[j_min]

        if child_i is None or child_j is None:
            raise RuntimeError(f"UPGMA failed: invalid child nodes {i_min}, {j_min}")
        h_i = max(0.0, float(height - node_heights[i_min]))
        h_j = max(0.0, float(height - node_heights[j_min]))
        if h_i < 0 or h_j < 0:
            raise RuntimeError(f"UPGMA failed: negative branch length {h_i}, {h_j}")

        # --- Merge nodes ---
        print(f"Merging nodes {i_min} and {j_min} with distance {dist_min}")
        new_node = TreeNode.internal([c for c in [child_i, child_j]], [h_i, h_j])

        nodes[i_min] = new_node
        node_heights[i_min] = height
        active[j_min] = False

        remaining -= 1

        # update distances: arithmetic mean weighted by cluster sizes
        for k in range(n0):
            if k == i_min or not active[k]:
                continue
            D[i_min, k] = float((D[i_min, k] * cluster_size[i_min] + D[j_min, k] * cluster_size[j_min]) / (cluster_size[i_min] + cluster_size[j_min]))
            D[k, i_min] = D[i_min, k]

        # update cluster size
        cluster_size[i_min] += cluster_size[j_min]

    # Get the last active node to be the root
    last_active_nodes = [nodes[i] for i, act in enumerate(active) if act]
    if len(last_active_nodes) != 1 or last_active_nodes[0] is None:
        raise RuntimeError(f"UPGMA failed: expected 1 active node, found {len(last_active_nodes)}")

    root_node = last_active_nodes[0]

    # If the root is a leaf (only 1 leaf in tree), wrap in a dummy root
    if root_node.is_leaf():
        root_node = TreeNode(children=[root_node], distances=[0.0])

    return Tree(root_node)


# --- NEIGHBOUR JOINING ---

def neighbor_joining(distances): 
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

    n0: int = distances.shape[0]
    if n0 < 3:
        raise ValueError("At least 3 nodes are required")

    D = pnp.array(distances, dtype=pnp.float64)

    #nodes: List[TreeNode] = [TreeNode(index=i) for i in range(n0)]
    nodes: List[Optional[TreeNode]] = []
    for i in range(n0):
        nodes.append(TreeNode.leaf(i))
    active: List[bool] = [True for _ in range(n0)]
    remaining: int = n0

    while remaining > 3:
        # total distances per row 
        total = pnp.zeros(D.shape[0], dtype=pnp.float64)
        for i in range(D.shape[0]):
            if active[i]:
                total[i] = sum((D[i, j]) for j in range(D.shape[0]) if active[j])

        # compute Q-matrix 
        Q = pnp.full(D.shape, MAX_FLOAT, dtype=pnp.float64)
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
            raise RuntimeError("Neighbor-Joining failed: no valid pair found")

        dist_ij = D[i_min, j_min]
        delta: float = 0.0
        if remaining > 2:
            delta = (total[i_min] - total[j_min]) / (remaining - 2)
        limb_i: float = 0.5 * (dist_ij + delta)
        limb_j: float = 0.5 * (dist_ij - delta)

        child_i = nodes[i_min]
        child_j = nodes[j_min]
        if child_i is None or child_j is None:
            raise RuntimeError(f"NJ failed: invalid nodes {i_min}, {j_min}")
        if limb_i < 0 or limb_j < 0:
            raise RuntimeError(f"NJ failed: negative limb lengths {limb_i}, {limb_j}")

        # create new node 
        print(f"Merging nodes {i_min} and {j_min} with distance {dist_ij}")
        new_node = TreeNode.internal([c for c in [child_i, child_j]], [float(limb_i), float(limb_j)])
        nodes[i_min] = new_node
        active[j_min] = False

        remaining -= 1

        # update distances for new cluster i_min 
        mask = [k for k in range(D.shape[0]) if active[k] and (k != i_min or not active[k])]
        for k in mask:
            D[i_min, k] = float(0.5 * (D[i_min, k] + D[j_min, k] - dist_ij))
            D[k, i_min] = D[i_min, k]

    # final combine remaining nodes into a root
    # collect active nodes
    final_nodes = [nodes[i] for i in range(D.shape[0]) if active[i]]
    final_idx = [i for i in range(D.shape[0]) if active[i]]

    if any(child is None for child in final_nodes):
        raise RuntimeError("NJ failed: some final nodes are None")

    if len(final_nodes) == 2:
        a, b = final_nodes
        ia, ib = final_idx
        d = D[ia, ib]
        root = TreeNode(children=[a, b], distances=[float(max(d/2.0, 0.0)), float(max(d/2.0, 0.0))])
    elif len(final_nodes) == 3:
        a, b, c = final_nodes
        ia, ib, ic = final_idx
        da = 0.5 * (D[ia, ic] + D[ia, ib] - D[ib, ic])
        db = 0.5 * (D[ib, ia] + D[ib, ic] - D[ia, ic])
        dc = 0.5 * (D[ic, ia] + D[ic, ib] - D[ia, ib])
        root = TreeNode(children=[a, b, c], distances=[float(max(da, 0.0)), float(max(db, 0.0)), float(max(dc, 0.0))])
    else:
        raise RuntimeError("Neighbor-Joining failed")

    return Tree(root)