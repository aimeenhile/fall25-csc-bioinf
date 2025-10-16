# ChatGPT Model GPT-5

```python
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
    assert active[i_min] and active[j_min], "Merging inactive nodes!"

    dist_min = float(D[i_min, j_min])
    height = dist_min / 2.0

    child_i = nodes[i_min]
    child_j = nodes[j_min]

    nodes[i_min] = TreeNode(
        children=[child_i, child_j],
        distances=[float(height - node_heights[i_min]),
                   float(height - node_heights[j_min])]
    )
    node_heights[i_min] = height
    active[j_min] = False

    # update distances
    for k in range(n0):
        if k == i_min or not active[k]:
            continue
        mean = ((D[i_min, k] * cluster_size[i_min] +
                 D[j_min, k] * cluster_size[j_min]) /
                (cluster_size[i_min] + cluster_size[j_min]))
        D[i_min, k] = D[k, i_min] = mean

    cluster_size[i_min] += cluster_size[j_min]
    remaining -= 1

# final root
for idx in range(n0 - 1, -1, -1):
    if active[idx]:
        return Tree(nodes[idx])

raise RuntimeError("UPGMA failed")

```

```python
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
        # total distances per active node
        total = np.zeros(D.shape[0], dtype=np.float64)
        for i in range(D.shape[0]):
            if not active[i]:
                continue
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

        # find minimum pair
        i_min, j_min = _find_min_pair_triangular(Q, active)
        if i_min == -1 or j_min == -1:
            break

        dist_ij = float(D[i_min, j_min])
        delta = 0.0
        if remaining > 2:
            delta = (total[i_min] - total[j_min]) / (remaining - 2)
        limb_i = 0.5 * (dist_ij + delta)
        limb_j = 0.5 * (dist_ij - delta)

        # CODON-SAFE: capture children first
        child_i = nodes[i_min]
        child_j = nodes[j_min]

        # create new node
        new_node = TreeNode(children=[child_i, child_j], distances=[float(limb_i), float(limb_j)])

        # update nodes and active array
        nodes[i_min] = new_node
        active[j_min] = False
        nodes[j_min] = None  # mark inactive safely

        # update distances: compute distances to new cluster
        mask = [k for k in range(D.shape[0]) if active[k] and k != i_min]
        for k in mask:
            D[i_min, k] = D[k, i_min] = 0.5 * (D[i_min, k] + D[j_min, k] - dist_ij)

        remaining -= 1

    # FINAL MERGE: safely combine remaining nodes
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
        raise RuntimeError("Neighbor-Joining failed: unexpected number of remaining nodes")

    return Tree(root)
```