from __future__ import annotations
import os, sys
import numpy as np

IS_CODON = os.environ.get("CODON_RUNTIME", "") == "1"

if IS_CODON:
    # Codon-compatible NumPy bridge
    from python import numpy as pnp
    import numpy.pybridge
    import phylo as phylo
    def loadtxt(path: str) -> np.ndarray:
        return pnp.loadtxt(path, dtype=pnp.float64)
else:
    import numpy as pnp
    import biotite.sequence.phylo as phylo
    def loadtxt(path: str) -> np.ndarray:
        return pnp.loadtxt(path, dtype=pnp.float64)

def compare_trees(tree1, tree2, tol=1e-3):
    return tree1 == tree2 or abs(tree1.get_distance(0, 1) - tree2.get_distance(0, 1)) < tol


def test_upgma():
    dist = loadtxt("/data/distances.txt")
    tree = phylo.upgma(dist)
    print("UPGMA tree:", tree.to_newick())
    assert len(tree.leaves) == dist.shape[0]


def test_neighbor_joining():
    dist = loadtxt("/data/distances.txt")
    tree = phylo.neighbor_joining(dist)
    print("NJ tree:", tree.to_newick())
    assert len(tree.leaves) == dist.shape[0]


def main():
    print("Running tests under", "Codon" if IS_CODON else "Python", "runtime...")
    try:
        test_upgma()
        print("✅ UPGMA passed")
    except Exception as e:
        print("❌ UPGMA failed:", e)
    try:
        test_neighbor_joining()
        print("✅ Neighbor Joining passed")
    except Exception as e:
        print("❌ Neighbor Joining failed:", e)


if __name__ == "__main__":
    main()