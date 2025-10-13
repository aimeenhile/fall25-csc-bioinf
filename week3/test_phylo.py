from __future__ import annotations
import time
#import numpy as np

IS_CODON = False

if IS_CODON:
    # Codon-compatible NumPy bridge
    from python import numpy as pnp
    import numpy.pybridge
    from bio_codon.codon_phylo import Tree, TreeNode, upgma, neighbor_joining
    def loadtxt(path: str) -> np.ndarray:
        return pnp.loadtxt(path, dtype=pnp.float64)
else:
    import numpy as pnp
    from biotite.sequence.phylo.tree import Tree, TreeNode
    from biotite.sequence.phylo.upgma import upgma
    from biotite.sequence.phylo.nj import neighbor_joining
    def loadtxt(path: str) -> np.ndarray:
        return pnp.loadtxt(path, dtype=pnp.float64)

def compare_trees(tree1, tree2, tol=1e-3):
    return tree1 == tree2 or abs(tree1.get_distance(0, 1) - tree2.get_distance(0, 1)) < tol


def test_upgma():
    dist = loadtxt("/data/distances.txt")
    tree = upgma(dist)
    print("UPGMA tree:", tree.to_newick())
    assert len(tree.leaves) == dist.shape[0]


def test_neighbor_joining():
    dist = loadtxt("/data/distances.txt")
    tree = neighbor_joining(dist)
    print("NJ tree:", tree.to_newick())
    assert len(tree.leaves) == dist.shape[0]


def main():
    print("Running tests under", "Codon" if IS_CODON else "Python", "runtime...")

    start_time = time.time()

    try:
        test_upgma()
        print("✅ UPGMA passed")
    except Exception as e:
        print("❌ UPGMA failed:", e)
    upgma_time = int((time.time() - start_time) * 1000)

    start_time = time.time()
    try:
        test_neighbor_joining()
        print("✅ Neighbor Joining passed")
    except Exception as e:
        print("❌ Neighbor Joining failed:", e)
    nj_time = int((time.time() - start_time) * 1000)

    print("----------------------------")
    print(f"Language    Runtime")
    print(f"-------------------")
    print(f"{'Codon' if IS_CODON else 'Python'}      {upgma_time + nj_time}ms")


if __name__ == "__main__":
    main()