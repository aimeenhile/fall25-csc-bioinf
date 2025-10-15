from __future__ import annotations
import time
#import numpy as np

IS_CODON = False

if IS_CODON:
    # Codon-compatible NumPy bridge
    from python import numpy as pnp
    import numpy.pybridge
    from bio_codon.codon_phylo import Tree, TreeNode, upgma, neighbor_joining

    data_dir = "fall25-csc-bioinf/week3/"

    distances: np.ndarray[int,2] = pnp.loadtxt("/data/distances.txt", dtype=pnp.int64)
    with open("/data/newick_upgma.txt", "r") as f:
        newick_upgma = f.read().strip()

    """
    
    def loadtxt(path: str) -> pnp.ndarray:
        with open(path, "r") as f:
            lines = f.readlines()
        mat = []
        for line in lines:
            row = [float(x) for x in line.strip().split()]
            mat.append(row)
        return pnp.array(mat, dtype=pnp.float64)
    """

else:
    import numpy as pnp
    import sys, os
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "biotite")))

    from biotite.tree import Tree, TreeNode
    from biotite.upgma import upgma
    from biotite.nj import neighbor_joining

    data_dir = os.path.join(os.path.dirname(__file__), "data")
    dist_path = os.path.join(data_dir, "distances.txt")
    newick_path = os.path.join(data_dir, "newick_upgma.txt")

    distances = pnp.loadtxt(dist_path, dtype=pnp.int64)
    with open(newick_path, "r") as f:
        newick_upgma = f.read().strip()

    """
        def loadtxt(path: str) -> np.ndarray:
        return pnp.loadtxt(path, dtype=pnp.float64)
    """

if IS_CODON:
    orig_from_newick = TreeNode.from_newick

    @staticmethod
    def safe_from_newick(newick: str) -> TreeNode:
        # replace labels like '29.0' -> '29'
        import re
        def fix_label(match):
            lbl = match.group(0)
            try:
                i = int(float(lbl))
                return str(i)
            except:
                return lbl
        newick_fixed = re.sub(r"[0-9]+\.[0-9]+", fix_label, newick)
        return orig_from_newick(newick_fixed)

    TreeNode.from_newick = safe_from_newic
    

def compare_trees(tree1, tree2, tol=1e-3):
    #return tree1 == tree2 or abs(tree1.get_distance(0, 1) - tree2.get_distance(0, 1)) < tol
    n = len(tree1.leaves)
    for i in range(n):
        for j in range(n):
            if abs(tree1.get_distance(i, j) - tree2.get_distance(i, j)) > tol:
                return False
            if tree1.get_distance(i, j, topological=True) != tree2.get_distance(i, j, topological=True):
                return False
    return True


def test_upgma():
    tree = upgma(distances)
    print("UPGMA tree:", tree.to_newick())
    assert len(tree.leaves) == distances.shape[0]

    ref_tree = Tree.from_newick(newick_upgma)
    assert compare_trees(tree, ref_tree)


def test_neighbor_joining():
    tree = neighbor_joining(distances)
    print("NJ tree:", tree.to_newick())
    assert len(tree.leaves) == distances.shape[0]


def main():
    print("Running tests under", "Codon" if IS_CODON else "Python", "runtime...")

    start_time1 = time.time()

    try:
        test_upgma()
        print("✅ UPGMA passed")
    except Exception as e:
        print("❌ UPGMA failed:", e)
    upgma_time = int((time.time() - start_time1) * 1000)

    start_time2 = time.time()
    try:
        test_neighbor_joining()
        print("✅ Neighbor Joining passed")
    except Exception as e:
        print("❌ Neighbor Joining failed:", e)
    nj_time = int((time.time() - start_time2) * 1000)

    print("----------------------------")
    print(f"Language    Runtime")
    print(f"-------------------")
    print(f"{'Codon' if IS_CODON else 'Python'}      {upgma_time + nj_time}ms")


if __name__ == "__main__":
    main()