import time
import numpy as np
from python import numpy as pnp
import numpy.pybridge

from bio_codon.codon_phylo import Tree, TreeNode, upgma, neighbor_joining

data_dir = "fall25-csc-bioinf/week3/"

distances: np.ndarray[int,2] = pnp.loadtxt("/data/distances.txt", dtype=pnp.int64)
with open("/data/newick_upgma.txt", "r") as f:
    newick_upgma = f.read().strip()
        


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
    print("Running tests under Codon")

    start_time = time.time()

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

    codon_runtime = int((time.time() - start_time) * 1000)

    print("----------------------------")
    print(f"Language    Runtime")
    print(f"-------------------")
    print(f"{'codon'}      {codon_runtime}ms")


if __name__ == "__main__":
    main()