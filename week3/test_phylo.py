import time
import numpy as np

import sys
import os

from biotite.tree import Tree, TreeNode
from biotite.upgma import upgma
from biotite.nj import neighbor_joining

data_dir = os.path.join(os.path.dirname(__file__), "data")
dist_path = os.path.join(data_dir, "distances.txt")
newick_path = os.path.join(data_dir, "newick_upgma.txt")

distances = np.loadtxt(dist_path, dtype=np.int64)
with open(newick_path, "r") as f:
    newick_upgma = f.read().strip()


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
    print("Running tests under Python")

    start_time = time.time()

    test_upgma()
    test_neighbor_joining()

    py_runtime = int((time.time() - start_time) * 1000)

    print("----------------------------")
    print(f"Language    Runtime")
    print(f"-------------------")
    print(f"{'python'}      {py_runtime}ms")


if __name__ == "__main__":
    main()