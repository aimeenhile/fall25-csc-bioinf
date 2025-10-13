from os.path import join
import numpy as np
import pytest
import biotite
import biotite.sequence.phylo as phylo
from tests.util import data_dir

import numpy as pnp
from week3.biotite.tree import Tree, TreeNode
from week3.biotite.upgma import upgma
from week3.biotite.nj import neighbor_joining

@pytest.fixture
def distances():
    # Distances are based on the example
    # "Dendrogram of the BLOSUM62 matrix"
    # with the small modification M[i,j] += i+j
    # to reduce ambiguity in the tree construction.
    return np.loadtxt(join(data_dir("sequence"), "distances.txt"), dtype=int)


@pytest.fixture
def upgma_newick():
    # Newick notation of the tree created from 'distances.txt',
    # created via DendroUPGMA
    with open(join(data_dir("sequence"), "newick_upgma.txt"), "r") as file:
        newick = file.read().strip()
    return newick


@pytest.fixture
def tree(distances):
    return phylo.upgma(distances)


def test_distances(tree):
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert leaf.distance_to(tree.root) == dist
    # Example topological distances
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10


def test_upgma(tree, upgma_newick):
    """
    Compare the results of `upgma()` with DendroUPGMA.
    """
    ref_tree = phylo.Tree.from_newick(upgma_newick)
    # Cannot apply direct tree equality assertion because the distance
    # might not be exactly equal due to floating point rounding errors
    for i in range(len(tree)):
        for j in range(len(tree)):
            # Check for equal distances and equal topologies
            assert tree.get_distance(i, j) == pytest.approx(
                ref_tree.get_distance(i, j), abs=1e-3
            )
            assert tree.get_distance(i, j, topological=True) == ref_tree.get_distance(
                i, j, topological=True
            )


def test_neighbor_joining():
    """
    Compare the results of `neighbor_join()` with a known tree.
    """
    dist = np.array([
        [ 0,  5,  4,  7,  6,  8],
        [ 5,  0,  7, 10,  9, 11],
        [ 4,  7,  0,  7,  6,  8],
        [ 7, 10,  7,  0,  5,  9],
        [ 6,  9,  6,  5,  0,  8],
        [ 8, 11,  8,  9,  8,  0],
    ])  # fmt: skip

    ref_tree = phylo.Tree(
        phylo.TreeNode(
            [
                phylo.TreeNode(
                    [
                        phylo.TreeNode(
                            [
                                phylo.TreeNode(index=0),
                                phylo.TreeNode(index=1),
                            ],
                            [1, 4],
                        ),
                        phylo.TreeNode(index=2),
                    ],
                    [1, 2],
                ),
                phylo.TreeNode(
                    [
                        phylo.TreeNode(index=3),
                        phylo.TreeNode(index=4),
                    ],
                    [3, 2],
                ),
                phylo.TreeNode(index=5),
            ],
            [1, 1, 5],
        )
    )

    test_tree = phylo.neighbor_joining(dist)

    assert test_tree == ref_tree