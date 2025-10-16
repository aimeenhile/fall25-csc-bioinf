Estimated time: 12h

# Steps:

## 1. Clone the source biotite file:
- Clone the four files: __init__.py, nj.pyx, tree.pyx, upgma.pyx from the source [[GitHub](https://github.com/biotite-dev/biotite/tree/v1.4.0/src/biotite/sequence/phylo)]

## 2. Set up test file:
- Find the relevant test under [[test_phylo.py](https://github.com/biotite-dev/biotite/blob/v1.4.0/tests/sequence/test_phylo.py)] in the source GitHub: def test_upgma(tree, upgma_newick), def test_neighbor_joining()
- Set up test file for Python tests and test files for Codon tests

## 3. Set up the Codon codes:
- Copy the source code to the Codon folders
- Set up working directory
- Modify the codon code, add typehints

## 4. Modify CI and evaluate.sh file 
### CI
- Add system dependencies so that CPython files are compiled 
- Add Python dependencies so that the source files are compiled and run correctly
### evaluate.sh
- Compile Cython files inside evaluate.sh and save the build files into the working directory

