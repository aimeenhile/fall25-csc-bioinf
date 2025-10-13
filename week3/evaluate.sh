#!/usr/bin/env bash

set -euo pipefail

echo "=== Compiling Cython files in week3/biotite ==="
python - <<'EOF'
import numpy, sys
from Cython.Build import cythonize
from setuptools import Extension

extensions = cythonize([
    Extension("biotite.nj", ["week3/biotite/nj.pyx"], include_dirs=[numpy.get_include()]),
    Extension("biotite.tree", ["week3/biotite/tree.pyx"], include_dirs=[numpy.get_include()]),
    Extension("biotite.upgma", ["week3/biotite/upgma.pyx"], include_dirs=[numpy.get_include()])
])
EOF

echo "============================="
echo "🐍 Running Python tests"
echo "============================="
python3 week3/test_phylo.py

echo "============================="
echo "🧬 Running Codon tests"
echo "============================="
codon run week3/bio_codon/test_phylo.py