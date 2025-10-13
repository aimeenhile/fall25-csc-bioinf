#!/usr/bin/env bash

set -euo pipefail

echo "=== Compiling Cython files in week3/biotite ==="
python - <<'EOF'
import numpy
from setuptools import setup, Extension
from Cython.Build import cythonize
import sys
import os

extensions = [
    Extension("biotite.tree", ["week3/biotite/tree.pyx"], include_dirs=[numpy.get_include()])
    Extension("biotite.nj", ["week3/biotite/nj.pyx"], include_dirs=[numpy.get_include()])
    Extension("biotite.upgma", ["week3/biotite/upgma.pyx"], include_dirs=[numpy.get_include()])
]

setup(
    name="biotite",
    ext_modules=cythonize(extensions, compiler_directives={'language_level': "3"}),
    script_args=["build_ext", "--inplace"]
)
EOF

echo "============================="
echo "🐍 Running Python tests"
echo "============================="
python3 week3/test_phylo.py

echo "============================="
echo "🧬 Running Codon tests"
echo "============================="
codon run -release -D IS_CODON=true week3/test_phylo.py