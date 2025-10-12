#!/usr/bin/env bash

set -euo pipefail

cd week3

echo "=== Compiling Cython files in week3/biotite ==="
python3 -m pip install --upgrade pip setuptools wheel cython > /dev/null 2>&1
cythonize -i /biotite/*.pyx

echo "============================="
echo "🐍 Running Python tests"
echo "============================="
python3 week3/test_phylo.py

echo "============================="
echo "🧬 Running Codon tests"
echo "============================="
export PATH=${HOME}/.codon/bin:$PATH
export CODON_RUNTIME=1
codon run week3/bio_codon/test_phylo.py