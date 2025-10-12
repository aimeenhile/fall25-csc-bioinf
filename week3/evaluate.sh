#!/usr/bin/env bash

set -euo pipefail

echo "=== Compiling Cython files in week3/biotite ==="
python -m Cython.Build.cythonize -i week3/biotite/*.pyx

echo "============================="
echo "🐍 Running Python tests"
echo "============================="
python3 week3/test_phylo.py

echo "============================="
echo "🧬 Running Codon tests"
echo "============================="
codon run week3/bio_codon/test_phylo.py