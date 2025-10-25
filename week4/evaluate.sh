#!/bin/bash

set -euo pipefail

cd week4

# Run Python version
echo f"{'Method':<25}{'Language':<12}{'Runtime'}"
echo "-------------------------------------------------"
codon run -release codon/main.py data | while read -r line; do
    echo -e "$line"
done
python python/main.py data | while read -r line; do
    echo -e "$line"
done
