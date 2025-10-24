#!/bin/bash

set -euo pipefail

cd week4

# Run Python version
echo -e "Method\t\tLanguage\tRuntime"
echo "--------------------------------------"
python python/main.py data >out.txt | while read -r line; do
    echo -e "$line"
done
codon run -release codon/main.py data >out.txt | while read -r line; do
    echo -e "$line"
done