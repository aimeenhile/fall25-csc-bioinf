#!/bin/bash

set -euo pipefail

cd week4

# Run Python version
printf "%-25s%-12s%s\n" "Method" "Language" "Runtime"
printf "%-25s%-12s%s\n" "-------------------------" "------------" "-------"

codon run -release codon/main.py data 

python python/main.py data 