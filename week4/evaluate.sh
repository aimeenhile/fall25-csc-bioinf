#!/bin/bash

set -euo pipefail

cd week4

# Run Python version
printf "%-25s%-12s%s\n" "Method" "Language" "Runtime"
printf "%-25s%-12s%s\n" "-------------------------" "------------" "-------"

codon run -release codon/main.py data | while read -r line; do
    method=$(echo "$line" | cut -f1)
    language=$(echo "$line" | cut -f2)
    runtime=$(echo "$line" | cut -f3)
    printf "%-25s%-12s%s\n" "$method" "$language" "$runtime"
done
python python/main.py data | while read -r line; do
    method=$(echo "$line" | cut -f1)
    language=$(echo "$line" | cut -f2)
    runtime=$(echo "$line" | cut -f3)
    printf "%-25s%-12s%s\n" "$method" "$language" "$runtime"
done
