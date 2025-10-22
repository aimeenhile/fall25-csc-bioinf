#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")/week4"

# Path to data folder
DATA_DIR="data"

# Run Python version
echo -e "Method\tLanguage\tRuntime"
echo "--------------------------------------"
python python/main.py "$DATA_DIR" | while read -r line; do
    echo -e "$line"
done

