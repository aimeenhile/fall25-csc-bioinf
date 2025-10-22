#!/bin/bash

set -euo pipefail

cd week4

# Path to data folder
DATA_DIR="./week4/data"

# Run Python version
echo -e "Method\tLanguage\tRuntime"
echo "--------------------------------------"
python main/python/main.py "$DATA_DIR" | while read -r line; do
    echo -e "$line"
done

