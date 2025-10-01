#!/bin/bash

# Add "-x" for debugging purposes; however, do not submit stuff with -x
set -euo pipefail

# The script starts at the repository root; move to week1
cd week2

# Run tests
echo "=========================================="
echo "Running Python tests..."
echo "=========================================="
python test.py > python_output.txt 2>&1
echo "Python tests completed."

echo ""
echo "=========================================="
echo "Running Codon tests..."
echo "=========================================="
codon run -release test.py > codon_output.txt 2>&1
echo "Codon tests completed."

# Compare results
echo "Comparing results..."
diff python_output.txt codon_output.txt || echo "Differences found between Python and Codon outputs."

echo "Evaluation completed."