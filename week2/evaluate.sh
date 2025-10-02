#!/bin/bash

# Exit immediately if a command exits with a non-zero status or on pipefail.
set -euo pipefail

# The script starts at the repository root 
# Add the repository root (PWD) to PYTHONPATH so 'bio_codon' package can be found.
export PYTHONPATH="${PWD}:${PYTHONPATH:-}"

# Change to the working directory 
cd week2

# Run tests
echo "=========================================="
echo "Running Python tests (BioPython)..."
echo "=========================================="
python test.py > python_output.txt 2>&1
echo "Python tests completed."

echo ""
echo "=========================================="
echo "Running Codon tests (bio_codon)..."
echo "=========================================="
# Codon uses the PYTHONPATH set above to locate bio_codon source
codon run -release test.py > codon_output.txt 2>&1
echo "Codon tests completed."

# --- Compare Results ---
echo "Comparing results..."
# Use diff and provide a clear success/failure message
if diff python_output.txt codon_output.txt; then
    echo "✅ SUCCESS: Python and Codon outputs are identical."
else
    echo "❌ FAILURE: Differences found between Python and Codon outputs."
    # Exit with error code 1 to fail the GitHub Action job
    exit 1
fi

echo "Evaluation completed."
