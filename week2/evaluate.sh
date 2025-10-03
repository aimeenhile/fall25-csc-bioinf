#!/bin/bash

# Add "--failfast" to python test.py for debugging purposes:
# python test.py --failfast > python_output.txt 2>&1

set -euo pipefail

# The script starts at the repository root; move to week2
cd week2

# Run tests
echo "=========================================="
echo "Running Python tests (BioPython)..."
echo "=========================================="

python test.py > python_output.txt 2>&1 || { 
    echo "------------------------------------------"
    echo "PYTHON TEST FAILED - TRACEBACK BELOW:"
    echo "------------------------------------------"
    cat python_output.txt
    echo "------------------------------------------"
    exit 1 
}
echo "Python tests completed successfully."

echo ""
echo "=========================================="
echo "Running Codon tests..."
echo "=========================================="

codon run --no-opt test.py > codon_output.txt 2>&1 || { 
    echo "------------------------------------------"
    echo "CODON TEST FAILED - TRACEBACK BELOW:"
    echo "------------------------------------------"
    cat codon_output.txt
    echo "------------------------------------------"
    exit 1 
}
echo "Codon tests completed."

# Compare results
echo "Comparing results..."
diff python_output.txt codon_output.txt || echo "Differences found between Python and Codon outputs."

echo "Evaluation completed."
