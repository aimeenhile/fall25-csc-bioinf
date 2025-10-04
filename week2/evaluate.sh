#!/bin/bash

# Add "--failfast" to python test.py for debugging purposes:
# python bio_codon/motifs/test.py --failfast > python_output.txt 2>&1

set -euo pipefail

# The script starts at the repository root; move to week2
cd week2

# Run tests
echo "=========================================="
echo "Running Python tests (BioPython)..."
echo "=========================================="

python bio_codon/motifs/test_python.py > python_output.txt 2>&1 || { 
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

codon run bio_codon/motifs/test_codon.py > codon_output.txt 2>&1 || { 
    echo "------------------------------------------"
    echo "CODON TEST FAILED - TRACEBACK BELOW:"
    echo "------------------------------------------"
    cat codon_output.txt
    echo "------------------------------------------"
    exit 1
}
echo "Codon tests completed."

# Compare results
echo ""
echo "Comparing results..."
diff python_output.txt codon_output.txt || echo "Differences found between Python and Codon outputs."

echo "Evaluation completed."
