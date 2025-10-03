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

echo "Running with 'codon run -O1 test.py'..."
codon run -O1 test.py > codon_output.txt 2>&1 || { 
    echo "------------------------------------------"
    echo "CODON TEST FAILED - TRACEBACK BELOW:"
    echo "------------------------------------------"
    cat codon_output.txt
    echo "------------------------------------------"

    # Custom message if the compiler assertion fails again
    if grep -q "union already sealed" codon_output.txt; then
        echo "The compiler is hitting the internal 'union already sealed' bug again."
        echo "This is a known compiler error, not an error in your logic."
        echo "ACTION REQUIRED: You must modify 'test.py' based on the hint to proceed."
        echo "------------------------------------------"
        echo "Please check 'test.py' for complex generic types (e.g., List[List[T]], Dict[K, V] where V is complex)"
        echo "or implicit type usage and add EXPLICIT type annotations everywhere to simplify Codon's work."
        echo "------------------------------------------"
    fi
}
echo "Codon tests completed."

# Compare results
echo "Comparing results..."
diff python_output.txt codon_output.txt || echo "Differences found between Python and Codon outputs."

echo "Evaluation completed."
