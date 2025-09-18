#!/usr/bin/env bash
set -euxo pipefail

# Path to your datasets
DATASETS=("data/data1")  # add more datasets if needed

# Output header
echo -e "Dataset\tLanguage\tRuntime\tN50"
echo "-------------------------------------------------------------------------------------------------------"

for DATA in "${DATASETS[@]}"; do
    # Python run
    START=$(date +%s)
    PYTHON_OUTPUT=$(python3 code/main.py "$DATA")
    END=$(date +%s)
    PYTHON_RUNTIME=$((END - START))
    
    # Extract N50 (assuming last line of output contains the N50)
    PYTHON_N50=$(echo "$PYTHON_OUTPUT" | tail -n1)
    
    # Format runtime as HH:MM:SS
    PYTHON_RUNTIME_FMT=$(printf '%02d:%02d:%02d' $((PYTHON_RUNTIME/3600)) $((PYTHON_RUNTIME%3600/60)) $((PYTHON_RUNTIME%60)))

    echo -e "$(basename "$DATA")\tpython\t$PYTHON_RUNTIME_FMT\t$PYTHON_N50"

    # Codon run
    START=$(date +%s)
    CODON_OUTPUT=$(codon run -release code/main.codon.py "$DATA")
    END=$(date +%s)
    CODON_RUNTIME=$((END - START))

    # Extract N50 similarly (last line of output)
    CODON_N50=$(echo "$CODON_OUTPUT" | tail -n1)
    CODON_RUNTIME_FMT=$(printf '%02d:%02d:%02d' $((CODON_RUNTIME/3600)) $((CODON_RUNTIME%3600/60)) $((CODON_RUNTIME%60)))

    echo -e "$(basename "$DATA")\tcodon\t$CODON_RUNTIME_FMT\t$CODON_N50"
done
Rank	Nickname	Submission Time	Submission Count	Genome_Fraction(%)	Duplication ratio	NGA50	Misassemblies	Mismatches per 100kbp
20	crayon	2019/06/18 11:27:47pm	22	99.886	1.982	9118.8	2.0	0
18	crayon	2019/06/18 11:29:13pm	7	99.942	2.0052	9129.2	3.0	0
20	crayon	2019/06/18 11:29:59pm	9	78.6	1.6	7859.2	0.0	0
5	crayon	2019/06/18 11:30:42pm	10	78.2948	1.607	55757.8	11.0