#!/usr/bin/env bash
set -euxo pipefail

# Root directories
CODE_DIR="./code"
DATA_DIR="./data"

# Column widths 
COL_DATASET=12
COL_LANG=10
COL_RUNTIME=10
COL_N50=8

# Output header
printf "%-${COL_DATASET}s %-${COL_LANG}s %-${COL_RUNTIME}s %-${COL_N50}s\n" "Dataset" "Language" "Runtime" "N50"
printf "%s\n" "-------------------------------------------------------------------------------------------------------"

# Helper function to run a command and capture N50 from its Markdown table output
run_and_capture_n50() {
    local lang=$1
    local cmd=$2
    local dataset=$3

    start=$(date +%s)

    # Capture the stdout of the script
    output=$($cmd "$DATA_DIR/$dataset")
    
    end=$(date +%s)
    runtime=$((end - start))

    # Extract N50 from Markdown table 
    N50=$(echo "$output" | grep -E "^\|.*$dataset.*\|" | awk -F'|' '{gsub(/ /,"",$5); print $5}')
    [[ -z "$N50" ]] && N50="NA"

    # Format runtime as h:mm:ss
    hh=$((runtime/3600))
    mm=$(( (runtime%3600)/60 ))
    ss=$((runtime%60))
    runtime_str=$(printf "%d:%02d:%02d" $hh $mm $ss)

    # Print nicely aligned output
    printf "%-${COL_DATASET}s %-${COL_LANG}s %-${COL_RUNTIME}s %-${COL_N50}s\n" "$dataset" "$lang" "$runtime_str" "$N50"
}

# Loop over datasets
for dataset in $(ls "$DATA_DIR" | sort); do
    if [ -d "$DATA_DIR/$dataset" ]; then
        # Run Python version
        run_and_capture_n50 "python" "python3 $CODE_DIR/main.py" "$dataset"

        # Run Codon version
        run_and_capture_n50 "codon" "codon run -release $CODE_DIR/main.codon.py" "$dataset"
    fi
done

