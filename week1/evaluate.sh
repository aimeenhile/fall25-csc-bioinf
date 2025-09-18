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

# Helper function to run a command and capture N50 
run_and_capture_n50() {
    local lang=$1
    local cmd=$2
    local dataset=$3

    # Run the script and capture Markdown table output
    output=$($cmd "$DATA_DIR")

    # Extract the line corresponding to the dataset
    dataset_line=$(echo "$output" | grep -E "^\| .* $dataset .* \|")

    if [[ -z "$dataset_line" ]]; then
        echo "Error: Dataset $dataset not found in output"
        N50="NA"
        runtime_str="NA"
    else
        # Extract N50 (7th column)
        N50=$(echo "$dataset_line" | awk -F'|' '{gsub(/ /,"",$7); print $7}')
        [[ -z "$N50" ]] && N50="NA"

        # Extract runtime (last column)
        runtime_str=$(echo "$dataset_line" | awk -F'|' '{gsub(/ /,"",$NF); print $NF}')
        [[ -z "$runtime_str" ]] && runtime_str="NA"
    fi


    # Print nicely aligned table
    printf "%-${COL_DATASET}s %-${COL_LANG}s %-${COL_RUNTIME}s %-${COL_N50}s\n" "$dataset" "$lang" "$runtime_str" "$N50"
}

# Loop over datasets
for dataset in $(ls "$DATA_DIR" | sort); do
    if [ -d "$DATA_DIR/$dataset" ]; then
        # Run Python and Codon in parallel
        run_and_capture_n50 "python" "python $CODE_DIR/main.py" "$dataset" &
        pid_python=$!

        run_and_capture_n50 "codon" "codon run -release $CODE_DIR/main.codon.py" "$dataset" &
        pid_codon=$!

        # Wait for both to finish
        wait $pid_python
        wait $pid_codon
    fi
done