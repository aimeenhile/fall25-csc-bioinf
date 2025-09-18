#!/usr/bin/env bash
set -euo pipefail

# Root directories
CODE_DIR="./code"
DATA_DIR="./data"

# Column widths
COL_DATASET=12
COL_LANG=10
COL_RUNTIME=10
COL_N50=8

# CPU detection
CPU_CORES=$(nproc)
MAX_JOBS=$(( CPU_CORES > 4 ? 4 : CPU_CORES < 2 ? 2 : CPU_CORES ))

# Temporary directory for outputs
TMPDIR=$(mktemp -d)
declare -A JOB_FILES

# Start all jobs (each dataset separately)
for dataset in $(ls "$DATA_DIR" | sort); do
    if [ -d "$DATA_DIR/$dataset" ]; then
        # Python
        py_out="$TMPDIR/$dataset.python.out"
        python "$CODE_DIR/main.py" "$DATA_DIR/$dataset" > "$py_out" 2>&1 &
        JOB_FILES["$dataset:python"]=$py_out

        # Codon
        codon_out="$TMPDIR/$dataset.codon.out"
        codon run -release "$CODE_DIR/main.codon.py" "$DATA_DIR/$dataset" > "$codon_out" 2>&1 &
        JOB_FILES["$dataset:codon"]=$codon_out

        # Throttle parallel jobs
        while [[ $(jobs -rp | wc -l) -ge $MAX_JOBS ]]; do
            sleep 1
        done
    fi
done

# Wait for all background jobs to finish
wait

# --- Helper to extract N50 and runtime ---
extract_metrics() {
    local file="$1"
    local dataset="$2"
    # Grab the first Markdown table line matching the dataset
    local line
    line=$(grep -E "^\| .* $dataset .* \|" "$file" | head -n1)
    if [[ -z "$line" ]]; then
        echo "NA NA"
    else
        # Column 7=N50, Column 10=Runtime
        local N50 runtime
        N50=$(echo "$line" | awk -F'|' '{gsub(/ /,"",$7); print $7}')
        runtime=$(echo "$line" | awk -F'|' '{gsub(/ /,"",$10); print $10}')
        [[ -z "$N50" ]] && N50="NA"
        [[ -z "$runtime" ]] && runtime="NA"
        echo "$N50 $runtime"
    fi
}

# --- Print Full Markdown table from all outputs ---
echo "=== Full Markdown Table ==="
for dataset in $(ls "$DATA_DIR" | sort); do
    for lang in python codon; do
        out_file="${JOB_FILES[$dataset:$lang]}"
        grep -E "^\|" "$out_file" || true
    done
done

# --- Print Summary Table ---
echo ""
echo "=== Summary Table (Dataset Language Runtime N50) ==="
printf "%-${COL_DATASET}s %-${COL_LANG}s %-${COL_RUNTIME}s %-${COL_N50}s\n" "Dataset" "Language" "Runtime" "N50"
printf "%s\n" "--------------------------------------------------------"

for dataset in $(ls "$DATA_DIR" | sort); do
    for lang in python codon; do
        out_file="${JOB_FILES[$dataset:$lang]}"
        read N50 runtime <<< $(extract_metrics "$out_file" "$dataset")
        printf "%-${COL_DATASET}s %-${COL_LANG}s %-${COL_RUNTIME}s %-${COL_N50}s\n" "$dataset" "$lang" "$runtime" "$N50"
    done
done

# Cleanup
rm -rf "$TMPDIR"
