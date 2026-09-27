#!/usr/bin/env bash

set -u
set -o pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

TOOL_DIR="$PROJECT/tools/UBCG_ver2"
META="$PROJECT/metadata/photobacterium_primary_cohort_2026-09-12.tsv"
DATA_DIR="$PROJECT/data/raw/photobacterium_primary_500_ncbi/ncbi_dataset/data"

OUT_DIR="$PROJECT/results/ubcg_500"
LOG_DIR="$PROJECT/logs/ubcg_500"

THREADS="${THREADS:-8}"
LIMIT="${LIMIT:-0}"

mkdir -p "$OUT_DIR" "$LOG_DIR"

echo "============================================================"
echo "Photobacterium UBCG2 production extraction"
echo "============================================================"
echo "Project : $PROJECT"
echo "Metadata: $META"
echo "Output  : $OUT_DIR"
echo "Threads : $THREADS"
echo "Limit   : $LIMIT (0 = all genomes)"
echo

# ------------------------------------------------------------
# Basic file checks
# ------------------------------------------------------------

for required in \
    "$META" \
    "$TOOL_DIR/UBCG2.jar" \
    "$TOOL_DIR/hmm/ubcg_v2.hmm"
do
    if [ ! -e "$required" ]; then
        echo "ERROR: required file not found:"
        echo "$required"
        exit 1
    fi
done

# ------------------------------------------------------------
# Read primary-cohort accessions from metadata
# ------------------------------------------------------------

mapfile -t ACCESSIONS < <(
python - "$META" <<'PY'
import sys
import pandas as pd

f = sys.argv[1]

df = pd.read_csv(f, sep="\t", dtype=str)

col = "Assembly Accession"

if col not in df.columns:
    raise SystemExit(f"ERROR: column '{col}' not found")

x = (
    df[col]
    .dropna()
    .astype(str)
    .str.strip()
)

for acc in x:
    print(acc)
PY
)

TOTAL=${#ACCESSIONS[@]}
UNIQUE=$(printf '%s\n' "${ACCESSIONS[@]}" | sort -u | wc -l)

echo "Metadata rows/accessions : $TOTAL"
echo "Unique accessions        : $UNIQUE"

if [ "$TOTAL" -ne 500 ] || [ "$UNIQUE" -ne 500 ]; then
    echo
    echo "ERROR: expected exactly 500 unique accessions."
    exit 1
fi

echo "Primary cohort check     : PASS"
echo

# ------------------------------------------------------------
# Check whether an existing UCG is valid JSON
# ------------------------------------------------------------

valid_ucg () {
    local file="$1"

    [ -s "$file" ] || return 1

    python - "$file" >/dev/null 2>&1 <<'PY'
import json
import sys

f = sys.argv[1]

with open(f) as handle:
    d = json.load(handle)

assert isinstance(d, dict)
assert isinstance(d.get("run_info"), dict)
assert isinstance(d.get("data"), dict)
assert len(d["data"]) > 0
PY
}

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

FAIL_FILE="$LOG_DIR/failures_latest.tsv"

printf "Assembly_Accession\tReason\n" > "$FAIL_FILE"

success=0
skipped=0
failed=0
attempted=0
processed_this_run=0

for i in "${!ACCESSIONS[@]}"
do
    ACC="${ACCESSIONS[$i]}"
    NUMBER=$((i + 1))

    if [ "$LIMIT" -gt 0 ] && [ "$processed_this_run" -ge "$LIMIT" ]; then
        echo
        echo "LIMIT=$LIMIT reached."
        break
    fi

    UCG="$OUT_DIR/${ACC}.ucg"
    LOG="$LOG_DIR/${ACC}.log"

    echo
    echo "============================================================"
    echo "[$NUMBER/$TOTAL] $ACC"
    echo "============================================================"

    # Resume support
    if valid_ucg "$UCG"; then
        echo "STATUS: already complete — skipping"
        skipped=$((skipped + 1))
        continue
    fi

    GENOME=$(find "$DATA_DIR/$ACC" \
        -maxdepth 1 \
        -type f \
        -name '*_genomic.fna' \
        -print \
        -quit)

    if [ -z "$GENOME" ]; then
        echo "STATUS: FAILED — genome FASTA not found"
        printf "%s\tGenome_FASTA_not_found\n" "$ACC" >> "$FAIL_FILE"
        failed=$((failed + 1))
        processed_this_run=$((processed_this_run + 1))
        continue
    fi

    echo "Genome: $GENOME"
    echo "Log   : $LOG"

    attempted=$((attempted + 1))
    processed_this_run=$((processed_this_run + 1))

    (
        cd "$TOOL_DIR" || exit 1

        java -jar UBCG2.jar \
            -i "$GENOME" \
            -ucg_dir "$OUT_DIR" \
            -label "$ACC" \
            -hmm hmm/ubcg_v2.hmm \
            -t "$THREADS"

    ) > "$LOG" 2>&1

    EXITCODE=$?

    if [ "$EXITCODE" -eq 0 ] && valid_ucg "$UCG"; then
        echo "STATUS: SUCCESS"
        success=$((success + 1))
    else
        echo "STATUS: FAILED"
        echo "See: $LOG"
        printf "%s\tUBCG2_failed_or_invalid_output\n" "$ACC" >> "$FAIL_FILE"
        failed=$((failed + 1))
    fi
done

# ------------------------------------------------------------
# Count all currently valid production UCG files
# ------------------------------------------------------------

VALID_TOTAL=0

for ACC in "${ACCESSIONS[@]}"
do
    if valid_ucg "$OUT_DIR/${ACC}.ucg"; then
        VALID_TOTAL=$((VALID_TOTAL + 1))
    fi
done

echo
echo "============================================================"
echo "RUN SUMMARY"
echo "============================================================"
echo "Expected cohort        : $TOTAL"
echo "Attempted this run     : $attempted"
echo "Successful this run    : $success"
echo "Skipped already valid  : $skipped"
echo "Failed this run        : $failed"
echo "Valid UCG files total  : $VALID_TOTAL"
echo
echo "Failure table:"
echo "$FAIL_FILE"
echo "============================================================"
