#!/bin/bash
set -euo pipefail

# -------- CONFIG --------
# P2Rank binary
P2RANK_BIN="${P2RANK_BIN:-$HOME/soft/p2rank_2.5/prank}"

# Dataset root (adjust per dataset, or override from env/CLI)
KINASE_FOLDER="${KINASE_FOLDER:-Kinase_Type_I}"

INPUT_DIR="${KINASE_FOLDER}/structures"
OUTPUT_DIR="${KINASE_FOLDER}/p2rank_predictions"
# ------------------------

mkdir -p "$OUTPUT_DIR"

echo "[i] Using P2Rank:     $P2RANK_BIN"
echo "[i] Input .cif dir:   $INPUT_DIR"
echo "[i] Output dir:       $OUTPUT_DIR"
echo

shopt -s nullglob
cifs=( "$INPUT_DIR"/*.cif )
if (( ${#cifs[@]} == 0 )); then
    echo "[!] No .cif files found in $INPUT_DIR"
    exit 1
fi

for cif_file in "${cifs[@]}"; do
    echo "[run] P2Rank on: $(basename "$cif_file")"
    "$P2RANK_BIN" predict -f "$cif_file" -o "$OUTPUT_DIR" -visualizations 0
done

echo
echo "All predictions finished. Output in: $OUTPUT_DIR"
