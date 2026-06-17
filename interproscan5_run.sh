#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.faa output_interproscan.tsv" >&2
    exit 1
fi

IN_FILE="$1"
OUT_FILE="$2"

INTERPROSCAN="$HOME/data/interproscan-5.78-109.0/interproscan.sh"
CPU=16

if [ ! -s "$IN_FILE" ]; then
    echo "エラー: 入力faaファイルが存在しない、または空です: $IN_FILE" >&2
    exit 1
fi

if [ -e "$OUT_FILE" ]; then
    echo "エラー: 出力ファイルが既に存在します: $OUT_FILE" >&2
    exit 1
fi

echo "Running InterProScan..."
echo "Input : $IN_FILE"
echo "Output: $OUT_FILE"

bash "$INTERPROSCAN" \
    -i "$IN_FILE" \
    -f tsv \
    -o "$OUT_FILE" \
    --goterms \
    --pathways \
    --iprlookup \
    --cpu "$CPU"

echo "InterProScan completed."
echo "Output saved to: $OUT_FILE"