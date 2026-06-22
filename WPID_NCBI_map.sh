#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 interproscan_output.tsv WP_NCBI_map.tsv" >&2
    exit 1
fi

IN="$1"
OUT="$2"

if [ ! -s "$IN" ]; then
    echo "エラー: 入力edge tableが存在しない、または空です: $IN" >&2
    exit 1
fi

if [ -e "$OUT" ]; then
    echo "エラー: 出力ファイルが既に存在します: $OUT" >&2
    exit 1
fi

if ! command -v efetch >/dev/null 2>&1; then
    echo "エラー: efetch が見つかりません。NCBI Entrez Directをインストールしてください。" >&2
    exit 1
fi

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

IDS="$TMP_DIR/wp_ids.txt"
FASTA="$TMP_DIR/wp_products.faa"

# 今回はすべてのIDがWP_で始まること、headerは存在しないことを確認しているので、必要に応じて条件を変更してください。
# WP_XXXXXXXXX.X形式のIDを抽出して一意にする。
awk -F'\t' '
$1 ~ /^WP_[0-9]+\.[0-9]+$/ {
    print $1
}
' "$IN" | sort -u > "$IDS"

N="$(wc -l < "$IDS")"

if [ "$N" -eq 0 ]; then
    echo "WP_で始まるsource IDが見つかりませんでした。" >&2
    exit 1
fi

echo "WP ID数: $N" >&2
echo "NCBI ProteinからFASTA headerを取得します..." >&2

: > "$FASTA"

# 100件ずつNCBI Proteinから取得
split -l 100 "$IDS" "$TMP_DIR/chunk_"

for chunk in "$TMP_DIR"/chunk_*; do
    ids_csv="$(paste -sd, "$chunk")"

    efetch -db protein -format fasta -id "$ids_csv" >> "$FASTA"

    # NCBIへのリクエストが多すぎると一時的にアクセスがブロックされる可能性があるため、適度に待機する。
    sleep 0.4
done

awk '
BEGIN {
    OFS = "\t"
    print "protein_id", "ncbi_product"
}

# >で始まる行のみを処理,headerからIDとproduct/functionを抽出
/^>/ {
    header = substr($0, 2)

    # IDのみを抽出
    id = header
    sub(/ .*/, "", id)
    
    # headerからIDを削除して、残りをproduct/functionとする
    product = header
    sub(/^[^ ]+ /, "", product)

    # 末尾の [organism] を削除
    sub(/ \[[^]]+\]$/, "", product)

    if (product == "" || product == id) {
        product = "-"
    }

    print id, product
}
' "$FASTA" > "$OUT"

echo "完了: $OUT"