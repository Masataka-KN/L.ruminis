#!/usr/bin/env bash
set -euo pipefail

# 今回はNCBI Datasets CLIから一括ダウンロードした.faaファイルを使用しました。
acidphilus_faa="$HOME/data/genome_data/acidphilus/ncbi_dataset/data/GCF_000011985.1/protein.faa"
ruminis_faa="$HOME/data/genome_data/ruminis/ncbi_dataset/data/GCF_000225845.1/protein.faa"

OUTPUT_DIR="$HOME/data/interproscan_output"

if [ -d "$OUTPUT_DIR" ]; then
    echo "エラー: ディレクトリが存在します: $OUTPUT_DIR" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

CPU=16

run_interproscan() {
    local input_faa="$1"
    local output_file="$2"

    echo "Running InterProScan on $input_faa..."

    bash $HOME/data/interproscan-5.78-109.0/interproscan.sh \
        -i "$input_faa" \
        -f tsv \
        -o "$output_file" \
        --goterms \
        --pathways \
        --iprlookup \
        --cpu "$CPU"
    
    echo "InterProScan completed for $input_faa. Output saved to $output_file."
}

# Run InterProScan for acidphilus　$1=*protein.faa, $2=出力ファイル名
run_interproscan "$acidphilus_faa" "$OUTPUT_DIR/acidphilus_interproscan.tsv"

# Run InterProScan for ruminis
run_interproscan "$ruminis_faa" "$OUTPUT_DIR/ruminis_interproscan.tsv"

echo "All InterProScan analyses completed."