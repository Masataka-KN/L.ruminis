#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_interproscan.tsv output_edges.tsv" >&2
    exit 1
fi

IN="$1"
OUT="$2"

if [ ! -s "$IN" ]; then
    echo "エラー: 入力TSVが存在しない、または空です: $IN" >&2
    exit 1
fi

if [ -e "$OUT" ]; then
    echo "エラー: 出力ファイルが既に存在します: $OUT" >&2
    exit 1
fi

TMP="${OUT}.tmp"

# Header

#TO DO : remove dups
printf "source\ttarget\tinteraction\tevidence_db\tsignature_acc\tsignature_desc\tinterpro_acc\tinterpro_desc\tprotein_length\tstart\tend\tscore\tstatus\n" > "$TMP"

# 今回はinterproscanのoutputをtsvにしている。output formatは状況に合わせて適宜変更してください。
awk -F'\t' '
BEGIN {
    OFS = "\t"
}

# For SIF format, nodes and interactions only.
# interproscan output format上, --goterms: $14, --pathways: $15に出力される。
# ここでは13列以上あればOKとする。 MD5_checksumなどの列は除いた。
# https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html

# pathway_annotationはターゲットメインが真核生物の上に、莫大になるので削除しました。

NF >= 13 {
    protein_accession = $1
    sequence_length = $3
    database = $4
    signature_accession = $5
    signature_description = $6
    start_location = $7
    stop_location = $8
    e_value = $9
    status_t_true = $10
    interpro_accession = $12
    interpro_description = $13
    go_annotation = (NF >= 14 ? $14 : "-")
    pathways_annotation = (NF >= 15 ? $15 : "-")

    # Cytoscapeでのnetwork edgeのsource, target, interactionを出力
    # 後にhas_interpro OR has_signatureなどCytoscape側でFilteringできるようにした。

    if (interpro_accession != "" && interpro_accession != "-") {
    print protein_accession, interpro_accession, "has_interpro", database, signature_accession, signature_description, interpro_accession, interpro_description, sequence_length, start_location, stop_location, e_value, status_t_true
    }

    if (signature_accession != "" && signature_accession != "-") {
        print protein_accession, signature_accession, "has_signature", database, signature_accession, signature_description, interpro_accession, interpro_description, sequence_length, start_location, stop_location, e_value, status_t_true
    }

    if (go_annotation != "" && go_annotation != "-") {
    n = split(go_annotation, go, /\|/)  # GO annotationは複数ある場合があるので、|で分割してそれぞれ出力する
    for (i = 1; i <= n; i++) {
        g = go[i]
        gsub(/^[ \t]+|[ \t]+$/, "", g)

        while (match(g, /GO:[0-9]+/)) {
            goid = substr(g, RSTART, RLENGTH)

            print protein_accession, goid, "has_GO", database, signature_accession, signature_description, interpro_accession, interpro_description, sequence_length, start_location, stop_location, e_value, status_t_true

            # 見つけたGO IDを削って、次を探す
            g = substr(g, RSTART + RLENGTH)
        }
    }
    }
}
' "$IN" | sort -u >> "$TMP"

mv "$TMP" "$OUT"

echo "完了: $OUT"