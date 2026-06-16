#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 edges.tsv WP_NCBI.tsv combined_node_labels.tsv" >&2
    exit 1
fi

EDGE_IN="$1"
MAP_IN="$2"
OUT="$3"

# Cytoscapeで長すぎるlabelになるのを防ぐ
MAX_LABEL_LEN=70

if [ ! -s "$EDGE_IN" ]; then
    echo "エラー: edge tableが存在しない、または空です: $EDGE_IN" >&2
    exit 1
fi

if [ ! -s "$MAP_IN" ]; then
    echo "エラー: WP-NCBI対応表が存在しない、または空です: $MAP_IN" >&2
    exit 1
fi

if [ -e "$OUT" ]; then
    echo "エラー: 出力ファイルが既に存在します: $OUT" >&2
    exit 1
fi

TMP="$(mktemp)"
trap 'rm -f "$TMP"' EXIT

awk -F'\t' -v MAX="$MAX_LABEL_LEN" '
BEGIN {
    OFS = "\t"
}

# 複数ファイルを読むときに、NRは全体の行数、FNRは各ファイルごとの行数になることを利用している。
# ここでは最初のファイルのみ処理
NR == FNR {
    # headerらしい行はskip
    if (FNR == 1 && $1 !~ /^WP_[0-9]+/) {
        next
    }

    id = clean_text($1)

    # 2列目以降を結合してproduct/functionとする
    product = ""
    for (i = 2; i <= NF; i++) {
        if (product == "") {
            product = $i
        } else {
            product = product " " $i
        }
    }

    product = clean_text(product)

　　 # WP_NCBI対応表にIDがあればproduct/functionを保存
    if (id != "" && id != "-" && product != "" && product != "-") {
        product = shorten(product, MAX)

        # versionあり
        product_map[id] = product

        # versionなしにも対応
        id_nover = id
        sub(/\.[0-9]+$/, "", id_nover)
        product_map[id_nover] = product
    }

    next
}


# 2つめのファイルの処理
FNR == 1 {
    next
}

{   
    source      = clean_text($1)
    target      = clean_text($2)
    interaction = clean_text($3)
    signature_accession = clean_text(getfield(5))
    signature_description = clean_text(getfield(6))
    interpro_accession = clean_text(getfield(7))
    interpro_description = clean_text(getfield(8))

    # source側: protein node
    if (source != "" && source != "-") {
        node_type[source] = "protein"
        source_nover = source
        sub(/\.[0-9]+$/, "", source_nover)

        product = "-"
        
        if (source in product_map) {
            product = product_map[source]
        }
        else if (source_nover in product_map) {
            product = product_map[source_nover]
        }

        if (product != "-") {
            product_short[source] = product
            node_label[source] = source " | " product
        }
        else {
            # WP_NCBI対応表にない場合のみInterProScanの説明でfallback
            if (!(source in node_label)) {
                if (signature_description != "" && signature_description != "-") {
                    node_label[source] = source " | " shorten(signature_description, MAX)
                }
                else if (interpro_description != "" && interpro_description != "-") {
                    node_label[source] = source " | " shorten(interpro_description , MAX)
                }
                else {
                    node_label[source] = source
                }
            }

            product_short[source] = "-"
        }

        # protein nodeに表示しない補助列として保持
        if (!(source in signature_short) && signature_description != "" && signature_description != "-") {
            signature_short[source] = shorten(signature_description, MAX)
        }

        if (!(source in interpro_short) && interpro_description != "" && interpro_description != "-") {
            interpro_short[source] = shorten(interpro_description, MAX)
        }
    }

    # target側: annotation node
    if (target != "" && target != "-") {
        if (!(target in node_type)) {
            if (interaction == "has_interpro") {
                node_type[target] = "interpro_domain"
            }
            else if (interaction == "has_signature") {
                node_type[target] = "signature_domain"
            }
            else if (interaction == "has_GO") {
                node_type[target] = "GO"
            }
            else if (interaction == "has_pathway") {
                node_type[target] = "pathway"
            }
            else {
                node_type[target] = "unknown"
            }
        }

        if (!(target in node_label)) {
            if (interaction == "has_interpro") {
                if (iprdesc != "" && iprdesc != "-") {
                    node_label[target] = target " | " shorten(iprdesc, MAX)
                } else {
                    node_label[target] = target
                }
            }
            else if (interaction == "has_signature") {
                if (sigdesc != "" && sigdesc != "-") {
                    node_label[target] = target " | " shorten(sigdesc, MAX)
                } else {
                    node_label[target] = target
                }
            }
            else {
                node_label[target] = target
            }
        }
    }
}

END {
    print "name", "node_label", "node_type", "product_short", "signature_short", "interpro_short"

    for (id in node_type) {
        ps = "-"
        ss = "-"
        is = "-"

        if (id in product_short) {
            ps = product_short[id]
        }

        if (id in signature_short) {
            ss = signature_short[id]
        }

        if (id in interpro_short) {
            is = interpro_short[id]
        }

        if (!(id in node_label) || node_label[id] == "") {
            node_label[id] = id
        }

        print id, node_label[id], node_type[id], ps, ss, is
    }
}

function getfield(i) {
    if (NF >= i) {
        return $i
    } else {
        return "-"
    }
}

function clean_text(s) {
    gsub(/\r/, "", s)
    gsub(/^[ \t]+/, "", s)
    gsub(/[ \t]+$/, "", s)
    gsub(/[ \t][ \t]+/, " ", s)

    sub(/^MULTISPECIES:[ \t]*/, "", s)
    sub(/[ \t]*\[[^]]+\]$/, "", s)

    if (s == "" || s == "-") {
        return "-"
    }

    return s
}

function shorten(s, maxlen) {
    if (s == "" || s == "-") {
        return "-"
    }

    if (length(s) > maxlen) {
        return substr(s, 1, maxlen - 3) "..."
    }

    return s
}
' "$MAP_IN" "$EDGE_IN" > "$TMP"

{
    head -n 1 "$TMP"
    tail -n +2 "$TMP" | sort -k3,3 -k1,1
} > "$OUT"

echo "完了: $OUT"
echo "edge table: $EDGE_IN"
echo "WP-NCBI map: $MAP_IN"
echo "最大node label長: $MAX_LABEL_LEN"