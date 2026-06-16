#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 cytoscape_edges.tsv ompr_phob_edges.tsv [strict|loose]" >&2
    exit 1
fi

IN="$1"
OUT="$2"
MODE="${3:-strict}"

if [ ! -s "$IN" ]; then
    echo "エラー: 入力edge tableが存在しない、または空です: $IN" >&2
    exit 1
fi

if [ -e "$OUT" ]; then
    echo "エラー: 出力ファイルが既に存在します: $OUT" >&2
    exit 1
fi

                                                                              
#各DB由来のIDはあるが、IPR_ID未割り当ての場合が多い。
#ここでは、Interpro_ID, Pfam_ID, Interpro_description, 各DB_descriptionを正規表現に含めているが、
#明らかに　awk '{print $4}' ruminis.tsv | sort -u の結果からDBが多すぎて除外されるものが多すぎる。
#各DBのannotationをどこまで含めるか検討する. 

case "$MODE" in
    strict)
        # OmpR/PhoB-type DNA-binding domainを中心に抽出
        REGEX='ipr001867|pf00486|ompr|phob|ompR/phob'
        ;;
    loose)
        # receiver domainも含める。response regulator全般まで広がるので注意
        REGEX='ipr001867|pf00486|ipr001789|pf00072|ompr|phob|response regulator|receiver domain|chey-like'
        ;;
    *)
        echo "エラー: modeは strict または loose を指定してください: $MODE" >&2
        exit 1
        ;;
esac

TMP_PROTEINS="$(mktemp)"
TMP_OUT="$(mktemp)"

trap 'rm -f "$TMP_PROTEINS" "$TMP_OUT"' EXIT

# OmpR/PhoB関連edgeを持つprotein IDを抽出
awk -F'\t' -v re="$REGEX" '
BEGIN {
    OFS = "\t"
}
# headerをスキップ
NR == 1 {
    next
}

{
    line = tolower($0)
    # ~: 正規表現にマッチするlineの$1: NCBI_protein_id, WP_... を出力
    if (line ~ re) {
        print $1
    }
}
' "$IN" | sort -u > "$TMP_PROTEINS"

N_PROTEINS="$(wc -l < "$TMP_PROTEINS")"

if [ "$N_PROTEINS" -eq 0 ]; then
    echo "該当proteinが見つかりませんでした。modeを loose にするか、検索語を追加してください。" >&2
    exit 1
fi

# 1つ目のファイルの1列目にあるIDだけを使って、2つ目のTSVをフィルタする
# headerは上記で排除されているので、1行目から処理.
awk -F'\t' '
BEGIN {
    OFS = "\t"
}
# 1つ目のファイル 1列目をキーにして配列keepに保存
NR == FNR {
    keep[$1] = 1 
}
# 2つ目のファイルはheaderを出力してから、1列目がkeepにある行だけを出力
FNR == 1 {
    print # headerを出力
    next
}

$1 in keep {
    print 
}
' "$TMP_PROTEINS" "$IN" > "$TMP_OUT"

mv "$TMP_OUT" "$OUT"

PROTEIN_LIST="${OUT%.tsv}.proteins.txt"
cp "$TMP_PROTEINS" "$PROTEIN_LIST"

N_EDGES="$(($(wc -l < "$OUT") - 1))"

echo "完了: $OUT"
echo "抽出protein数: $N_PROTEINS"
echo "抽出edge数: $N_EDGES"
echo "protein list: $PROTEIN_LIST"
echo "mode: $MODE"