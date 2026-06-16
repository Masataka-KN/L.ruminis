#!/usr/bin/env python3
import sys
from pathlib import Path

MAX_LABEL_LEN = 70


def clean_text(s):
    s = str(s).replace("\r", "").strip()
    while "  " in s:
        s = s.replace("  ", " ")

    if s.startswith("MULTISPECIES:"):
        s = s.replace("MULTISPECIES:", "", 1).strip()

    if s.endswith("]") and "[" in s:
        s = s[:s.rfind("[")].strip()

    return s if s else "-"


def shorten(s, maxlen=MAX_LABEL_LEN):
    if not s or s == "-":
        return "-"
    return s[:maxlen - 3] + "..." if len(s) > maxlen else s


def remove_version(wp_id):
    if "." in wp_id:
        return wp_id.rsplit(".", 1)[0]
    return wp_id


if len(sys.argv) != 4:
    print("Usage: script.py edges.tsv WP_NCBI.tsv combined_node_labels.tsv", file=sys.stderr)
    sys.exit(1)

edge_file = Path(sys.argv[1])
map_file = Path(sys.argv[2])
out_file = Path(sys.argv[3])

if not edge_file.is_file() or edge_file.stat().st_size == 0:
    sys.exit(f"エラー: edge tableが存在しない、または空です: {edge_file}")

if not map_file.is_file() or map_file.stat().st_size == 0:
    sys.exit(f"エラー: WP-NCBI対応表が存在しない、または空です: {map_file}")

if out_file.exists():
    sys.exit(f"エラー: 出力ファイルが既に存在します: {out_file}")


# 1. WP_ID: product の辞書を作る
product_map = {}

with map_file.open() as f:
    for line_no, line in enumerate(f, start=1):
        cols = line.rstrip("\n").split("\t")

        # header行はスキップ
        if line_no == 1 and not cols[0].startswith("WP_"):
            continue

        wp_id = clean_text(cols[0])
        # 複数列ある場合は結合してproductとする
        product = clean_text(" ".join(cols[1:]))

        if wp_id != "-" and product != "-":
            product = shorten(product)
            product_map[wp_id] = product
            product_map[remove_version(wp_id)] = product


node_info = {}

def add_node(name, node_type, label=None, product="-", signature="-", interpro="-"):
    if name == "-":
        return

    if name not in node_info:
        node_info[name] = {
            "node_label": label or name,
            "node_type": node_type,
            "product_short": product,
            "signature_short": signature,
            "interpro_short": interpro,
        }

    else:
        if label and node_info[name]["node_label"] == name:
            node_info[name]["node_label"] = label

        if product != "-":
            node_info[name]["product_short"] = product
        if signature != "-":
            node_info[name]["signature_short"] = signature
        if interpro != "-":
            node_info[name]["interpro_short"] = interpro


# 2. edge_tableを読んでnode tableを作る
with edge_file.open() as f:
    header = next(f, None)

    for line in f:
        cols = line.rstrip("\n").split("\t")

        source = clean_text(cols[0]) if len(cols) > 0 else "-"
        target = clean_text(cols[1]) if len(cols) > 1 else "-"
        interaction = clean_text(cols[2]) if len(cols) > 2 else "-"

        signature_description = clean_text(cols[5]) if len(cols) > 5 else "-"
        interpro_description = clean_text(cols[7]) if len(cols) > 7 else "-"

        # source protein
        product = product_map.get(source, product_map.get(remove_version(source), "-"))

        if product != "-":
            label = f"{source} | {product}"
        elif signature_description != "-":
            label = f"{source} | {shorten(signature_description)}"
        elif interpro_description != "-":
            label = f"{source} | {shorten(interpro_description)}"
        else:
            label = source

        add_node(
            source,
            "protein",
            label=label,
            product=product,
            signature=shorten(signature_description),
            interpro=shorten(interpro_description),
        )

        # target annotation
        if interaction == "has_interpro":
            node_type = "interpro_domain"
            desc = interpro_description
        elif interaction == "has_signature":
            node_type = "signature_domain"
            desc = signature_description
        elif interaction == "has_GO":
            node_type = "GO"
            desc = "has_GO"
        elif interaction == "has_pathway":
            node_type = "pathway"
            desc = "has_pathway"
        else:
            node_type = "unknown"
            desc = "-"

        if desc != "-":
            label = f"{target} | {shorten(desc)}"
        else:
            label = target

        add_node(target, node_type, label=label)


# 3. node_infoを出力する
with out_file.open("w") as out:
    out.write("name\tnode_label\tnode_type\tproduct_short\tsignature_short\tinterpro_short\n")

    # node_infoをnode_typeとnameの両方でソートして出力
    for name in sorted(node_info, key=lambda x: (node_info[x]["node_type"], x)):
        d = node_info[name]
        # 1行にまとめて出力
        out.write(
            f"{name}\t{d['node_label']}\t{d['node_type']}\t"
            f"{d['product_short']}\t{d['signature_short']}\t{d['interpro_short']}\n"
        )

print(f"完了: {out_file}")