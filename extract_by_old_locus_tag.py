#!/usr/bin/env python3
import csv
import re
import sys
from urllib.parse import unquote
from collections import defaultdict

# csvが "...", "..."の形式になっているので、"などを除去
def clean_cell(x):
    if x is None:
        return ""
    return (
        x.strip()
         .strip('"')
         .strip("'")
         .strip("“")
         .strip("”")
         .strip()
    )


# gff file 9列目にNCBI IDのWP_...、locus_tag, old_lous_tagが";"区切りで"..."="..."でまとまっているので、
# それらを辞書形式で保存.
def parse_gff_attributes(attr_text):
    attrs = {}

    for item in attr_text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)

            # GFFではURL encodingされている場合があるためunquoteを追加.
            key = clean_cell(unquote(key))
            value = clean_cell(unquote(value))

            attrs[key] = value

    return attrs


# locus_tagが複数ある場合に対応させました.
def split_tags(x):
    if not x:
        return []

    tags = []
    for v in x.split(","):
        v = clean_cell(v)
        if v:
            tags.append(v)

    return tags


# csvからold_locus_tagを抽出.
def load_target_locus_tags(csv_file):
    
    target_locus = set()
    
    # Excelで作成したcsvはBOM付きUTF-8の可能性があるので、encoding="utf-8-sig"を指定しました.
    with open(csv_file, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, skipinitialspace=True)

        if reader.fieldnames is None:
            sys.exit("エラー: CSVが空です")

        reader.fieldnames = [clean_cell(h) for h in reader.fieldnames]

        if "Locus" not in reader.fieldnames:
            sys.exit(
                'エラー: CSVに "Locus" 列がありません。列名を確認してください。\n'
                f"見つかった列名: {reader.fieldnames}"
            )

        for row in reader:
            locus = clean_cell(row.get("Locus", ""))

            if locus:
                target_locus.add(locus)

    return target_locus


def extract_wp_id_from_attrs(attrs):
    
    # gff 9列目 Name="..."にはそのままWP_...があるが、今回はID=cds-WP_...に合わせてIDからも拾えるようにしました。
    # 他のKeysでも書き換えれば対応できます。
    candidate_keys = ["ID", "protein_id", "Name", "Dbxref"]

    for key in candidate_keys:
        if key in attrs:
            value = attrs[key]

            if value.startswith("cds-WP_"):
                return value.replace("cds-", "", 1)


            m = re.search(r"WP_\d+\.\d+", value)
            if m:
                return m.group(0)

    joined = ";".join(f"{k}={v}" for k, v in attrs.items())
    m = re.search(r"WP_\d+\.\d+", joined)
    if m:
        return m.group(0)

    return None


# gffのProtein=gene-NQ...などからfeature prefixを除く
def normalize_gff_id(x):

    x = clean_cell(x)

    for prefix in ["gene-", "cds-"]:
        if x.startswith(prefix):
            return x.replace(prefix, "", 1)

    return x


# protein.faaと対応させるため、今回はfeature_type=CDSのみ使いますが
# 他のfeatureの結果も見れるように対応させました。
def extract_wp_ids_from_gff(gff_file, target_locus):
 
    # feature_id -> locus_tag の対応. 複数対応の可能性があるのでsetで保存しました.
    locus_feature_ids = defaultdict(set)

    locus_feature_rows = set()
    wp_ids = set()
    locus_to_wp = defaultdict(set)
    matched_locus_in_gff = set()

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            feature_type = cols[2]
            attrs = parse_gff_attributes(cols[8])
            
            # 今回はold_locus_tagを使いますが、再現性の高さからlocus_tagでもヒットするようにしています。
            tags = []
            tags.extend(split_tags(attrs.get("old_locus_tag", "")))
            tags.extend(split_tags(attrs.get("locus_tag", "")))

            matched_tags = [t for t in tags if t in target_locus]

            if not matched_tags:
                continue

            matched_locus_in_gff.update(matched_tags)

            if "ID" in attrs:
                raw_feature_id = attrs["ID"]
                normalized_feature_id = normalize_gff_id(raw_feature_id)

                for tag in matched_tags:
                    locus_feature_ids[normalized_feature_id].add(tag)
                    locus_feature_rows.add(
                        (
                            feature_type,
                            raw_feature_id,
                            normalized_feature_id,
                            tag
                        )
                    )

    # CDSのみ
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"): # gff file 最初のコメント行をスキップ.
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            feature_type = cols[2]
            if feature_type != "CDS":
                continue

            attrs = parse_gff_attributes(cols[8])

            wp_id = extract_wp_id_from_attrs(attrs)
            if not wp_id:
                continue

            # 今回はold_locus_tagを使いますが、一応再現性の高さからlocus_tagでもヒットするようにしています。
            # locus_tagとold_locus_tagの重複がないことを確認してください。
            matched_tags = set()

            cds_tags = []
            cds_tags.extend(split_tags(attrs.get("old_locus_tag", "")))
            cds_tags.extend(split_tags(attrs.get("locus_tag", "")))

            for tag in cds_tags:
                if tag in target_locus:
                    matched_tags.add(tag)

            if "ID" in attrs:
                cds_id = normalize_gff_id(attrs["ID"])
                matched_tags.update(locus_feature_ids.get(cds_id, set()))

            # CDSの行でParent=...はあるが、old_locus_tag, locus_tagが無いものはParentのIDに従うようにしました.
            #sed -n '9p' genomic.gff | grep 'Parent' | grep -v 'old_locus_tag' | wc -l
            # 0  より今回は不要ですが、再現性の高さのためにこのようにしました。
            parents = split_tags(attrs.get("Parent", ""))

            for parent in parents:
                parent = normalize_gff_id(parent)
                
                if parent in target_locus:
                    matched_tags.add(parent)
                
                # ParentのIDがすでに拾われていた場合に対応させました
                matched_tags.update(locus_feature_ids.get(parent, set()))

            # 重複除去
            matched_tags = sorted(matched_tags)

            if not matched_tags:
                continue

            for tag in matched_tags:
                wp_ids.add(wp_id)
                # 複数のlocus_tagがある場合に対応させるため、setで保存しました。
                locus_to_wp[tag].add(wp_id)

    return wp_ids, locus_to_wp, matched_locus_in_gff, locus_feature_ids, locus_feature_rows



def get_fasta_id(header):

    header = header.strip()

    m = re.search(r"WP_\d+\.\d+", header)
    if m:
        return m.group(0)

    # >WP...
    first = header[1:].split()[0]
    if first.startswith("WP_"):
        return first

    return None


def extract_fasta(faa_file, wp_ids, output_faa):
    extracted = set()

    with open(faa_file) as inp, open(output_faa, "w") as out:
        write_seq = False

        for line in inp:
            if line.startswith(">"):
                current_id = get_fasta_id(line)

                if current_id in wp_ids:
                    write_seq = True
                    extracted.add(current_id)
                    out.write(line)
                else:
                    write_seq = False
            else:
                if write_seq:
                    out.write(line)

    return extracted


def write_mapping(locus_to_wp, mapping_out):
    with open(mapping_out, "w") as out:
        out.write("locus\twp_id\n")

        for locus in sorted(locus_to_wp):
            for wp in sorted(locus_to_wp[locus]):
                out.write(f"{locus}\t{wp}\n")


def write_locus_feature_ids(locus_feature_rows, locus_feature_ids_out):
    
    # feature != CDSのものも出力するようにしました。
    with open(locus_feature_ids_out, "w") as out:
        out.write("feature_type\traw_feature_id\tnormalized_feature_id\tlocus\n")

        for feature_type, raw_feature_id, normalized_feature_id, locus in sorted(locus_feature_rows):
            out.write(
                f"{feature_type}\t{raw_feature_id}\t{normalized_feature_id}\t{locus}\n"
            )


def write_list(items, output_file):
    with open(output_file, "w") as out:
        for x in sorted(items):
            out.write(x + "\n")


def main():
    if len(sys.argv) != 5:
        sys.exit(
            "Usage:\n"
            "  python extract_faa_by_old_locus_tag.py input.csv genome.gff genome.faa output.faa\n"
        )

    csv_file = sys.argv[1]
    gff_file = sys.argv[2]
    faa_file = sys.argv[3]
    output_faa = sys.argv[4]

    mapping_out = output_faa + ".mapping.tsv"
    missing_wp_out = output_faa + ".missing_wp.txt"
    unmatched_locus_out = output_faa + ".unmatched_locus.txt"
    locus_feature_ids_out = output_faa + ".locus_all_feature_ids.tsv"

    target_locus = load_target_locus_tags(csv_file)

    wp_ids, locus_to_wp, matched_locus_in_gff, locus_feature_ids, locus_feature_rows = extract_wp_ids_from_gff(
        gff_file,
        target_locus
    )

    extracted = extract_fasta(faa_file, wp_ids, output_faa)

    write_mapping(locus_to_wp, mapping_out)
    write_locus_feature_ids(locus_feature_rows, locus_feature_ids_out)

    missing_wp = wp_ids - extracted
    unmatched_locus = target_locus - set(locus_to_wp.keys())

    write_list(missing_wp, missing_wp_out)
    write_list(unmatched_locus, unmatched_locus_out)

    print(f"CSV locus数: {len(target_locus)}")
    print(f"GFF内で見つかったlocus数: {len(matched_locus_in_gff)}")
    print(f"old_locus_tag/locus_tagに対応したfeature ID数: {len(locus_feature_ids)}")
    print(f"GFFで対応したWP数: {len(wp_ids)}")
    print(f"FAAから抽出した配列数: {len(extracted)}")
    print(f"出力FAA: {output_faa}")
    print(f"対応表: {mapping_out}")
    print(f"feature ID対応表: {locus_feature_ids_out}")
    print(f"FAA内で見つからなかったWP_ID: {missing_wp_out}")
    print(f"GFFでWPに対応できなかったlocus: {unmatched_locus_out}")


if __name__ == "__main__":
    main()