#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse the KEGG BRITE hierarchy (ko00001) into a static pathway -> category map.

Input : br:ko00001 text dump (https://rest.kegg.jp/get/br:ko00001)
Output: TSV with columns: pathway_id, level1, level2, pathway_name

Each HOG-level KEGG_Pathway entry (e.g. "00010,00030,...") can then be mapped
to a level-1 / level-2 functional category for composition analysis.

Usage:
    python3 build_kegg_category_map.py -i ko00001.txt -o kegg_pathway_categories.tsv
"""

import argparse
import re


def parse_arguments():
    parser = argparse.ArgumentParser(description="Build KEGG pathway -> category map from BRITE ko00001")
    parser.add_argument("-i", "--input", required=True, help="BRITE ko00001 text dump")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    return parser.parse_args()


def main():
    args = parse_arguments()
    level1 = level2 = ""
    rows = []
    c_re = re.compile(r"^C\s+(\d{5})\s+(.*)$")
    with open(args.input) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("A "):
                level1 = line[2:].strip()
            elif line.startswith("B "):
                level2 = line[2:].strip()
            elif line.startswith("C "):
                m = c_re.match(line)
                if m:
                    rows.append((m.group(1), level1, level2, m.group(2)))
    with open(args.output, "w") as fh:
        fh.write("pathway_id\tlevel1\tlevel2\tpathway_name\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print(f"{len(rows)} pathways mapped -> {args.output}")


if __name__ == "__main__":
    main()
