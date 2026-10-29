#!/usr/bin/env python3
"""
统计染色体上滑窗内的重复序列比例
输入: RepeatMasker *.fasta.out 文件 + 基因组 fasta 文件
输出: TSV (chr_id, window_start, window_end, repeat_family, proportion)
"""

import argparse
import os
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="统计染色体滑窗内重复序列比例")
    parser.add_argument("-r", "--repeatmasker", required=True, help="RepeatMasker 输出的 *.fasta.out 文件")
    parser.add_argument("-f", "--fasta", required=True, help="基因组 fasta 文件，用于获取染色体长度")
    parser.add_argument("-w", "--window", type=int, default=10000, help="窗口大小 (默认 10000)")
    parser.add_argument("-s", "--step", type=int, default=10000, help="滑动步长 (默认 10000)")
    parser.add_argument("-o", "--output", default=None, help="输出文件 (默认 stdout)")
    return parser.parse_args()


def get_chrom_lengths(fasta_path):
    """从 fasta 文件获取每条染色体的长度"""
    chrom_lengths = {}
    current_chrom = None
    length = 0
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_chrom is not None:
                    chrom_lengths[current_chrom] = length
                current_chrom = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
        if current_chrom is not None:
            chrom_lengths[current_chrom] = length
    return chrom_lengths


def parse_repeatmasker(rm_path):
    """
    解析 RepeatMasker *.fasta.out 文件
    跳过前 3 行标题，从第 4 行开始解析
    返回: list of dict, 每个 dict 包含 chr_id, begin, end, repeat_family, span
    """
    repeats = []
    with open(rm_path, "r") as f:
        # 跳过前 3 行标题
        for _ in range(3):
            next(f)
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 10:
                continue
            try:
                entry = {
                    "score": float(fields[0]),
                    "perc_div": float(fields[1]),
                    "perc_del": float(fields[2]),
                    "perc_ins": float(fields[3]),
                    "chr_id": fields[4],
                    "begin": int(fields[5]),
                    "end": int(fields[6]),
                    "strand": fields[8],
                    "repeat_match": fields[9],
                    "repeat_family": fields[10],
                    "span": int(fields[6]) - int(fields[5]),
                }
                repeats.append(entry)
            except (ValueError, IndexError):
                continue
    return repeats


def compute_sliding_window(repeats, chrom_lengths, window_size, step_size):
    """
    计算每条染色体上滑窗内各 repeat_family 的比例
    返回: list of dict, 每个 dict 包含 chr_id, start, end, repeat_family, proportion
    """
    # 按染色体分组 repeats
    repeats_by_chr = defaultdict(list)
    for r in repeats:
        repeats_by_chr[r["chr_id"]].append(r)

    results = []

    for chrom_id, length in sorted(chrom_lengths.items()):
        chrom_repeats = repeats_by_chr.get(chrom_id, [])
        if not chrom_repeats:
            continue

        # 遍历该染色体的所有窗口
        start = 0
        while start + window_size <= length:
            end = start + window_size
            # 统计窗口内各 repeat_family 覆盖的 bp 数
            family_bp = defaultdict(int)
            for r in chrom_repeats:
                r_start = r["begin"]
                r_end = r["end"]
                # 计算 overlap
                overlap_start = max(start, r_start)
                overlap_end = min(end, r_end)
                if overlap_end > overlap_start:
                    family_bp[r["repeat_family"]] += overlap_end - overlap_start

            # 输出每个有 repeat 的 family
            for family, bp in sorted(family_bp.items()):
                proportion = bp / window_size
                results.append({
                    "chr_id": chrom_id,
                    "start": start,
                    "end": end,
                    "repeat_family": family,
                    "proportion": proportion,
                })

            start += step_size

    return results


def main():
    args = parse_args()

    # 获取染色体长度
    chrom_lengths = get_chrom_lengths(args.fasta)
    print(f"Loaded {len(chrom_lengths)} chromosomes from {args.fasta}", flush=True)

    # 解析 RepeatMasker 文件
    repeats = parse_repeatmasker(args.repeatmasker)
    print(f"Parsed {len(repeats)} repeat annotations from {args.repeatmasker}", flush=True)

    # 计算滑窗统计
    results = compute_sliding_window(repeats, chrom_lengths, args.window, args.step)

    # 输出结果
    header = "chr_id\tstart\tend\trepeat_family\tproportion"
    if args.output:
        with open(args.output, "w") as f:
            f.write(header + "\n")
            for r in results:
                f.write(f"{r['chr_id']}\t{r['start']}\t{r['end']}\t{r['repeat_family']}\t{r['proportion']:.6f}\n")
        print(f"Output saved to {args.output}", flush=True)
    else:
        print(header)
        for r in results:
            print(f"{r['chr_id']}\t{r['start']}\t{r['end']}\t{r['repeat_family']}\t{r['proportion']:.6f}")


if __name__ == "__main__":
    main()
