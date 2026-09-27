#!/usr/bin/env python3
"""Streaming QC for standard four-line, Phred+33 Nanopore FASTQ files."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path


def analyse(path, sample):
    lengths = Counter()
    quality_hist = Counter()
    reads = bases = 0
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rb') as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip(b'\r\n')
            plus = handle.readline()
            qual_line = handle.readline()
            qual = qual_line.rstrip(b'\r\n')
            if (not header.startswith(b'@') or not plus.startswith(b'+')
                    or not seq or not qual_line or len(seq) != len(qual)):
                raise ValueError(f'Malformed FASTQ record {reads + 1}')
            if min(qual) < 33 or max(qual) > 126:
                raise ValueError(f'Invalid Phred+33 quality at record {reads + 1}')
            reads += 1
            bases += len(seq)
            lengths[len(seq)] += 1
            quality_hist.update(qual)
    if not reads:
        raise ValueError('Empty FASTQ')
    midpoints = [(reads + 1) // 2, (reads + 2) // 2]
    cumulative = 0
    medians = []
    for length, count in sorted(lengths.items()):
        previous = cumulative
        cumulative += count
        medians.extend(length for rank in midpoints if previous < rank <= cumulative)
    covered = 0
    for length, count in sorted(lengths.items(), reverse=True):
        covered += length * count
        if covered >= bases / 2:
            n50 = length
            break
    return dict(sample=sample, reads=reads, bases=bases,
                mean_length=bases / reads, median_length=sum(medians) / 2,
                n50=n50, min_length=min(lengths), max_length=max(lengths),
                mean_base_phred=sum((q - 33) * n for q, n in quality_hist.items()) / bases)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastq', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    Path(args.output).write_text(json.dumps(analyse(args.fastq, args.sample), indent=2) + '\n')


if __name__ == '__main__':
    main()
