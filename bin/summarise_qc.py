#!/usr/bin/env python3
"""Combine raw metrics and primary-alignment QC; do not exclude samples."""
import argparse
import csv
import json
from pathlib import Path


def summarise(raw, flagstat, mapq20, min_reads, min_mapped_pct):
    # Include both SAM QC categories. These flags are NOT our sample QC verdict.
    categories = [flagstat['QC-passed reads'], flagstat['QC-failed reads']]
    primary = sum(c['primary'] for c in categories)
    mapped = sum(c['primary mapped'] for c in categories)
    if primary != raw['reads']:
        raise ValueError(f"Primary records ({primary}) != FASTQ reads ({raw['reads']})")
    if not 0 <= mapq20 <= mapped <= primary:
        raise ValueError('Inconsistent mapping counts')
    pct = 100 * mapped / primary
    reasons = []
    if raw['reads'] < min_reads:
        reasons.append(f"read_count<{min_reads}")
    if pct < min_mapped_pct:
        reasons.append(f"primary_mapped_pct<{min_mapped_pct:g}")
    return dict(raw, primary_records=primary, primary_mapped=mapped,
                primary_mapped_pct=pct, primary_mapq20=mapq20,
                primary_mapq20_pct=100 * mapq20 / primary,
                min_reads=min_reads, min_mapped_pct=min_mapped_pct,
                qc_status='FAIL' if reasons else 'PASS',
                qc_reasons=';'.join(reasons) or 'meets_screening_thresholds')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['raw', 'flagstat', 'mapq', 'output']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--min-reads', type=int, required=True)
    parser.add_argument('--min-mapped-pct', type=float, required=True)
    args = parser.parse_args()
    result = summarise(json.loads(Path(args.raw).read_text()),
                       json.loads(Path(args.flagstat).read_text()),
                       int(Path(args.mapq).read_text()), args.min_reads, args.min_mapped_pct)
    with open(args.output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result), delimiter='\t')
        writer.writeheader()
        writer.writerow(result)


if __name__ == '__main__':
    main()
