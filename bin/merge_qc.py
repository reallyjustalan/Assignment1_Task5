#!/usr/bin/env python3
"""Merge per-sample TSV files into a stable, single-header QC table."""
import csv
import sys

rows = []
fields = None
for name in sys.argv[1:]:
    with open(name, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if fields is not None and fields != reader.fieldnames:
            raise ValueError('Inconsistent QC columns')
        fields = reader.fieldnames
        rows.extend(reader)
if not rows:
    raise ValueError('No QC rows')
writer = csv.DictWriter(sys.stdout, fieldnames=fields, delimiter='\t')
writer.writeheader()
writer.writerows(sorted(rows, key=lambda row: row['sample']))
