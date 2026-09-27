#!/usr/bin/env bash
# Real end-to-end tests using Bambu's bundled chr9 fixture, NOT assignment data.
set -euo pipefail
project="$(cd "$(dirname "$0")/.." && pwd)"
source "$project/env.sh"
run="$project/tests/runs/$(date +%Y%m%d_%H%M%S)"
mkdir -p "$run"
cd "$run"
extdata=$(Rscript -e 'cat(system.file("extdata", package="bambu"))')
reference="$extdata/Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa"
gtf="$extdata/Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf"
samtools fastq -F 2304 "$extdata/SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam" \
    | gzip > fixture.fastq.gz
# Same fixture twice ONLY to exercise both alignment branches and BAM collection.
printf 'sample,protocol,fastq\nfixture_directRNA,directRNA,%s\nfixture_cDNA,cDNA,%s\n' \
    "$run/fixture.fastq.gz" "$run/fixture.fastq.gz" > samples.csv
common=(run "$project/main.nf" --samples "$run/samples.csv" --genome "$reference"
        --gtf "$gtf" --outdir "$run/results" -work-dir "$run/work" -ansi-log false)
nextflow "${common[@]}" --mode annotated -with-report annotated.report.html \
    -with-trace annotated.trace.tsv 2>&1 | tee annotated.log
annotated_run=$(nextflow log -q | tail -n 1)
nextflow "${common[@]}" --mode unannotated -resume "$annotated_run" -with-report unannotated.report.html \
    -with-trace unannotated.trace.tsv 2>&1 | tee unannotated.log
python3 - <<'PY'
import csv
from pathlib import Path
for mode in ['annotated', 'unannotated']:
    with open(f'{mode}.trace.tsv') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert len(rows) == 9, rows
    for row in rows:
        expected = 'CACHED' if mode == 'unannotated' and not row['name'].startswith('BAMBU') else 'COMPLETED'
        assert row['status'] == expected, row
    root = Path('results/bambu') / mode / 'bambu_output'
    assert (root / 'se.rds').is_file()
    assert list(root.glob('*.gtf'))
    assert (root / 'gene_counts.tsv').is_file()
print('PASS: both real smoke runs completed; upstream tasks cached; Bambu reran.')
PY
printf 'Smoke evidence: %s\n' "$run"
