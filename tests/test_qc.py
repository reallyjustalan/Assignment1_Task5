import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from fastq_qc import analyse
from summarise_qc import summarise


class QCTests(unittest.TestCase):
    def test_lengths_quality_and_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'reads.fastq.gz'
            with gzip.open(path, 'wt') as out:
                out.write('@a\nAC\n+\nII\n@b\nACGT\n+\n5555\n')
            result = analyse(path, 'test')
        self.assertEqual(result['reads'], 2)
        self.assertEqual(result['bases'], 6)
        self.assertEqual(result['median_length'], 3)
        self.assertEqual(result['n50'], 4)
        self.assertAlmostEqual(result['mean_base_phred'], 160 / 6)

    def test_rejects_truncated_record_and_empty_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'reads.fastq'
            for text in ['', '@a\nAC\n+\nI\n', '@a\nAC\n+\n']:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    analyse(path, 'test')

    def test_primary_not_all_alignment_records(self):
        flags = {'QC-passed reads': {'primary': 100, 'primary mapped': 70, 'mapped': 140},
                 'QC-failed reads': {'primary': 0, 'primary mapped': 0}}
        result = summarise({'sample': 'x', 'reads': 100}, flags, 60, 100, 70)
        self.assertEqual(result['primary_mapped_pct'], 70)
        self.assertEqual(result['qc_status'], 'PASS')
        result = summarise({'sample': 'x', 'reads': 100}, flags, 60, 101, 71)
        self.assertEqual(result['qc_status'], 'FAIL')
        self.assertEqual(len(result['qc_reasons'].split(';')), 2)
        with self.assertRaises(ValueError):
            summarise({'sample': 'x', 'reads': 99}, flags, 60, 100, 70)

    def test_merge_single_header_sorted(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp) / f'{s}.tsv' for s in ['z', 'a']]
            for path in paths:
                path.write_text(f'sample\tqc_status\n{path.stem}\tPASS\n')
            script = Path(__file__).resolve().parents[1] / 'bin/merge_qc.py'
            output = subprocess.check_output([sys.executable, str(script), *map(str, paths)], text=True)
            self.assertEqual(output.splitlines(), ['sample\tqc_status', 'a\tPASS', 'z\tPASS'])


if __name__ == '__main__':
    unittest.main()
