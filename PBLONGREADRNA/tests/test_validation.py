"""Regression checks for failures that must stop runs and sample-complete QC."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('validate_samples', ROOT/'bin/validate_samples.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'reads.bam').write_bytes(b'test')
        self.sheet = self.root/'samples.csv'

    def write_sheet(self, rows):
        with self.sheet.open('w', newline='', encoding='utf-8-sig') as handle:
            writer = csv.writer(handle)
            writer.writerow(['sample_id', 'flnc_bam', 'tissue'])
            writer.writerows(rows)

    def test_csv_quotes_bom_and_relative_paths(self):
        self.write_sheet([['SampleA', 'reads.bam', 'kidney, "left"']])
        result = module.validate(self.sheet)
        self.assertEqual(result[0]['tissue'], 'kidney, "left"')
        self.assertEqual(result[0]['flnc_bam'], str(self.root/'reads.bam'))

    def test_start_step_discovery_uses_aligned_bam(self):
        (self.root/'aligned.bam').write_bytes(b'test'); (self.root/'aligned.bam.bai').write_bytes(b'i')
        with self.sheet.open('w', newline='') as handle:
            w = csv.writer(handle)
            w.writerow(['sample_id', 'flnc_bam', 'tissue', 'start_step', 'aligned_bam'])
            w.writerows([['A', 'reads.bam', 'x', '', ''], ['B', '', 'y', 'discovery', 'aligned.bam']])
        result = module.validate(self.sheet)
        self.assertEqual([r['start_step'] for r in result], ['pbmm2', 'discovery'])
        self.assertEqual(result[1]['aligned_bam'], str(self.root/'aligned.bam'))

    def test_discovery_without_index_or_bad_step_stops_run(self):
        (self.root/'aligned.bam').write_bytes(b'test')
        for step, bam, msg in [('discovery', 'aligned.bam', 'index not found'), ('sqanti3', 'aligned.bam', 'start_step must be')]:
            with self.sheet.open('w', newline='') as handle:
                w = csv.writer(handle)
                w.writerow(['sample_id', 'flnc_bam', 'tissue', 'start_step', 'aligned_bam'])
                w.writerow(['A', '', 'x', step, bam])
            with self.assertRaisesRegex(ValueError, msg): module.validate(self.sheet)

    def test_duplicate_id_stops_run(self):
        self.write_sheet([['A','reads.bam','x'], ['A','reads.bam','y']])
        with self.assertRaisesRegex(ValueError, 'Duplicate'): module.validate(self.sheet)

    def test_missing_bam_stops_run(self):
        self.write_sheet([['A','missing.bam','x']])
        with self.assertRaisesRegex(ValueError, 'BAM not found'): module.validate(self.sheet)

    def test_empty_sheet_stops_run(self):
        self.write_sheet([])
        with self.assertRaisesRegex(ValueError, 'no samples'): module.validate(self.sheet)

    def test_invalid_sample_id_stops_run(self):
        self.write_sheet([['../escape','reads.bam','x']])
        with self.assertRaisesRegex(ValueError, 'invalid sample_id'): module.validate(self.sheet)

    def test_qc_counts_failed_reads_and_checks_missing_samples(self):
        flagstat = self.root/'flagstat.txt'
        flagstat.write_text('90 + 10 in total (QC-passed reads + QC-failed reads)\n80 + 5 mapped (88.89% : 50.00%)\n')
        classification = self.root/'classification.txt'
        classification.write_text('isoform\tstructural_category\na\tfull-splice_match\nb\tnovel_in_catalog\n')
        output = self.root/'A.qc.csv'
        cmd = [sys.executable,str(ROOT/'bin/qc_summary.py')]
        subprocess.run(cmd+['collect','--sample','A','--tissue','x','--flagstat',str(flagstat),
                            '--classification',str(classification),'--output',str(output)],check=True)
        with output.open() as handle: row = next(csv.DictReader(handle))
        self.assertEqual((row['total_reads'],row['mapped_reads'],float(row['mapping_rate'])),('100','85',85.0))
        result = subprocess.run(cmd+['merge','--input-dir',str(self.root),'--expected','A,B',
                                     '--output',str(self.root/'summary.csv')],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('differ from input',result.stderr)
        self.assertFalse((self.root/'summary.csv').exists())

    def test_qc_merge_requires_each_method_for_same_sample(self):
        columns=['scope','sample_id','method','tissue']
        for method in ['isoquant','isocall','flair']:
            with (self.root/f'{method}.qc.csv').open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=columns)
                writer.writeheader();writer.writerow(dict(scope='sample',sample_id='A',method=method,tissue='x'))
        cmd=[sys.executable,str(ROOT/'bin/qc_summary.py'),'merge','--input-dir',str(self.root),
             '--expected','sample|A|isoquant,sample|A|isocall,sample|A|flair','--output',str(self.root/'summary.csv')]
        subprocess.run(cmd,check=True)
        (self.root/'flair.qc.csv').unlink()
        result=subprocess.run(cmd,capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('differ from input',result.stderr)

    def test_malformed_classification_stops_qc(self):
        flagstat = self.root/'flagstat.txt'
        flagstat.write_text('10 + 0 in total\n10 + 0 mapped (100.00% : N/A)\n')
        classification = self.root/'classification.txt'
        classification.write_text('isoform\twrong_column\na\tFSM\n')
        result = subprocess.run([sys.executable,str(ROOT/'bin/qc_summary.py'),'collect','--sample','A',
                                 '--tissue','x','--flagstat',str(flagstat),'--classification',str(classification),
                                 '--output',str(self.root/'A.qc.csv')],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('structural_category',result.stderr)


if __name__ == '__main__': unittest.main()
