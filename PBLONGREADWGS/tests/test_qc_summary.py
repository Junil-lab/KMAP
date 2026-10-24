import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest
spec=importlib.util.spec_from_file_location('qc',Path(__file__).resolve().parents[1]/'scripts/qc_summary.py')
qc=importlib.util.module_from_spec(spec);spec.loader.exec_module(qc)
class QCTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.p=Path(self.tmp.name)
        (self.p/'flag').write_text('90 + 10 in total (QC-passed reads + QC-failed reads)\n80 + 5 mapped (88.89% : 50.00%)\n70 + 0 primary mapped (90%)\n')
        (self.p/'sv').write_text('3\n')
        (self.p/'phase').write_text('sample_name\tchromosome\tnum_phased\tnum_variants\tnum_blocks\tblock_ng50\nA\tchr1\t4\t8\t2\t100\nA\tchr2\t1\t3\t1\t200\n')
    def collect(self):return qc.collect('A',self.p/'flag',self.p/'sv',self.p/'phase')
    def test_failed_qc_reads_counted_and_metric_named(self):
        r=self.collect();self.assertEqual(r['total_reads'],100);self.assertEqual(r['mapped_reads'],85);self.assertEqual(r['mapping_rate'],'85.00%');self.assertEqual(r['max_chromosome_block_ng50'],200);self.assertNotIn('block_n50',r)
    def test_bad_stats_fail(self):
        (self.p/'flag').write_text('bad');self.assertRaises(ValueError,self.collect)
    def test_wrong_sample_fails(self):
        (self.p/'phase').write_text((self.p/'phase').read_text().replace('A\t','B\t'));self.assertRaises(ValueError,self.collect)
    def test_merge_rejects_missing_duplicate_or_failed_sample(self):
        r=self.collect();f=self.p/'qc.csv'
        def write():
            with f.open('w') as h:
                w=csv.DictWriter(h,fieldnames=qc.FIELDS);w.writeheader();w.writerow(r)
        write();self.assertEqual(len(qc.merge([f],['A'])),1)
        self.assertRaises(ValueError,qc.merge,[f],['A','B'])
        self.assertRaises(ValueError,qc.merge,[f,f],['A'])
        r['error']='failed';write();self.assertRaises(ValueError,qc.merge,[f],['A'])
if __name__=='__main__':unittest.main()
