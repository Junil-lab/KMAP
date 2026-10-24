import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest
spec=importlib.util.spec_from_file_location('merge_qc',Path(__file__).resolve().parents[1]/'scripts/merge_qc.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class MergeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.p=Path(self.tmp.name)
        self.q=self.p/'q.csv';self.s=self.p/'s.csv'
        self.q.write_text('sample_id,n_cells_after_qc\nA,10\n')
        self.s.write_text('sample_id,n_cells_after_qc,n_solo_doublets,n_solo_singlets\nA,10,2,8\n')
    def test_complete(self):self.assertEqual(m.merge([self.q],[self.s],['A'])[0]['n_solo_singlets'],'8')
    def test_missing_samples(self):self.assertRaises(ValueError,m.merge,[self.q],[self.s],['A','B'])
    def test_duplicate_rows(self):self.assertRaises(ValueError,m.merge,[self.q,self.q],[self.s],['A'])
    def test_missing_file(self):self.assertRaises(FileNotFoundError,m.merge,[self.p/'absent'],[self.s],['A'])
    def test_lost_cells(self):
        self.s.write_text('sample_id,n_cells_after_qc,n_solo_doublets,n_solo_singlets\nA,9,2,7\n');self.assertRaises(ValueError,m.merge,[self.q],[self.s],['A'])
    def test_counts_do_not_add_up(self):
        self.s.write_text('sample_id,n_cells_after_qc,n_solo_doublets,n_solo_singlets\nA,10,2,7\n');self.assertRaises(ValueError,m.merge,[self.q],[self.s],['A'])
if __name__=='__main__':unittest.main()
