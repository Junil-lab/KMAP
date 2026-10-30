import importlib.util
from pathlib import Path
import sys,tempfile,unittest
from unittest.mock import patch
import anndata as ad
import numpy as np
import pandas as pd
spec=importlib.util.spec_from_file_location('filter_solo',Path(__file__).resolve().parents[1]/'scripts/filter_solo_singlets.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class BarcodeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.p=Path(self.tmp.name)
    def run_case(self,names,labels):
        full=ad.AnnData(np.ones((2,3)),obs=pd.DataFrame(index=['A_c1','A_c2']));full.raw=full.copy()
        prob=np.linspace(0.2,0.8,len(names))  # corrected SOLO probabilities (correct_solo_scores.py output)
        solo=ad.AnnData(np.ones((len(names),1)),obs=pd.DataFrame({'is_doublet':labels,'logit_scores':prob,'softmax_scores':prob,'double_softmax_scores_original':1/(1+np.exp(1-2*prob))},index=names))
        full.write_h5ad(self.p/'full.h5ad');solo.write_h5ad(self.p/'solo.h5ad')
        argv=['filter','--full-h5ad',str(self.p/'full.h5ad'),'--solo-h5ad',str(self.p/'solo.h5ad'),'--sample-id','A','--annotated-output',str(self.p/'ann.h5ad'),'--output',str(self.p/'out.h5ad'),'--summary',str(self.p/'summary.csv'),'--doublet-barcodes',str(self.p/'doublets.csv')]
        with patch.object(sys,'argv',argv):m.main()
    def test_reordered_barcodes_preserve_alignment(self):
        self.run_case(['A_c2','A_c1'],[True,False]);a=ad.read_h5ad(self.p/'out.h5ad');self.assertEqual(list(a.obs_names),['A_c1']);self.assertEqual(a.n_vars,3)
    def test_missing_barcode_fails_without_output(self):
        self.assertRaises(RuntimeError,self.run_case,['A_c1'],[False]);self.assertFalse((self.p/'out.h5ad').exists())
    def test_extra_barcode_fails(self):self.assertRaises(RuntimeError,self.run_case,['A_c1','A_c2','extra'],[False,False,False])
    def test_duplicate_barcode_fails(self):self.assertRaises(RuntimeError,self.run_case,['A_c1','A_c1'],[False,False])
    def test_string_boolean_fails(self):self.assertRaises(RuntimeError,self.run_case,['A_c1','A_c2'],['False','True'])
if __name__=='__main__':unittest.main()
