"""Small end-to-end checks for cell QC, gene preservation, and SOLO label transfer."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
if not (SCRIPTS / 'qc_for_solo.py').exists(): SCRIPTS = Path(__file__).parent

def module(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f'{name}.py')
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj

qc, filt, scores = module('qc_for_solo'), module('filter_solo_singlets'), module('correct_solo_scores')

class RawGenesTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.x = np.array([[0,4,3,2,1,0], [0,3,4,2,0,0], [0,2,2,3,0,0],
                           [0,3,2,1,0,0], [0,2,4,2,0,0], [0,2,3,4,0,0],
                           [9,1,1,0,0,0], [0,1,0,0,0,0]], dtype=np.float32)
        self.a = ad.AnnData(sparse.csr_matrix(self.x),obs=pd.DataFrame(index=[f'c{i}' for i in range(8)]),
                           var=pd.DataFrame(index=['MT-A','G1','G2','G3','RARE','ZERO']))
        self.a.write_h5ad(self.root/'input.h5ad')
        args = ['qc','--input',str(self.root/'input.h5ad'),'--sample-id','S',
                '--output-full',str(self.root/'full.h5ad'),'--output-hvg',str(self.root/'hvg.h5ad'),
                '--metrics',str(self.root/'qc.csv'),'--min-genes','2','--min-counts','3',
                '--min-cells-gene','3','--hvg-n-top','3','--mt-threshold','20']
        with patch.object(sys,'argv',args): qc.main()
        self.full = ad.read_h5ad(self.root/'full.h5ad')
        self.hvg = ad.read_h5ad(self.root/'hvg.h5ad')
        probabilities = np.array([0.1,0.8,0.2,0.3,0.9,0.4])
        self.solo = ad.AnnData(obs=pd.DataFrame({
            'is_doublet': probabilities > 0.5,
            'logit_scores': probabilities,
            'softmax_scores': probabilities,
            'double_softmax_scores_original': scores.expected_compressed(probabilities),
        }, index=self.full.obs_names))

    def run_filter(self):
        self.solo.write_h5ad(self.root/'solo.h5ad')
        args=['filter','--full-h5ad',str(self.root/'full.h5ad'),'--solo-h5ad',str(self.root/'solo.h5ad'),
              '--sample-id','S','--annotated-output',str(self.root/'annotated.h5ad'),
              '--output',str(self.root/'singlets.h5ad'),'--summary',str(self.root/'summary.csv'),
              '--doublet-barcodes',str(self.root/'doublets.csv')]
        with patch.object(sys,'argv',args): filt.main()

    def test_mitochondrial_filter_is_opt_in(self):
        import csv
        args = ['qc', '--input', str(self.root/'input.h5ad'), '--sample-id', 'S',
                '--output-full', str(self.root/'default_full.h5ad'),
                '--output-hvg', str(self.root/'default_hvg.h5ad'),
                '--metrics', str(self.root/'default_qc.csv'),
                '--min-genes', '2', '--min-counts', '3',
                '--min-cells-gene', '3', '--hvg-n-top', '3']
        with patch.object(sys, 'argv', args): qc.main()
        default_full = ad.read_h5ad(self.root/'default_full.h5ad')
        self.assertEqual(default_full.n_obs, 7)
        self.assertIn('S_c6', default_full.obs_names)
        self.assertNotIn('S_c6', self.full.obs_names)
        self.assertEqual(default_full.raw.n_vars, 6)
        with (self.root/'default_qc.csv').open() as f:
            metrics = next(csv.DictReader(f))
        self.assertEqual(metrics['mt_filter_enabled'], 'False')
        self.assertEqual(metrics['mt_threshold'], '')
        self.assertEqual(metrics['cells_mt_ge_threshold'], '')

    def test_full_and_hvg_raw_retain_rare_and_zero_genes(self):
        self.assertEqual(self.full.shape,(6,6));self.assertEqual(self.hvg.n_vars,3)
        for a in [self.full,self.hvg]:
            self.assertEqual(a.raw.shape,(6,6))
            np.testing.assert_array_equal(a.raw.X.toarray(),self.x[:6])
            self.assertTrue(a.raw.var_names.equals(self.a.var_names))
        np.testing.assert_array_equal(self.full.X.toarray(),self.x[:6])
        np.testing.assert_array_equal(self.full.layers['counts'].toarray(),self.x[:6])
        self.assertNotIn('RARE',self.hvg.var_names)
        self.assertNotIn('ZERO',self.hvg.var_names)

    def test_label_transfer_retains_all_genes_and_correct_raw_cells(self):
        self.solo=self.solo[::-1].copy()  # verify barcode matching, not row positions
        self.run_filter()
        a=ad.read_h5ad(self.root/'annotated.h5ad');b=ad.read_h5ad(self.root/'singlets.h5ad')
        self.assertEqual(a.shape,(6,6));self.assertEqual(b.shape,(4,6))
        np.testing.assert_array_equal(a.raw.X.toarray(),self.x[:6])
        np.testing.assert_array_equal(b.raw.X.toarray(),self.x[[0,2,3,5]])
        self.assertFalse(b.obs['solo_is_doublet'].any())
        np.testing.assert_allclose(a.obs['solo_softmax_scores'], a.obs['solo_doublet_probability'])
        np.testing.assert_allclose(a.obs['solo_double_softmax_scores_original'], scores.expected_compressed(a.obs['solo_doublet_probability']))

    def test_solo_double_softmax_repair(self):
        d = self.root / 'sample.solo'; d.mkdir()
        probabilities = np.array([0.1, 0.8, 0.2, 0.9], dtype=np.float32)
        compressed = scores.expected_compressed(probabilities).astype(np.float32)
        ad.AnnData(obs=pd.DataFrame({
            'is_doublet': probabilities > 0.5,
            'logit_scores': probabilities,
            'softmax_scores': compressed,
        }, index=[f'c{i}' for i in range(4)])).write_h5ad(d / 'soloed.h5ad')
        np.save(d / 'logit_scores.npy', probabilities)
        np.save(d / 'is_doublet.npy', probabilities > 0.5)
        for stem in ['softmax_scores', 'no_updates_softmax_scores']:
            np.save(d / (stem + '.npy'), compressed)
            np.savetxt(d / (stem + '.csv'), compressed)
        self.assertEqual(scores.solo_dir(d, repair=True), 4)
        self.assertEqual(scores.solo_dir(d, repair=True), 4)  # idempotent
        out = ad.read_h5ad(d / 'soloed.h5ad')
        np.testing.assert_allclose(out.obs['softmax_scores'], probabilities)
        np.testing.assert_allclose(out.obs['double_softmax_scores_original'], compressed)
        np.testing.assert_allclose(out.obs['doublet_probability'], probabilities)
        np.testing.assert_allclose(np.load(d / 'softmax_scores.npy'), probabilities)

    def test_missing_solo_barcode_fails_without_silent_cell_loss(self):
        self.solo=self.solo[:-1].copy()
        with self.assertRaisesRegex(RuntimeError,'missing'):self.run_filter()
        self.assertFalse((self.root/'annotated.h5ad').exists())

    def test_legacy_gene_filtered_input_rejected(self):
        self.full.raw=None;self.full.write_h5ad(self.root/'full.h5ad')
        with self.assertRaisesRegex(RuntimeError,'raw'):self.run_filter()

if __name__=='__main__':unittest.main()
