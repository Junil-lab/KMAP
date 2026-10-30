import csv
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest

MODULE = Path(__file__).resolve().parent.parent / 'scripts/cleanup_work.py'
if not MODULE.exists():
    MODULE = Path(__file__).with_name('cleanup_work.py')
spec = importlib.util.spec_from_file_location('cleanup_work', MODULE)
cleanup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cleanup)

class CleanupTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / 'results'
        self.out.mkdir()
        self.rows = []
        self.work = {}
        self.external = self.root / 'input.h5'
        self.external.write_text('original data')
        for i, process in enumerate(['CELLRANGER_COUNT', *cleanup.STEPS]):
            work = self.root / 'work' / f'{i:02x}' / f'{i:030x}'
            work.mkdir(parents=True)
            (work / '.command.sh').write_text('#!/bin/bash\n')
            (work / 'staged_input').symlink_to(self.external)
            self.work[process] = work
            self.rows.append(dict(name=process if process == 'QC_MERGE' else process + ' (sample)', status='COMPLETED', exit='0', workdir=str(work)))
            if process == 'CELLRANGER_COUNT':
                (work / 'outs').mkdir()
                (work / 'outs/raw.h5').write_text('preserve cellranger')
                continue
            subdir, required, optional = cleanup.STEPS[process]
            dest = self.out / 'sample' / subdir if subdir else self.out
            dest.mkdir(parents=True, exist_ok=True)
            for pattern in required + optional:
                name = pattern.format(s='sample').replace('*', '2000')
                source = work / name
                if name.endswith('.solo'):
                    source.mkdir()
                    (source / 'soloed.h5ad').write_text('solo result')
                    shutil.copytree(source, dest / name)
                else:
                    source.write_text('result contents')
                    shutil.copy2(source, dest / name)
            (work / 'ckpt.tar.gz').write_text('temporary checkpoint')
        self.trace = self.root / 'trace.csv'
        self.write_trace()

    def write_trace(self):
        with self.trace.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['name', 'status', 'exit', 'workdir'])
            writer.writeheader(); writer.writerows(self.rows)

    def test_cleanup_preserves_cellranger_results_and_external_inputs(self):
        targets = cleanup.plan_cleanup(self.trace, self.out)
        self.assertEqual(len(targets), 5)
        for target in targets:
            shutil.rmtree(target)
        self.assertTrue((self.work['CELLRANGER_COUNT'] / 'outs/raw.h5').is_file())
        self.assertEqual(self.external.read_text(), 'original data')
        self.assertTrue((self.out / 'scrnaseq_qc_summary.csv').is_file())
        self.assertTrue((self.out / 'sample/04_solo/sample.solo/soloed.h5ad').is_file())

    def test_bad_copy_prevents_all_cleanup(self):
        (self.out / 'scrnaseq_qc_summary.csv').write_text('corrupt')
        with self.assertRaises(ValueError):
            cleanup.plan_cleanup(self.trace, self.out)
        self.assertTrue(all(w.exists() for w in self.work.values()))

    def test_output_symlink_is_rejected(self):
        output = self.out / 'scrnaseq_qc_summary.csv'
        output.unlink()
        output.symlink_to(self.work['QC_MERGE'] / output.name)
        with self.assertRaises(ValueError):
            cleanup.plan_cleanup(self.trace, self.out)

    def test_failed_run_is_rejected(self):
        self.rows[-1]['status'] = 'FAILED'
        self.write_trace()
        with self.assertRaises(ValueError):
            cleanup.plan_cleanup(self.trace, self.out)

    def test_cached_tasks_are_checked(self):
        for row in self.rows: row['status'] = 'CACHED'
        self.write_trace()
        self.assertEqual(len(cleanup.plan_cleanup(self.trace, self.out)), 5)

if __name__ == '__main__':
    unittest.main()
