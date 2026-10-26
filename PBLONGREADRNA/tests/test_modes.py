"""Reject invalid mode requests before launching biological tools."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('check_config', ROOT/'bin/check_config.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

class ModeValidation(unittest.TestCase):
    def test_invalid_requests(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = root/'config'
            config.write_text("params.methods = 'isocall,flair'\n")
            sheet = root/'samples.csv'
            sheet.write_text('sample_id,flnc_bam,tissue\nA,/test.bam,x\n')
            cases = [(['--analysis_mode','typo'], 'Invalid --analysis_mode'),
                     (['--analysis_mode','multi','--methods','isoquant'], 'requires isocall or flair'),
                     (['--analysis_mode','multi','--samplesheet',str(sheet)], 'at least two'),
                     (['--analysis_mode','both','--samplesheet',str(sheet)], 'at least two')]
            for args, message in cases:
                with self.subTest(args=args), self.assertRaisesRegex(ValueError, message):
                    module.validate(config, root/'images', args)

if __name__ == '__main__':
    unittest.main()
