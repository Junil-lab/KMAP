"""CLI selection errors must stop before launching Nextflow or containers."""
from pathlib import Path
import subprocess
import unittest
ROOT=Path(__file__).resolve().parents[1]
class CliSelection(unittest.TestCase):
    def test_help_describes_default_and_modes(self):
        r=subprocess.run(['bash',str(ROOT/'run.sh'),'--help'],text=True,capture_output=True)
        self.assertEqual(r.returncode,0)
        for token in ('--methods','default: all three','--analysis_mode','IsoQuant always'):
            self.assertIn(token,r.stdout)
    def test_missing_selection_values(self):
        for option in ('--methods','--analysis_mode'):
            for suffix in ([],['--check-only']):
                with self.subTest(option=option,suffix=suffix):
                    r=subprocess.run(['bash',str(ROOT/'run.sh'),option,*suffix],text=True,capture_output=True)
                    self.assertEqual(r.returncode,2)
                    self.assertIn('missing value for '+option,r.stderr)
if __name__=='__main__':unittest.main()
