"""Check atomic installation and existing-image preservation in the downloader."""
import json,os,shutil,subprocess,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class FetchOfficial(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
  self.root=Path(self.tmp.name);(self.root/'bin').mkdir();(self.root/'containers').mkdir()
  for name in ('bin/fetch_sqanti3.sh','containers/sqanti3_source.json'):shutil.copy2(ROOT/name,self.root/name)
  mock=self.root/'bin/singularity'
  mock.write_text('#!/usr/bin/env bash\nset -eu\nif [[ "$1" == pull ]]; then\n echo "$3" >> "$TEST_PULL_LOG"\n if [[ "${TEST_FAIL:-0}" == 1 ]]; then exit 42; fi\n printf test-image > "$2"\nfi\n')
  mock.chmod(0o755)
  self.env=dict(os.environ,PATH=str(self.root/'bin')+':'+os.environ['PATH'],TEST_PULL_LOG=str(self.root/'pull.log'))
  self.dest=self.root/'containers/sqanti3_v6.0.1.sif'
 def run_fetch(self):
  return subprocess.run(['bash',str(self.root/'bin/fetch_sqanti3.sh')],env=self.env,capture_output=True,text=True)
 def test_pinned_source_and_metadata(self):
  r=self.run_fetch();self.assertEqual(r.returncode,0,r.stderr)
  source=json.loads((self.root/'containers/sqanti3_source.json').read_text())
  self.assertEqual((self.root/'pull.log').read_text().strip(),source['uri'])
  self.assertIn('@sha256:',source['uri']);self.assertEqual(self.dest.read_bytes(),b'test-image')
  self.assertEqual(json.loads(Path(str(self.dest)+'.source.json').read_text())['installed_bytes'],10)
 def test_existing_image_is_not_overwritten(self):
  self.dest.write_bytes(b'existing');self.assertEqual(self.run_fetch().returncode,0)
  self.assertEqual(self.dest.read_bytes(),b'existing');self.assertFalse((self.root/'pull.log').exists())
 def test_failed_download_leaves_no_sif(self):
  self.env['TEST_FAIL']='1';self.assertNotEqual(self.run_fetch().returncode,0)
  self.assertFalse(self.dest.exists());self.assertFalse(list((self.root/'containers').glob('.sqanti3-download.*')))
if __name__=='__main__':unittest.main()
