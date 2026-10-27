"""Release packaging, assembly, corruption and existing-image protection."""
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('release',ROOT/'bin/container_release.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Releases(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.src=self.root/'source';self.src.mkdir()
        for i,name in enumerate(m.IMAGES): (self.src/name).write_bytes(bytes([i+1])*(7+i*21))
        self.assets=self.root/'assets';self.dest=self.root/'installed'
        m.package(self.src,self.assets,16)
    def test_roundtrip_and_repeat(self):
        m.install(self.dest,self.assets);m.install(self.dest,self.assets)
        for name in m.IMAGES:self.assertEqual((self.src/name).read_bytes(),(self.dest/name).read_bytes())
    def test_corrupt_part_never_installed(self):
        (self.assets/m.IMAGES[0]).write_bytes(b'bad')
        with self.assertRaisesRegex(ValueError,'checksum/size'):m.install(self.dest,self.assets)
        self.assertFalse((self.dest/m.IMAGES[0]).exists())
        self.assertFalse(list(self.dest.glob('.release-*')))
    def test_existing_image_preserved(self):
        self.dest.mkdir();p=self.dest/m.IMAGES[0];p.write_bytes(b'old')
        with self.assertRaisesRegex(ValueError,'Existing image differs'):m.install(self.dest,self.assets)
        self.assertEqual(p.read_bytes(),b'old')
    def test_whole_hash_checked(self):
        p=self.assets/m.MANIFEST;data=json.loads(p.read_text());data['images'][0]['sha256']='0'*64;p.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,'Reassembled'):m.install(self.dest,self.assets)
        self.assertFalse((self.dest/m.IMAGES[0]).exists())
    def test_path_traversal_rejected(self):
        p=self.assets/m.MANIFEST;data=json.loads(p.read_text());data['images'][0]['parts'][0]['name']='../escape';p.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,'asset name'):m.install(self.dest,self.assets)
    def test_missing_part_never_installed(self):
        (self.assets/(m.IMAGES[1]+'.part001')).unlink()
        with self.assertRaises(OSError):m.install(self.dest,self.assets)
        self.assertFalse((self.dest/m.IMAGES[1]).exists())
    def test_package_rejects_existing_output(self):
        with self.assertRaises(FileExistsError):m.package(self.src,self.assets,16)
    def test_setup_requires_explicit_tag(self):
        p=subprocess.run(['bash',str(ROOT/'setup.sh'),'--download-containers'],capture_output=True,text=True)
        self.assertEqual(p.returncode,2);self.assertIn('requires --release-tag',p.stderr)
if __name__=='__main__':unittest.main()
