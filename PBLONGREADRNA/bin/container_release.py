#!/usr/bin/env python3
"""Package split release assets or install and verify them without replacing images."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
from urllib.parse import quote
from urllib.request import urlopen

IMAGES = tuple('kmap_isoseq_'+x+'.sif' for x in ('main','discovery'))
MANIFEST = 'isoseq-containers.json'
CHUNK = 8 * 1024 * 1024


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for data in iter(lambda: f.read(CHUNK), b''):
            h.update(data)
    return h.hexdigest()


def package(source, output, part_bytes=1900*1024*1024):
    if not 0 < part_bytes < 2*1024**3:
        raise ValueError('Part size must be positive and below 2 GiB')
    for name in IMAGES:
        if not (source/name).is_file() or not (source/name).stat().st_size:
            raise ValueError('Missing/empty image: '+str(source/name))
    output.mkdir(parents=True, exist_ok=False)
    images = []
    for name in IMAGES:
        src = source/name
        total = src.stat().st_size
        parts = []; whole = hashlib.sha256(); index = 0
        with src.open('rb') as f:
            remaining = total
            while remaining:
                asset = name if total <= part_bytes else f'{name}.part{index:03d}'
                size = min(part_bytes, remaining); h = hashlib.sha256()
                with (output/asset).open('xb') as out:
                    left = size
                    while left:
                        data = f.read(min(CHUNK, left))
                        if not data: raise ValueError('Input image changed during packaging')
                        out.write(data); h.update(data); whole.update(data); left -= len(data)
                parts.append(dict(name=asset, bytes=size, sha256=h.hexdigest()))
                remaining -= size; index += 1
            if f.read(1): raise ValueError('Input image changed during packaging')
        images.append(dict(name=name, bytes=total, sha256=whole.hexdigest(), parts=parts))
        print('Packaged', name, 'parts='+str(len(parts)), flush=True)
    manifest = dict(schema_version=1, images=images)
    # Manifest is written last: incomplete packaging cannot look ready to install.
    (output/MANIFEST).write_text(json.dumps(manifest, indent=2)+'\n')
    (output/'SHA256SUMS').write_text(''.join(f"{p['sha256']}  {p['name']}\n" for i in images for p in i['parts'])+
                                  f'{digest(output/MANIFEST)}  {MANIFEST}\n')
    return manifest


def validate_manifest(manifest):
    if manifest.get('schema_version') != 1: raise ValueError('Unsupported manifest version')
    images = manifest.get('images', [])
    if len(images) != len(IMAGES) or {i['name'] for i in images} != set(IMAGES):
        raise ValueError('Manifest must contain exactly the two KMAP main/discovery images')
    assets = set()
    for item in images:
        parts = item.get('parts', [])
        if not parts or item['bytes'] <= 0 or sum(p['bytes'] for p in parts) != item['bytes']:
            raise ValueError('Invalid image/part sizes')
        for record in [item, *parts]:
            if not re.fullmatch(r'[0-9a-f]{64}', record['sha256']): raise ValueError('Invalid SHA256')
        for n, part in enumerate(parts):
            expected = item['name'] if len(parts) == 1 and part['name'] == item['name'] else f"{item['name']}.part{n:03d}"
            if part['name'] != expected or part['name'] in assets: raise ValueError('Invalid asset name/order')
            if not 0 < part['bytes'] < 2*1024**3: raise ValueError('Invalid part size')
            assets.add(part['name'])
    return images


def install(destination, assets_dir=None, repo=None, tag=None):
    if assets_dir is None:
        if not repo or not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+', repo) or not tag:
            raise ValueError('A repository and explicit release tag are required')
        base = f'https://github.com/{repo}/releases/download/{quote(tag, safe="")}/'
    def stream(name):
        return (assets_dir/name).open('rb') if assets_dir is not None else urlopen(base+quote(name), timeout=120)
    with stream(MANIFEST) as f:
        data = f.read(1024*1024+1)
    if len(data)>1024*1024: raise ValueError('Manifest too large')
    images = validate_manifest(json.loads(data))
    destination.mkdir(parents=True, exist_ok=True)
    # Validate all existing files before installing anything. Never overwrite a SIF.
    for item in images:
        target = destination/item['name']
        if target.exists():
            if target.stat().st_size != item['bytes'] or digest(target) != item['sha256']:
                raise ValueError(f'Existing image differs: {target}; move it to a backup before retrying')
    for item in images:
        target = destination/item['name']
        if target.exists():
            print('Verified existing', target.name, flush=True); continue
        with tempfile.TemporaryDirectory(prefix='.release-', dir=destination) as temp:
            assembled = Path(temp)/item['name']; whole = hashlib.sha256(); total = 0
            with assembled.open('wb') as out:
                for part in item['parts']:
                    h = hashlib.sha256(); size = 0
                    with stream(part['name']) as f:
                        for block in iter(lambda: f.read(CHUNK), b''):
                            size += len(block)
                            if size > part['bytes']: raise ValueError('Oversized asset: '+part['name'])
                            h.update(block); whole.update(block); out.write(block)
                    if size != part['bytes'] or h.hexdigest() != part['sha256']:
                        raise ValueError('Asset checksum/size mismatch: '+part['name'])
                    total += size
            if total != item['bytes'] or whole.hexdigest() != item['sha256']:
                raise ValueError('Reassembled SIF checksum mismatch: '+item['name'])
            # Atomic creation without replacing a concurrently created target.
            os.link(assembled, target)
        print('Installed and verified', target.name, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('package')
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True, help='New directory; existing directories are rejected')
    p.add_argument('--part-mib', type=int, default=1900)
    p = sub.add_parser('install')
    p.add_argument('--destination', type=Path, required=True)
    p.add_argument('--assets-dir', type=Path)
    p.add_argument('--repo', default='hi6167579-droid/KMAP_Pipeline')
    p.add_argument('--tag')
    args = parser.parse_args()
    try:
        if args.action == 'package': package(args.source, args.output, args.part_mib*1024*1024)
        else: install(args.destination, args.assets_dir, args.repo, args.tag)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(1, 'ERROR: '+str(exc)+'\n')

if __name__ == '__main__': main()
