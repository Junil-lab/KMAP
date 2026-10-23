#!/usr/bin/env python3
"""Remove only successful non-Cell-Ranger task directories after publication.

Called by run.sh only after Nextflow exits successfully. Dry-run unless --apply.
Uses only the Python standard library; never traverses staged input symlinks.
"""
import argparse
import csv
import hashlib
from pathlib import Path
import re
import shutil

STEPS = {
    'CELLBENDER_REMOVE_BACKGROUND': ('02_cellbender', ['{s}.cellbender_cleaned.h5', '{s}.cellbender.log'], ['{s}.cellbender_cleaned*.csv', '{s}.cellbender_cleaned*.pdf']),
    'SCANPY_QC_FOR_SOLO': ('03_scanpy_qc', ['{s}_qc_full.h5ad', '{s}_qc_hvg*.h5ad', '{s}.qc_metrics.csv'], []),
    'SOLO_DOUBLETS': ('04_solo', ['{s}.solo'], []),
    'FILTER_SOLO_SINGLETS': ('05_singlets', ['{s}_qc_solo_annotated_full.h5ad', '{s}_qc_solo_singlets_full.h5ad', '{s}.solo_summary.csv', '{s}.solo_doublet_barcodes.csv'], []),
    'QC_MERGE': ('', ['scrnaseq_qc_summary.csv'], []),
}

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.digest()

def verify_copy(src, dst, work):
    if src.is_symlink():
        raise ValueError(f'Output is a symlink; preserving work: {src}')
    if dst.is_symlink() or dst.resolve().is_relative_to(work):
        raise ValueError(f'Published output is not an independent copy: {dst}')
    if src.is_dir():
        if not dst.is_dir():
            raise ValueError(f'Missing published directory: {dst}')
        for item in src.iterdir():
            verify_copy(item, dst / item.name, work)
    elif not dst.is_file() or src.stat().st_size != dst.stat().st_size or digest(src) != digest(dst):
        raise ValueError(f'Missing or different published file: {dst}')

def plan_cleanup(trace, outdir):
    with trace.open(newline='') as f:
        reader = csv.DictReader(f)
        if not {'name', 'status', 'exit', 'workdir'}.issubset(reader.fieldnames or []):
            raise ValueError('Trace must contain name, status, exit and workdir; no files removed.')
        rows = list(reader)
    if not rows or any(r['status'] not in ('COMPLETED', 'CACHED') or r['exit'] != '0' for r in rows):
        raise ValueError('Trace is empty or contains unsuccessful tasks; no files removed.')
    if not any(r['name'] == 'QC_MERGE' for r in rows):
        raise ValueError('No successful QC_MERGE in trace; no files removed.')
    targets = []
    protected = []
    for row in rows:
        name = row['name']
        match = re.fullmatch(r'([A-Z_]+)(?: \(([^/\\]+)\))?', name)
        if not match:
            raise ValueError(f'Unexpected process name: {name}')
        process, sample = match.groups()
        work = Path(row['workdir']).absolute()
        if process == 'CELLRANGER_COUNT':
            protected.append(work.resolve())
            continue
        if process not in STEPS:
            raise ValueError(f'Unexpected process: {process}')
        # Nextflow local work paths have two hex characters and a 30-character hash.
        if not re.fullmatch(r'[0-9a-f]{2}', work.parent.name) or not re.fullmatch(r'[0-9a-f]{30}', work.name):
            raise ValueError(f'Unexpected work directory: {work}')
        if work.is_symlink() or work.resolve() != work or not (work / '.command.sh').is_file():
            raise ValueError(f'Unsafe or missing task directory: {work}')
        if outdir.is_relative_to(work) or work.is_relative_to(outdir):
            raise ValueError(f'Results and work overlap: {work}')
        subdir, required, optional = STEPS[process]
        dest = outdir / sample / subdir if subdir else outdir
        if subdir and (not sample or sample in ('.', '..')):
            raise ValueError(f'Invalid sample: {sample}')
        for pattern in required + optional:
            hits = list(work.glob(pattern.format(s=sample)))
            if pattern in required and not hits:
                raise ValueError(f'Missing task output: {work / pattern.format(s=sample)}')
            for src in hits:
                verify_copy(src, dest / src.name, work)
        targets.append(work)
    for target in targets:
        if any(p.is_relative_to(target) or target.is_relative_to(p) for p in protected):
            raise ValueError(f'Cell Ranger path overlaps cleanup target: {target}')
    return list(dict.fromkeys(targets))

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--trace', required=True, type=Path)
    p.add_argument('--outdir', required=True, type=Path)
    p.add_argument('--apply', action='store_true')
    args = p.parse_args()
    # Validate every output before deleting any task directory.
    targets = plan_cleanup(args.trace, args.outdir.resolve())
    for work in targets:
        if args.apply:
            shutil.rmtree(work)  # unlink staged symlinks; do not follow them
        print(('REMOVED' if args.apply else 'WOULD_REMOVE') + '\t' + str(work), flush=True)
    print(f'Validated {len(targets)} non-Cell-Ranger tasks. Cell Ranger work retained.')

if __name__ == '__main__':
    main()
