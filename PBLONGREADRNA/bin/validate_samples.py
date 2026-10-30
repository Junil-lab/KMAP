#!/usr/bin/env python3
"""Validate CSV and resolve BAM paths relative to the samplesheet directory.

Required columns: sample_id, flnc_bam, tissue.
Optional columns: start_step (pbmm2 = default, or discovery) and aligned_bam.
start_step=discovery skips pbmm2 and starts transcript discovery from aligned_bam
(coordinate-sorted, with a .bai index next to it); flnc_bam may then be empty.
"""
import argparse
import csv
import json
from pathlib import Path
import re


START_STEPS = ('pbmm2', 'discovery')
FIELDS = ['sample_id', 'flnc_bam', 'tissue', 'start_step', 'aligned_bam']


def resolve_bam(sheet, sample, value, column='flnc_bam'):
    value = value.strip()
    if not value:
        raise ValueError(f'Empty BAM path: {sample} ({column})')
    bam = Path(value)
    if not bam.is_absolute():
        bam = sheet.parent / bam
    if not bam.is_file():
        raise ValueError(f'BAM not found: {sample}: {bam}')
    return str(bam.resolve())


def validate(sheet):
    sheet = Path(sheet).resolve(strict=True)
    with sheet.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, skipinitialspace=True)
        required = {'sample_id', 'flnc_bam', 'tissue'}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError('Samplesheet requires sample_id,flnc_bam,tissue columns')
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError('Duplicate CSV column names')
        samples, seen = [], set()
        for line, row in enumerate(reader, 2):
            if None in row or any(row.get(key) is None for key in required):
                raise ValueError(f'CSV row {line}: inconsistent number of columns')
            sample = row['sample_id'].strip()
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample):
                raise ValueError(f'CSV row {line}: invalid sample_id {sample!r}')
            if sample in seen:
                raise ValueError(f'Duplicate sample_id: {sample}')
            seen.add(sample)
            step = (row.get('start_step') or '').strip().lower() or 'pbmm2'
            if step not in START_STEPS:
                raise ValueError(f'CSV row {line}: start_step must be one of {", ".join(START_STEPS)}, not {step!r}')
            flnc, aligned = '', ''
            if step == 'pbmm2':
                flnc = resolve_bam(sheet, sample, row['flnc_bam'])
            else:
                aligned = resolve_bam(sheet, sample, row.get('aligned_bam') or '', 'aligned_bam')
                if not (Path(aligned + '.bai').is_file() or Path(aligned[:-4] + '.bai').is_file()):
                    raise ValueError(f'BAM index not found: {sample}: {aligned}.bai (run samtools index)')
            samples.append(dict(sample_id=sample, flnc_bam=flnc,
                                tissue=row['tissue'].strip() or 'NA',
                                start_step=step, aligned_bam=aligned))
    if not samples:
        raise ValueError('Samplesheet has no samples')
    return samples


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('samplesheet')
    parser.add_argument('--output')
    parser.add_argument('--json', action='store_true')
    args = parser.parse_args()
    try:
        samples = validate(args.samplesheet)
        if args.output:
            with open(args.output, 'w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS)
                writer.writeheader()
                writer.writerows(samples)
        if args.json:
            print(json.dumps(samples))
        else:
            print(f'Validated {len(samples)} samples')
    except (OSError, ValueError, csv.Error) as exc:
        parser.exit(1, f'ERROR: {exc}\n')
