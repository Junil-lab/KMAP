#!/usr/bin/env python3
"""Validate CSV and resolve BAM paths relative to the samplesheet directory."""
import argparse
import csv
import json
from pathlib import Path
import re


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
            value = row['flnc_bam'].strip()
            if not value:
                raise ValueError(f'Empty BAM path: {sample}')
            bam = Path(value)
            if not bam.is_absolute():
                bam = sheet.parent / bam
            if not bam.is_file():
                raise ValueError(f'BAM not found: {sample}: {bam}')
            samples.append(dict(sample_id=sample, flnc_bam=str(bam.resolve()),
                                tissue=row['tissue'].strip() or 'NA'))
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
                writer = csv.DictWriter(handle, fieldnames=['sample_id', 'flnc_bam', 'tissue'])
                writer.writeheader()
                writer.writerows(samples)
        if args.json:
            print(json.dumps(samples))
        else:
            print(f'Validated {len(samples)} samples')
    except (OSError, ValueError, csv.Error) as exc:
        parser.exit(1, f'ERROR: {exc}\n')
