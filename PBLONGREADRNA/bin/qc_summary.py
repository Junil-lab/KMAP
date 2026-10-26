#!/usr/bin/env python3
"""Collect strict sample QC and require all expected samples at merge time."""
import argparse
import csv
from pathlib import Path
import re


def collect(args):
    total_reads = mapped_reads = 0
    for flagstat in args.flagstat:
        text = Path(flagstat).read_text()
        total = re.search(r'^(\d+) \+ (\d+) in total', text, re.M)
        mapped = re.search(r'^(\d+) \+ (\d+) mapped \(', text, re.M)
        if not total or not mapped:
            raise ValueError(f'Unrecognized samtools flagstat output: {flagstat}')
        total_reads += sum(map(int, total.groups()))
        mapped_reads += sum(map(int, mapped.groups()))
    categories = {}
    with open(args.classification) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or 'structural_category' not in reader.fieldnames:
            raise ValueError('Missing structural_category in SQANTI3 classification')
        for row in reader:
            category = row['structural_category']
            if not category:
                raise ValueError('Empty structural_category in SQANTI3 classification')
            categories[category] = categories.get(category, 0) + 1
    result = dict(scope=args.scope, sample_id=args.sample, method=args.method, tissue=args.tissue,
                  total_reads=total_reads, mapped_reads=mapped_reads,
                  mapping_rate=round(100 * mapped_reads / total_reads, 4) if total_reads else 0,
                  total_isoforms=sum(categories.values()),
                  FSM=categories.get('full-splice_match', 0),
                  ISM=categories.get('incomplete-splice_match', 0),
                  NIC=categories.get('novel_in_catalog', 0),
                  NNC=categories.get('novel_not_in_catalog', 0))
    with open(args.output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result))
        writer.writeheader()
        writer.writerow(result)


def merge(args):
    rows, header = [], None
    for path in sorted(Path(args.input_dir).glob('*.qc.csv')):
        with path.open() as handle:
            reader = csv.DictReader(handle)
            if header is None:
                header = reader.fieldnames
            if reader.fieldnames != header:
                raise ValueError(f'Inconsistent QC columns: {path}')
            sample_rows = list(reader)
            if len(sample_rows) != 1:
                raise ValueError(f'Expected one QC row: {path}')
            rows.extend(sample_rows)
    expected = args.expected.split(',')
    key = (lambda row: f"{row['scope']}|{row['sample_id']}|{row['method']}") if any('|' in item for item in expected) else (lambda row: row['sample_id'])
    actual = [key(row) for row in rows]
    if len(actual) != len(set(actual)) or sorted(actual) != sorted(expected):
        raise ValueError(f'QC samples differ from input: expected={expected}, actual={actual}')
    rows.sort(key=lambda row: expected.index(key(row)))
    with open(args.output, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest='command', required=True)
    collect_parser = subs.add_parser('collect')
    for option in ('sample', 'tissue', 'classification', 'output'):
        collect_parser.add_argument('--' + option, required=True)
    collect_parser.add_argument('--flagstat', nargs='+', required=True)
    collect_parser.add_argument('--scope', default='sample', choices=['sample','cohort'])
    collect_parser.add_argument('--method', default='isoquant', choices=['isoquant','isocall','flair'])
    merge_parser = subs.add_parser('merge')
    for option in ('input-dir', 'expected', 'output'):
        merge_parser.add_argument('--' + option, required=True)
    options = parser.parse_args()
    try:
        (collect if options.command == 'collect' else merge)(options)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')
