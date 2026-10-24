#!/usr/bin/env python3
"""Strict WGS QC collection and sample-complete merging (standard library only)."""
import argparse
import csv
from pathlib import Path
import re

FIELDS = ['sample_id', 'status', 'error', 'mapping_rate', 'total_reads', 'mapped_reads',
          'total_sv', 'phased_variants', 'total_variants', 'phase_blocks', 'max_chromosome_block_ng50']

def collect(sample, flagstat, sv_count, summary):
    text = Path(flagstat).read_text()
    def metric(label):
        hits = re.findall(r'^(\d+) \+ (\d+) ' + label, text, flags=re.M)
        if len(hits) != 1:
            raise ValueError('Missing/ambiguous flagstat metric: ' + label)
        return sum(map(int, hits[0]))
    total = metric(r'in total\b')
    mapped = metric(r'mapped \(')
    if total <= 0 or not 0 <= mapped <= total:
        raise ValueError('Invalid alignment counts')
    sv = int(Path(sv_count).read_text().strip())
    if sv < 0:
        raise ValueError('Negative SV count')
    with Path(summary).open() as h:
        reader = csv.DictReader(h, delimiter='\t')
        required = {'sample_name', 'chromosome', 'num_phased', 'num_variants', 'num_blocks', 'block_ng50'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Missing HiPhase summary columns')
        rows = list(reader)
    if not rows or any(r['sample_name'] != sample for r in rows):
        raise ValueError('Empty HiPhase summary or incorrect sample identity')
    if len({r['chromosome'] for r in rows}) != len(rows):
        raise ValueError('Duplicate chromosome in HiPhase summary')
    values = [[int(r[k]) for k in ['num_phased', 'num_variants', 'num_blocks', 'block_ng50']] for r in rows]
    if any(min(v) < 0 or v[0] > v[1] for v in values):
        raise ValueError('Invalid HiPhase summary values')
    return dict(zip(FIELDS, [sample, 'success', '', f'{100*mapped/total:.2f}%', total, mapped, sv,
                           sum(v[0] for v in values), sum(v[1] for v in values),
                           sum(v[2] for v in values), max(v[3] for v in values)]))

def merge(paths, expected):
    if not expected or len(set(expected)) != len(expected):
        raise ValueError('Expected samples must be nonempty and unique')
    rows = []
    for path in paths:
        with Path(path).open() as h:
            reader = csv.DictReader(h)
            if reader.fieldnames != FIELDS:
                raise ValueError('Unexpected QC schema: ' + str(path))
            part = list(reader)
        if len(part) != 1 or part[0]['status'] != 'success' or part[0]['error']:
            raise ValueError('Invalid per-sample QC: ' + str(path))
        rows.extend(part)
    ids = [r['sample_id'] for r in rows]
    if len(ids) != len(set(ids)) or set(ids) != set(expected):
        raise ValueError(f'QC sample mismatch: expected={expected}, observed={ids}')
    return sorted(rows, key=lambda r: r['sample_id'])

def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='mode', required=True)
    c = sub.add_parser('collect')
    for arg in ['sample', 'flagstat', 'sv-count', 'summary', 'output']:
        c.add_argument('--' + arg, required=True)
    m = sub.add_parser('merge')
    m.add_argument('--files', nargs='+', required=True)
    m.add_argument('--expected-samples', nargs='+', required=True)
    m.add_argument('--output', required=True)
    a = p.parse_args()
    rows = [collect(a.sample, a.flagstat, a.sv_count, a.summary)] if a.mode == 'collect' else merge(a.files, a.expected_samples)
    with Path(a.output).open('w', newline='') as h:
        w = csv.DictWriter(h, fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
if __name__ == '__main__':
    main()
