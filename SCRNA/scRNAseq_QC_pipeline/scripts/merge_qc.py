#!/usr/bin/env python3
"""Merge only complete, unique, mutually consistent sample QC summaries."""
import argparse
import csv
from pathlib import Path

def read_many(paths, expected, required):
    if not expected or len(set(expected)) != len(expected):
        raise ValueError('Expected samples must be nonempty and unique')
    rows = {}
    for path in paths:
        with Path(path).open(newline='') as h:
            reader = csv.DictReader(h)
            if not set(required).issubset(reader.fieldnames or []):
                raise ValueError('Missing required QC columns: ' + str(path))
            part = list(reader)
        if len(part) != 1:
            raise ValueError('Expected one sample per summary: ' + str(path))
        r = part[0]; sample = r['sample_id']
        if sample in rows:
            raise ValueError('Duplicate summary for sample: ' + sample)
        rows[sample] = r
    if set(rows) != set(expected):
        raise ValueError(f'QC sample mismatch: expected={expected}, observed={list(rows)}')
    return rows

def merge(qc_paths, solo_paths, expected):
    qc = read_many(qc_paths, expected, ['sample_id', 'n_cells_after_qc'])
    solo = read_many(solo_paths, expected, ['sample_id', 'n_cells_after_qc', 'n_solo_doublets', 'n_solo_singlets'])
    merged = []
    for sample in sorted(expected):
        q, s = qc[sample], solo[sample]
        n, ns, nd, nkeep = [int(x) for x in [q['n_cells_after_qc'], s['n_cells_after_qc'], s['n_solo_doublets'], s['n_solo_singlets']]]
        if min(n, ns, nd, nkeep) < 0 or n != ns or nd + nkeep != n:
            raise ValueError('Inconsistent QC/SOLO cell counts for ' + sample)
        row = dict(q)
        row.update({k + '_solo' if k in row else k: v for k, v in s.items() if k != 'sample_id'})
        merged.append(row)
    return merged

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--qc-files', nargs='+', required=True)
    p.add_argument('--solo-files', nargs='+', required=True)
    p.add_argument('--expected-samples', nargs='+', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    rows = merge(a.qc_files, a.solo_files, a.expected_samples)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with Path(a.output).open('w', newline='') as h:
        w = csv.DictWriter(h, fieldnames=fields); w.writeheader(); w.writerows(rows)
if __name__ == '__main__':
    main()
