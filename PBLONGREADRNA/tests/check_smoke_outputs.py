#!/usr/bin/env python3
"""Validate all three real discovery branches, their cohort outputs, and complete QC."""
import csv
from pathlib import Path
import re
import sys

root = Path(sys.argv[1])
with (root/'qc_summary/qc_summary.csv').open() as handle:
    rows = list(csv.DictReader(handle))
expected = {('sample',s,m) for s in ('SmokeA','SmokeB') for m in ('isoquant','isocall','flair')}
expected |= {('cohort','cohort',m) for m in ('isocall','flair')}
actual = {(r['scope'],r['sample_id'],r['method']) for r in rows}
assert len(rows) == 8 and actual == expected, (actual,expected)
for row in rows:
    scope,sample,method = row['scope'],row['sample_id'],row['method']
    assert int(row['total_reads']) == (816 if scope == 'cohort' else 408), row
    assert float(row['mapping_rate']) == 100.0, row
    assert int(row['total_isoforms']) > 0, row
    if scope == 'sample':
        assert row['tissue'] == {'SmokeA':'synthetic, tissue','SmokeB':'synthetic "quoted" tissue'}[sample],row
    directory = root/('cohort' if scope == 'cohort' else 'samples/'+sample)/method
    prefix = f'{scope}.{sample}.{method}'
    sq = directory/'sqanti3'
    for suffix in ('_classification.txt','_junctions.txt','_corrected.gtf','_corrected.fasta','_SQANTI3_report.html'):
        path = sq/(prefix+suffix)
        assert path.is_file() and path.stat().st_size > 0,path
    proteins = sq/(prefix+'_corrected.faa')
    assert proteins.is_file() and proteins.read_text().startswith('>'),proteins
    with (sq/(prefix+'_classification.txt')).open() as handle:
        classes = list(csv.DictReader(handle,delimiter='\t'))
    assert any(r['within_CAGE_peak'] == 'TRUE' for r in classes), 'CAGE annotation not applied'
    assert any(r['within_polyA_site'] == 'TRUE' for r in classes), 'polyA annotation not applied'
    if method == 'isoquant':
        assert int(row['FSM']) == 24 and int(row['NIC']) == 12,row
        iq=directory/'discovery'/sample
        model_ids=set(re.findall(r'transcript_id "([^"]+)"',(iq/f'{sample}.transcript_models.gtf').read_text()))
        with (iq/f'{sample}.discovered_transcript_counts.tsv').open() as handle:
            count_ids={line.split('\t')[0] for line in handle if not line.startswith('#')}
        assert model_ids.issubset(count_ids),model_ids-count_ids
    elif method == 'flair':
        path = directory/'discovery'/(prefix+'.counts.tsv' if scope == 'cohort' else sample+'.isoform.counts.txt')
        assert path.is_file() and path.stat().st_size > 0,path
    elif method == 'isocall':
        path=directory/'discovery'/(prefix+'.count_matrix.txt')
        assert path.is_file() and path.stat().st_size > 0,path
print('PASS: IsoQuant, isocall, FLAIR; 6 sample + 2 cohort QC rows; native counts, ORFs, HTML and optional annotations')
