"""Verify caller selection, final scopes and real biological outputs."""
import csv
from pathlib import Path
import sys
root, mode, methods = Path(sys.argv[1]), sys.argv[2], set(sys.argv[3].split(','))
samples = tuple(sys.argv[4].split(',')) if len(sys.argv) > 4 else ('SmokeA','SmokeB')
do_cohort = mode in ('multi','both') or (mode == 'auto' and len(samples) > 1)
expected = {('sample',sample,method) for sample in samples for method in methods
            if method == 'isoquant' or mode in ('single','both','auto')}
if do_cohort:
    expected |= {('cohort','cohort',method) for method in methods & {'isocall','flair'}}
with (root/'qc_summary/qc_summary.csv').open() as f:
    rows = list(csv.DictReader(f))
assert {(r['scope'],r['sample_id'],r['method']) for r in rows} == expected
assert len(rows) == len(expected)
for row in rows:
    scope,sample,method = row['scope'],row['sample_id'],row['method']
    assert int(row['total_reads']) == (408 * len(samples) if scope == 'cohort' else 408), row
    assert float(row['mapping_rate']) == 100.0 and int(row['total_isoforms']) > 0, row
    folder = root/('cohort' if scope == 'cohort' else 'samples/'+sample)/method
    prefix = f'{scope}.{sample}.{method}'
    for suffix in ('_classification.txt','_corrected.faa','_SQANTI3_report.html'):
        path = folder/'sqanti3'/(prefix+suffix)
        assert path.is_file() and path.stat().st_size, path
        if suffix == '_corrected.faa': assert path.read_text().startswith('>')
    if scope == 'cohort':
        suffix = '.counts.tsv' if method == 'flair' else '.count_matrix.txt'
        header = (folder/'discovery'/(prefix+suffix)).read_text().splitlines()[0]
        assert all(sample in header for sample in samples), header
with (root/'pipeline_info/trace.csv').open() as f:
    trace = list(csv.DictReader(f))
assert all(r['status'] in ('COMPLETED','CACHED') and r['exit'] == '0' for r in trace), trace
names = [r['name'] for r in trace]
for method,process in [('isoquant','ISOQUANT'),('flair','FLAIR_TRANSCRIPTOME'),('isocall','ISOCALL_PROFILE')]:
    assert any(process in n for n in names) == (method in methods), (process, names)
for process,enabled in [
    ('ISOCALL_SAMPLE_CALL', 'isocall' in methods and mode != 'multi'),
    ('ISOCALL_COHORT_CALL', 'isocall' in methods and do_cohort),
    ('FLAIR_COMBINE', 'flair' in methods and do_cohort)]:
    assert any(process in n for n in names) == enabled, (process, names)
print(f'PASS: methods={",".join(sorted(methods))}, mode={mode}, QC rows={len(rows)}; trace, ORFs, reports, counts')
