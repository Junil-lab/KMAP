"""Check actual workflow scope, reports, counts and process selection."""
import csv
from pathlib import Path
import sys
root, mode = Path(sys.argv[1]), sys.argv[2]
with (root/'qc_summary/qc_summary.csv').open() as handle:
    rows = list(csv.DictReader(handle))
expected = set()
if mode in ('single','both'):
    expected |= {('sample',s,m) for s in ('SmokeA','SmokeB') for m in ('isocall','flair')}
if mode in ('multi','both'):
    expected |= {('cohort','cohort',m) for m in ('isocall','flair')}
actual = {(r['scope'],r['sample_id'],r['method']) for r in rows}
assert actual == expected and len(rows) == len(expected), (actual, expected)
for scope,sample,method in expected:
    directory = root/('cohort' if scope == 'cohort' else 'samples/'+sample)/method
    prefix = f'{scope}.{sample}.{method}'
    for suffix in ('_classification.txt','_corrected.faa','_SQANTI3_report.html'):
        path = directory/'sqanti3'/(prefix+suffix)
        assert path.is_file() and path.stat().st_size, path
    if scope == 'cohort':
        suffix = '.counts.tsv' if method == 'flair' else '.count_matrix.txt'
        path = directory/'discovery'/(prefix+suffix)
        header = path.read_text().splitlines()[0]
        assert 'SmokeA' in header and 'SmokeB' in header, header
with (root/'pipeline_info/trace.csv').open() as handle:
    names = [r['name'] for r in csv.DictReader(handle)]
assert any('ISOCALL_SAMPLE_CALL' in n for n in names) == (mode != 'multi')
assert any('ISOCALL_COHORT_CALL' in n for n in names) == (mode != 'single')
assert any('FLAIR_COMBINE' in n for n in names) == (mode != 'single')
print(f'PASS: {mode}, {len(rows)} QC rows, ORF/report outputs, process routing and cohort sample counts')
