"""Check that SQANTI3 ORF output preserves colon-containing transcript IDs."""
import csv
from pathlib import Path
import sys
root=Path(sys.argv[1]);prefix=sys.argv[2]
with (root/(prefix+'_classification.txt')).open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==36,len(rows)
assert all(r['isoform'].startswith('id:') for r in rows)
expected={r['isoform'] for r in rows if r['coding']=='coding'}
proteins={line[1:].split()[0] for line in (root/(prefix+'_corrected.faa')).read_text().splitlines() if line.startswith('>')}
assert expected and proteins==expected,(expected,proteins)
print(f'PASS: {len(rows)} colon IDs retained; all {len(expected)} coding IDs match protein FASTA')
