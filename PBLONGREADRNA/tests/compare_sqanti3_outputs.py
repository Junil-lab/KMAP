"""Compare structural/coding classifications and protein sequences across images."""
import csv
import json
from pathlib import Path
import sys
baseline,candidate=map(Path,sys.argv[1:3])
fields=('structural_category','coding','ORF_length','CDS_start','CDS_end')
def rows(path):
    with path.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        assert set(fields).issubset(reader.fieldnames),reader.fieldnames
        return {r['isoform']:tuple(r[k] for k in fields) for r in reader}
def fasta(path):
    out={};name=None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name=line[1:].split()[0];assert name not in out;out[name]=''
        elif line.strip():
            assert name is not None;out[name]+=line.strip()
    return out
results=[]
paths=sorted(baseline.glob('**/sqanti3/*_classification.txt'))
assert len(paths)==8,len(paths)
for a in paths:
    rel=a.relative_to(baseline);b=candidate/rel
    left,right=rows(a),rows(b)
    assert left==right, 'Classification differs: '+str(rel)
    prefix=a.name.removesuffix('_classification.txt')
    proteins=prefix+'_corrected.faa'
    pa,pb=fasta(a.parent/proteins),fasta(b.parent/proteins)
    assert pa==pb,'Protein sequences differ: '+str(rel)
    results.append(dict(result=str(rel),isoforms=len(left),proteins=len(pa),status='identical'))
summary={'status':'passed','compared_classification_fields':fields,'cases':results}
(candidate/'sqanti3_image_comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
