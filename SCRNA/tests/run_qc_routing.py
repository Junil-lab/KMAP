#!/usr/bin/env python3
"""Exercise three scRNA entry routes with stubbed tools and real summary merging."""
import argparse,csv,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--nextflow',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[1];out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
for d in ['scripts','configs']:shutil.copytree(root/d,out/d)
s=(root/'main.nf').read_text()
stubs={
'CELLRANGER_COUNT':'mkdir -p ${meta.id}/outs; touch ${meta.id}/outs/raw_feature_bc_matrix.h5 ${meta.id}.cellranger.log',
'CELLBENDER_REMOVE_BACKGROUND':'touch ${meta.id}.cellbender_cleaned.h5 ${meta.id}.cellbender.log',
'SCANPY_QC_FOR_SOLO':'''touch ${meta.id}_qc_full.h5ad ${meta.id}_qc_hvg${params.hvg_n_top}.h5ad
    printf 'sample_id,n_cells_after_qc\\n${meta.id},10\\n' > ${meta.id}.qc_metrics.csv''',
'SOLO_DOUBLETS':'mkdir -p ${meta.id}.solo; touch ${meta.id}.solo/soloed.h5ad',
'FILTER_SOLO_SINGLETS':'''touch ${meta.id}_qc_solo_annotated_full.h5ad ${meta.id}_qc_solo_singlets_full.h5ad ${meta.id}.solo_doublet_barcodes.csv
    printf 'sample_id,n_cells_after_qc,n_solo_doublets,n_solo_singlets\\n${meta.id},10,2,8\\n' > ${meta.id}.solo_summary.csv'''}
for name,body in stubs.items():
 start=s.index('process '+name+' {');pos=s.index('\n}\n',start)
 s=s[:pos]+'\n    stub:\n    """\n    '+body+'\n    """\n'+s[pos:]
(out/'main.nf').write_text(s)
(out/'nextflow.config').write_text((root/'nextflow.config').read_text()+f'''\nparams.outdir = '{out}/results'
params.samplesheet = '{out}/samples.csv'
singularity.enabled = false
process {{ withName: '.*' {{ cpus = 1; memory = '256 MB'; container = null }} }}
executor.queueSize = 4
''')
(out/'input').mkdir();(out/'ref').mkdir();(out/'raw.h5').touch();(out/'clean.h5').touch()
with (out/'samples.csv').open('w') as h:
 w=csv.writer(h);w.writerow(['sample_id','fastqs','transcriptome','raw_h5','cellbender_h5','start_step'])
 for step in ['cellranger','cellbender','qc']:w.writerow(['sample_'+step,str(out/'input'),str(out/'ref'),str(out/'raw.h5'),str(out/'clean.h5'),step])
with (out/'run.log').open('w') as h:subprocess.run([a.nextflow,'run',str(out/'main.nf'),'-stub-run','-work-dir',str(out/'work')],cwd=out,stdout=h,stderr=subprocess.STDOUT,check=True)
with (out/'results/scrnaseq_qc_summary.csv').open() as h:rows=list(csv.DictReader(h))
assert {r['sample_id'] for r in rows}=={'sample_cellranger','sample_cellbender','sample_qc'}
assert all(r['n_solo_singlets']=='8' for r in rows)
print('PASS: all 3 entry routes produce sample-complete summaries; biological tools stubbed.')
