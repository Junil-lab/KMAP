#!/usr/bin/env python3
"""Run real Nextflow routing/QC with stubbed biological tools for all entry routes.
Requires Java 17+ and Nextflow. No BAM analysis is performed.
"""
import argparse,csv,os,re,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--nextflow',required=True);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
root=Path(__file__).resolve().parents[1];out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
shutil.copytree(root/'scripts',out/'scripts')
s=(root/'main.nf').read_text()
stubs={
'JASMINE':'touch ${meta.id}.methyl.bam ${meta.id}.methyl.bam.bai ${meta.id}.jasmine.log',
'PBMM2':'touch ${meta.id}.aligned.sorted.bam ${meta.id}.aligned.sorted.bam.bai ${meta.id}.pbmm2.log ${meta.id}.flagstat.txt',
'DEEPVARIANT':'touch ${meta.id}.deepvariant.vcf.gz ${meta.id}.deepvariant.vcf.gz.tbi ${meta.id}.deepvariant.g.vcf.gz ${meta.id}.deepvariant.g.vcf.gz.tbi',
'SAWFISH':'touch ${meta.id}.sv.vcf.gz ${meta.id}.sv.vcf.gz.tbi ${meta.id}.sv.run.stats.json',
'TRGT':'touch ${meta.id}.trgt.vcf.gz ${meta.id}.trgt.vcf.gz.tbi ${meta.id}.trgt.spanning.sorted.bam ${meta.id}.trgt.spanning.sorted.bam.bai',
'HIPHASE':'''touch ${meta.id}.haplotagged.bam ${meta.id}.haplotagged.bam.bai ${meta.id}.snv.phased.vcf.gz ${meta.id}.snv.phased.vcf.gz.tbi ${meta.id}.sv.phased.vcf.gz ${meta.id}.sv.phased.vcf.gz.tbi ${meta.id}.tr.phased.vcf.gz ${meta.id}.tr.phased.vcf.gz.tbi ${meta.id}.blocks.tsv ${meta.id}.stats.csv
    printf 'sample_name\\tchromosome\\tnum_phased\\tnum_variants\\tnum_blocks\\tblock_ng50\\n${meta.id}\\tchr1\\t2\\t3\\t1\\t100\\n' > ${meta.id}.summary.tsv''',
'CPG_TOOLS':'''sleep 1
    touch ${meta.id}.combined.bed.gz ${meta.id}.combined.bed.gz.tbi'''}
for name,body in stubs.items():
 start=s.index('process '+name+' {');pos=s.index('\n}\n',start)
 s=s[:pos]+'\n    stub:\n    """\n    '+body+'\n    """\n\n'+s[pos:]
s=s.replace('samtools flagstat "${bam}" > alignment.flagstat.txt',"printf '10 + 0 in total (QC-passed reads + QC-failed reads)\\n9 + 0 mapped (90.0%% : N/A)\\n' > alignment.flagstat.txt")
s=s.replace('bcftools view -H "${sv_vcf}" | wc -l > sv_count.txt','echo 3 > sv_count.txt')
(out/'main.nf').write_text(s)
config=(root/'nextflow.config').read_text()+f'''\nparams.outdir = '{out}/results'
params.samplesheet = '{out}/samples.csv'
params.vcf_dir = '{out}/vcfs'
singularity.enabled = false
process {{ withName: '.*' {{ cpus = 1; memory = '256 MB'; container = null }} }}
executor.queueSize = 4
'''
(out/'nextflow.config').write_text(config)
rows=[]
for step in ['jasmine','pbmm2','deepvariant','hiphase']:
 sample='sample_'+step;bam=out/(sample+'.bam');bam.touch();Path(str(bam)+'.bai').touch();rows.append([sample,str(bam),'XY',step])
 if step=='hiphase':
  for folder,stem in [('03_small_variants','deepvariant'),('04_structural_var','sv'),('05_tandem_repeats','trgt')]:
   d=out/'vcfs'/sample/folder;d.mkdir(parents=True)
   for ext in ['.vcf.gz','.vcf.gz.tbi']:(d/(sample+'.'+stem+ext)).touch()
with (out/'samples.csv').open('w') as h:
 w=csv.writer(h);w.writerow(['sample_id','bam','sex','start_step']);w.writerows(rows)
with (out/'run.log').open('w') as h:subprocess.run([a.nextflow,'run',str(out/'main.nf'),'-stub-run','-work-dir',str(out/'work')],cwd=out,stdout=h,stderr=subprocess.STDOUT,check=True)
with (out/'results/qc_summary.csv').open() as h:result=list(csv.DictReader(h))
assert {r['sample_id'] for r in result}=={r[0] for r in rows}
assert all(r['status']=='success' and r['total_reads']=='10' for r in result)
with (out/'results/pipeline_info/trace.csv').open() as h:trace=list(csv.DictReader(h))
for sample,*_ in rows:
 q=next(r for r in trace if r['name']=='QC_COLLECT ('+sample+')')
 c=next(r for r in trace if r['name']=='CPG_TOOLS ('+sample+')')
 assert q['start']>=c['complete'],(q,c)
print('PASS: 4 entry routes produce complete QC after CpG; biological tools stubbed.')
