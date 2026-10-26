#!/usr/bin/env python3
"""Benchmark all methods on one full sample, without cached tasks or concurrent analysis jobs."""
import argparse
import csv
from datetime import datetime, timezone
import importlib.util
import hashlib
import json
from pathlib import Path
import subprocess
import time
from zoneinfo import ZoneInfo

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('samples',ROOT/'bin/validate_samples.py')
samples_module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(samples_module)


def seconds(text):
    if not text or text == '-': return 0.0
    units={'ms':0.001,'s':1,'m':60,'h':3600,'d':86400}
    import re
    return sum(float(n)*units[u] for n,u in re.findall(r'([\d.]+)\s*(ms|[smhd])',text))


def summarize(outdir, elapsed, started, ended):
    traces=list(csv.DictReader((outdir/'results/pipeline_info/trace.csv').open()))
    if not traces or any(r['status']!='COMPLETED' or r['exit']!='0' for r in traces):
        raise ValueError('Benchmark contains failed, ignored, incomplete or cached tasks')
    results={method:{'discovery_seconds':0.,'alignment_and_input_seconds':0.,'sqanti3_seconds':0.,'qc_seconds':0.}
             for method in ['isoquant','isocall','flair']}
    for row in traces:
        name=row['name']; value=seconds(row['realtime'])
        if name.startswith('PBMM2 '):
            for method in results: results[method]['alignment_and_input_seconds']+=value
        elif name.startswith('ISOQUANT '): results['isoquant']['discovery_seconds']+=value
        elif name.startswith('ISOCALL_'): results['isocall']['discovery_seconds']+=value
        elif name.startswith('FLAIR_TRANSCRIPTOME '):
            results['flair']['discovery_seconds']+=value
        elif name.startswith(('SQANTI3 ','QC_COLLECT ')):
            for method in results:
                if name.endswith('.'+method+')'):
                    results[method]['sqanti3_seconds' if name.startswith('SQANTI3 ') else 'qc_seconds']+=value
    for values in results.values():
        values['total_active_seconds']=sum(values.values())
    report={'started_at_kst':started.astimezone(ZoneInfo('Asia/Seoul')).isoformat(),
            'ended_at_kst':ended.astimezone(ZoneInfo('Asia/Seoul')).isoformat(),
            'whole_run_wall_seconds':elapsed,'task_count':len(traces),'methods':results,
            'timing_definition':'Sum of actual per-stage elapsed times, excluding scheduler waits; One pbmm2 alignment is shared by all three methods and charged to each independent-route total. FLAIR uses native transcriptome counts without read realignment. Whole-run wall includes all methods, preflight and reporting.',
            'task_trace':traces}
    with (outdir/'results/qc_summary/qc_summary.csv').open() as handle: report['qc']=list(csv.DictReader(handle))
    (outdir/'benchmark_results.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samplesheet',required=True)
    parser.add_argument('--sample',required=True)
    parser.add_argument('--outdir',required=True)
    parser.add_argument('--threads',type=int,default=16)
    parser.add_argument('--config',action='append',default=[],help='Additional Nextflow reference/settings config; benchmark resource limits take priority')
    args=parser.parse_args()
    if args.threads<1: parser.error('--threads must be positive')
    samples=samples_module.validate(args.samplesheet)
    matches=[r for r in samples if r['sample_id']==args.sample]
    if len(matches)!=1: parser.error('Sample not found uniquely in samplesheet')
    outdir=Path(args.outdir).resolve()
    outdir.mkdir(parents=True,exist_ok=False)
    with (outdir/'samples.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['sample_id','flnc_bam','tissue'])
        writer.writeheader();writer.writerows(matches)
    config=outdir/'benchmark.config'
    config.write_text(f"executor {{ cpus = {args.threads}; memory = '128 GB'; queueSize = 1 }}\nprocess.maxForks = 1\nparams.run_cohort = false\n")
    command=['bash',str(ROOT/'run.sh'),'--samplesheet',str(outdir/'samples.csv'),
             '--outdir',str(outdir/'results'),'--threads',str(args.threads),
             '--methods','isoquant,isocall,flair']
    for user_config in args.config: command.extend(['-c',str(Path(user_config).resolve(strict=True))])
    command.extend(['-c',str(config),'-work-dir',str(outdir/'work')])
    started=datetime.now(timezone.utc); tick=time.monotonic()
    source_files=[ROOT/'main.nf',ROOT/'nextflow.config',ROOT/'run.sh']+list((ROOT/'modules').glob('*.nf'))
    manifest={'source_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in source_files},
              'container_checksums':(ROOT/'containers/SHA256SUMS').read_text(),
              'sample':matches[0],'input_bam_bytes':Path(matches[0]['flnc_bam']).stat().st_size,
              'threads_per_task':args.threads,'max_concurrent_tasks':1,'command':command,
              'started_at':started.isoformat(),'status':'running'}
    (outdir/'benchmark_status.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Benchmark started: {outdir}',flush=True)
    with (outdir/'benchmark.log').open('w') as log:
        result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    elapsed=time.monotonic()-tick; ended=datetime.now(timezone.utc)
    manifest.update(status='complete' if result.returncode==0 else 'failed',exit_code=result.returncode,
                    ended_at=ended.isoformat(),wall_seconds=elapsed)
    (outdir/'benchmark_status.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if result.returncode: raise SystemExit(result.returncode)
    report=summarize(outdir,elapsed,started,ended)
    print(json.dumps(report['methods'],indent=2),flush=True)
