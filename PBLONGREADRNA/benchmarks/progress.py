#!/usr/bin/env python3
"""Read full-sample benchmark progress; --watch maintains pipeline_status.json."""
import argparse
import csv
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import time
from zoneinfo import ZoneInfo


def inspect(root):
    status=json.loads((root/'benchmark_status.json').read_text())
    trace=root/'results/pipeline_info/trace.csv'
    tasks=list(csv.DictReader(trace.open())) if trace.is_file() else []
    log=root/'benchmark.log'
    lines=log.read_text(errors='replace').splitlines() if log.exists() else []
    events=[line for line in lines if 'process >' in line]
    progress={'updated_at_kst':datetime.now(ZoneInfo('Asia/Seoul')).isoformat(),
              'status':status['status'],'sample_id':status['sample']['sample_id'],
              'input_bam_gib':round(status['input_bam_bytes']/1024**3,3),
              'completed_tasks':sum(t['status']=='COMPLETED' for t in tasks),
              'expected_tasks':14,'latest_event':events[-1] if events else 'preflight',
              'completed_stage_times':{t['name']:t['realtime'] for t in tasks if t['status']=='COMPLETED'}}
    # Byte position is input consumed, not percent of alignment work or an ETA.
    source=Path(status['sample']['flnc_bam']).resolve()
    work=str(root/'work')
    for line in subprocess.check_output(['ps','-eo','pid=,args='],text=True).splitlines():
        parts=line.strip().split(None,1)
        if len(parts)!=2: continue
        pid,args=parts
        first=os.path.basename(args.split()[0])
        if first not in ('pbmm2','samtools') or (' align ' not in args and ' fastq ' not in args): continue
        try:
            if not str(Path(f'/proc/{pid}/cwd').resolve()).startswith(work+'/'): continue
            for fd in Path(f'/proc/{pid}/fd').iterdir():
                if fd.resolve()!=source: continue
                info=Path(f'/proc/{pid}/fdinfo/{fd.name}').read_text()
                pos=int(next(x.split(':')[1] for x in info.splitlines() if x.startswith('pos:')))
                progress['active_input_bytes_consumed_percent']=round(100*pos/status['input_bam_bytes'],1)
                progress['active_input_process']=first
        except (OSError,StopIteration): continue
    return progress


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('benchmark_dir')
    parser.add_argument('--watch',action='store_true')
    args=parser.parse_args();root=Path(args.benchmark_dir).resolve()
    while True:
        result=inspect(root)
        if args.watch:
            temporary=root/'pipeline_status.json.tmp'
            temporary.write_text(json.dumps(result,indent=2)+'\n')
            temporary.replace(root/'pipeline_status.json')
        else: print(json.dumps(result,indent=2))
        if not args.watch or result['status']!='running': break
        time.sleep(20)
