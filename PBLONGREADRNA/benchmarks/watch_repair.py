#!/usr/bin/env python3
"""Track a parser-fix resume separately from the original cold benchmark."""
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time

root = Path(sys.argv[1]).resolve()
while True:
    log = (root/'repair_resume.log').read_text(errors='replace')
    trace = root/'results/pipeline_info/trace.csv'
    tasks = list(csv.DictReader(trace.open())) if trace.exists() else []
    status = 'running'
    if 'Pipeline completed:' in log:
        status = 'complete'
    elif 'ERROR ~' in log or 'Execution cancelled' in log:
        status = 'failed'
    progress = {'updated_at': datetime.now(timezone.utc).isoformat(), 'status': status,
                'phase': 'SQANTI3 parser-fix resume', 'expected_tasks': 14,
                'completed_tasks': sum(t['status'] in ('COMPLETED','CACHED') for t in tasks),
                'tasks': tasks, 'log': 'repair_resume.log'}
    if status == 'complete':
        sample = json.loads((root/'benchmark_status.json').read_text())['sample']['sample_id']
        validation = {}
        for method in ('isoquant','isocall','flair'):
            directory = root/'results/samples'/sample/method/'sqanti3'
            prefix = f'sample.{sample}.{method}'
            with (directory/(prefix+'_classification.txt')).open() as handle:
                ids = {row['isoform'] for row in csv.DictReader(handle, delimiter='\t')}
            faa = directory/(prefix+'_corrected.faa')
            protein_ids = [line[1:].split()[0] for line in faa.open() if line.startswith('>')]
            assert set(protein_ids).issubset(ids), method
            assert (directory/(prefix+'_SQANTI3_report.html')).stat().st_size > 0
            validation[method] = {'transcripts': len(ids), 'protein_records': len(protein_ids),
                                  'colon_protein_ids': sum('::' in i for i in protein_ids)}
        progress['validation'] = validation
        (root/'repair_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    tmp = root/'pipeline_status.json.tmp'
    tmp.write_text(json.dumps(progress,indent=2)+'\n');tmp.replace(root/'pipeline_status.json')
    if status != 'running':
        manifest_path = root/'benchmark_status.json'
        manifest = json.loads(manifest_path.read_text())
        manifest['status'] = status
        manifest['repair'].update(status=status, ended_at=datetime.now(timezone.utc).isoformat())
        manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
        break
    time.sleep(20)
