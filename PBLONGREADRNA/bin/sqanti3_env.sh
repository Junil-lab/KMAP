#!/usr/bin/env bash
# Sourced inside the SQANTI3 image; official Docker ENTRYPOINT is skipped by exec.
if [[ -x /conda/miniconda3/envs/sqanti3/bin/python ]]; then
    export PATH="/conda/miniconda3/envs/sqanti3/bin:$PATH"
fi
if [[ "${SQANTI3_QC_SCRIPT:-auto}" == auto ]]; then
    SQANTI3_QC_SCRIPT=''
    for candidate in /opt/SQANTI3/sqanti3_qc.py /opt2/sqanti3/6.0.1/SQANTI3-6.0.1/sqanti3_qc.py; do
        if [[ -f "$candidate" ]]; then SQANTI3_QC_SCRIPT="$candidate"; break; fi
    done
fi
if [[ -z "$SQANTI3_QC_SCRIPT" || ! -f "$SQANTI3_QC_SCRIPT" ]]; then
    echo 'ERROR: SQANTI3 QC script not found; configure --sqanti3_qc with its container path' >&2
    return 1
fi
export SQANTI3_QC_SCRIPT
