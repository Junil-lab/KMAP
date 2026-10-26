#!/usr/bin/env python3
"""Run SQANTI3 with its TD2 parser fixed to preserve colon-containing IDs.

The installed SQANTI3 source and SIF remain unchanged. Only the in-process
parser is patched, before SQANTI3 imports it.
"""
import inspect
from pathlib import Path
import runpy
import sys


def patch_parser(qc_script):
    sys.path.insert(0, str(Path(qc_script).resolve().parent))
    from src import parsers
    source = inspect.getsource(parsers.extract_variables)
    broken = r'(?P<id>[^\s:]+)'
    fixed = r'(?P<id>[^\s]+)'
    if broken in source:
        if source.count(broken) != 1:
            raise RuntimeError('Unexpected SQANTI3 ORF parser; refusing ambiguous patch')
        exec(compile(source.replace(broken, fixed), '<SQANTI3 full-ID parser fix>', 'exec'), parsers.__dict__)
    elif fixed not in source:
        raise RuntimeError('Unrecognized SQANTI3 ORF parser; review compatibility')
    return parsers


if __name__ == '__main__':
    qc_script = sys.argv.pop(1)
    patch_parser(qc_script)
    sys.argv[0] = qc_script
    runpy.run_path(qc_script, run_name='__main__')
