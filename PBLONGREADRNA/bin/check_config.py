#!/usr/bin/env python3
"""Validate effective Nextflow file parameters, including CLI overrides."""
import argparse
import ast
import csv
import json
from pathlib import Path
import sys

FILE_PARAMS = ('genome_fa', 'gencode_gtf', 'main_sif', 'sqanti3_sif', 'discovery_sif',
               'polya_motif_list', 'polya_peak_bed', 'cage_peak_bed')
OPTIONAL = {'polya_motif_list', 'polya_peak_bed', 'cage_peak_bed'}


def validate(config_path, images_path, arguments):
    params = {}
    for line in Path(config_path).read_text().splitlines():
        key, sep, raw = line.partition(' = ')
        if not sep or not key.startswith('params.'):
            continue
        raw = raw.strip()
        try:
            params[key[7:]] = ast.literal_eval(raw)
        except (ValueError, SyntaxError):
            params[key[7:]] = {'null': None, 'false': False, 'true': True}.get(raw, raw)
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument('-params-file')
    parser.add_argument('--samplesheet', default=argparse.SUPPRESS)
    for name in FILE_PARAMS:
        parser.add_argument('--' + name, default=argparse.SUPPRESS)
    parser.add_argument('--methods', default=argparse.SUPPRESS)
    parser.add_argument('--analysis_mode', default=argparse.SUPPRESS)
    overrides, _ = parser.parse_known_args(arguments)
    if overrides.params_file:
        path = Path(overrides.params_file)
        if path.suffix.lower() == '.json':
            supplied = json.loads(path.read_text())
        else:
            try:
                import yaml
            except ImportError as exc:
                raise ValueError('YAML params files require PyYAML; alternatively use JSON') from exc
            supplied = yaml.safe_load(path.read_text())
        if not isinstance(supplied, dict):
            raise ValueError('Params file must contain a mapping/object')
        params.update(supplied)
    for key, value in vars(overrides).items():
        if key in FILE_PARAMS or key in ('methods', 'analysis_mode', 'samplesheet'):
            params[key] = value
    methods = str(params.get('methods', 'isoquant,isocall,flair')).split(',')
    methods = [x.strip().lower() for x in methods]
    if not methods or any(x not in ('isoquant','isocall','flair') for x in methods):
        raise ValueError('Invalid --methods; choose isoquant,isocall,flair')
    if str(params.get('analysis_mode', 'auto')).lower() not in ('auto','single','multi','both'):
        raise ValueError('Invalid --analysis_mode; choose auto,single,multi,both')
    discovery_required = any(x in ('isocall','flair') for x in methods)
    if str(params.get('analysis_mode', 'auto')).lower() in ('multi','both'):
        if not discovery_required:
            raise ValueError('--analysis_mode multi/both requires isocall or flair')
        if params.get('samplesheet'):
            with Path(params['samplesheet']).open(encoding='utf-8-sig', newline='') as handle:
                if len(list(csv.DictReader(handle))) < 2:
                    raise ValueError('--analysis_mode multi/both requires at least two samples')
    for name in FILE_PARAMS:
        if name == 'discovery_sif' and not discovery_required:
            continue
        value = params.get(name)
        if name in OPTIONAL and (value is False or value is None or str(value).lower() in ('false', 'null', '')):
            params[name] = None
            continue
        if not isinstance(value, str) or not value:
            raise ValueError(f'Missing required parameter: {name}')
        path = Path(value).absolute()
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f'Missing/empty file for {name}: {path}')
        params[name] = str(path)
        print(f'[OK] {name}: {path}')
    fai = Path(params['genome_fa'] + '.fai')
    if not fai.is_file() or fai.stat().st_size == 0:
        raise ValueError(f'Missing/empty FASTA index: {fai}; run samtools faidx')
    names = [Path(params[name]).name for name in OPTIONAL if params.get(name)]
    if len(set(names)) != len(names):
        raise ValueError('Optional annotation filenames must be unique')
    images = [params['main_sif'], params['sqanti3_sif']]
    if discovery_required: images.append(params['discovery_sif'])
    Path(images_path).write_text('\n'.join(images) + '\n')


if __name__ == '__main__':
    try:
        validate(sys.argv[1], sys.argv[2], sys.argv[3:])
    except (OSError, ValueError, IndexError) as exc:
        sys.exit(f'ERROR: {exc}')
