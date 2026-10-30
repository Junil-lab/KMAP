#!/usr/bin/env python3
"""Repair SOLO CLI's extra softmax while retaining the historical values.

The installed SOLO CLI labels `predict(soft=True)` probabilities as logit_scores,
then applies scipy.softmax to those probabilities. This script is specific to
that output contract and verifies it numerically before changing any score.
"""
import argparse
import csv
import json
import os
from pathlib import Path

import h5py
import numpy as np

ATOL = 2e-6


def expected_compressed(prob):
    return 1 / (1 + np.exp(1 - 2 * prob))


def read_scores(path, prob_key, compressed_key, archive_key):
    with h5py.File(path, 'r') as h:
        obs = h['obs']
        prob = obs[prob_key][:].astype(np.float64)
        compressed = obs[archive_key][:].astype(np.float64) if archive_key in obs else obs[compressed_key][:].astype(np.float64)
        current = obs[compressed_key][:].astype(np.float64)
        labels = obs['is_doublet'][:] if 'is_doublet' in obs else obs['solo_is_doublet'][:]
    if not len(prob) or not (len(prob) == len(compressed) == len(labels)):
        raise ValueError(f'Score/label lengths do not match: {path}')
    if not np.isfinite(prob).all() or (prob < -ATOL).any() or (prob > 1 + ATOL).any():
        raise ValueError(f'Invalid original probabilities in {path}')
    if not np.allclose(compressed, expected_compressed(prob), atol=ATOL, rtol=0):
        raise ValueError(f'Historical score does not match double-softmax formula: {path}')
    if not (np.allclose(current, compressed, atol=ATOL, rtol=0) or np.allclose(current, prob, atol=ATOL, rtol=0)):
        raise ValueError(f'Unknown current score values: {path}')
    # SOLO can use an expected-doublet count and choose a threshold other than 0.5.
    # Leave the original classifications untouched.
    return prob, compressed


def write_obs(path, prob_key, compressed_key, archive_key, probability_key):
    prob, compressed = read_scores(path, prob_key, compressed_key, archive_key)
    with h5py.File(path, 'r+') as h:
        obs = h['obs']
        if archive_key not in obs:
            archive = obs.create_dataset(archive_key, data=compressed)
            archive.attrs['encoding-type'] = 'array'
            archive.attrs['encoding-version'] = '0.2.0'
        if probability_key not in obs:
            ds = obs.create_dataset(probability_key, data=prob)
            ds.attrs['encoding-type'] = 'array'
            ds.attrs['encoding-version'] = '0.2.0'
        else:
            obs[probability_key][:] = prob
        obs[compressed_key][:] = prob
        order = [x.decode() if isinstance(x, bytes) else str(x) for x in obs.attrs['column-order']]
        for name in [archive_key, probability_key]:
            if name not in order:
                order.append(name)
        obs.attrs['column-order'] = np.asarray(order, dtype=h5py.string_dtype('utf-8'))
    return len(prob)


def write_array(path, values):
    temp = path.with_name(path.name + '.repair_tmp')
    if path.suffix == '.npy':
        with temp.open('wb') as f:
            np.save(f, values)
    else:
        np.savetxt(temp, values, fmt='%.18e')
    os.replace(temp, path)


def check_array(path, expected):
    vals = np.load(path) if path.suffix == '.npy' else np.loadtxt(path)
    if vals.shape != expected.shape or not np.allclose(vals, expected, atol=ATOL, rtol=0):
        raise ValueError(f'Inconsistent SOLO score file: {path}')
    return vals


def solo_arrays(directory, repair=False):
    prob = np.load(directory / 'logit_scores.npy')
    compressed = expected_compressed(prob)
    for stem, original in [('softmax_scores', prob), ('no_updates_softmax_scores', prob)]:
        for suffix in ('.npy', '.csv'):
            path = directory / (stem + suffix)
            if not path.exists():
                continue
            archive = directory / (stem + '_double_softmax_original' + suffix)
            old = check_array(archive, compressed) if archive.exists() else check_array(path, compressed)
            if repair:
                if not archive.exists():
                    write_array(archive, old)
                write_array(path, original)
    sim = directory / 'logit_scores_sim.npy'
    if sim.exists():
        p_sim = np.load(sim)
        stem = 'no_updates_softmax_scores_sim'
        path = directory / (stem + '.npy')
        if path.exists():
            archive = directory / (stem + '_double_softmax_original.npy')
            old = check_array(archive, expected_compressed(p_sim)) if archive.exists() else check_array(path, expected_compressed(p_sim))
            if repair:
                if not archive.exists():
                    write_array(archive, old)
                write_array(path, p_sim)


def solo_dir(directory, repair=False):
    h5ad = directory / 'soloed.h5ad'
    prob, _ = read_scores(h5ad, 'logit_scores', 'softmax_scores', 'double_softmax_scores_original')
    array_prob = np.load(directory / 'logit_scores.npy')
    if not np.allclose(prob, array_prob, atol=ATOL, rtol=0):
        raise ValueError(f'H5AD/NPY probability mismatch: {directory}')
    with h5py.File(h5ad, 'r') as h:
        labels = h['obs']['is_doublet'][:].astype(bool)
    if not np.array_equal(labels, np.load(directory / 'is_doublet.npy').astype(bool)):
        raise ValueError(f'H5AD/NPY doublet labels mismatch: {directory}')
    solo_arrays(directory, repair)
    if repair:
        write_obs(h5ad, 'logit_scores', 'softmax_scores', 'double_softmax_scores_original', 'doublet_probability')
        marker = directory / 'score_correction.json'
        marker.write_text(json.dumps({
            'reason': 'SOLO CLI applied an extra softmax to predict(soft=True) probabilities',
            'probability_source': 'logit_scores (legacy misnomer; already a probability)',
            'corrected_score': 'softmax_scores',
            'original_transformed_score': 'double_softmax_scores_original',
            'n_cells': len(prob),
            'plots': 'Original SOLO PDF plots were generated before correction and still show transformed scores.',
        }, indent=2) + '\n')
    return len(prob)


def repair_final_h5ad(path, repair=False):
    prob, _ = read_scores(path, 'solo_logit_scores', 'solo_softmax_scores', 'solo_double_softmax_scores_original')
    if repair:
        write_obs(path, 'solo_logit_scores', 'solo_softmax_scores', 'solo_double_softmax_scores_original', 'solo_doublet_probability')
    return len(prob)


def repair_barcodes(path, repair=False):
    with path.open(newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return 0
    for row in rows:
        prob = float(row['solo_logit_scores'])
        original = float(row.get('solo_double_softmax_scores_original', row['solo_softmax_scores']))
        if abs(original - expected_compressed(prob)) > ATOL:
            raise ValueError(f'Barcode CSV score does not match: {path}')
        if repair:
            row['solo_double_softmax_scores_original'] = original
            row['solo_softmax_scores'] = prob
            row['solo_doublet_probability'] = prob
    if repair:
        tmp = path.with_name(path.name + '.repair_tmp')
        with tmp.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp, path)
    return len(rows)


def results_root(root):
    samples = sorted(d for d in root.iterdir() if d.is_dir() and (d / '04_solo').is_dir())
    if not samples:
        raise ValueError(f'No SOLO samples found in {root}')
    jobs = []
    # Validate every file before changing any of them.
    for sample in samples:
        d = sample / '04_solo' / (sample.name + '.solo')
        full = sample / '05_singlets' / (sample.name + '_qc_solo_annotated_full.h5ad')
        singlet = sample / '05_singlets' / (sample.name + '_qc_solo_singlets_full.h5ad')
        barcodes = sample / '05_singlets' / (sample.name + '_solo_doublet_barcodes.csv')
        if not barcodes.exists():
            matches = list((sample / '05_singlets').glob('*.solo_doublet_barcodes.csv'))
            if len(matches) != 1:
                raise ValueError(f'Cannot identify barcode CSV: {sample}')
            barcodes = matches[0]
        n = solo_dir(d)
        if repair_final_h5ad(full) != n:
            raise ValueError(f'Full cell count differs from SOLO: {sample}')
        ns = repair_final_h5ad(singlet)
        nd = repair_barcodes(barcodes)
        if ns + nd != n:
            raise ValueError(f'Singlet/doublet count differs: {sample}')
        jobs.append((d, full, singlet, barcodes, n))
    for d, full, singlet, barcodes, n in jobs:
        solo_dir(d, repair=True)
        repair_final_h5ad(full, repair=True)
        repair_final_h5ad(singlet, repair=True)
        repair_barcodes(barcodes, repair=True)
        print(f'REPAIRED {d.parent.parent.name}: {n} QC cells', flush=True)
    print(f'REPAIRED {len(jobs)} samples', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--solo-dir', type=Path)
    group.add_argument('--results-root', type=Path)
    args = parser.parse_args()
    if args.solo_dir:
        print(f'REPAIRED SOLO scores: {solo_dir(args.solo_dir, repair=True)} cells')
    else:
        results_root(args.results_root)


if __name__ == '__main__':
    main()
