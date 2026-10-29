#!/usr/bin/env python3
import argparse
import csv

import anndata as ad
import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Transfer SOLO scores and write annotated full and singlet-only H5AD files."
    )
    parser.add_argument("--full-h5ad", required=True)
    parser.add_argument("--solo-h5ad", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--annotated-output", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--doublet-barcodes", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    full = ad.read_h5ad(args.full_h5ad)
    solo = ad.read_h5ad(args.solo_h5ad)

    required = ["is_doublet", "softmax_scores", "logit_scores", "double_softmax_scores_original"]
    missing = [col for col in required if col not in solo.obs.columns]
    if missing:
        raise RuntimeError(f"SOLO H5AD is missing required obs columns: {missing}")

    if not full.obs_names.is_unique or not solo.obs_names.is_unique:
        raise RuntimeError("Duplicate cell barcodes in QC or SOLO output.")
    if full.raw is None or not full.var_names.equals(full.raw.var_names):
        raise RuntimeError("Full QC input must preserve all input genes in both X and raw; rerun Scanpy QC.")
    missing = full.obs_names.difference(solo.obs_names)
    extra = solo.obs_names.difference(full.obs_names)
    if len(missing) or len(extra):
        raise RuntimeError(
            f"QC/SOLO barcode mismatch: {len(missing)} missing, {len(extra)} unexpected. "
            "No cells were silently removed."
        )
    if full.n_obs == 0:
        raise RuntimeError("QC matrix contains no cells.")
    full = full.copy()
    solo_obs = solo.obs.loc[full.obs_names]
    labels = solo_obs["is_doublet"]
    if labels.isna().any() or not labels.isin([True, False, 0, 1]).all():
        raise RuntimeError("SOLO is_doublet must contain nonmissing boolean or 0/1 values.")

    probability = pd.to_numeric(solo_obs["logit_scores"], errors="raise").to_numpy()
    corrected = pd.to_numeric(solo_obs["softmax_scores"], errors="raise").to_numpy()
    historical = pd.to_numeric(solo_obs["double_softmax_scores_original"], errors="raise").to_numpy()
    if not np.isfinite(probability).all() or not np.allclose(corrected, probability, atol=2e-6, rtol=0):
        raise RuntimeError("SOLO probabilities were not corrected before singlet filtering.")
    if not np.allclose(historical, 1 / (1 + np.exp(1 - 2 * probability)), atol=2e-6, rtol=0):
        raise RuntimeError("Archived SOLO scores do not match the historical double softmax.")
    full.obs["solo_is_doublet"] = labels.astype(bool).to_numpy()
    for col in ("softmax_scores", "logit_scores"):
        if col in solo_obs.columns:
            full.obs[f"solo_{col}"] = pd.to_numeric(solo_obs[col], errors="coerce").to_numpy()

    full.obs["solo_doublet_probability"] = full.obs["solo_softmax_scores"].to_numpy()
    full.obs["solo_double_softmax_scores_original"] = pd.to_numeric(
        solo_obs["double_softmax_scores_original"], errors="raise"
    ).to_numpy()

    full.write_h5ad(args.annotated_output, compression="gzip")

    singlets = full[~full.obs["solo_is_doublet"].to_numpy()].copy()
    singlets.write_h5ad(args.output, compression="gzip")

    doublet_obs = full.obs.loc[full.obs["solo_is_doublet"], :].copy()
    barcode_table = pd.DataFrame({"barcode": doublet_obs.index})
    if "solo_softmax_scores" in doublet_obs.columns:
        barcode_table["solo_softmax_scores"] = doublet_obs["solo_softmax_scores"].to_numpy()
    if "solo_logit_scores" in doublet_obs.columns:
        barcode_table["solo_logit_scores"] = doublet_obs["solo_logit_scores"].to_numpy()
    barcode_table["solo_doublet_probability"] = doublet_obs["solo_doublet_probability"].to_numpy()
    barcode_table["solo_double_softmax_scores_original"] = doublet_obs["solo_double_softmax_scores_original"].to_numpy()
    barcode_table.to_csv(args.doublet_barcodes, index=False)

    n_cells_after_qc = int(full.n_obs)
    n_doublets = int(full.obs["solo_is_doublet"].sum())
    n_singlets = int(singlets.n_obs)
    summary = {
        "sample_id": args.sample_id,
        "n_cells_after_qc": n_cells_after_qc,
        "n_solo_doublets": n_doublets,
        "n_solo_singlets": n_singlets,
        "n_genes_full": int(full.n_vars),
        "n_genes_raw": int(full.raw.n_vars),
        "solo_doublet_pct": round(n_doublets / max(n_cells_after_qc, 1) * 100, 3),
        "annotated_h5ad": args.annotated_output,
        "output_h5ad": args.output,
    }

    with open(args.summary, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary.keys()))
        writer.writeheader()
        writer.writerow(summary)


if __name__ == "__main__":
    main()
