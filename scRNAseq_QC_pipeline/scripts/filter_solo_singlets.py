#!/usr/bin/env python3
import argparse
import csv

import anndata as ad
import pandas as pd


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

    required = ["is_doublet"]
    missing = [col for col in required if col not in solo.obs.columns]
    if missing:
        raise RuntimeError(f"SOLO H5AD is missing required obs columns: {missing}")

    shared = full.obs_names.intersection(solo.obs_names)
    if len(shared) == 0:
        raise RuntimeError("No overlapping barcodes between full QC H5AD and SOLO H5AD.")

    if len(shared) != full.n_obs:
        full = full[shared].copy()
    else:
        full = full.copy()

    solo_obs = solo.obs.loc[full.obs_names]
    full.obs["solo_is_doublet"] = solo_obs["is_doublet"].astype(bool).to_numpy()
    for col in ("softmax_scores", "logit_scores"):
        if col in solo_obs.columns:
            full.obs[f"solo_{col}"] = pd.to_numeric(solo_obs[col], errors="coerce").to_numpy()

    full.write_h5ad(args.annotated_output, compression="gzip")

    singlets = full[~full.obs["solo_is_doublet"].to_numpy()].copy()
    singlets.write_h5ad(args.output, compression="gzip")

    doublet_obs = full.obs.loc[full.obs["solo_is_doublet"], :].copy()
    barcode_table = pd.DataFrame({"barcode": doublet_obs.index})
    if "solo_softmax_scores" in doublet_obs.columns:
        barcode_table["solo_softmax_scores"] = doublet_obs["solo_softmax_scores"].to_numpy()
    if "solo_logit_scores" in doublet_obs.columns:
        barcode_table["solo_logit_scores"] = doublet_obs["solo_logit_scores"].to_numpy()
    barcode_table.to_csv(args.doublet_barcodes, index=False)

    n_cells_after_qc = int(full.n_obs)
    n_doublets = int(full.obs["solo_is_doublet"].sum())
    n_singlets = int(singlets.n_obs)
    summary = {
        "sample_id": args.sample_id,
        "n_cells_after_qc": n_cells_after_qc,
        "n_solo_doublets": n_doublets,
        "n_solo_singlets": n_singlets,
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
