#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser(description="Run basic scRNA-seq QC and create HVG input for SOLO.")
    parser.add_argument("--input", required=True, help="CellBender/10x H5 or H5AD input")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-full", required=True, help="Full post-QC H5AD output")
    parser.add_argument("--output-hvg", required=True, help="HVG subset H5AD output for SOLO")
    parser.add_argument("--metrics", required=True, help="Per-sample QC metrics CSV")
    parser.add_argument("--mt-threshold", type=float, default=None, help="Optional mitochondrial percentage cutoff; disabled by default")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-counts", type=int, default=500)
    parser.add_argument("--min-cells-gene", type=int, default=5)
    parser.add_argument("--hvg-n-top", type=int, default=4000)
    return parser.parse_args()


def read_input(path):
    suffixes = "".join(Path(path).suffixes)
    if suffixes.endswith(".h5ad"):
        return ad.read_h5ad(path)
    return sc.read_10x_h5(path, gex_only=True)


def detect_mt_genes(adata):
    names = pd.Index(adata.var_names.astype(str))
    upper_names = names.str.upper()
    mask = (
        names.str.startswith(("MT-", "mt-", "Mt-", "MT_", "mt_"))
        | upper_names.str.startswith("MTRNR")
    )
    if "gene_ids" in adata.var:
        gene_ids = pd.Index(adata.var["gene_ids"].astype(str))
        gene_upper = gene_ids.str.upper()
        mask = mask | gene_ids.str.startswith(("MT-", "mt-", "Mt-", "MT_", "mt_")) | gene_upper.str.startswith("MTRNR")
    return np.asarray(mask, dtype=bool)


def ensure_counts_layer(adata):
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()


def sparse_var(matrix, mean):
    squared = matrix.copy()
    squared.data **= 2
    return np.asarray(squared.mean(axis=0)).ravel() - np.square(mean)


def set_top_hvg_by_dispersion(adata, n_top):
    counts = adata.layers["counts"]
    mean = np.asarray(counts.mean(axis=0)).ravel()
    if sparse.issparse(counts):
        var = sparse_var(counts, mean)
    else:
        var = np.asarray(counts.var(axis=0)).ravel()
    score = var / (mean + 1e-8)
    score[~np.isfinite(score)] = -np.inf
    top = np.argsort(score)[-n_top:]
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[top] = True
    adata.var["highly_variable"] = mask
    return "dispersion_fallback"


def compute_hvg_stable(adata, n_top):
    ensure_counts_layer(adata)
    n_top = min(int(n_top), adata.n_vars)
    if n_top <= 0:
        raise ValueError("No genes remain after QC filtering.")
    if adata.n_vars <= n_top:
        adata.var["highly_variable"] = True
        return "all_genes"

    errors = []
    for span in (0.3, 0.5, 0.8):
        try:
            sc.pp.highly_variable_genes(
                adata,
                flavor="seurat_v3",
                n_top_genes=n_top,
                layer="counts",
                subset=False,
                span=span,
            )
            if int(adata.var["highly_variable"].sum()) > 0:
                return f"seurat_v3_span_{span}"
        except Exception as exc:
            errors.append(f"seurat_v3 span={span}: {type(exc).__name__}: {exc}")

    try:
        tmp = adata.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        sc.pp.highly_variable_genes(
            tmp,
            flavor="seurat",
            n_top_genes=n_top,
            subset=False,
        )
        adata.var["highly_variable"] = tmp.var["highly_variable"].to_numpy(dtype=bool)
        if int(adata.var["highly_variable"].sum()) > 0:
            return "seurat_log1p_fallback"
    except Exception as exc:
        errors.append(f"seurat fallback: {type(exc).__name__}: {exc}")

    method = set_top_hvg_by_dispersion(adata, n_top)
    if errors:
        adata.uns["hvg_fallback_errors"] = errors
    return method


def summarize_vector(values, fn, default="NA"):
    if len(values) == 0:
        return default
    return float(fn(values))


def main():
    args = parse_args()
    adata = read_input(args.input)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs_names = [f"{args.sample_id}_{barcode}" for barcode in adata.obs_names.astype(str)]
    adata.obs["sample_id"] = args.sample_id

    ensure_counts_layer(adata)
    adata.var["mt"] = detect_mt_genes(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)

    n_cells_input = int(adata.n_obs)
    n_genes_input = int(adata.n_vars)
    total_counts_input_values = adata.obs["total_counts"].to_numpy()
    n_genes_input_values = adata.obs["n_genes_by_counts"].to_numpy()
    mt_values = adata.obs["pct_counts_mt"].to_numpy()
    cells_mt_ge_threshold = (
        int((adata.obs["pct_counts_mt"] >= args.mt_threshold).sum())
        if args.mt_threshold is not None else ""
    )

    keep_cells = (
        (adata.obs["n_genes_by_counts"] >= args.min_genes)
        & (adata.obs["total_counts"] >= args.min_counts)
    )
    if args.mt_threshold is not None:
        keep_cells &= adata.obs["pct_counts_mt"] < args.mt_threshold
    adata = adata[keep_cells].copy()
    if adata.n_obs == 0:
        raise RuntimeError(
            f"{args.sample_id}: no cells remain after QC "
            f"(mt filter: {args.mt_threshold if args.mt_threshold is not None else 'off'}, "
            f"genes >= {args.min_genes}, counts >= {args.min_counts})"
        )

    # Snapshot corrected counts after CELL QC, before any GENE filtering.
    # AnnData.raw preserves all variables when the training copy is subset.
    adata.raw = adata.copy()
    adata.uns["counts_provenance"] = "input counts after cell QC; no gene filtering or normalization"
    solo_adata = adata.copy()
    sc.pp.filter_genes(solo_adata, min_cells=args.min_cells_gene)
    if solo_adata.n_vars == 0:
        raise RuntimeError(f"{args.sample_id}: no training genes remain after min_cells={args.min_cells_gene}.")

    hvg_method = compute_hvg_stable(solo_adata, args.hvg_n_top)
    hvg_mask = solo_adata.var["highly_variable"].to_numpy(dtype=bool)
    hvg_adata = solo_adata[:, hvg_mask].copy()
    adata.var["solo_gene_filter_pass"] = adata.var_names.isin(solo_adata.var_names)
    adata.var["highly_variable"] = adata.var_names.isin(hvg_adata.var_names)
    if "hvg_fallback_errors" in solo_adata.uns:
        adata.uns["hvg_fallback_errors"] = solo_adata.uns["hvg_fallback_errors"]
    assert adata.n_vars == n_genes_input == adata.raw.n_vars == hvg_adata.raw.n_vars

    adata.write_h5ad(args.output_full, compression="gzip")
    hvg_adata.write_h5ad(args.output_hvg, compression="gzip")

    metrics = {
        "sample_id": args.sample_id,
        "input_h5": args.input,
        "n_cells_input": n_cells_input,
        "n_genes_input": n_genes_input,
        "n_cells_after_qc": int(adata.n_obs),
        "n_genes_after_qc": int(adata.n_vars),
        "n_genes_solo_prefilter": int(solo_adata.n_vars),
        "n_genes_raw": int(adata.raw.n_vars),
        "n_hvg": int(hvg_adata.n_vars),
        "removed_cells_total": int(n_cells_input - adata.n_obs),
        "removed_cells_total_pct": round((n_cells_input - adata.n_obs) / max(n_cells_input, 1) * 100, 3),
        "cells_mt_ge_threshold": cells_mt_ge_threshold,
        "mt_threshold": args.mt_threshold if args.mt_threshold is not None else "",
        "mt_filter_enabled": args.mt_threshold is not None,
        "min_genes": args.min_genes,
        "min_counts": args.min_counts,
        "min_cells_gene": args.min_cells_gene,
        "hvg_n_top": args.hvg_n_top,
        "hvg_method": hvg_method,
        "median_total_counts_input": summarize_vector(total_counts_input_values, np.median),
        "median_n_genes_input": summarize_vector(n_genes_input_values, np.median),
        "median_pct_mt_input": summarize_vector(mt_values, np.median),
        "median_total_counts_after_qc": summarize_vector(adata.obs["total_counts"].to_numpy(), np.median),
        "median_n_genes_after_qc": summarize_vector(adata.obs["n_genes_by_counts"].to_numpy(), np.median),
        "median_pct_mt_after_qc": summarize_vector(adata.obs["pct_counts_mt"].to_numpy(), np.median),
    }

    with open(args.metrics, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics.keys()))
        writer.writeheader()
        writer.writerow(metrics)


if __name__ == "__main__":
    main()
