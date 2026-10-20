#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Merge Scanpy QC and SOLO summary CSV files.")
    parser.add_argument("--qc-files", nargs="*", default=[])
    parser.add_argument("--solo-files", nargs="*", default=[])
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_many(paths):
    frames = []
    for path in paths:
        if Path(path).is_file() and Path(path).stat().st_size > 0:
            frames.append(pd.read_csv(path))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def main():
    args = parse_args()
    qc = read_many(args.qc_files)
    solo = read_many(args.solo_files)

    if qc.empty and solo.empty:
        pd.DataFrame({"status": ["no_samples"]}).to_csv(args.output, index=False)
        return
    if qc.empty:
        solo.to_csv(args.output, index=False)
        return
    if solo.empty:
        qc.to_csv(args.output, index=False)
        return

    merged = qc.merge(solo, on="sample_id", how="outer", suffixes=("", "_solo"))
    merged.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
