## NAS validation (2026-10-26)

All checks below ran from the NAS copy with NAS-hosted containers and data.
Machine-readable results: [nas_validation_20261026.json](nas_validation_20261026.json).

- 11 Python regression tests passed, including CLI missing values and help defaults.
- Fresh default run: all three tools, 33/33 real tasks completed with exit 0;
  6 sample-level and 2 cohort QC rows, native counts, ORF FASTA, HTML and optional
  CAGE/polyA annotations verified.
- 11 selection cases passed: all 7 nonempty caller subsets in single mode;
  isocall-only, FLAIR-only and all-three multi mode; all-three both mode.
- One-row CSV with default methods/auto mode passed: 3 sample QC rows and no
  cohort calls. Selection and one-row runs reused matching smoke tasks with
  `-resume`; trace checks verified unselected processes were absent.
- IsoQuant remained per-sample in multi/both; isocall/FLAIR cohort count headers
  retained both sample names.
- The NAS study sample-10 CSV and FLNC BAM passed path checks and samtools
  quickcheck; the complete human sample was not rerun for this release.

Reproduce with `tests/run_smoke.sh`, `tests/run_selection_matrix.sh` and
`tests/run_single_example.sh` in that order (see the main README).
These are functional tests on synthetic reads, not biological accuracy validation.

## Shared alignment validation (2026-10-21)

Current workflow: one pbmm2 ISOSEQ alignment per sample shared by IsoQuant, isocall and FLAIR transcriptome. No input FASTQ conversion or quantify realignment.

`test_runs/shared_alignment_smoke`: 33/33 tasks completed, exit 0; exactly two PBMM2 tasks for two full synthetic samples. All eight sample/cohort QC rows passed native count, protein FASTA, HTML and optional annotation checks. Eight Python unit tests passed.

Full human sample `25_CMC-R_1` launched separately in `test_runs/shared_alignment_full_20261021_25_CMC-R_1`; completion and timing must be read from that run, not inferred from this smoke validation.

Earlier validation below describes the previous alignment architecture.

# Local functional validation — 2026-10-21

- Three actual Singularity images: pbmm2/IsoQuant; FLAIR/isocall; SQANTI3/TD2/PSAURON/R.
- Two synthetic samples, each containing 408 FLNC reads against a 100 kb synthetic genome.
- All 40 analysis/QC tasks completed with exit 0, without cache reuse.
- All 40 tasks reused by explicit session-UUID resume.
- All 8 expected QC rows present: 2 samples × 3 methods, plus isocall and FLAIR cohort analyses.
- Every method/sample and cohort produced 36 classified isoforms (24 FSM, 12 NIC), protein predictions, and HTML reports.
- FLAIR sample/cohort counts and native TPM outputs checked; isocall call count matrices checked.
- IsoQuant model IDs matched discovered-transcript count IDs.
- PolyA/CAGE annotations checked; quoted CSV metadata preserved.
- Separate single-sample benchmark-driver test completed all 17 tasks with no cohort tasks and wrote per-method timing JSON.
- Eight regression tests passed, including missing/duplicate samples, malformed classification, failed-read counting, and missing method rows.

Evidence is stored locally under `test_runs/multitool_sif/results/validation_summary.json` and `test_runs/benchmark_smoke_checked/benchmark_results.json`. These are functional checks on synthetic data, not estimates of human full-sample runtime or biological accuracy. Actual full-sample timing is a separate benchmark.
