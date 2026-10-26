# KMAP PacBio IsoSeq: IsoQuant, isocall and FLAIR

Nextflow DSL2 + Singularity workflow for full-length non-concatemer (FLNC) PacBio reads. Three discovery methods run from the same complete input samples. Each method's transcript models receive SQANTI3 structural QC, TD2/PSAURON ORF prediction and an HTML report. CCS generation, demultiplexing and FLNC preparation are upstream steps.

```text
FLNC BAM → pbmm2 ISOSEQ (once per sample) ┬→ IsoQuant → SQANTI3 + ORF + QC
                                         ├→ isocall profile → merge → call → SQANTI3 + ORF + QC
                                         └→ FLAIR transcriptome -b → SQANTI3 + ORF + QC
```

No input FASTQ/FASTA conversion or additional read alignment is performed. Transcript FASTA and predicted protein FASTA are retained as analysis outputs.

For multiple input samples, additional cohort analyses run by default:

- FLAIR combines sample transcriptomes with their native read maps/counts and runs SQANTI3 on the combined GTF.
- isocall jointly merges all sample profiles, calls a cohort transcriptome/count matrix, and runs SQANTI3 on the cohort GTF.

Sample-level and cohort-level outputs are separate. Cohort QC rows describe the joint transcriptome, not the transcriptome of an individual sample. The isocall `prep-isoforms` annotation is prepared once and shared. Its sample-level calls use one profile each; cohort calls use all profiles.

## Requirements and setup

Linux, Bash, Python 3, Singularity, Java 17+, and these images:

| Image in `containers/` | Contents |
|---|---|
| `kmap_isoseq_main.sif` | pbmm2, samtools, IsoQuant, QC scripts' Python |
| `kmap_isoseq_discovery.sif` | isocall 1.3.0, FLAIR 3.0.0+master at pinned commit, minimap2/samtools |
| `kmap_isoseq_sqanti3.sif` | SQANTI3 6.0.1, TD2/PSAURON, R reporting |

```bash
bash setup.sh
# Installs Java 17 locally if needed, uses pinned Nextflow 26.04.3,
# and checks references and actual tool/dependency execution.
```

References default to hg38 FASTA + `.fai`, GENCODE v45 GTF, and optional human polyA/CAGE annotations. See `resources/ref/README.md`. Paths can be overridden using `--genome_fa`, `--gencode_gtf`, `--main_sif`, `--discovery_sif`, `--sqanti3_sif`, or a Nextflow `-c` config. Explicitly configured missing files stop validation. Missing default optional annotations are omitted; disable one with `--cage_peak_bed false`.

For a fresh Git checkout, obtain the three images separately or build them with
the recipes below before `setup.sh`. The images are not committed to Git; the
existing WGS release image does not contain this RNA workflow.

Container building/rebuilding:

```bash
bash def/build_containers.sh
bash def/build_containers.sh discovery --rebuild
```

A successful rebuild backs up the prior image. FLAIR's upstream commit and the local patch required for empty/no-SQ SAM streams are recorded in `def/kmap_isoseq_discovery.def` and `def/patches/`. The patch reproduces the existing lab installation's fix. Code/reference paths inside tasks are staged or configurable; no lab server path is required by the distribution source.

## Input and running

Copy `input/samples.example.csv` to a local CSV:

```csv
sample_id,flnc_bam,tissue
SampleA,/data/SampleA/flnc.bam,lung
SampleB,/data/SampleB/flnc.bam,kidney
```

All three columns are required. IDs must be unique, start with a letter/digit, and contain only letters/digits/`_`/`-`/`.`. Blank tissue becomes `NA`. Relative BAM paths resolve from the CSV directory. Standard CSV quoting, UTF-8 BOM and Windows line endings are supported.

```bash
# All three methods, all input reads, plus cohort analysis when there are >1 samples:
bash run.sh --samplesheet input/samples.csv --threads 16 --check-only
bash run.sh --samplesheet input/samples.csv --threads 16
bash run.sh --samplesheet input/samples.csv --threads 16 -resume

# Select a method; disable additional cohort analyses:
bash run.sh --samplesheet input/samples.csv --methods isoquant
bash run.sh --samplesheet input/samples.csv --methods isocall,flair --run_cohort false

# Cluster-specific config; default SLURM queue is cpu:
bash run.sh --samplesheet input/samples.csv -profile slurm -c cluster.config
```

`--threads` is a per-task CPU request, not a cap on the entire workflow. Set `executor.cpus`, `executor.memory`, and/or `process.maxForks` in a local config to control total concurrent use. `NXF_WORK` or `-work-dir` sets intermediates; `--outdir` sets published results. Resume requires the original working directory's `.nextflow/` cache and retained work files. JSON params files are supported; YAML additionally requires PyYAML on the launch host.

## Outputs

```text
<outdir>/samples/<sample_id>/
├── alignment/pbmm2/                  # shared alignment for all three methods
├── isoquant/{discovery,sqanti3}/
├── isocall/{profile,discovery,sqanti3}/
└── flair/{discovery,sqanti3}/
<outdir>/cohort/
├── isocall/{discovery,sqanti3}/
└── flair/{discovery,sqanti3}/
<outdir>/reference/isocall/known_isoforms.gz
<outdir>/qc_summary/                   # per-method CSVs and complete qc_summary.csv
<outdir>/pipeline_info/                # trace, timeline, report, run logs
```

Output layout changed from the prior IsoQuant-only version; that code is preserved under local `backups/`.

IsoQuant's `*.discovered_transcript_counts.tsv` / `*.discovered_transcript_tpm.tsv` correspond to `*.transcript_models.gtf` (known + novel models). `*.transcript_counts.tsv` / `*.transcript_tpm.tsv` are reference-only. The sample subdirectory produced by IsoQuant is retained inside `discovery/`.

FLAIR retains native `transcriptome` read counts (`*.isoform.counts.txt`) and cohort combine counts. Separate `quantify` realignment is not run. isocall produces its native count matrix. Counts across different discovery methods are not directly interchangeable.

QC contains `scope`, `sample_id`, `method`, tissue, alignment totals/rate and SQANTI3 categories. Cohort alignment counts sum all participating sample flagstats. Mapping rate counts mapped alignment records, including QC-passed and QC-failed records; it is not a deduplicated molecule/primary-read metric.

Core task failures retry up to twice, then terminate. QC errors terminate immediately. Final QC requires every expected sample × method row and enabled cohort row. A protein FASTA is optional when no coding isoforms are predicted; HTML reports are required.

## Verification

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
bash tests/run_smoke.sh
bash tests/run_selection_matrix.sh test_runs/smoke
bash tests/run_single_example.sh test_runs/smoke
```

The smoke test generates a synthetic 100 kb reference and two samples with 408 FLNC reads each. It exercises all three discovery methods, sample-level quantification, cohort combine/joint calling, SQANTI3/ORF/HTML, optional polyA/CAGE annotations and sample-complete QC. Results are isolated under `test_runs/`.

This checks functional execution, not biological accuracy or human-dataset runtime/memory requirements. SIFs, study inputs, reference data, results, local runtimes, work files and backups are excluded from Git. The real study sample list is local; `input/samples.example.csv` is the public template.

## Full-sample runtime benchmark

```bash
python3 benchmarks/run_sample.py --samplesheet input/samples.csv \
    --sample 25_CMC-R_1 --threads 16 --outdir test_runs/full_sample_run
```

This selects one complete sample and creates a separate samplesheet and fresh output/work directories. It runs all three methods without `-resume` and limits the local executor to one task at a time. Timing results are written to `benchmark_results.json`, with live state in `benchmark_status.json` and logs in `benchmark.log`.

The report separates input/alignment, discovery/quantification, and SQANTI3/ORF/report generation. The single pbmm2 alignment is shared by all three methods; that cost is included in each independent route total, and the full workflow wall time is also reported. Per-route active times exclude scheduler waits. Tool installation/container build time is excluded. The output directory must be new, to prevent inadvertent cache/output reuse.

## Single-sample and multi-sample analysis modes

Use `--analysis_mode` to select the FLAIR/isocall analysis scope:

| Mode | Independent sample calls and SQANTI3 | Cohort calls and SQANTI3 |
| --- | --- | --- |
| `single` | Yes, for every row in the CSV | No |
| `multi` | No final sample-level SQANTI3 | Yes; requires at least two samples |
| `both` | Yes | Yes; requires at least two samples |
| `auto` (default) | Yes | Enabled for multiple samples when `--run_cohort true` |

```bash
bash run.sh --samplesheet input/samples.csv --methods isocall,flair --analysis_mode single --outdir results/single
bash run.sh --samplesheet input/samples.csv --methods isocall,flair --analysis_mode multi --outdir results/multi
bash run.sh --samplesheet input/samples.csv --methods isocall,flair --analysis_mode both --outdir results/both
```

`single` can process several samples independently; to run just one sample, supply a one-row CSV. Explicit modes override the legacy `run_cohort` setting. IsoQuant remains sample-level in every mode if selected; multi/both require at least one of FLAIR/isocall.

In `multi`, FLAIR still builds each sample transcriptome/read map before `combine`; these intermediate discovery outputs remain published. isocall builds each sample profile, merges all profiles and calls only the joint transcriptome. Both methods retain sample counts in their cohort count outputs. FLAIR combine counts use native read maps; no separate quantify realignment is added. Only selected final scopes enter SQANTI3 and the merged QC summary.

Use a separate output directory for each mode so that files from a previous mode are not mistaken for current outputs. Keep the Nextflow work/cache directories and add `-resume` when resuming the same run.

## Tool selection quick reference

The default `--methods isoquant,isocall,flair` runs all three tools. Supply one
name or any comma-separated subset to run only those tools. Each selected final
transcriptome goes through SQANTI3, ORF prediction and QC. No unselected caller
is launched. Shared pbmm2 alignment runs once per input sample.

| Request | Options |
| --- | --- |
| All three, independent samples | `--analysis_mode single` |
| IsoQuant only | `--methods isoquant --analysis_mode single` |
| isocall only, joint calling | `--methods isocall --analysis_mode multi` |
| FLAIR only, independent samples | `--methods flair --analysis_mode single` |
| isocall + FLAIR, joint calling | `--methods isocall,flair --analysis_mode multi` |
| All three, sample + cohort outputs | `--analysis_mode both` |

The scope option applies jointly to selected isocall/FLAIR branches. IsoQuant
always remains per sample, including when selected alongside a multi-mode tool.
`multi` and `both` require at least two CSV rows and at least one of isocall/FLAIR.
The NAS-local `input/sample1.csv` selects study sample 10 and should use `single`
or `auto`. Public checkouts contain the generic `input/samples.example.csv`;
study sample lists and absolute NAS paths are not published.

The selection matrix reuses the synthetic smoke test's retained cache and checks
all seven nonempty caller subsets in single mode, each cohort-capable caller in
multi mode, and joint multi/both modes. Trace checks reject unselected callers;
QC rows, protein FASTA, HTML reports and cohort count headers are verified.
