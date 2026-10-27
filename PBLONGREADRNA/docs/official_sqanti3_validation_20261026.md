# Official SQANTI3 v6.0.1 validation

The default SQANTI3 container is converted from the authors' Linux amd64 image:

```text
docker://anaconesalab/sqanti3@sha256:e770fb7a136e4c67abb6b50d7a118bd8724d430dee49d4ceb2ef306fae7767f4
```

The tested SIF was 7,170,994,176 bytes. `setup.sh` pins this OCI digest and
downloads the image through `bin/fetch_sqanti3.sh`. The KMAP code selects the
image's Conda environment when it is run with `singularity exec`; it does not
modify the upstream SIF.

## Functional compatibility

Validation on 2026-10-26 produced these results:

- Container dependency checks passed for SQANTI3, TD2/PSAURON, gffread,
  pandoc and the required R packages.
- A fresh two-sample workflow completed all 33 tasks. IsoQuant, isocall and
  FLAIR generated six sample-level and two cohort SQANTI3/QC results, including
  ORF FASTA and HTML reports.
- The official image was compared with the previously validated image in eight
  sample/cohort caller cases. All 36 isoform records per case matched for
  structural category, coding status, ORF length and CDS coordinates; all ten
  predicted protein sequences per case were identical.
- A dedicated isocall fixture retained all 36 transcript IDs containing a
  colon. All ten coding IDs matched their protein FASTA IDs.
- Twenty-two Python regression tests passed, including release installation,
  official image fetching, caller/mode selection and input/QC validation.

SQANTI3 6.0.1's TD2 result parser originally stopped transcript IDs at a colon.
`bin/run_sqanti3.py` applies a narrow in-process parser correction from
`[^\\s:]+` to `[^\\s]+`, preserving the complete isocall ID. It refuses to run
if the installed parser has an unrecognized shape, so a future incompatible
upstream change cannot silently alter IDs.

## Twelve-sample isocall cohort check

An additional real-data check completed on 2026-10-27 using this exact route:

```text
12 alignments -> 12 isocall profiles -> one merge -> one call
              -> one SQANTI3 QC/ORF run
```

The merged call contained 211,617 isoforms, including 15,611 IDs with a colon.
SQANTI3 classified all 211,617 IDs and predicted 164,611 coding isoforms. The
164,611 protein FASTA identifiers exactly matched the coding identifiers. The
isocall count matrix retained all 12 sample columns, and the automated check
reported no errors.

These checks establish workflow and identifier compatibility. They do not
measure biological accuracy or validate every possible reference and dataset.
