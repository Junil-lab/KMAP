# Containers: KMAP release + official SQANTI3

SQANTI3 is downloaded directly from its authors' Docker Hub image and converted
to SIF. KMAP's release only needs to distribute the main and discovery images.
No SQANTI3 split-file upload is needed.

## Users

Supply the references described in `resources/ref/README.md`, then run:

```bash
bash setup.sh
```

Setup checks host tools, prepares Java/Nextflow, automatically fetches a missing
`containers/sqanti3_v6.0.1.sif`, checks references and runs dependency checks.
Main/discovery SIFs must already be present, or can be fetched from a KMAP RNA
release **once that release is published**:

```bash
bash setup.sh --download-containers --release-tag isoseq-v1.0.0
```

`isoseq-v1.0.0` is the dedicated KMAP RNA container release tag.
The existing WGS `v1.0.0` images do not contain this RNA workflow. An explicit
release tag is required to avoid picking another module's latest release.
The default repository is `Junil-lab/KMAP` (`--release-repo` can
change it). Automatic KMAP release downloads support public releases.

For an already downloaded main/discovery asset folder:

```bash
bash setup.sh --assets-dir /path/to/downloaded-assets
```

This still downloads official SQANTI3 if missing. For offline checks, first
place all three SIFs in their expected locations and use:

```bash
bash setup.sh --no-download-sqanti3
```

Java/Nextflow must also be installed for fully offline setup. Reference FASTA,
GTF and sample BAM data are not downloaded by this script.

## Official SQANTI3 provenance and compatibility

The tested release is `anaconesalab/sqanti3:v6.0.1`, pinned by immutable OCI digest:

```text
sha256:e770fb7a136e4c67abb6b50d7a118bd8724d430dee49d4ceb2ef306fae7767f4
```

`containers/sqanti3_source.json` records the source. To fetch only this image:

```bash
bash bin/fetch_sqanti3.sh
```

Singularity performs the conversion; Docker and a host Conda environment are not
required. Downloads use a temporary file and only install a successfully inspected
SIF. Existing files are preserved. The tested image is Linux amd64 and is large:
about 7.48 GB of compressed OCI layers; the tested SIF was 7,170,994,176 bytes.
Allow additional temporary/cache disk space for extraction and SIF creation.
Network speed can dominate installation time. Downloading and conversion are not
biological analysis time.

The official image's Conda PATH is selected inside the container because Nextflow
uses `singularity exec`, bypassing Docker's Conda entrypoint. `sqanti3_qc=auto`
detects its script path. The existing in-process ORF parser fix preserves colon-
containing transcript IDs without modifying the SIF. The legacy KMAP image remains
supported via `--sqanti3_sif ...`; the default uses a distinct official filename so
old task cache entries for the legacy container are not accidentally reused.

Use the OCI digest for reproducible source selection; generated SIF file hashes
can differ with conversion metadata. `containers/SHA256SUMS` lists the two KMAP
release images, not a universal checksum for every locally converted SQANTI3 SIF.
Official instructions: https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation

## Maintainers: prepare a KMAP RNA release

```bash
python3 bin/container_release.py package \
  --source containers --output release_assets/isoseq-v1.0.0
```

Use a new output directory. With the current main/discovery images, upload all
four output files to the same dedicated RNA release:

```text
kmap_isoseq_main.sif
kmap_isoseq_discovery.sif
isoseq-containers.json
SHA256SUMS
```

The installer verifies asset sizes/hashes and final SIF hashes, preserves matching
existing images and rejects differing existing images. Main/discovery images
larger than 1900 MiB would automatically be split and restored by the same helper.
GitHub requires each asset to be under 2 GiB:
https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases

The current publisher does not package SQANTI3, references or study data.
`def/build_containers.sh` builds KMAP main/discovery images and fetches the official
SQANTI3 image. `--rebuild` applies to KMAP-built images; it does not overwrite an
existing upstream SIF. The legacy SQANTI3 definition remains for audit/custom use
but is no longer the default build route. Packaging does not publish a release.
