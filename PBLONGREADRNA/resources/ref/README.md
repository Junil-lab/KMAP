# Reference inputs

Default paths in `nextflow.config`:

- `hg38.fa` and `hg38.fa.fai`: uncompressed GRCh38/hg38 genome and its samtools FASTA index.
- `gencode.v45.annotation.gtf`: matching GENCODE v45 human gene annotation.
- Optional: `polyA_motifs.human.txt`, `atlas.clusters.2.0.GRCh38.96.chr.bed`, and `hg38_fair+new_CAGE_peaks_phase1and2.bed`.

Genome and GTF must use matching contig identifiers and assembly coordinates. Build the FASTA index with `samtools faidx hg38.fa`. References may be local files or symlinks; the workflow stages references and their index as process inputs. Override the paths using `--genome_fa` and `--gencode_gtf` or a Nextflow config file.

Reference data are distributed separately from the code. Keep a record of their source URLs and checksums when preparing a study. This directory's current local symlinks are not included in Git. Optional annotation files are omitted when absent by default; explicitly setting a missing path is an error.
