# Existing-script review and integration decisions

Reviewed the lab's existing `run_isoquant.sh`, `run_isoquant_novel.sh`, `run_isocall.sh`, `run_flair.sh`, `run_sqanti3_isocall.sh` and installed CLI/source before integration.

- Existing isocall ran prep/profile/merge but only printed the call command. The workflow executes call and checks its GTF/count outputs.
- All three callers now share one pbmm2 ISOSEQ alignment. FLAIR starts with transcriptome -b and retains its native read-assignment counts; no input FASTQ/FASTA conversion or quantify realignment is performed. Transcript and protein FASTA outputs remain necessary for downstream SQANTI3 ORF prediction.
- Existing FLAIR had per-sample transcriptome discovery, cohort combine and cohort quantification. Discovery and combine remain with separate sample/cohort QC. Native transcriptome counts replace the additional quantify realignment in the shared-BAM workflow.
- Installed FLAIR is `3.0.0+master`, commit `c17913a2cad4ae35e5903f23061e17ba0ef5631b`. `quantify -i` expects FASTA; BED goes to `--isoform_bed`. Confirmed by CLI help and `flair_quantify.py` argument parsing/minimap2 call. See https://flair.readthedocs.io/en/latest/modules.html#flair-quantify.
- Existing combine manifest referenced `.read.map.txt`; the installed transcriptome command writes `.isoform.read.map.txt`. The workflow stages the actual output.
- Existing FLAIR script expected `.counts.tpm.tsv`; this version actually writes `.tpm.tsv`. Output declarations use actual filenames.
- Existing installed FLAIR modified `filter_transcriptome_align.py` to accept empty/no-SQ SAM headers and iterate `until_eof=True`. The patch is retained explicitly in the image build, rather than depending on editable host source.
- Fixed missing-sample skips and cross-directory output assumptions by declaring all inputs/outputs in Nextflow. Method/cohort metadata is retained in QC; incomplete runs cannot produce a successful complete summary.
- The original source scripts are not modified. The previous validated IsoQuant-only workflow is backed up locally.

- SQANTI3 TD2 parser truncated colon-containing IDs, causing collisions and out-of-range CDS coordinates on full isocall data. `bin/run_sqanti3.py` patches the parser in memory to preserve full IDs before launching SQANTI3. The SIF remains unchanged. Validated all 56,287 TD2 ORFs against exact corrected FASTA IDs and bounds, including 71 IDs containing `::`.
