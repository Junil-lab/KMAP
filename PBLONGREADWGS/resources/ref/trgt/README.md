# Download and merge trgt repeat catalogs

```bash
# download pathogenic TR catalog and add tag to include fail reads for GA-rich motifs
wget -q https://raw.githubusercontent.com/dashnowlab/STRchive/refs/tags/v2.4.4/data/STRchive-disease-loci.hg38.TRGT.bed
awk -v OFS="\t" \
  '/CANVAS_RFC1|FRDA_FXN|SCA27B_FGF14/ {print $0 ";INCLUDE_FAIL_READS"} !/CANVAS_RFC1|FRDA_FXN|SCA27B_FGF14/ {print}' \
  STRchive-disease-loci.hg38.TRGT.bed \
  > STRchive-disease-loci.hg38.TRGT.TAG_INCLUDE_FAIL_READS.bed

wget -q https://zenodo.org/records/8329210/files/adotto_repeats.hg38.bed.gz
gunzip adotto_repeats.hg38.bed.gz

# merge
cat adotto_repeats.hg38.bed STRchive-disease-loci.hg38.TRGT.TAG_INCLUDE_FAIL_READS.bed \
| bedtools sort \
  -faidx ./human_GRCh38_no_alt_analysis_set.fasta.fai \
  -i /dev/stdin \
| bgzip -c \
> adotto_strchive_20250827.hg38.bed.gz

# clean up
rm adotto_repeats.hg38.bed STRchive-disease-loci.hg38.TRGT.TAG_INCLUDE_FAIL_READS.bed

# check for duplicate TRIDs
zcat adotto_strchive_20250827.hg38.bed.gz | awk '{split($4,name,";");print name[1]}' | sort | uniq -d
```

## File used by this pipeline

`adotto_strchive_20250827.hg38.noFailTag.bed.gz` (937,195 loci): the catalog above from the
PacBio HiFi-human-WGS-WDL v4.0.0 reference bundle, with the `;INCLUDE_FAIL_READS` tag removed,
the same way the WDL's `process_trgt_catalog` task does it:

```bash
zcat adotto_strchive_20250827.hg38.bed.gz | sed 's/;INCLUDE_FAIL_READS//' | bgzip -c \
  > adotto_strchive_20250827.hg38.noFailTag.bed.gz
```
