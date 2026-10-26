#!/usr/bin/env bash
set -euo pipefail
export PATH=/opt/conda/bin:$PATH
mamba install -y -n base --override-channels -c conda-forge -c bioconda \
  'python=3.11' 'numpy>=1.26,<2' 'scipy>=1.11,<1.12' 'pandas>=2.2,<2.3' \
  bx-python bcbio-gff pybedtools gffutils gffread seqtk bedtools samtools minimap2 \
  biopython pysam scikit-learn seaborn matplotlib psutil pyyaml jinja2 \
  pandoc r-base=4.4 r-ggplot2 r-ggplotify r-scales r-reshape r-gridextra \
  bioconductor-noiseq r-rmarkdown r-htmltools r-dt r-plyr r-plotly r-dplyr \
  r-tidyr r-stringr r-optparse r-ggbeeswarm
python -m pip install --no-cache-dir 'psauron==1.0.8' 'TD2==1.0.8' 'gtfparse==2.7.1' 'pdf2image==1.17.0'
conda list --explicit > /opt/sqanti3-conda-explicit.txt
python -m pip freeze > /opt/sqanti3-pip-freeze.txt
conda clean -afy
