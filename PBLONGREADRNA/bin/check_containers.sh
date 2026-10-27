#!/usr/bin/env bash
set -euo pipefail
[[ $# == 2 || $# == 3 ]] || { echo 'usage: check_containers.sh MAIN.sif SQANTI3.sif [DISCOVERY.sif]' >&2; exit 2; }
for image in "$@"; do
    [[ -s "$image" ]] || { echo "ERROR: missing/empty container: $image" >&2; exit 1; }
done
singularity exec --cleanenv "$1" bash -euo pipefail -c '
    pbmm2 --version
    isoquant --version
    samtools --version > /dev/null
    python3 -c "import csv, pathlib"
'
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
singularity exec --cleanenv -B "$SCRIPT_DIR" "$2" bash -euo pipefail -c '
    source "$1/sqanti3_env.sh"
    python "$SQANTI3_QC_SCRIPT" --help > /dev/null
    TD2.LongOrfs --help > /dev/null
    TD2.Predict --help > /dev/null
    command -v gffread > /dev/null
    command -v pandoc > /dev/null
    Rscript -e '\''for (p in c("ggplot2","ggplotify","scales","reshape","gridExtra","NOISeq","rmarkdown","htmltools","DT","plyr","plotly","dplyr")) stopifnot(requireNamespace(p, quietly=TRUE))'\''
' _ "$SCRIPT_DIR"
if [[ $# == 3 ]]; then
    singularity exec --cleanenv "$3" bash -euo pipefail -c '
        flair --version
        flair transcriptome --help >/dev/null
        flair quantify --help >/dev/null
        isocall --version
        samtools --version >/dev/null
        minimap2 --version
    '
fi
echo 'Container checks passed'
