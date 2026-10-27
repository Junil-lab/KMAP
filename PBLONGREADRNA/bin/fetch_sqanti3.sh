#!/usr/bin/env bash
# Download the tested upstream SQANTI3 image; never replace an existing SIF.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ "${1:-}" == --help || "${1:-}" == -h ]]; then
    echo 'Usage: bash bin/fetch_sqanti3.sh [destination.sif]'; exit 0
fi
(( $# <= 1 )) || { echo 'ERROR: at most one destination argument is allowed' >&2; exit 2; }
DEST="${1:-$PIPELINE_DIR/containers/sqanti3_v6.0.1.sif}"
if [[ -e "$DEST" ]]; then
    [[ -s "$DEST" ]] || { echo "ERROR: existing file is empty: $DEST" >&2; exit 1; }
    echo "Using existing SIF: $DEST"; exit 0
fi
[[ "$(uname -m)" == x86_64 ]] || { echo 'ERROR: tested upstream image requires Linux x86_64' >&2; exit 1; }
SOURCE="$PIPELINE_DIR/containers/sqanti3_source.json"
URI="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["uri"])' "$SOURCE")"
mkdir -p "$(dirname "$DEST")"
FETCH_TMP="$(mktemp -d "$(dirname "$DEST")/.sqanti3-download.XXXXXX")"
trap 'rm -rf "$FETCH_TMP"' EXIT
singularity pull "$FETCH_TMP/image.sif" "$URI"
singularity inspect "$FETCH_TMP/image.sif" > /dev/null
# A hard link atomically installs the completed image and refuses an existing target.
ln "$FETCH_TMP/image.sif" "$DEST"
python3 - "$SOURCE" "$DEST" <<'PY'
import json,sys
from pathlib import Path
source,destination=map(Path,sys.argv[1:])
record=json.loads(source.read_text());record['installed_bytes']=destination.stat().st_size
Path(str(destination)+'.source.json').write_text(json.dumps(record,indent=2)+'\n')
PY
echo "Installed official SQANTI3: $DEST"
