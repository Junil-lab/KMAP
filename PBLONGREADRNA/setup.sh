#!/usr/bin/env bash
# Bootstrap Java/Nextflow locally, then validate references and container tools.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNLOAD=0
DOWNLOAD_SQANTI3=1
RELEASE_TAG=''
RELEASE_REPO='Junil-lab/KMAP'
ASSETS_DIR=''
usage() {
    cat <<'HELP'
Usage: bash setup.sh [options]
  (no options)                 Prepare Java/Nextflow, fetch missing official SQANTI3, check local files
  --download-containers        Download KMAP main/discovery SIFs from a published release
  --release-tag TAG            Required with --download-containers (no implicit latest)
  --release-repo OWNER/REPO    Default: Junil-lab/KMAP
  --assets-dir DIR            Install from an already downloaded release asset folder
  --no-download-sqanti3        Check local files without fetching missing SQANTI3
  -h, --help                   Show this help
Split SIF assets are reassembled and SHA256-checked automatically.
Existing matching images are reused; differing images are never overwritten.
Reference genome, annotation and sample BAMs must be provided separately.
HELP
}
while (( $# )); do
    case "$1" in
        --no-download-sqanti3) DOWNLOAD_SQANTI3=0; shift ;;
        --download-containers) DOWNLOAD=1; shift ;;
        --release-tag|--release-repo|--assets-dir)
            (( $# >= 2 )) && [[ -n "$2" && "$2" != -* ]] || { echo "ERROR: missing value for $1" >&2; exit 2; }
            case "$1" in
                --release-tag) RELEASE_TAG="$2" ;;
                --release-repo) RELEASE_REPO="$2" ;;
                --assets-dir) ASSETS_DIR="$(realpath -e "$2")" ;;
            esac
            shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done
if (( DOWNLOAD )) && [[ -n "$ASSETS_DIR" ]]; then
    echo 'ERROR: choose --download-containers or --assets-dir, not both' >&2; exit 2
fi
if (( DOWNLOAD )) && [[ -z "$RELEASE_TAG" ]]; then
    echo 'ERROR: --download-containers requires --release-tag TAG' >&2; exit 2
fi
echo 'KMAP IsoSeq setup — Java / Nextflow / containers / references'
echo '[1/6] Required host tools'
for tool in python3 singularity; do
    command -v "$tool" >/dev/null || { echo "[FAIL] $tool is required" >&2; exit 1; }
    echo "[OK] $tool"
done
echo '[2/6] Java 17+'
source "$PIPELINE_DIR/bin/select_java.sh"
if ! select_java; then
    command -v curl >/dev/null || { echo 'ERROR: curl required to install Java 17' >&2; exit 1; }
    case "$(uname -m)" in x86_64) ARCH=x64 ;; aarch64) ARCH=aarch64 ;; *) echo 'ERROR: install Java 17+ and set JAVA_HOME' >&2; exit 1 ;; esac
    JAVA_TMP="$(mktemp -d "$PIPELINE_DIR/.java17-install.XXXXXX")"
    trap 'rm -rf "$JAVA_TMP"' EXIT
    curl --fail --location --retry 3 -o "$JAVA_TMP/jdk.tar.gz" \
        "https://api.adoptium.net/v3/binary/latest/17/ga/linux/${ARCH}/jdk/hotspot/normal/eclipse"
    mkdir "$JAVA_TMP/jdk"
    tar -xzf "$JAVA_TMP/jdk.tar.gz" -C "$JAVA_TMP/jdk" --strip-components=1
    java_ok "$JAVA_TMP/jdk/bin/java" || { echo 'ERROR: Java installation failed' >&2; exit 1; }
    if [[ -e "$PIPELINE_DIR/.java17" ]]; then mv "$PIPELINE_DIR/.java17" "$PIPELINE_DIR/.java17.backup.$(date +%s)"; fi
    mv "$JAVA_TMP/jdk" "$PIPELINE_DIR/.java17"
    select_java
fi
java -version
echo '[3/6] Nextflow 26.04.3'
export NXF_VER=26.04.3
if [[ ! -x "$PIPELINE_DIR/nextflow" ]]; then
    curl --fail --location --retry 3 -o "$PIPELINE_DIR/nextflow.download" \
        "https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow"
    chmod +x "$PIPELINE_DIR/nextflow.download"
    mv "$PIPELINE_DIR/nextflow.download" "$PIPELINE_DIR/nextflow"
fi
"$PIPELINE_DIR/nextflow" -version
echo '[4/6] Singularity containers'
if (( DOWNLOAD )); then
    python3 "$PIPELINE_DIR/bin/container_release.py" install --destination "$PIPELINE_DIR/containers" \
        --repo "$RELEASE_REPO" --tag "$RELEASE_TAG"
elif [[ -n "$ASSETS_DIR" ]]; then
    python3 "$PIPELINE_DIR/bin/container_release.py" install --destination "$PIPELINE_DIR/containers" \
        --assets-dir "$ASSETS_DIR"
fi
if (( DOWNLOAD_SQANTI3 )); then
    bash "$PIPELINE_DIR/bin/fetch_sqanti3.sh"
fi
MISSING=0
for name in kmap_isoseq_main.sif kmap_isoseq_discovery.sif sqanti3_v6.0.1.sif; do
    image="$PIPELINE_DIR/containers/$name"
    if [[ -s "$image" ]]; then echo "[OK] $image";
    else echo "[FAIL] missing/empty $image"; MISSING=1; fi
done
if (( MISSING )); then
    echo 'Install the main/discovery release images, or build with: bash def/build_containers.sh'
    echo 'After the RNA release is published: bash setup.sh --download-containers --release-tag TAG'
    echo 'See docs/container_releases.md for exact asset names and split-file installation.'
fi
echo '[5/6] References and effective configuration'
CHECK_TMP="$(mktemp -d)"
trap 'rm -rf "$CHECK_TMP" "${JAVA_TMP:-}"' EXIT
"$PIPELINE_DIR/nextflow" config "$PIPELINE_DIR" -profile local -flat > "$CHECK_TMP/config"
if ! python3 "$PIPELINE_DIR/bin/check_config.py" "$CHECK_TMP/config" "$CHECK_TMP/images"; then
    echo '[FAIL] Supply missing files; see resources/ref/README.md. Reference data are not downloaded by setup.' >&2
    exit 1
fi
(( MISSING == 0 )) || exit 1
mapfile -t IMAGES < "$CHECK_TMP/images"
echo '[6/6] Execute container tool/dependency checks'
bash "$PIPELINE_DIR/bin/check_containers.sh" "${IMAGES[@]}"
echo '[OK] Environment checks passed. No sample analysis was started.'
echo 'Run: bash run.sh --samplesheet input/samples.csv'
