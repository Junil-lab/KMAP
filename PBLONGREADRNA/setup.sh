#!/usr/bin/env bash
# Bootstrap Java/Nextflow locally, then validate references and container tools.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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
export NXF_VER=26.04.3
if [[ ! -x "$PIPELINE_DIR/nextflow" ]]; then
    curl --fail --location --retry 3 -o "$PIPELINE_DIR/nextflow.download" \
        "https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow"
    chmod +x "$PIPELINE_DIR/nextflow.download"
    mv "$PIPELINE_DIR/nextflow.download" "$PIPELINE_DIR/nextflow"
fi
command -v python3 >/dev/null || { echo 'ERROR: python3 required' >&2; exit 1; }
command -v singularity >/dev/null || { echo 'ERROR: singularity required' >&2; exit 1; }
"$PIPELINE_DIR/nextflow" -version
CHECK_TMP="$(mktemp -d)"
trap 'rm -rf "$CHECK_TMP" "${JAVA_TMP:-}"' EXIT
"$PIPELINE_DIR/nextflow" config "$PIPELINE_DIR" -profile local -flat > "$CHECK_TMP/config"
python3 "$PIPELINE_DIR/bin/check_config.py" "$CHECK_TMP/config" "$CHECK_TMP/images"
mapfile -t IMAGES < "$CHECK_TMP/images"
bash "$PIPELINE_DIR/bin/check_containers.sh" "${IMAGES[@]}"
echo 'Environment checks passed. Run: bash run.sh --samplesheet input/samples.csv'
