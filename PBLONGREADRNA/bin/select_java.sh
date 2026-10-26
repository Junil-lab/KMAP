#!/usr/bin/env bash
# Source from run.sh/setup.sh after defining PIPELINE_DIR.
java_ok() {
    local version
    [[ -x "$1" ]] || return 1
    version="$("$1" -version 2>&1 | awk -F '"' 'NR==1 {print $2}' | cut -d. -f1)"
    [[ "$version" =~ ^[0-9]+$ ]] && (( version >= 17 ))
}
select_java() {
    local candidate
    for candidate in "${JAVA_HOME:-/nonexistent}/bin/java" "$PIPELINE_DIR/.java17/bin/java" "$(command -v java || true)"; do
        if java_ok "$candidate"; then
            export JAVA_HOME="$(dirname "$(dirname "$(readlink -f "$candidate")")")"
            export PATH="$JAVA_HOME/bin:$PATH"
            return 0
        fi
    done
    return 1
}
