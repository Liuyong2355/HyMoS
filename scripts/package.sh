#!/usr/bin/env bash
set -euo pipefail

SOURCE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
build_dir="$SOURCE_ROOT/build-release"
output_dir="$SOURCE_ROOT/dist"
make_binary=true
make_source=false
epoch="${SOURCE_DATE_EPOCH:-0}"
version="$(sed -nE 's/^project\(HyMoS VERSION ([0-9]+\.[0-9]+\.[0-9]+).*/\1/p' "$SOURCE_ROOT/CMakeLists.txt")"

usage() {
    cat <<'EOF'
Usage: scripts/package.sh [options]
  --build-dir PATH  Configured build directory (default: build-release)
  --output-dir PATH  Archive directory (default: dist)
  --source          Also create the source archive
  --source-only     Create only the source archive; no build needed
  -h, --help        Show help

Run scripts/build.sh before packaging a binary release.
SOURCE_DATE_EPOCH controls archive timestamps; the default is 0.
The binary archive contains HyMoS/ at its top level.
EOF
}
while (($#)); do
    case "$1" in
        --build-dir|--output-dir)
            (($# >= 2)) || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --build-dir) build_dir="$2" ;;
                --output-dir) output_dir="$2" ;;
            esac
            shift 2 ;;
        --source) make_source=true; shift ;;
        --source-only) make_source=true; make_binary=false; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ "$epoch" =~ ^[0-9]+$ && -n "$version" ]] ||
    { echo "Invalid SOURCE_DATE_EPOCH or project version." >&2; exit 2; }
build_dir="$(realpath -m -- "$build_dir")"
mkdir -p -- "$output_dir"
output_dir="$(realpath -e -- "$output_dir")"
temporary=""
trap 'if [[ -n "$temporary" && -f "$temporary" ]]; then rm -f -- "$temporary"; fi' EXIT

archive() {
    local parent="$1" filename="$2"
    shift 2
    temporary="$(mktemp "$output_dir/.HyMoS-archive.XXXXXX")"
    tar --sort=name --mtime="@$epoch" --owner=0 --group=0 --numeric-owner \
        --format=posix --pax-option=exthdr.name=%d/PaxHeaders/%f,delete=atime,delete=ctime \
        --mode='u+rwX,go+rX,go-w' -C "$parent" -cf - "$@" |
        gzip -n >"$temporary"
    chmod 644 "$temporary"
    mv -f -- "$temporary" "$output_dir/$filename"
    temporary=""
    (cd -- "$output_dir" && sha256sum -- "$filename" >"$filename.sha256")
    printf '%s\n' "$output_dir/$filename"
}

if "$make_binary"; then
    [[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] ||
        { echo "This binary release targets Linux x86_64." >&2; exit 1; }
    [[ -f "$build_dir/CMakeCache.txt" ]] ||
        { echo "Build first with scripts/build.sh." >&2; exit 1; }
    # Reassemble so the archive includes the current documentation and licenses.
    cmake --build "$build_dir" --target hymos_baltamatica_bundle --parallel 4
    bundle="$build_dir/plugin/HyMoS"
    [[ -f "$bundle/main.so" && -f "$bundle/lib/libhymos.so.5" ]] ||
        { echo "The HyMoS bundle is incomplete." >&2; exit 1; }
    (cd -- "$bundle" && sha256sum --check --quiet SHA256SUMS)
    # Package only the solver. The SDK and OS supply these remaining libraries.
    for library in "$bundle/main.so" "$bundle/lib/libhymos.so.5"; do
        dependency_text="$(LC_ALL=C readelf -d "$library")"
        while IFS= read -r dependency; do
            case "$dependency" in
                libhymos.so.5|libbex.so|libstdc++.so.6|libgcc_s.so.1|libc.so.6|libm.so.6|libgomp.so.1) ;;
                *) echo "Unaccounted runtime dependency: $dependency ($library)" >&2; exit 1 ;;
            esac
        done < <(printf '%s\n' "$dependency_text" | sed -n 's/.*Shared library: \[\(.*\)\]/\1/p')
    done
    archive "$build_dir/plugin" "HyMoS-$version-linux-x86_64.tar.gz" HyMoS
fi

if "$make_source"; then
    source_files=(CMakeLists.txt)
    for entry in .gitattributes .gitignore README.md LICENSE LICENSE.md LICENSE.txt COPYING NOTICE NOTICE.md \
                 AUTHORS AUTHORS.md CITATION.cff AI_USAGE.md THIRD_PARTY_NOTICES.md \
                 apps cases docs examples licenses NRxx plugins scripts src tests submission; do
        [[ ! -e "$SOURCE_ROOT/$entry" ]] || source_files+=("$entry")
    done
    # An explicit allowlist excludes local builds, runs, backups and old packages.
    archive "$SOURCE_ROOT" "HyMoS-$version-source.tar.gz" \
        --transform='s,^,HyMoS/,' \
        --exclude='__pycache__' --exclude='*.pyc' --exclude='*.swp' \
        --exclude='.DS_Store' --exclude='Thumbs.db' \
        --exclude='*.aux' --exclude='*.log' --exclude='*.out' --exclude='*.toc' \
        --exclude='*.bbl' --exclude='*.blg' --exclude='*.fls' \
        --exclude='*.fdb_latexmk' --exclude='*.synctex.gz' --exclude='*.xdv' \
        --exclude='tests/output' \
        --exclude='examples/*/.control.*' --exclude='examples/*/.hymos_control_runtime.sh' \
        --exclude='examples/*/hymos_control.sh' --exclude='examples/*/run_*.sh' \
        --exclude='examples/*/*.dat' --exclude='examples/*/*.bin' --exclude='examples/*/res_step' \
        --exclude='examples/*/runs' \
        "${source_files[@]}"
fi
