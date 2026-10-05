#!/usr/bin/env bash
set -euo pipefail

SOURCE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
build_dir="$SOURCE_ROOT/build-release"
host_prefix="/opt/Baltamatica"
jobs=4
cxx="${CXX:-}"
tests=ON

usage() {
    cat <<'EOF'
Usage: scripts/build.sh [options]
  --build-dir PATH  Build directory (default: build-release)
  --prefix PATH     Baltamatica installation (default: /opt/Baltamatica)
  --cxx PATH        C++ compiler (default: CXX, g++-9, then g++)
  --jobs N          Parallel build jobs (default: 4)
  --no-tests        Build without running the library/CLI regression check
  -h, --help        Show help
EOF
}
while (($#)); do
    case "$1" in
        --build-dir|--prefix|--cxx|--jobs)
            (($# >= 2)) || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --build-dir) build_dir="$2" ;;
                --prefix) host_prefix="$2" ;;
                --cxx) cxx="$2" ;;
                --jobs) jobs="$2" ;;
            esac
            shift 2 ;;
        --no-tests) tests=OFF; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo "--jobs must be positive" >&2; exit 2; }
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] ||
    { echo "This plugin release targets Linux x86_64." >&2; exit 1; }
if [[ -z "$cxx" ]]; then
    cxx="$(command -v g++-9 || command -v g++)"
fi
host_prefix="$(realpath -e -- "$host_prefix")"
[[ -f "$host_prefix/lib/cmake/Baltamatica/BaltamaticaConfig.cmake" ]] ||
    { echo "Baltamatica SDK was not found under $host_prefix" >&2; exit 1; }
build_dir="$(realpath -m -- "$build_dir")"
[[ "$build_dir" != "$SOURCE_ROOT" && "$build_dir" != / ]] ||
    { echo "Use a separate build directory." >&2; exit 1; }
cmake -S "$SOURCE_ROOT" -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    "-DCMAKE_CXX_COMPILER=$cxx" \
    "-DBaltamatica_DIR=$host_prefix/lib/cmake/Baltamatica" \
    -DHYMOS_BUILD_BALTAMATICA_PLUGIN=ON \
    "-DHYMOS_BUILD_TESTS=$tests"
cmake --build "$build_dir" --parallel "$jobs"
if [[ "$tests" == ON ]]; then
    ctest --test-dir "$build_dir" --output-on-failure
fi
printf '\nPlugin: %s/plugin/HyMoS\n' "$build_dir"
