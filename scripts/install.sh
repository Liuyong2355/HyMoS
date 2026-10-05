#!/usr/bin/env bash
set -euo pipefail

SCRIPT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
host_prefix="/opt/Baltamatica"
plugin_dir=""
replace=false
dry_run=false
prefix_given=false
if [[ -f "$SCRIPT_ROOT/main.so" ]]; then
    bundle="$SCRIPT_ROOT"
else
    bundle="$SCRIPT_ROOT/build-release/plugin/HyMoS-Baltamatica"
fi

usage() {
    cat <<'EOF'
Usage: scripts/install.sh [options]
  --prefix PATH      Baltamatica installation (default: /opt/Baltamatica)
  --plugin-dir PATH  Plugin search directory; HyMoS-Baltamatica is created inside it
  --bundle PATH      Built or extracted HyMoS-Baltamatica plugin directory
  --replace          Replace a marked HyMoS-Baltamatica package, retaining a backup
  --dry-run          Check the package and destination without installing
  -h, --help         Show help

The default destination is <prefix>/plugins/HyMoS-Baltamatica.
An existing unmarked directory is never replaced.
libbex.so is provided by the host installation, not by this package.
EOF
}
while (($#)); do
    case "$1" in
        --prefix|--plugin-dir|--bundle)
            (($# >= 2)) || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --prefix) host_prefix="$2"; prefix_given=true ;;
                --plugin-dir) plugin_dir="$2" ;;
                --bundle) bundle="$2" ;;
            esac
            shift 2 ;;
        --replace) replace=true; shift ;;
        --dry-run) dry_run=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done
if [[ -z "$plugin_dir" ]] || "$prefix_given"; then
    host_prefix="$(realpath -e -- "$host_prefix")"
    [[ -f "$host_prefix/lib/libbex.so" ]] ||
        { echo "Host library not found: $host_prefix/lib/libbex.so" >&2; exit 1; }
fi
[[ -n "$plugin_dir" ]] || plugin_dir="$host_prefix/plugins"
bundle="$(realpath -e -- "$bundle")"
[[ -f "$bundle/.hymos-plugin" && -f "$bundle/main.so" &&
   -f "$bundle/lib/libhymos.so.5" && -f "$bundle/SHA256SUMS" ]] ||
    { echo "This is not a complete HyMoS-Baltamatica plugin package: $bundle" >&2; exit 1; }
grep -qx 'name=HyMoS-Baltamatica' "$bundle/.hymos-plugin" ||
    { echo "The package identity is not HyMoS-Baltamatica." >&2; exit 1; }
(cd -- "$bundle" && sha256sum --check --quiet SHA256SUMS)
plugin_dir="$(realpath -m -- "$plugin_dir")"
destination="$plugin_dir/HyMoS-Baltamatica"
[[ "$bundle" != "$destination" ]] ||
    { echo "The source package is already at the destination." >&2; exit 1; }
[[ ! -L "$destination" ]] ||
    { echo "Refusing to replace a symbolic link: $destination" >&2; exit 1; }
if [[ -e "$destination" ]]; then
    [[ -d "$destination" && -f "$destination/.hymos-plugin" ]] &&
        grep -qx 'name=HyMoS-Baltamatica' "$destination/.hymos-plugin" ||
        { echo "Refusing to replace an unmarked plugin: $destination" >&2; exit 1; }
    "$replace" ||
        { echo "HyMoS-Baltamatica is already installed; use --replace to keep a backup and update." >&2; exit 1; }
fi
if "$dry_run"; then
    printf 'Package verified. Destination: %s\n' "$destination"
    exit 0
fi

mkdir -p -- "$plugin_dir"
stage="$(mktemp -d "$plugin_dir/.HyMoS-Baltamatica-install.XXXXXX")"
backup=""
cleanup() {
    # stage is created by mktemp inside the resolved plugin directory.
    if [[ -d "$stage" && "$stage" == "$plugin_dir"/.HyMoS-Baltamatica-install.* ]]; then
        rm -rf -- "$stage"
    fi
}
trap cleanup EXIT
cp -a -- "$bundle" "$stage/HyMoS-Baltamatica"
(cd -- "$stage/HyMoS-Baltamatica" && sha256sum --check --quiet SHA256SUMS)
if [[ -e "$destination" ]]; then
    backup="$(mktemp -d "$plugin_dir/.HyMoS-Baltamatica-backup.XXXXXX")"
    mv -- "$destination" "$backup/HyMoS-Baltamatica"
fi
if ! mv -- "$stage/HyMoS-Baltamatica" "$destination"; then
    if [[ -n "$backup" && ! -e "$destination" ]]; then
        mv -- "$backup/HyMoS-Baltamatica" "$destination"
    fi
    echo "Installation failed; the previous package was restored when possible." >&2
    exit 1
fi
printf 'Installed: %s\n' "$destination"
[[ -z "$backup" ]] || printf 'Previous package: %s/HyMoS-Baltamatica\n' "$backup"
printf "In Baltamatica: load_plugin('HyMoS-Baltamatica'); hymos_setup('1D');\n"
