#!/usr/bin/env bash
# Clean builds of the production comparator and the causal controls.
set -euo pipefail
work=$(realpath -m "$1")
shift
repo=$(git rev-parse --show-toplevel)
d="$repo/benchmarks/archive_overlap"
base=21714700baecee21e617f0c472b1136e4be33b4a
fixed=40b25a5a0aff981a06bd8222d11ec075a4405fb8
for arm in parent fixed observed-parent observed-fixed save merge hold counted-hold; do
    source="$work/$arm/src-tree"
    mkdir -p "$work/$arm"
    mkdir "$source"
    revision=$fixed
    if [[ $arm == *parent ]]; then revision=$base; fi
    git -C "$repo" archive "$revision" src tools/build.py | tar -x -C "$source"
    case "$arm" in
        observed-*) python3 "$d/development/instrument.py" "$source" ;;
        save|merge) (cd "$repo" && python3 "$d/development/ablate.py" "$source" "$arm") ;;
        hold|counted-hold)
            (cd "$source" && patch -p1 -i "$d/refill/once_full.patch")
            if [[ $arm == counted-hold ]]; then (cd "$source" && patch -p1 -i "$d/refill/count.patch"); fi ;;
    esac
    python3 "$source/tools/build.py" --output-dir "$work/$arm" "$@"
done
