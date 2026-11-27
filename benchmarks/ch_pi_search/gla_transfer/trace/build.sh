#!/bin/sh
# Build from the exact baseline used to generate the frozen GLA pools.
set -eu
out=$(realpath -m "$1")
shift
repo=$(git rev-parse --show-toplevel)
trace=$repo/benchmarks/ch_pi_search/gla_transfer/trace
mkdir -p "$out"
base=$(mktemp -d)
trap 'rm -rf "$base"' EXIT
git -C "$repo" archive 3f21e02acabf023965d8160e8d375e1efd8edc43 src tools/build.py \
    benchmarks/archive_overlap/development/observer.h \
    benchmarks/ch_pi_search/search.patch benchmarks/ch_pi_search/search_grid.h > "$base/source.tar"
tar xf "$base/source.tar" -C "$base"
cd "$base"
cp benchmarks/ch_pi_search/search_grid.h benchmarks/archive_overlap/development/observer.h "$trace/checkpoints.h" src/lib/
patch -p1 -i benchmarks/ch_pi_search/search.patch
patch -p1 -i "$trace/checkpoints.patch"
python tools/build.py --output-dir "$out" "$@"
