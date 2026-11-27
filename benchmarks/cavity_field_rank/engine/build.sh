#!/bin/sh
# From repo root. Builds the opt-in field reranker on the coordinate/ID observer.
set -eu
out=$(realpath -m "$1")
shift
base=$PWD/benchmarks/cavity_triplet
rank=$PWD/benchmarks/cavity_field_rank/engine
loops=$PWD/benchmarks/receptor_linkage/aurora/internal_loops
mkdir -p "$out/src" "$out/tools" "$out/tmp"
cp -R src/lib src/main "$out/src/"
cp tools/build.py "$out/tools/"
cp "$loops/proposal_crossover/cross_hunt.h" "$loops/rigid_capture/inject.h" "$base/source_snapshot.h" "$rank/native_rank.h" "$out/src/lib/"
(
    cd "$out"
    for p in "$loops/proposal_crossover/calls.patch" "$loops/rigid_capture/diagnostic.patch" "$loops/rigid_capture/residency.patch" "$base/trace.patch" "$base/field_count.patch" "$rank/observer.patch"; do
        patch -p1 < "$p"
    done
)
TMPDIR="$out/tmp" python "$out/tools/build.py" --output-dir "$out/bin" "$@"
