#!/bin/sh
# Assemble the existing MC and AD4 observers plus a deferred seven-value input.
set -eu
src=$(pwd)
out=$1
shift
mkdir -p "$out/src" "$out/tools" "$out/tmp"
cp -R src/lib src/main "$out/src/"
cp tools/build.py "$out/tools/"
triplet="$src/benchmarks/cavity_triplet"
loops="$src/benchmarks/receptor_linkage/aurora/internal_loops"
cp "$loops/proposal_crossover/cross_hunt.h" "$loops/rigid_capture/inject.h" "$triplet/source_snapshot.h" "$out/src/lib/"
(cd "$out" && patch -p1 < "$loops/proposal_crossover/calls.patch" && patch -p1 < "$loops/rigid_capture/diagnostic.patch" && patch -p1 < "$loops/rigid_capture/residency.patch" && patch -p1 < "$triplet/trace.patch" && patch -p1 < "$triplet/field_count.patch" && patch -p1 < "$src/benchmarks/fragment_anchor/fkbp/live_ad4.patch")
TMPDIR="$out/tmp" python "$out/tools/build.py" --output-dir "$out/bin" "$@"
