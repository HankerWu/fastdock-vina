#!/bin/sh
set -eu
# Run at the repository root. Argument is an empty absolute build directory.
root=$(pwd)
p=$1
d="$root/benchmarks/cavity_triplet"
old="$root/benchmarks/receptor_linkage/aurora/internal_loops"
mkdir -p "$p/src" "$p/tmp"
cp -R src tools "$p/src/"
cp "$old/proposal_crossover/cross_hunt.h" "$old/rigid_capture/inject.h" "$d/source_snapshot.h" "$p/src/src/lib/"
(cd "$p/src" && patch -p1 < "$old/proposal_crossover/calls.patch" && patch -p1 < "$old/rigid_capture/diagnostic.patch" && patch -p1 < "$old/rigid_capture/residency.patch" && patch -p1 < "$d/trace.patch")
shift
TMPDIR="$p/tmp" python "$p/src/tools/build.py" --output-dir "$p/bin" --jobs 1 "$@"
