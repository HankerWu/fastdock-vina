#!/bin/sh
set -eu
root=$(pwd)
out=$1
shift
mkdir -p "$out/src" "$out/tools" "$out/tmp"
cp -R "$root/src/lib" "$root/src/main" "$out/src/"
cp "$root/tools/build.py" "$out/tools/"
d="$root/benchmarks/cavity_triplet/online"
cp "$d/live_match.h" "$out/src/lib/"
(cd "$out" && patch -p1 < "$d/mc.patch")
TMPDIR="$out/tmp" python "$out/tools/build.py" --output-dir "$out/bin" "$@"
