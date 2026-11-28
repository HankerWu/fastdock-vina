#!/usr/bin/env bash
set -eu
root=$(git rev-parse --show-toplevel)
out=$(realpath -m "$1")
boost=$(realpath "$2")
mkdir "$out"
for spec in stable:666f49d2570c047470d3c281df30d18c5ebd3ce1 candidate:0bcae6ef09857f97656291bb4821d20b3b619bf5; do
    arm=${spec%%:*}
    dst="$out/$arm"
    mkdir "$dst"
    git -C "$root" archive "${spec#*:}" src tools/build.py | tar -xC "$dst"
    python "$root/benchmarks/bfgs_source_replay/instrument.py" "$dst"
    python "$dst/tools/build.py" --boost-prefix "$boost" --output-dir "$dst/build" --jobs "${JOBS:-4}"
done
