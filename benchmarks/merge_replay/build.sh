#!/bin/sh
set -eu
root=$(git rev-parse --show-toplevel)
work=$(realpath -m "$1")
boost=$(realpath "$2")
mkdir "$work"
for arm in stable candidate; do
    rev=666f49d2570c047470d3c281df30d18c5ebd3ce1
    if [ "$arm" = candidate ]; then rev=0bcae6ef09857f97656291bb4821d20b3b619bf5; fi
    src=$work/$arm
    mkdir "$src"
    git -C "$root" archive "$rev" src tools/build.py | tar -xC "$src"
    python "$root/benchmarks/merge_replay/observe.py" "$src"
    python "$src/tools/build.py" --boost-prefix "$boost" --jobs "${JOBS:-3}" --output-dir "$src/build"
done
