#!/usr/bin/env bash
# Recompute the primary endpoints from archived inputs and raw outputs.
set -eu
repo=$(git rev-parse --show-toplevel)
d="$repo/benchmarks/bond_distance"
w=$(mktemp -d)
trap 'rm -rf "$w"' EXIT
py=${PYTHON:-python3}
export OPENBLAS_NUM_THREADS=1
mkdir "$w/macroinput" "$w/kdinput"
tar xJf "$repo/benchmarks/macrocycle_labels/observations.tar.xz" -C "$w/macroinput"
tar xJf "$d/macro/extra_inputs.tar.xz" -C "$w/macroinput"
tar xJf "$repo/benchmarks/search_backtracking/assay/inputs.tar.xz" -C "$w/kdinput"
for study in macro tail tail/cost kd; do
    out="$w/results/$study"
    mkdir -p "$out"
    cp "$d/$study/runs.csv" "$out"
    tar xJf "$d/$study/raw.tar.xz" -C "$out"
    if [ "$study" = kd ]; then
        "$py" "$d/analyze.py" kd "$w/kdinput" "$out"
        cmp "$out/kd_summary.csv" "$d/kd/kd_summary.csv"
        cmp "$out/kd_pairs.csv" "$d/kd/kd_pairs.csv"
    else
        "$py" "$d/analyze.py" macro "$w/macroinput" "$out"
        cmp "$out/models.csv" "$d/$study/models.csv"
        cmp "$out/structures.csv" "$d/$study/structures.csv"
    fi
done
"$py" "$d/chem_graph.py" "$w/macroinput" "$d/macro/parent_pairs.csv" "$w/chemical_pairs.csv"
cmp "$w/chemical_pairs.csv" "$d/chemical_pairs.csv"
mkdir "$w/ordinary"
tar xJf "$d/pose/observations.tar.xz" -C "$w/ordinary"
"$py" "$d/pose/audit.py" "$w/ordinary"
echo 'Archived poses, Kd pairs, closure chemistry and phase counters verified.'
