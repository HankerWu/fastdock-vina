#!/bin/bash
# Compare the official parser's own typed adjacency (including order) with explicit Pro.
set -eu
up=$1
boost=$2
work=$3
pro_probe=$4
export TMPDIR="$work/tmp"
mkdir -p "$TMPDIR" "$work/official_graph"
reader=$(dirname "$0")/../../receptor_linkage/aurora/internal_loops/rigid_graph_control/read_graph.cpp
# The upstream header has only the ordinary parser. Reuse the identical adjacency dumper.
sed '/^    model m = argc == 2 ? /,+1c\    model m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);' "$reader" > "$work/official_graph/reader.cpp"
objects=()
for obj in "$up"/build/obj/*.o; do
    if [[ $(basename "$obj") != main.o ]]; then objects+=("$obj"); fi
done
probe="$work/official_graph/read_graph"
c++ -std=c++14 -O3 -DNDEBUG -pthread "-I$up/src/lib" "-I$boost/include" \
    "$work/official_graph/reader.cpp" "${objects[@]}" "-L$boost/lib" "-Wl,-rpath,$boost/lib" \
    -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options -o "$probe"
for ligand in "$work"/inputs/*-ligand.pdbqt; do
    id=$(basename "$ligand" -ligand.pdbqt)
    "$probe" "$ligand" > "$work/official_graph/$id.csv"
    "$pro_probe" "$ligand" "${ligand%.pdbqt}.bonds" > "$work/official_graph/$id-pro.csv"
    cmp "$work/official_graph/$id.csv" "$work/official_graph/$id-pro.csv"
done
