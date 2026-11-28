#!/bin/bash
# SRC exact candidate export, BUILD its tools/build.py output, ORIGINAL round-start CLI, BOOST Boost prefix.
set -eu
root=$(git rev-parse --show-toplevel)
D=$root/benchmarks/fork_round/source_audit
mkdir "$WORK"
tar xf "$root/benchmarks/flex_covalent_transfer/inputs.tar.xz" -C "$WORK" prepared/3FMZ/{partition_flex.pdbqt,partition_rigid.pdbqt,ligand.pdbqt} prepared/3FR2/ligand.pdbqt
map=$WORK/wall.C_H.map
cat > "$map" <<'EOF'
GRID_PARAMETER_FILE NULL
GRID_DATA_FILE NULL
MACROMOLECULE NULL
SPACING 1
NELEMENTS 2 2 2
CENTER 3000 0 0
EOF
for ((i=0;i<27;i++)); do echo 0 >> "$map"; done
objects=()
for o in "$BUILD"/obj/*.o; do [[ ${o##*/} == main.o ]] || objects+=("$o"); done
c++ -O3 -DNDEBUG -std=c++14 -I"$SRC/src/lib" -I"$BOOST/include" "$D/check.cpp" "${objects[@]}" -pthread \
  -L"$BOOST/lib" -Wl,-rpath,"$BOOST/lib" -lboost_program_options -lboost_filesystem -lboost_serialization -lboost_system -lboost_thread -o "$WORK/check"
"$WORK/check" "$WORK/prepared/3FMZ/ligand.pdbqt" "$WORK/prepared/3FR2/ligand.pdbqt" "$WORK/prepared/3FMZ/partition_rigid.pdbqt" "$WORK/prepared/3FMZ/partition_flex.pdbqt" "$WORK/wall" > "$WORK/results.tsv"
for flag in version help help_advanced; do
  "$BUILD/vina" --"$flag" > "$WORK/candidate-$flag"
  "$ORIGINAL" --"$flag" > "$WORK/original-$flag"
  cmp "$WORK/candidate-$flag" "$WORK/original-$flag"
done
