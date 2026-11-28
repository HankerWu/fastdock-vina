#!/bin/bash
# Arguments: candidate's source, build, Boost prefix, prepared assay directory and output directory.
set -eu
src=$1
build=$2
boost=$3
prepared=$4
out=$5
mkdir -p "$out"
objects=()
for obj in "$build"/obj/*.o; do
    if [[ ${obj##*/} != main.o ]]; then objects+=("$obj"); fi
done
c++ -std=c++14 -O3 -DNDEBUG -I"$src/src/lib" -I"$boost/include" \
    benchmarks/fork_round/check_tangent.cpp "${objects[@]}" -pthread -L"$boost/lib" -Wl,-rpath,"$boost/lib" \
    -lboost_program_options -lboost_filesystem -lboost_serialization -lboost_system -lboost_thread -o "$out/tangent"
ligands=()
for p in "$prepared"/*.pdbqt; do
    if [[ ! -f "${p%.pdbqt}.box" ]]; then ligands+=("$p"); fi
done
"$out/tangent" "${ligands[@]}"
c++ -std=c++14 -O3 -I"$src/src/lib" -I"$boost/include" benchmarks/check_armijo_numeric.cpp \
    "$build/obj/quaternion.o" "$build/obj/random.o" -o "$out/armijo"
"$out/armijo"
