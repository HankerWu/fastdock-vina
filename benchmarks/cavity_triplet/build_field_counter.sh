#!/bin/bash
# Link the single-worker counter to the same objects as build_observer.sh.
set -eu
src=$1
obj=$2
target=$3
boost=$4
mkdir -p "$target"
mkdir -p "$target/tmp"
export TMPDIR="$(cd "$target/tmp" && pwd)"
cp "$src/src/lib/grid.cpp" "$target/grid.cpp"
patch "$target/grid.cpp" < "${5:-benchmarks/cavity_triplet/field_count.patch}"
c++ -O3 -DNDEBUG -std=c++14 -pthread -fvisibility=hidden "-I$src/src/lib" "-I$boost/include" -c "$target/grid.cpp" -o "$target/grid.o"
objects=()
for f in "$obj"/*.o; do
    if [[ $(basename "$f") != grid.o ]]; then objects+=("$f"); fi
done
c++ -O3 -std=c++14 -pthread "$target/grid.o" "${objects[@]}" "-L$boost/lib" "-Wl,-rpath,$boost/lib" -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options -o "$target/vina"
