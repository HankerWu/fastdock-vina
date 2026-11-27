#!/usr/bin/env bash
# BUILD is an unmodified tools/build.py product from exactly the same source tree.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
build=$1
boost=$2
out=$3
mkdir -p "$out"
patch -s -o "$out/quasi_newton.cpp" "$root/src/lib/quasi_newton.cpp" "$here/refine_observer.patch"
flags=(-std=c++14 -O3 -DNDEBUG -pthread -fvisibility=hidden "-I$root/src/lib" "-I$boost/include")
libs=("-L$boost/lib" "-Wl,-rpath,$boost/lib" -lboost_thread -lboost_filesystem -lboost_serialization -lboost_program_options)
c++ "${flags[@]}" -c "$out/quasi_newton.cpp" -o "$out/quasi_newton.o"
native=()
instrumented=()
for obj in "$build"/obj/*.o; do
   name=$(basename "$obj")
   if [[ $name != main.o && $name != split.o ]]; then native+=("$obj"); fi
   if [[ $name == quasi_newton.o ]]; then instrumented+=("$out/quasi_newton.o");
   elif [[ $name != split.o ]]; then instrumented+=("$obj"); fi
done
c++ "${flags[@]}" "${instrumented[@]}" "${libs[@]}" -o "$out/vina"
c++ "${flags[@]}" "$here/refine_graph.cpp" "${native[@]}" "${libs[@]}" -o "$out/graph"
