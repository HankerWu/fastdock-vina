#!/usr/bin/env bash
# Same objects as stable Pro, replace two translation units with passive observers.
set -eu
build=$1
boost=$2
out=$3
root=$PWD
here=$root/benchmarks/ache_receptor_state/trace
basin=$root/benchmarks/receptor_linkage/aurora/internal_loops/basin_trace
mkdir -p "$out/src/lib"
cp src/lib/{vina,monte_carlo}.cpp "$out/src/lib/"
(cd "$out" && patch -s -p1 < "$here/observer.patch")
flags=(-std=c++14 -DNDEBUG -O3 -pthread -fvisibility=hidden "-I$root/src/lib" "-I$boost/include" "-I$here" "-I$basin")
c++ "${flags[@]}" -c "$out/src/lib/monte_carlo.cpp" -o "$out/monte_carlo.o" &
worker=$!
c++ "${flags[@]}" -c "$out/src/lib/vina.cpp" -o "$out/vina.o"
wait "$worker"
objects=()
for o in "$build"/obj/*.o; do
    case ${o##*/} in vina.o|monte_carlo.o) objects+=("$out/${o##*/}");; split.o) ;; *) objects+=("$o");; esac
done
c++ "${flags[@]}" "${objects[@]}" -L"$boost/lib" -Wl,-rpath,"$boost/lib" -lboost_thread -lboost_filesystem -lboost_program_options -lboost_serialization -o "$out/vina"
