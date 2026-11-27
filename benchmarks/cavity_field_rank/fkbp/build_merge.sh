#!/usr/bin/env bash
# Relink the existing one-worker cavity_triplet observation build with an AD4 merge snapshot.
set -eu
source_tree=$1
original_objects=$2
boost_prefix=$3
work=$4
mkdir -p "$work/tmp"
work=$(cd "$work" && pwd)
export TMPDIR="$work/tmp"
python - "$source_tree/src/lib/vina.cpp" "$work/vina.cpp" <<'PY'
from pathlib import Path
import sys
source, target = map(Path, sys.argv[1:])
s = source.read_text()
anchor = 'parallelmc(m_model, poses, m_precalculated_byatom, m_ad4grid, m_ad4grid.corner1(), m_ad4grid.corner2(), generator, m_progress_callback, nullptr, false, counters);'
assert s.count(anchor) == 1
target.write_text(s.replace(anchor, anchor + '''
        if (std::getenv("FD_RIGID_COORDS")) {
            model snapshot(m_model);
            unsigned rank = 0;
            for (auto& pose : poses) {
                snapshot.set(pose.c);
                rigid_stage_xyz("merged_ad4", 0, ++rank, snapshot);
            }
        }'''))
PY
c++ -O3 -DNDEBUG -std=c++14 -pthread -fvisibility=hidden "-I$source_tree/src/lib" "-I$boost_prefix/include" -c "$work/vina.cpp" -o "$work/vina.o"
objects=()
for object in "$original_objects"/*.o; do
    if [[ $(basename "$object") != vina.o ]]; then objects+=("$object"); fi
done
c++ -O3 -std=c++14 -pthread "$work/vina.o" "${objects[@]}" "-L$boost_prefix/lib" "-Wl,-rpath,$boost_prefix/lib" -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options -o "$work/vina"
