#!/usr/bin/env bash
# Usage: bash check_geometry.sh ccd|reset|lock patched-source-root build-dir macro-input-dir output.csv [additional-prepared-case-directory ...]
# macro-input-dir contains the two expanded archives documented in bond_distance/README.md.
set -euo pipefail
arm=$1; source_root=$2; build_dir=$3; inputs=$4; output=$5; shift 5
test "$arm" = ccd || test "$arm" = reset || test "$arm" = lock
here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
boost_args=()
if test -n "${BOOST_PREFIX:-}"; then
    boost_args+=("-I$BOOST_PREFIX/include" "-L$BOOST_PREFIX/lib" "-Wl,-rpath,$BOOST_PREFIX/lib")
fi
boost_version=$("${CXX:-c++}" -x c++ "${boost_args[@]}" -dM -E -include boost/version.hpp /dev/null | awk '$2 == "BOOST_VERSION" { print $3 }')
boost_system=()
if test "$boost_version" -lt 108900; then boost_system+=(-lboost_system); fi
objs=()
for obj in "$build_dir"/obj/*.o; do
    if test "$(basename "$obj")" != main.o; then objs+=("$obj"); fi
done
"${CXX:-c++}" -std=c++14 -O3 -DNDEBUG -pthread "-I$source_root/src/lib" "${boost_args[@]}" \
    "$here/${arm}_check.cpp" "${objs[@]}" "${boost_system[@]}" -lboost_thread -lboost_serialization \
    -lboost_filesystem -lboost_program_options -o "$work/check"
if test "$arm" = lock; then
    echo 'structure,seed,locked_angles,allowed_moves,mc_proposals,null_proposals,bfgs_evals,max_gap_after_move_A,max_gap_after_bfgs_A,max_projected_directional_derivative_error' > "$work/results.csv"
elif test "$arm" = ccd; then
    echo 'structure,seed,glue_pairs,identifiable_angles,flat_angles,initial_sum_squared_A2,projected_sum_squared_A2,reference_sum_squared_A2,max_angle_error_rad,max_atom_error_A,max_rigid_equivariance_error_A' > "$work/results.csv"
else
    echo 'structure,seed,glue_pairs,reset_angles,untouched_angles,angles_affecting_distance,max_reference_gap,max_initial_gap,max_corrected_gap,max_gap_difference_from_prepared,max_change_from_untouched_rotor' > "$work/results.csv"
fi
{
    "$work/check" 5QCT "$inputs/tutorial/solution/BACE_1_ligand.pdbqt"
    "$work/check" 1YET "$inputs/hsp/1YET/ligand.pdbqt"
    "$work/check" 1FKB "$inputs/extra/1FKB/ligand.pdbqt"
    "$work/check" 1FKF "$inputs/extra/1FKF/ligand.pdbqt"
    "$work/check" cyclododecane "$source_root/benchmarks/fixtures/cyclododecane.pdbqt"
    "$work/check" no_closure_control "$source_root/benchmarks/bond_distance/baricitinib.pdbqt"
    "$work/check" mixed_same_tag "$inputs/extra/1FKB/ligand.pdbqt" "$inputs/hsp/1YET/ligand.pdbqt" \
        "$source_root/benchmarks/bond_distance/baricitinib.pdbqt"
    for prepared in "$@"; do
        "$work/check" "$(basename "$prepared")" "$prepared/ligand.pdbqt"
    done
} >> "$work/results.csv"
cp "$work/results.csv" "$output"
