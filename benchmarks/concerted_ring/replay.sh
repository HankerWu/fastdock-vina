#!/usr/bin/env bash
# Call from a full checkout. BOOST_PREFIX is optional; PYTHON must import RDKit, NumPy and SciPy.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work=$(realpath -m "$1")
py=${PYTHON:-python3}
jobs=${JOBS:-2}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$work/inputs" "$work/mer" "$work/new"
build=(--jobs "$jobs")
if [[ -n ${BOOST_PREFIX:-} ]]; then build+=(--boost-prefix "$BOOST_PREFIX"); fi
for arm in parent lock null concerted; do
    case "$arm" in
        parent) rev=2546cd5 ;;
        lock) rev=fda3793 ;;
        *) rev=fed1222 ;;
    esac
    git -C "$repo" worktree add --quiet --detach "$work/src-$arm" "$rev"
    if [[ $arm == concerted ]]; then
        git -C "$work/src-$arm" apply --unidiff-zero "$here/implementation.patch"
        cp "$here/concerted_ring.h" "$work/src-$arm/src/lib/"
    fi
    "$py" "$work/src-$arm/tools/build.py" --output-dir "$work/bin-$arm" "${build[@]}"
done
"$py" "$repo/benchmarks/pair_gradient/build_arms.py" "$work/up" --arms upstream upstream-counts "${build[@]}"
tar -xJf "$repo/benchmarks/macrocycle_proposal/external/inputs.tar.xz" -C "$work/inputs"
tar -xJf "$repo/benchmarks/macrocycle_proposal/holdout/observations.tar.xz" -C "$work/mer" inputs/5K0X inputs/raw/5K0X-bound.sdf
cp -r "$work/mer/inputs/5K0X" "$work/inputs/"
cp "$work/mer/inputs/raw/5K0X-bound.sdf" "$work/inputs/raw/"
tar -xJf "$here/new_inputs.tar.xz" -C "$work/new"
engines=("$work/up/upstream/build/vina" "$work/up/upstream-counts/build/vina" "$work/bin-parent/vina"
         "$work/bin-lock/vina" "$work/bin-null/vina" "$work/bin-concerted/vina")
bash "$here/check.sh" "$work/src-concerted" "$work/bin-concerted" "$work/inputs" "$work/check" > "$work/frozen.csv"
"$py" "$here/run.py" "$work/inputs" "$work/panel" "${engines[@]}"
"$py" "$here/run.py" "$work/new" "$work/transfer" "${engines[@]}" --ids 3VHA 8H6P
"$py" "$here/run.py" "$work/new" "$work/catd" "${engines[@]}" --ids 6QBH
"$py" "$here/run.py" "$work/new" "$work/lower" "${engines[@]}" --ids 6QBH 3VHA 8H6P --lower
"$py" "$here/run.py" "$work/new" "$work/fresh" "${engines[@]}" --ids 3VHA 8H6P 6QBH --fresh
obs=()
if [[ -n ${BOOST_PREFIX:-} ]]; then obs+=(--boost-prefix "$BOOST_PREFIX"); fi
for arm in null concerted; do
    "$py" "$here/stages.py" build "$work/src-$arm" "$work/bin-$arm/obj" "$work/obs-$arm" "${obs[@]}"
done
"$py" "$repo/benchmarks/macrocycle_proposal/observer.py" "$work/src-lock" "$work/bin-lock/obj" "$work/obs-lock" "${obs[@]}"
for pair in 'inputs panel stages' 'new transfer held_stages'; do
    read -r inputs panel stage <<< "$pair"
    "$py" "$here/stages.py" replay "$work/$inputs" "$work/$panel" "$work/$stage" "$work/obs-null/vina" "$work/obs-concerted/vina"
    "$py" "$here/diagnose.py" "$work/$inputs" "$work/$stage" "$work/$stage"
done
"$py" "$repo/benchmarks/macrocycle_proposal/external/external_stages.py" "$work/new" "$work/transfer" "$work/held_lock_stages" --lock "$work/obs-lock/vina"
