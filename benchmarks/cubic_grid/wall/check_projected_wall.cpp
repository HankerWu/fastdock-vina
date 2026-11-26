// Exact normal and tangential chain-rule checks against the native explicit evaluator.
// Run with benchmarks/fixtures/flex-rigid.pdbqt (a single C at a nonzero origin).
#include "non_cache.h"
#include "parse_pdbqt.h"

#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    model pair = parse_receptor_pdbqt(argv[1]);
    if (pair.grid_atoms().size() != 1) return 2;
    const vec center = pair.grid_atoms()[0].coords;
    pair.append(parse_ligand_pdbqt_from_string(
        "ROOT\nATOM      1    C UNL     1       0.000   0.000   0.000  1.00  0.00     0.000 C\nENDROOT\nTORSDOF 0\n",
        atom_type::XS));
    const auto sf = std::make_shared<ScoringFunction>(SF_VINA,
        flv{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50, 5*0.05846/0.1-1});
    const precalculate_byatom tables(*sf, pair);
    const auto bounds = [&](fl half) {
        grid_dims box;
        for (sz d = 0; d < 3; ++d) {
            box[d].begin = center[d] - half;
            box[d].end = center[d] + half;
            box[d].n_voxels = 8;
        }
        return box;
    };
    non_cache clipped(pair, bounds(2), sf, 0);
    non_cache unbounded(pair, bounds(8), sf, 0);
    clipped.prepare_for_model(pair);
    unbounded.prepare_for_model(pair);
    const vec cap(10, 10, 10);
    const conf initial = pair.get_initial_conf();
    auto evaluate = [&](const vec& pos, const non_cache& field, flv& gradient) {
        conf pose = initial;
        pose.ligands[0].rigid.position = center + pos;
        pair.set(pose);
        return pair.eval_deriv(tables, field, cap, gradient);
    };
    // No evaluation is on a non-differentiable box face.
    const vec offsets[] = {vec(2.4, 0.8, 0.6), vec(-2.7, 0.8, 0.6),
                           vec(2.6, -2.35, 0.7), vec(-2.6, 2.35, -2.2),
                           vec(1.7, 0.55, -0.9)};
    flv reported(pair.get_size().num_degrees_of_freedom()), reference(reported.size()), scratch(reported.size());
    unsigned failed = 0;
    std::cout << std::setprecision(12) << "slope,case,axis,gradient,projected_reference,finite_difference\n";
    for (fl slope : {0., 1.7, 100.}) {
        clipped.slope = slope;
        for (sz id = 0; id < sizeof(offsets)/sizeof(offsets[0]); ++id) {
            const vec offset = offsets[id];
            vec projection = offset;
            fl penalty = 0;
            for (sz d = 0; d < 3; ++d) {
                projection[d] = std::min(2., std::max(-2., offset[d]));
                penalty += slope * std::abs(offset[d]-projection[d]);
            }
            const fl base = evaluate(projection, unbounded, reference);
            const fl value = evaluate(offset, clipped, reported);
            if (std::abs(value - (base + penalty)) > 1e-10) ++failed;
            for (sz d = 0; d < 3; ++d) {
                const bool outside = offset[d] != projection[d];
                const fl expected = outside ? std::copysign(slope, offset[d]) : reference[d];
                const fl step = 1e-5;
                vec displacement(0,0,0);
                displacement[d] = step;
                const fl high = evaluate(offset+displacement, clipped, scratch);
                const fl low = evaluate(offset-displacement, clipped, scratch);
                const fl numeric = (high-low)/(2*step);
                std::cout << slope << ',' << id << ',' << d << ',' << reported[d] << ',' << expected << ',' << numeric << '\n';
                // Native radial energy and radial derivatives interpolate independently.
                // The smooth finite-difference reference is exact on a clamped normal.
                if (std::abs(expected - reported[d]) > 1e-10 ||
                    std::abs(numeric - reported[d]) > (outside ? 1e-7 : .02)) ++failed;
            }
        }
    }
    if (failed) std::cerr << "Failed projected-wall checks: " << failed << '\n';
    return failed ? 1 : 0;
}
