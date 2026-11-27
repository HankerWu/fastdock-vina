// Link against the prepared-angle ablation, including its model.o. Cartesian
// single-angle perturbations provide an independent check of branch indexing.
#include "mutate.h"
#include "parse_pdbqt.h"

#include <algorithm>
#include <iomanip>
#include <iostream>

namespace {
flv closure_distances(model& m, const conf& c, const interacting_pairs& pairs) {
    m.set(c);
    flv distances;
    for (const auto& pair : pairs)
        distances.push_back(std::sqrt(vec_distance_sqr(m.get_coords(pair.a), m.get_coords(pair.b))));
    return distances;
}

fl maximum_difference(const flv& a, const flv& b) {
    fl result = 0;
    for (sz i = 0; i < a.size(); ++i) result = std::max(result, std::abs(a[i] - b[i]));
    return result;
}

fl largest(const flv& values) {
    return values.empty() ? 0 : *std::max_element(values.begin(), values.end());
}
}

int main(int argc, char** argv) {
    if (argc < 3) return 2;
    model input = parse_ligand_pdbqt_from_file(argv[2], atom_type::XS);
    for (int i = 3; i < argc; ++i)
        input.append(parse_ligand_pdbqt_from_file(argv[i], atom_type::XS));
    const auto pairs = input.get_glue_pairs();
    const conf prepared = input.get_initial_conf();
    const flv reference = closure_distances(input, prepared, pairs);
    for (unsigned seed : {191009u, 458921u, 740317u, 384529u, 834721u}) {
        rng stream(seed);
        rng generator(static_cast<rng::result_type>(random_int(0, 1000000, stream)));
        model m = input;
        conf first = prepared;
        first.randomize(vec(-12,-12,-12), vec(12,12,12), generator);
        mutate_conf(first, m, 2, generator, false);
        const auto initial = closure_distances(m, first, pairs);
        conf reset = first;
        m.restore_closure_path_torsions(reset);
        const auto corrected = closure_distances(m, reset, pairs);
        sz n_reset = 0, n_retained = 0, n_cartesian = 0;
        fl max_unmarked_change = 0;
        for (sz i = 0; i < first.ligands.size(); ++i) {
            const auto& a = first.ligands[i];
            const auto& b = reset.ligands[i];
            if (a.rigid.orientation != b.rigid.orientation ||
                vec_distance_sqr(a.rigid.position,b.rigid.position) != 0) return 3;
            for (sz j = 0; j < a.torsions.size(); ++j) {
                const bool changed = a.torsions[j] != b.torsions[j];
                if (changed && b.torsions[j] != 0) return 4;
                n_reset += changed;
                n_retained += !changed;
                conf trial = first;
                trial.ligands[i].torsions[j] += 0.731;
                const fl variation = maximum_difference(initial, closure_distances(m, trial, pairs));
                if (!changed) max_unmarked_change = std::max(max_unmarked_change, variation);
                n_cartesian += variation > 1e-7;
            }
        }
        const fl residual = maximum_difference(reference, corrected);
        std::cout << std::setprecision(12) << argv[1] << ',' << seed << ',' << pairs.size()
                  << ',' << n_reset << ',' << n_retained << ',' << n_cartesian
                  << ',' << largest(reference) << ',' << largest(initial) << ',' << largest(corrected)
                  << ',' << residual << ',' << max_unmarked_change << '\n';
        if (residual > 1e-7 || max_unmarked_change > 1e-7) return 5;
    }
}
