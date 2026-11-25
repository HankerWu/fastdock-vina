#include <iomanip>
#include <iostream>

#include "non_cache.h"
#include "parse_pdbqt.h"

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    model m = parse_receptor_pdbqt(argv[1]);
    m.append(parse_ligand_pdbqt_from_file(argv[2], atom_type::XS));
    if (m.get_glue_pairs().empty()) return 2;

    grid_dims bounds;
    const vec center(15.19, 53.903, 16.917);
    VINA_FOR_IN(i, bounds) {
        bounds[i].begin = center[i] - 10;
        bounds[i].end = center[i] + 10;
        bounds[i].n_voxels = 54;
    }

    const auto sf = std::make_shared<ScoringFunction>(
        SF_VINA, flv{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50, 0.05846});
    non_cache receptor(m, bounds, sf, 100);
    receptor.prepare_for_model(m);
    const precalculate_byatom ligand_tables(*sf, m);
    const auto size = m.get_size();
    conf c = m.get_initial_conf();
    c.ligands.front().torsions.back() += 0.5;
    m.set(c);

    flv first(size.num_degrees_of_freedom()), second(first.size());
    const vec cap(1000, 1000, 1000);
    const fl first_energy = m.eval_deriv(ligand_tables, receptor, cap, first);
    const fl second_energy = m.eval_deriv(ligand_tables, receptor, cap, second);
    if (first_energy != second_energy) return 2;
    for (sz i = 0; i < first.size(); ++i) {
        if (first[i] != second[i]) {
            std::cerr << std::setprecision(17) << "Unchanged pose: gradient " << i << " changed from "
                      << first[i] << " to " << second[i] << '\n';
            return 1;
        }
    }
    std::cout << "Explicit receptor gradient is repeatable\n";
}
