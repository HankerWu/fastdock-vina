// Compare the score descriptor with longest finite internal-distance change.
// Different ligands are intentionally permitted to share an atom serial.
#include "parse_pdbqt.h"
#include "conf_independent.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    auto m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
    for (int i=2; i<argc; ++i)
        m.append(parse_ligand_pdbqt_from_file(argv[i], atom_type::XS));
    const auto ref = m.get_initial_conf();
    m.set(ref);
    const auto initial = m.get_chemical_heavy_atom_movable_coords();
    std::cout << std::setprecision(12) << conf_independent_inputs(m).num_tors << " "
              << m.num_ligands();
    for (size_t lig=0; lig<m.num_ligands(); ++lig) {
        // Pair-distance observability is an intraligand property.
        if (m.num_ligands()>1) break;
        for (size_t k=0; k<ref.ligands[lig].torsions.size(); ++k) {
            auto updated=ref;
            updated.ligands[lig].torsions[k] += 1.234;
            m.set(updated);
            const auto end = m.get_chemical_heavy_atom_movable_coords();
            double max_change=0;
            for (size_t i=0; i<end.size(); ++i)
                for (size_t j=0; j<i; ++j)
                    max_change=std::max(max_change,std::abs(
                        std::sqrt(vec_distance_sqr(initial[i],initial[j])) -
                        std::sqrt(vec_distance_sqr(end[i],end[j]))));
            std::cout << " " << max_change;
        }
    }
    std::cout << "\n";
}
