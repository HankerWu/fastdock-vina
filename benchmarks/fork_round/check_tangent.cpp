// Check the whole coupled Jp, including translation, lab-frame rotation and branched torsions.
#include "parse_pdbqt.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    fl largest_error = 0;
    unsigned checked = 0;
    for (int file = 1; file < argc; ++file) {
        model m = parse_ligand_pdbqt_from_file(argv[file], atom_type::XS);
        const conf initial = m.get_initial_conf();
        for (unsigned pose = 0; pose < 6; ++pose) {
            conf base = initial;
            auto& ligand = base.ligands.front();
            quaternion_increment(ligand.rigid.orientation, vec(0.21*pose, -0.34*pose, 0.12*pose));
            for (sz i = 0; i < ligand.torsions.size(); ++i)
                ligand.torsions[i] = 2.7*std::sin((i+1.)*pose);
            flv direction(m.get_size().num_degrees_of_freedom());
            for (sz i = 0; i < direction.size(); ++i)
                direction[i] = std::cos((pose+1.)*(i+1.));
            const fl h = 1e-5;
            conf trial = base;
            trial.set_incremented(base, direction, h);
            m.set(trial);
            const vecv plus = m.get_chemical_heavy_atom_movable_coords();
            trial.set_incremented(base, direction, -h);
            m.set(trial);
            const vecv minus = m.get_chemical_heavy_atom_movable_coords();
            fl fd = 0;
            for (sz i = 0; i < plus.size(); ++i)
                fd = std::max(fd, std::sqrt(vec_distance_sqr(plus[i], minus[i]))/(2*h));
            const fl analytic = m.max_heavy_atom_speed(base, direction);
            const fl error = std::abs(fd - analytic);
            if (!std::isfinite(error) || error > 1e-6)
                throw std::runtime_error(std::string("Coupled Cartesian speed mismatch: ") + argv[file]);
            largest_error = std::max(largest_error, error);
            ++checked;
        }
    }
    std::cout << std::setprecision(17) << "Ligands " << argc - 1 << ", poses " << checked
              << ", max abs speed error " << largest_error << '\n';
}
