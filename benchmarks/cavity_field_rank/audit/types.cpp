// Export topology-derived types only. The independent calculation reads raw coordinates.
#include "parse_pdbqt.h"
#include "precalculate.h"
#include "curl.h"
#include <iomanip>
#include <iostream>

// Independently check the lattice-free Python result with FastDock's native
// distance table. Input rows carry already clipped heavy positions and the
// separately computed wall penalty. All receiver atoms are enumerated here.
void check_pairs(const model& receptor, const model& ligand) {
    const flv weights = {-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50., 1.923};
    ScoringFunction sf(SF_VINA, weights);
    precalculate table(sf.get_cutoff());
    std::bitset<XS_TYPE_SIZE> fixed, mobile;
    for (const auto& a : receptor.grid_atoms())
        if (a.xs < XS_TYPE_SIZE) fixed.set(a.xs);
    for (auto t : ligand.get_movable_atom_types(atom_type::XS))
        if (t < XS_TYPE_SIZE) mobile.set(t);
    table.populate(sf, fixed, mobile);
    size_t count;
    double total;
    while (std::cin >> count >> total) {
        for (size_t i = 0; i < count; ++i) {
            size_t t;
            vec x;
            if (!(std::cin >> t >> x[0] >> x[1] >> x[2]) || !mobile.test(t))
                throw std::runtime_error("Invalid direct-space sample");
            double local = 0;
            for (const auto& a : receptor.grid_atoms()) {
                if (a.xs >= XS_TYPE_SIZE) continue;
                const fl r2 = vec_distance_sqr(a.coords, x);
                if (r2 < table.cutoff_sqr())
                    local += table.eval_fast(triangular_matrix_index_permissive(XS_TYPE_SIZE, t, a.xs), r2);
            }
            curl(local, 10.);
            total += local;
        }
        std::cout << total << "\n";
    }
}

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) throw std::runtime_error("types receptor.pdbqt ligand.pdbqt ligand.bonds");
    const model receptor = parse_receptor_pdbqt(argv[1]);
    const model ligand = parse_ligand_pdbqt_with_bonds(argv[2], atom_type::XS, argv[3]);
    std::cout << std::setprecision(17);
    if (argc == 5) {
        if (std::string(argv[4]) != "--pairs") throw std::runtime_error("unknown option");
        check_pairs(receptor, ligand);
        return 0;
    }
    for (const auto& a : receptor.grid_atoms()) {
        if (a.xs < XS_TYPE_SIZE)
            std::cout << "R " << a.xs << " " << a.coords[0] << " " << a.coords[1] << " " << a.coords[2] << "\n";
    }
    for (auto entry : ligand_heavy_indices_by_serial(ligand))
        std::cout << "L " << entry.first << " " << ligand.get_atom(entry.second).xs << "\n";
}
