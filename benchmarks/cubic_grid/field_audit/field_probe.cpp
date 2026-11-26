// Sample the same uncapped midpoint pair tables used for native Vina map construction.
// Read the native XS typing from an intact prepared ligand and a rigid receptor.
#include "curl.h"
#include "grid.h"
#include "parse_pdbqt.h"
#include "precalculate.h"

#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    const auto receptor = parse_receptor_pdbqt(argv[1]);
    const auto ligand = parse_ligand_pdbqt_from_file(argv[2], atom_type::XS);
    const auto by_serial = ligand_heavy_indices_by_serial(ligand);
    const flv weights{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439,
                      50., 5 * 0.05846 / 0.1 - 1};
    ScoringFunction sf(SF_VINA, weights, 8.0);
    precalculate pairs(8.0);
    std::bitset<XS_TYPE_SIZE> receptor_types, ligand_types;
    for (const auto& a : receptor.grid_atoms())
        if (a.xs < XS_TYPE_SIZE) receptor_types.set(a.xs);
    for (const auto& entry : by_serial) ligand_types.set(ligand.movable_atom(entry.second).xs);
    pairs.populate(sf, receptor_types, ligand_types);
    std::cout << std::setprecision(17);
    vec origin;
    double h;
    if (!(std::cin >> origin[0] >> origin[1] >> origin[2] >> h) || h <= 0) return 2;
    unsigned serial;
    vec point;
    while (std::cin >> serial >> point[0] >> point[1] >> point[2]) {
        const auto& probe = ligand.movable_atom(by_serial.at(serial));
        grid_dims stencil;
        for (sz axis = 0; axis < 3; ++axis) {
            auto& dim = stencil[axis];
            dim.n_voxels = 3;
            dim.begin = origin[axis] + (std::floor((point[axis]-origin[axis])/h) - 1)*h;
            dim.end = dim.begin + 3*h;
        }
        grid patch(stencil);
        double explicit_fast = 0, explicit_smooth = 0;
        vec derivative(0, 0, 0);
        int near_cutoff = 0;
        for (const auto& a : receptor.grid_atoms()) {
            if (a.xs >= XS_TYPE_SIZE) continue;
            const vec d = point - a.coords;
            const fl d2 = sqr(d);
            // Every node lies less than 2 sqrt(3) grid spacings from the query.
            if (d2 > sqr(8.0 + 2*std::sqrt(3.0)*h)) continue;
            const sz index = get_type_pair_index(atom_type::XS, probe, a);
            if (d2 < 64) {
                explicit_fast += pairs.eval_fast(index, d2);
                const auto e_g = pairs.eval_deriv(index, d2);
                explicit_smooth += e_g.first;
                derivative += e_g.second * d;
            }
            if (std::abs(std::sqrt(d2) - 8.0) < 0.01) ++near_cutoff;
            for (sz x = 0; x < 4; ++x)
                for (sz y = 0; y < 4; ++y)
                    for (sz z = 0; z < 4; ++z) {
                        const fl r2 = sqr(patch.index_to_argument(x,y,z) - a.coords);
                        if (r2 <= 64) patch.m_data(x,y,z) += pairs.eval_fast(index,r2);
                    }
        }
        vec grid_gradient(0, 0, 0);
        const fl native_linear = patch.evaluate(point, 1e6, 10., grid_gradient);
        vec cubic_gradient(0, 0, 0);
        const fl native_cubic = patch.evaluate_cubic(point, 1e6, 10., &cubic_gradient);
        curl(explicit_fast, 10.);
        curl(explicit_smooth, derivative, 10.);
        std::cout << probe.xs << ' ' << near_cutoff << ' ' << explicit_fast << ' ' << explicit_smooth
                  << ' ' << derivative[0] << ' ' << derivative[1] << ' ' << derivative[2]
                  << ' ' << native_linear << ' ' << grid_gradient[0] << ' ' << grid_gradient[1] << ' ' << grid_gradient[2]
                  << ' ' << native_cubic << ' ' << cubic_gradient[0] << ' ' << cubic_gradient[1] << ' ' << cubic_gradient[2];
        for (sz x = 0; x < 4; ++x)
            for (sz y = 0; y < 4; ++y)
                for (sz z = 0; z < 4; ++z) std::cout << ' ' << patch.m_data(x,y,z);
        std::cout << '\n';
    }
    return std::cin.eof() ? 0 : 2;
}
