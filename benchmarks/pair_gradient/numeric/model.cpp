// Compare native model::eval_deriv with central differences of its own returned
// energy on a supplied experimental complex. Lengths are Angstrom, angles radians.
#include "cache.h"
#include "non_cache.h"
#include "parse_pdbqt.h"

#include <iomanip>
#include <iostream>

namespace {
void audit(model& m, const precalculate_byatom& pairs, const igrid& receptor,
           const conf& c, const vec& cap, const char* stage) {
    flv gradient(m.get_size().num_degrees_of_freedom()), scratch(gradient.size());
    m.set(c);
    const double energy = m.eval_deriv(pairs, receptor, cap, gradient);
    for (sz i = 0; i < gradient.size(); ++i) {
        flv basis(gradient.size(), 0);
        basis[i] = 1;
        conf displaced = c;
        const auto value = [&](double offset) {
            displaced.set_incremented(c, basis, offset);
            m.set(displaced);
            return m.eval_deriv(pairs, receptor, cap, scratch);
        };
        for (double h : {1e-5, 1e-6, 1e-7}) {
            const double plus = value(h), minus = value(-h);
            const double difference = (plus - minus) / (2 * h);
            std::cout << stage << ',' << i << ',' << h << ',' << energy << ','
                      << gradient[i] << ',' << difference << '\n';
        }
    }
    m.set(c);
}
// Trace the largest whole-model residual back to physical atom pairs. Rebuild
// just four required nodes from the weighted potential rather than looking up
// any stored energies or forces. Rigid rotation preserves all internal terms.
void trace_y_torque(const model& m, const conf& c, const ScoringFunction& sf) {
    const vec y_axis(0, 1, 0);
    std::cout << std::setprecision(16)
              << "ligand_i,receptor_j,xs1,xs2,r,ligand_x,ligand_y,ligand_z,rec_x,rec_y,rec_z,"
                 "native_torque,linear_torque\n";
    struct pair_torque { sz receptor; fl r, native, linear; };
    VINA_FOR(i, m.num_movable_atoms()) {
        const auto& a = m.movable_atom(i);
        if (a.xs >= XS_TYPE_SIZE) continue;
        const auto& x = m.movable_coords(i);
        const vec displacement = cross_product(y_axis, x - c.ligands.front().rigid.position);
        fl energy = 0;
        std::vector<pair_torque> terms;
        const auto& receptor = m.grid_atoms();
        VINA_FOR(j, receptor.size()) {
            const auto& b = receptor[j];
            const vec d = x - b.coords;
            const fl r2 = sqr(d);
            if (r2 >= 64 || b.xs >= XS_TYPE_SIZE) continue;
            const fl q = 32 * r2;
            const int k = static_cast<int>(q);
            const fl t = q - k;
            const auto rnode = [](int n) { return std::sqrt(n / 32.0); };
            const auto e = [&](int n) { return sf.eval(a.xs, b.xs, rnode(n)); };
            const auto slope = [&](int n) {
                return n == 0 ? 0 : (e(n+1)-e(n-1))/((rnode(n+1)-rnode(n-1))*rnode(n));
            };
            const fl delta = e(k+1) - e(k);
            energy += e(k) + t * delta;
            const fl lever = displacement * d;
            terms.push_back({j, std::sqrt(r2), lever*(slope(k) + t*(slope(k+1)-slope(k))),
                              lever*64*delta});
        }
        const fl cap_derivative = energy > 0 ? sqr(1000/(1000+energy)) : 1;
        for (const auto& term : terms) {
            const auto& b = receptor[term.receptor];
            std::cout << i << ',' << term.receptor << ',' << a.xs << ',' << b.xs << ',' << term.r;
            for (sz axis = 0; axis < 3; ++axis) std::cout << ',' << x[axis];
            for (sz axis = 0; axis < 3; ++axis) std::cout << ',' << b.coords[axis];
            std::cout << ',' << term.native * cap_derivative << ','
                      << term.linear * cap_derivative << '\n';
        }
    }
}
} // namespace

int main(int argc, char** argv) {
    const bool trace = argc == 4 && std::string(argv[3]) == "--trace-y-torque";
    if (argc != 3 && !trace) return 2;
    model m = parse_receptor_pdbqt(argv[1]);
    m.append(parse_ligand_pdbqt_from_file(argv[2], atom_type::XS));
    const auto c = m.get_initial_conf();
    m.set(c);
    const vec center = m.ligand_heavy_centroid(0);
    grid_dims box;
    VINA_FOR_IN(i, box) {
        // Fixed off-grid center avoids forcing a particular atom onto a node.
        box[i].begin = center[i] - 12.031;
        box[i].end = box[i].begin + 24;
        box[i].n_voxels = 64;
    }
    const auto sf = std::make_shared<ScoringFunction>(SF_VINA,
        flv{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50, 5*0.05846/0.1-1});
    if (trace) {
        trace_y_torque(m, c, *sf);
        return 0;
    }
    const precalculate_byatom pairs(*sf, m);
    non_cache explicit_receptor(m, box, sf, 1e6);
    explicit_receptor.prepare_for_model(m);
    cache maps(box);
    explicit_receptor.populate_maps(maps, m, cache::map_types(m.get_movable_atom_types(atom_type::XS)), true);
    std::cout << std::setprecision(16)
              << "stage,dof,h,energy,returned,central_difference\n";
    audit(m, pairs, maps, c, vec(10, 10, 10), "hunt");
    audit(m, pairs, explicit_receptor, c, vec(1000, 1000, 1000), "explicit");
}
