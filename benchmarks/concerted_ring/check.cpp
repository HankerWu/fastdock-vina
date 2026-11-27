// One-step, common-random-number tests plus a central difference against
// independent full forward kinematics on real PDBQT ligand trees.
#include "concerted_ring.h"
#include "fixed_ring.h"
#include "mutate.h"
#include "parse_pdbqt.h"
#include <iomanip>
#include <iostream>

static fl gap(const model& m, const concerted_ring::ends& ends) {
    return std::max((m.get_coords(ends[0])-m.get_coords(ends[1])).norm(),
                    (m.get_coords(ends[2])-m.get_coords(ends[3])).norm());
}

static fl finite_difference(conf& c, model& m, const concerted_ring::ends& ends,
                            const std::vector<concerted_ring::axis>& axes) {
    using namespace concerted_ring;
    m.set(c);
    std::vector<six> jac(axes.size());
    geometry(m, ends, axes, jac);
    fl largest = 0;
    const fl h = 1e-6;
    for (sz k = 0; k < axes.size(); ++k) {
        fl& torsion = c.ligands[0].torsions[axes[k].torsion];
        const fl saved = torsion;
        torsion = saved + h;
        m.set(c);
        const vec a = m.get_coords(ends[0])-m.get_coords(ends[1]);
        const vec b = m.get_coords(ends[2])-m.get_coords(ends[3]);
        torsion = saved - h;
        m.set(c);
        const vec diff[] = {a - (m.get_coords(ends[0])-m.get_coords(ends[1])),
                            b - (m.get_coords(ends[2])-m.get_coords(ends[3]))};
        for (sz p = 0; p < 2; ++p)
            for (sz d = 0; d < 3; ++d)
                largest = std::max(largest, std::abs(diff[p][d]/(2*h)-jac[k][3*p+d]));
        torsion = saved;
    }
    m.set(c);
    return largest;
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    std::cout << "pdb,seed,trial,driver,n_ring_axes,draw_rad,applied,max_torsion_step_rad,max_endpoint_A,max_jacobian_error_A_per_rad\n" << std::setprecision(12);
    for (int f = 1; f < argc; ++f) {
        const path filename(argv[f]);
        model m = parse_ligand_pdbqt_from_file(filename.string(), atom_type::XS);
        const fixed_ring ring(m);
        const auto size = m.get_size();
        for (const unsigned seed : {191009u, 458921u, 740317u}) {
            rng init(seed), generator(seed+17);
            conf c(size);
            c.randomize(vec(-5,-5,-5),vec(5,5,5),init);
            ring.prepare(c);
            m.set(c);
            concerted_ring::ends ends;
            std::vector<concerted_ring::axis> axes;
            if (!concerted_ring::topology(m,0,ends,axes)) return 1;
            for (unsigned trial = 0; trial < 128; ++trial) {
                rng cloned = generator;
                const sz which = random_int(0, int(1+size.ligands[0]), cloned);
                const bool is_ring = !std::binary_search(ring.allowed_moves().begin(),ring.allowed_moves().end(),which);
                if (!is_ring) {
                    mutate_conf(c,m,2,generator,false,false,false,false,&ring.allowed_moves());
                    m.set(c);
                    continue;
                }
                const fl draw = random_fl(-pi, pi, cloned);
                const sz driver = which - 2;
                const fl error = finite_difference(c,m,ends,axes);
                const conf before = c;
                conf direct = c;
                const bool applied = concerted_ring::propose(direct,m,0,driver,draw,2);
                m.set(c);
                mutate_conf(c,m,2,generator,false,false,false,false,&ring.allowed_moves());
                if (generator != cloned || c.ligands[0].torsions != direct.ligands[0].torsions ||
                    c.ligands[0].rigid.orientation != before.ligands[0].rigid.orientation ||
                    vec_distance_sqr(c.ligands[0].rigid.position,before.ligands[0].rigid.position) != 0)
                    return 1;
                m.set(c);
                fl max_step = 0;
                for (sz i = 0; i < size.ligands[0]; ++i)
                    max_step = std::max(max_step,std::abs(c.ligands[0].torsions[i]-before.ligands[0].torsions[i]));
                const fl residual = gap(m,ends);
                if (error > 1e-6 || residual > 1e-6 || (!applied && max_step != 0)) return 1;
                std::cout << filename.parent_path().filename().string() << ',' << seed << ',' << trial << ',' << driver << ','
                          << axes.size() << ',' << draw << ',' << applied << ',' << max_step << ',' << residual << ',' << error << '\n';
            }
        }
    }
}
