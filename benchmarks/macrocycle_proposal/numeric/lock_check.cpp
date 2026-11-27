// Compare core MC and BFGS updates to the fixed internal-coordinate subspace.
#include "fixed_ring.h"
#include "mutate.h"
#include "non_cache.h"
#include "parse_pdbqt.h"
#include "quasi_newton.h"

#include <iomanip>
#include <iostream>

namespace {
// Smooth periodic test potential with known derivatives. Translation is
// linear so this also checks offsets across multiple ligand rigid blocks.
fl smooth_energy(const conf& c, flv* gradient = nullptr) {
    sz k=0;
    fl result=0;
    for (const auto& lig : c.ligands) {
        for (sz i=0;i<3;++i) {
            result+=(i+1)*lig.rigid.position[i];
            if (gradient) (*gradient)[k+i]=i+1;
        }
        k+=6;
        for (fl angle : lig.torsions) {
            result+=std::sin(angle)+0.25*std::cos(2*angle);
            if (gradient) (*gradient)[k]=std::cos(angle)-0.5*std::sin(2*angle);
            ++k;
        }
    }
    return result;
}

fl check_derivative(const conf& c, const fixed_ring& ring, const conf_size& size) {
    flv g(size.num_degrees_of_freedom()), d(g.size());
    smooth_energy(c,&g);
    ring.project(g);
    for (sz i=0;i<d.size();++i) d[i]=std::cos(0.83*(i+1));
    const fl analytical=scalar_product(g,d);
    ring.project(d);
    const fl h=1e-5;
    conf plus(size),minus(size);
    plus.set_incremented(c,d,h);
    minus.set_incremented(c,d,-h);
    return std::abs(analytical-(smooth_energy(plus)-smooth_energy(minus))/(2*h));
}

fl max_gap(const model& m, const interacting_pairs& pairs) {
    fl gap = 0;
    for (const auto& pair : pairs)
        gap = std::max(gap, std::sqrt(vec_distance_sqr(m.get_coords(pair.a),m.get_coords(pair.b))));
    return gap;
}
}

int main(int argc, char** argv) {
    if (argc < 3) return 2;
    model m = parse_ligand_pdbqt_from_file(argv[2], atom_type::XS);
    for (int i = 3; i < argc; ++i) m.append(parse_ligand_pdbqt_from_file(argv[i], atom_type::XS));
    const auto pairs = m.get_glue_pairs();
    const fixed_ring ring(m);
    const auto allowed = ring.active() ? &ring.allowed_moves() : nullptr;
    const auto size = m.get_size();
    const auto sf = std::make_shared<ScoringFunction>(SF_VINA,
        flv{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50, 0.05846});
    const precalculate_byatom table(*sf, m);
    grid_dims bounds;
    for (auto& axis : bounds) { axis.begin=-16; axis.end=16; axis.n_voxels=64; }
    // No protein is required to stress these coordinate constraints: force
    // comes from the actual intramolecular tables and box, with native BFGS.
    non_cache grid(m,bounds,sf,100);
    grid.prepare_for_model(m);
    flv diagonal(size.num_degrees_of_freedom(),1);
    ring.project(diagonal);
    const sz locks = std::count(diagonal.begin(),diagonal.end(),0);
    for (unsigned seed : {191009u,458921u,740317u,384529u,834721u}) {
        rng generator(seed);
        output_type pose(size,0);
        pose.c.randomize(vec(-3,-3,-3),vec(3,3,3),generator);
        ring.prepare(pose.c);
        quasi_newton minimize(size,16);
        minimize.ring = ring.active() ? &ring : nullptr;
        const int steps=128;
        fl mc_gap=0, bfgs_gap=0, fd_error=0;
        int evals=0, null_moves=0;
        for (int move=0; move<steps; ++move) {
            m.set(pose.c);
            const conf previous=pose.c;
            mutate_conf(pose.c,m,2,generator,false,false,false,false,allowed);
            bool changed=false;
            for (sz i=0;i<previous.ligands.size();++i) {
                const auto& a=previous.ligands[i]; const auto& b=pose.c.ligands[i];
                changed |= a.rigid.orientation != b.rigid.orientation ||
                           vec_distance_sqr(a.rigid.position,b.rigid.position)>0 || a.torsions != b.torsions;
            }
            null_moves += !changed;
            m.set(pose.c);
            mc_gap=std::max(mc_gap,max_gap(m,pairs));
            if ((move+1)%16==0) {
                fd_error=std::max(fd_error,check_derivative(pose.c,ring,size));
                minimize(m,table,grid,pose,vec(10,10,10),evals);
                m.set(pose.c);
                bfgs_gap=std::max(bfgs_gap,max_gap(m,pairs));
            }
            sz k=0;
            for (const auto& ligand : pose.c.ligands) {
                k+=6;
                for (fl angle : ligand.torsions) {
                    if (diagonal[k++] == 0 && angle != 0) return 1;
                }
            }
        }
        const sz moves = allowed ? allowed->size() : 2*size.ligands.size()+sum(size.ligands)+sum(size.flex);
        std::cout << std::setprecision(12) << argv[1] << ',' << seed << ',' << locks << ',' << moves << ','
                  << steps << ',' << null_moves << ',' << evals << ',' << mc_gap << ',' << bfgs_gap << ',' << fd_error << '\n';
        if (mc_gap>1e-7 || bfgs_gap>1e-7 || fd_error>1e-7 || null_moves) return 1;
    }
}
