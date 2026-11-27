// Compile against the archived placement core, omitting monte_carlo.o.
#include "monte_carlo.cpp"
#include "parse_pdbqt.h"
#include <iostream>
#include <stdexcept>

static void check(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    model m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
    const vec edge(12,12,12), low(-12,-12,-12);
    rng gen(829153);
    unsigned outside = 0, inside = 0;
    for (unsigned i=0; i<256; ++i) {
        conf source=m.get_initial_conf();
        source.randomize(low,edge,gen);
        source.ligands.front().rigid.position[0] += i%3 == 0 ? 24 : 0;
        m.set(source);
        const vecv atoms=m.get_chemical_heavy_atom_movable_coords();
        vec a=edge, b=low;
        bool external=false;
        for (const vec& x : atoms) for (sz k=0; k<3; ++k) {
            // Offsets rather than world extrema, independently of the implementation.
            const fl offset=x[k]-source.ligands.front().rigid.position[k];
            a[k]=std::min(a[k],offset);
            b[k]=std::max(b[k],offset);
            external |= x[k]<low[k] || x[k]>edge[k];
        }
        external ? ++outside : ++inside;
        for (auto mode : {initial_placement::fit, initial_placement::invalid_only}) {
            conf posed=source;
            fit_first_proposal(m,posed,low,edge,mode);
            const auto& rigid=posed.ligands.front().rigid;
            check(posed.ligands.front().torsions == source.ligands.front().torsions &&
                  rigid.orientation == source.ligands.front().rigid.orientation,"changed internal or rotational DOFs");
            for (sz k=0; k<3; ++k) {
                const fl old=source.ligands.front().rigid.position[k];
                const fl u=(old+12)/24;
                const fl expected=mode==initial_placement::invalid_only && !external ? old :
                    (-12-a[k])*(1-(u-std::floor(u))) + (12-b[k])*(u-std::floor(u));
                check(std::abs(rigid.position[k]-expected)<1e-11,"conditional quantile mismatch");
            }
            m.set(posed);
            for (const vec& x : m.get_chemical_heavy_atom_movable_coords())
                for(sz k=0;k<3;++k) check(std::abs(x[k])<=12+1e-11,"placed outside");
            conf oversized=source;
            fit_first_proposal(m,oversized,vec(0,0,0),vec(0.01,24,24),mode);
            check(vec_distance_sqr(oversized.ligands.front().rigid.position,source.ligands.front().rigid.position)==0,
                  "changed part of an infeasible proposal");
        }
    }
    check(inside && outside,"missing branch coverage");
    std::cout<<"passed: "<<inside<<" inside, "<<outside<<" outside, both modes and oversized fallback\n";
}
