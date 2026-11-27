#include "search_grid.h"
#include "non_cache.h"
#include "precalculate.h"
#include <iostream>

// Compare lattice samples to the analytic potential, then differentiate a
// translated/rotated/torsional test ligand through the actual search interface.
int main(int argc, char** argv) {
    if (argc != 2) throw std::runtime_error("Expected prepared ligand PDBQT");
    auto m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
    const auto indices = ligand_heavy_indices_by_serial(m);
    const auto sf = std::make_shared<ScoringFunction>(SF_VINA,
        flv{-0.035579,-0.005156,0.840245,-0.035069,-0.587439,50,5*0.05846/0.1-1});
    grid_dims dims;
    for (auto& axis : dims) { axis.begin = -12; axis.end = 12; axis.n_voxels = 64; }
    non_cache base(m, dims, sf, 1e6);
    base.prepare_for_model(m);
    std::stringstream spec;
    spec << "0.078 4.49 0.75 8\n";
    for (const auto& item : indices)
        if (m.movable_atom(item.second).el == EL_TYPE_C) spec << item.first << ' ';
    spec << "\n0 0 0\n";
    ch_pi_search_grid field(base, dims, m, spec);
    for (int i = 0; i <= 32; ++i) {
        const fl r = i*0.375;
        const fl direct = r < 8 ? -0.078*std::exp(-0.5*sqr((r-4.49)/0.75)) : 0;
        if (std::abs(field.at(vec(r,0,0))-direct)>1e-12) throw std::runtime_error("lattice sample mismatch");
    }
    const fl h = 1e-5;
    const vec sample(4.11,1.19,-0.14);
    vec grad;
    const fl scalar = field.at(sample);
    if (std::abs(scalar-field.at(sample,&grad))>1e-12) throw std::runtime_error("derivative energy mismatch");
    fl cart_error = 0, conf_error = 0;
    for (sz i = 0; i < 3; ++i) {
        vec plus=sample, minus=sample;
        plus[i]+=h; minus[i]-=h;
        cart_error=std::max(cart_error,std::abs((field.at(plus)-field.at(minus))/(2*h)-grad[i]));
    }
    const precalculate_byatom pairs(*sf,m);
    const conf initial=m.get_initial_conf();
    flv increments(m.get_size().num_degrees_of_freedom(), 0.0137);
    conf origin=initial;
    origin.set_incremented(initial,increments,1);
    m.set(origin);
    flv g0(increments.size()), g1(increments.size());
    const vec cap(1000,1000,1000);
    const fl e0=m.eval_deriv(pairs,base,cap,g0), e1=m.eval_deriv(pairs,field,cap,g1);
    const auto correction=[&]{return field.eval(m,1000)-base.eval(m,1000);};
    if (std::abs(e1-e0-correction())>1e-9 || std::abs(correction())<1e-8)
        throw std::runtime_error("integrated energy mismatch or inactive field");
    for (sz i=0;i<g0.size();++i) {
        flv axis(g0.size(),0); axis[i]=1;
        conf perturbed=origin;
        perturbed.set_incremented(origin,axis,h); m.set(perturbed); const fl plus=correction();
        perturbed.set_incremented(origin,axis,-h); m.set(perturbed); const fl minus=correction();
        conf_error=std::max(conf_error,std::abs((plus-minus)/(2*h)-(g1[i]-g0[i])));
    }
    if (cart_error>1e-8 || conf_error>1e-6) throw std::runtime_error("finite difference mismatch");
    std::cout << "node_checks=33 cartesian_max_error=" << cart_error
              << " conformation_dof=" << g0.size() << " conformation_max_error=" << conf_error << '\n';
}
