#include "cache.h"
#include "parse_pdbqt.h"
#include "precalculate.h"
#include "quasi_newton.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

fl compare(model& m, const conf& state, const flv& p) {
    const fl h = 1e-5;
    conf tmp = state;
    m.set(state);
    const vecv base = m.get_chemical_heavy_atom_movable_coords();
    tmp.set_incremented(state, p, h);
    m.set(tmp);
    const vecv plus = m.get_chemical_heavy_atom_movable_coords();
    tmp.set_incremented(state, p, -h);
    m.set(tmp);
    const vecv minus = m.get_chemical_heavy_atom_movable_coords();
    fl finite_difference = 0;
    VINA_FOR_IN(i, plus)
        finite_difference = std::max(finite_difference, (plus[i] - minus[i]).norm() / (2*h));
    const fl analytic = m.max_heavy_atom_speed(state, p);
    const vecv restored = m.get_chemical_heavy_atom_movable_coords();
    VINA_FOR_IN(i, restored)
        require(vec_distance_sqr(restored[i],base[i]) == 0, "speed helper did not restore base");
    const fl error = std::abs(analytic - finite_difference);
    require(std::isfinite(error) && error < 1e-6, "Cartesian tangent mismatch");
    return error;
}

void layout(model m, const char* name) {
    const conf initial = m.get_initial_conf();
    const auto layout = m.get_size();
    const sz n = layout.num_degrees_of_freedom();
    unsigned checked = 0;
    fl max_error = 0;
    for (sz pose=0; pose<6; ++pose) {
        flv offset(n);
        VINA_FOR_IN(i,offset) offset[i] = std::sin((i+2)*.81);
        conf state = initial;
        state.set_incremented(initial, offset, pose*.71);
        for (sz k=0; k<n+4; ++k) {
            flv p(n,0);
            if (k<n) p[k] = 1;
            else if (k<n+3) {
                const sz subset = n - sum(layout.flex);
                VINA_FOR_IN(i,p) {
                    const bool selected = k<n+2 || (layout.flex.empty() ? i>=6+layout.ligands.front() : i>=subset);
                    if (selected) p[i] = std::cos((.29+k-n)*(i+1)+.17*pose);
                }
            }
            max_error=std::max(max_error, compare(m,state,p));
            ++checked;
        }
    }
    std::cout << name << "\tligands=" << m.num_ligands() << "\tflex=" << m.num_flex() << "\tdof=" << n
              << "\tdirection_cases=" << checked << "\tmax_abs_error=" << max_error << '\n';
}

// Give the real single-atom model the derivative of a native wall at x=2999,
// while scripting just the energy acceptance/rejection sequence.
class scripted_grid : public igrid {
public:
    cache underlying;
    mutable flv x;
    bool accept_last;
    scripted_grid(const std::string& prefix, bool accept) : underlying(1024), accept_last(accept) {
        underlying.read(prefix);
    }
    fl eval(const model&, fl) const override { throw std::logic_error("unexpected score-only"); }
    fl eval_intra(model&, fl) const override { throw std::logic_error("unexpected intra"); }
    fl eval_deriv(model& m, fl v) const override {
        underlying.eval_deriv(m,v);
        x.push_back(m.movable_coords(0)[0]);
        if (x.size()==1) return 0;
        if (x.size()==2) return -1000;
        if (x.size()==12 && accept_last) return -2000;
        return 100;
    }
};

void failure(const std::string& prefix, bool accept_last, fl spacing) {
    const std::string pdbqt = "ROOT\nATOM      1  C   UNL A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\nENDROOT\nTORSDOF 0\n";
    model m = parse_ligand_pdbqt_from_string(pdbqt,atom_type::XS);
    output_type pose(m.get_initial_conf(),0);
    scripted_grid grid(prefix,accept_last);
    quasi_newton local(m.get_size(),2);
    precalculate_byatom p;
    int evals=0;
    local(m,p,grid,pose,vec(10,10,10),evals,nullptr,0,spacing);
    if (evals!=12) { std::cerr << "trials"; for (fl x:grid.x) std::cerr << " " << x; std::cerr << "\n"; }
    require(evals==12 && grid.x.size()==12,"expected initial, one accept, ten trials");
    require(grid.x[0]==0 && grid.x[1]==1024,"incorrect accepted start");
    for(sz j=0;j<9;++j)
        require(grid.x[j+2]==1024 + 1024.0/(1u<<j),"changed early trial");
    require(grid.x[11]==1024 + (spacing>0?spacing:2),"incorrect final trial");
    const fl answer=accept_last?grid.x[11]:grid.x[1];
    require(pose.e==(accept_last?-2000:-1000) && pose.c.ligands[0].rigid.position[0]==answer
            && m.movable_coords(0)[0]==answer,"returned energy, conformation or model mismatch");
    std::cout << "scripted\tspacing=" << spacing << "\taccept_tenth=" << accept_last
              << "\tevals=" << evals << "\ttenth_x=" << grid.x[11]
              << "\treturned_x=" << answer << "\treturned_e=" << pose.e << '\n';
}

int main(int argc,char** argv) {
    if(argc!=6) return 2;
    std::cout << std::unitbuf << std::setprecision(17);
    model a=parse_ligand_pdbqt_from_file(argv[1],atom_type::XS);
    const model b=parse_ligand_pdbqt_from_file(argv[2],atom_type::XS);
    model both=a;
    both.append(b);
    layout(both,"two_ligands");
    model flex=parse_receptor_pdbqt(argv[3],argv[4],atom_type::XS);
    flex.append(a);
    layout(flex,"ligand_flex");
    flex.append(b);
    layout(flex,"two_ligands_flex");
    for(const fl spacing: {0.,.375})
        for(const bool pass:{false,true}) failure(argv[5],pass,spacing);
}
