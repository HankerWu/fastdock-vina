// Analytic and central-difference oracle for AD4 grid assembly.
#include "vina.h"
#include "ad4cache.h"
#include "parse_pdbqt.h"
#include <iostream>
#include <stdexcept>

namespace {
void close(fl actual, fl expected, fl tol = 1e-9) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tol)
        throw std::runtime_error("AD4 field mismatch: " + std::to_string(actual) + " vs " + std::to_string(expected));
}

void field_grid() {
    grid_dims box;
    for (auto& axis : box) {
        axis.begin = -1;
        axis.end = 1;
        axis.n_voxels = 4;
    }
    grid a(box), e(box), d(box), transferred(box), no_e(box);
    for (fl q : {-0.45, 0., 0.45}) {
        for (sz k = 0; k < 5; ++k)
            for (sz j = 0; j < 5; ++j)
                for (sz i = 0; i < 5; ++i) {
                    const vec p = a.index_to_argument(i,j,k);
                    const fl affinity = 3 + 1.5*p[0] + p[1]*p[2];
                    const fl electric = 28 - 2*p[0] + 0.15*p[1];
                    const fl desolv = 1 + 2*p[2];
                    a.m_data(i,j,k) = affinity;
                    e.m_data(i,j,k) = electric;
                    d.m_data(i,j,k) = desolv;
                    transferred.m_data(i,j,k) = affinity + q*electric;
                    no_e.m_data(i,j,k) = 0;
                }
        for (fl cap : {0., 10., 1000., max_fl})
            for (const vec p : {vec(0.41,-0.21,0.12), vec(1.32,0.13,-1.19), vec(-1.23,0.16,0.33)}) {
                vec c = p;
                fl miss = 0;
                for (sz j=0; j<3; ++j) {
                    c[j] = std::min(1., std::max(-1., p[j]));
                    miss += std::abs(c[j]-p[j]);
                }
                const fl u = 3+1.5*c[0]+c[1]*c[2]+q*(28-2*c[0]+0.15*c[1])+std::abs(q)*(1+2*c[2]);
                const fl capped = (u>0 && cap!=max_fl) ? cap*u/(cap+u) : u;
                vec grad, transferred_grad;
                const fl slope = 17;
                close(a.evaluate_ad4(p,slope,cap,q,e,d,&grad), capped+slope*miss);
                close(a.evaluate_ad4(p,slope,cap,q,e,d,nullptr), capped+slope*miss);
                close(transferred.evaluate_ad4(p,slope,cap,q,no_e,d,&transferred_grad), capped+slope*miss);
                for (sz j=0; j<3; ++j) {
                    vec lo = p, hi = p;
                    lo[j] -= 1e-5;
                    hi[j] += 1e-5;
                    const fl fd = (a.evaluate_ad4(hi,slope,cap,q,e,d,nullptr) -
                                   a.evaluate_ad4(lo,slope,cap,q,e,d,nullptr)) / 2e-5;
                    close(grad[j], fd, 1e-7);
                    close(grad[j], transferred_grad[j]);
                    if (c[j]!=p[j]) close(grad[j], p[j]>c[j] ? slope : -slope);
                }
            }
    }
}

// Exercise the production decomposition with a mobile sidechain and a ligand.
void partitions(const std::string& prefix, const std::string& rigid, const std::string& flexible) {
    for (const auto& item : std::vector<std::pair<std::string,int>>{{"C",7},{"e",13},{"d",2}}) {
        ofile f(path(prefix+"."+item.first+".map"));
        f << "GRID_PARAMETER_FILE none\nGRID_DATA_FILE none\nMACROMOLECULE none\n"
          << "SPACING 4\nNELEMENTS 2 2 2\nCENTER 12 50 14\n";
        for (sz i=0; i<27; ++i) f << item.second << '\n';
    }
    model m = parse_receptor_pdbqt(rigid, flexible, atom_type::AD);
    const std::string ligand = "ROOT\nATOM      1   C1 LIG A   1       9.500  51.000  13.000  1.00  0.00    -0.450 C\nENDROOT\nTORSDOF 0\n";
    m.append(parse_ligand_pdbqt_from_string(ligand, atom_type::AD));
    m.set(m.get_initial_conf());
    ad4cache field(17);
    field.read(prefix);
    fl lig = 0, flex = 0;
    const fl cap = 10;
    for (sz i = 0; i < m.num_movable_atoms(); ++i) {
        const fl q = m.get_atom(i).charge;
        const fl raw = 7+q*13+std::abs(q)*2;
        (m.is_atom_in_ligand(i) ? lig : flex) += cap*raw/(cap+raw);
    }
    if (lig == 0 || flex == 0) throw std::runtime_error("missing ligand or movable sidechain");
    close(field.eval(m,cap), lig);
    close(field.eval_intra(m,cap), flex);
    close(field.eval_deriv(m,cap), lig+flex);
}

// An uninitialized client must not clear a caller's file. A completed
// no-solution search owns its requested output, and changing inputs resets it.
void output_lifecycle(const std::string& prefix) {
    Vina vina("ad4", 1, 331, 0);
    const std::string output = prefix+".pdbqt";
    const auto sentinel = [&]() { std::ofstream f(output); f << "existing"; };
    const auto content = [&]() {
        std::ifstream f(output);
        return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
    };
    sentinel();
    vina.write_poses(output, 5, 100);
    if (content()!="existing") throw std::runtime_error("unused Vina instance modified output");
    const std::string ligand =
        "ROOT\nATOM      1   C1 LIG A   1       9.500  51.000  13.000  1.00  0.00     0.000 C\n"
        "ATOM      2   C2 LIG A   1      29.500  51.000  13.000  1.00  0.00     0.000 C\nENDROOT\nTORSDOF 0\n";
    vina.load_maps(prefix);
    vina.set_ligand_from_string(ligand);
    vina.global_search(1, 5, 1, 200);
    vina.write_poses(output, 5, 100);
    if (!content().empty() || !vina.get_poses_energies(5,100).empty())
        throw std::runtime_error("infeasible completed search returned stale results");
    for (const auto& invalid : {std::make_pair(-1, 100.0), std::make_pair(5, -1.0)}) {
        sentinel();
        bool rejected = false;
        try { vina.write_poses(output, invalid.first, invalid.second); }
        catch (const vina_runtime_error&) { rejected = true; }
        if (!rejected || content() != "existing")
            throw std::runtime_error("invalid output arguments modified a previous file");
    }
    sentinel();
    vina.set_ligand_from_string(ligand);
    vina.write_poses(output, 5, 100);
    if (content()!="existing") throw std::runtime_error("new unsearched ligand modified output");
}
}

int main(int argc,char** argv) {
    if (argc != 4) return 2;
    field_grid();
    partitions(argv[1],argv[2],argv[3]);
    output_lifecycle(argv[1]);
    std::cout << "AD4 signed maps, redistribution, wall gradient and ligand/flex decomposition passed\n";
}
