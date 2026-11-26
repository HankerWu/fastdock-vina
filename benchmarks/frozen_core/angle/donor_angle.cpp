// Evaluate an opt-in donor-angle correction on frozen native Vina poses.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include "parse_pdbqt.h"
#include "scoring_function.h"

struct PolarSite {
    sz index;
    sz xs;
    szv hydrogens;
};

std::vector<PolarSite> polar_sites(const model& m, bool rigid) {
    std::vector<PolarSite> result;
    sz n = rigid ? m.grid_atoms().size() : m.num_atoms();
    for (sz i = 0; i < n; ++i) {
        const atom& a = m.get_atom(atom_index(i, rigid));
        // Only N/O can supply the directional donor under test. Metal
        // coordination cannot be assigned an explicit hydrogen orientation.
        if ((a.el != EL_TYPE_N && a.el != EL_TYPE_O) ||
            (!xs_is_donor(a.xs) && !xs_is_acceptor(a.xs))) continue;
        PolarSite site{i, a.xs, {}};
        if (xs_is_donor(a.xs)) {
            for (const bond& edge : a.bonds) {
                const atom_index h = edge.connected_atom_index;
                if (h.in_grid == rigid && m.get_atom(h).ad == AD_TYPE_HD) {
                    // Reject a spuriously multibonded polar H instead of
                    // assigning the same H to two different donors.
                    if (m.get_atom(h).bonds.size() != 1)
                        throw std::runtime_error("Ambiguous covalent owner of donor H");
                    site.hydrogens.push_back(h.i);
                }
            }
            if (site.hydrogens.empty())
                throw std::runtime_error("Vina N/O donor has no bonded HD atom");
        }
        result.push_back(site);
    }
    return result;
}

fl donor_weight(const vec& d, const vec& h, const vec& a) {
    const vec dh = h - d, ha = a - h;
    const fl norm = std::sqrt(dh.norm_sqr() * ha.norm_sqr());
    if (norm <= 1e-12) return 0;
    const fl c = std::max(fl(0), std::min(fl(1), (dh * ha) / norm));
    return c * c;
}

struct Terms {
    sz contacts = 0;
    sz dual_role = 0;
    fl radial = 0;
    fl penalty = 0;
};

Terms evaluate(const model& ligand, const model& receptor,
               const std::vector<PolarSite>& ls, const std::vector<PolarSite>& rs,
               const ScoringFunction& sf) {
    Terms sum;
    for (const PolarSite& l : ls) {
        const vec lc = ligand.get_coords(l.index);
        for (const PolarSite& r : rs) {
            const vec rc = receptor.grid_atoms()[r.index].coords;
            const fl r2 = vec_distance_sqr(lc, rc);
            if (r2 >= 16) continue;  // all N/O Vina H-bond radii sums are below 4 A
            const fl e = sf.vina_hbond(l.xs, r.xs, std::sqrt(r2));
            if (e == 0) continue;
            const bool ligand_donates = !l.hydrogens.empty() && xs_is_acceptor(r.xs);
            const bool receptor_donates = !r.hydrogens.empty() && xs_is_acceptor(l.xs);
            if (!ligand_donates && !receptor_donates)
                throw std::runtime_error("Cannot orient a native N/O H-bond pair");
            fl g = 0;
            if (ligand_donates)
                for (sz h : l.hydrogens)
                    g = std::max(g, donor_weight(lc, ligand.get_coords(h), rc));
            if (receptor_donates)
                for (sz h : r.hydrogens)
                    g = std::max(g, donor_weight(rc, receptor.grid_atoms()[h].coords, lc));
            // One radial contribution per heavy-atom pair even when both
            // polar directions or multiple donor hydrogens are available.
            ++sum.contacts;
            sum.dual_role += ligand_donates && receptor_donates;
            sum.radial += e;
            sum.penalty += -e * (1 - g);
        }
    }
    return sum;
}

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    try {
        if (donor_weight(vec(0,0,0), vec(1,0,0), vec(3,0,0)) != 1 ||
            donor_weight(vec(0,0,0), vec(1,0,0), vec(1,2,0)) != 0 ||
            donor_weight(vec(0,0,0), vec(1,0,0), vec(-2,0,0)) != 0)
            throw std::runtime_error("D-H-A angle convention failed");
        model ligand = parse_ligand_pdbqt_from_file(argv[2], atom_type::XS);
        const model receptor = parse_receptor_pdbqt(argv[1]);
        const auto ls = polar_sites(ligand, false), rs = polar_sites(receptor, true);
        // Default weights from Vina::set_vina_weights. Use the scorer's
        // radial H-bond term and conformation-independent reporting transform.
        ScoringFunction sf(SF_VINA, {-0.035579,-0.005156,0.840245,-0.035069,
                                     -0.587439,50.0,5 * 0.05846 / 0.1 - 1});
        const fl scale = sf.conf_independent(ligand, 1);
        std::ifstream input(argv[3]);
        if (!input) throw std::runtime_error("Cannot open docking pool");
        std::ostringstream block;
        std::string line;
        bool inside = false;
        sz rank = 0;
        fl native = 0;
        bool score_seen = false;
        std::cout << "rank,contacts,dual_role,raw_radial,raw_penalty,scale,score_delta,native_score,corrected_score\n" << std::setprecision(12);
        while (std::getline(input,line)) {
            if (line.compare(0,6,"MODEL ") == 0) {
                if (inside) throw std::runtime_error("Nested MODEL");
                inside = true;
                block.str("");
                score_seen = false;
            } else if (line == "ENDMDL") {
                if (!inside || !score_seen) throw std::runtime_error("Incomplete scored MODEL");
                replace_fixed_pose_coordinates(ligand, block.str());
                const Terms t = evaluate(ligand, receptor, ls, rs, sf);
                std::cout << ++rank << ',' << t.contacts << ',' << t.dual_role << ',' << t.radial
                          << ',' << t.penalty << ',' << scale << ',' << t.penalty * scale << ','
                          << native << ',' << native + t.penalty * scale << '\n';
                inside = false;
            } else if (inside) {
                const std::string prefix = "REMARK VINA RESULT:";
                if (line.compare(0, prefix.size(), prefix) == 0) {
                    if (score_seen) throw std::runtime_error("Duplicate native score");
                    native = std::stod(line.substr(prefix.size()));
                    score_seen = true;
                }
                block << line << '\n';
            }
        }
        if (!rank || inside) throw std::runtime_error("No complete pose pool");
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
