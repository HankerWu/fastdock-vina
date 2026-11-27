#include "vina.h"
#include "coords.h"
#include "parse_pdbqt.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

// The first geometric representative has worse native energy but higher custom
// rank. Removing it by the requested energy window must release its neighbour.
void check_windows() {
    output_container pool;
    auto append = [&](fl x, fl e) {
        auto* pose = new output_type(conf{}, e);
        pose->coords.push_back(vec(x, 0, 0));
        pool.push_back(pose);
    };
    append(0, 2);
    append(0.5, -3);
    append(1.0, -3);
    append(4.0, -3);
    require(ranked_output_indices(pool, 9, 1, boost::none) == pose_indices({0, 2, 3}), "existing rank order / exact radius");
    require(ranked_output_indices(pool, 9, 1, -3.) == pose_indices({1, 3}), "energy eligibility before geometry");
    require(ranked_output_indices(pool, 1, 1, -3.) == pose_indices({1}), "maximum returned count follows filtering");
    require(ranked_output_indices(pool, 2, 0, -3.) == pose_indices({1, 2}), "zero radius / equal scores");
    require(ranked_output_indices(pool, 0, 1, boost::none).empty(), "zero count");
    for (unsigned i=0; i<20; ++i) append(20+3*i, -3);
    require(ranked_output_indices(pool, 30, 1, -3.).size()==22, "more than nine eligible poses");
}

int main(int argc, char** argv) {
    require(argc==3, "pass Gal3 3T1L receptor and ligand 10");
    check_windows();
    Vina vina("vina", 1, 391781, 0);
    vina.set_receptor(argv[1]);
    vina.set_ligand_from_file(argv[2]);
    vina.compute_vina_maps(17.36344, -8.90152, -1.04760, 24, 24, 24, 0.375);
    const int search_limit = 9;
    const double spacing = 1.0;
    const double wide_range = std::numeric_limits<double>::max();
    vina.global_search(4, search_limit, spacing, 3000);
    const auto all = vina.get_poses_energies(search_limit + 1, wide_range);
    const auto xyz = vina.get_poses_coordinates(search_limit + 1, wide_range);
    require(all.size() > 1 && all.size() <= sz(search_limit) && xyz.size() == all.size(),
            "real fixture must return a nontrivial pool within the search cap");
    const auto ligand = parse_ligand_pdbqt_from_file(argv[2], atom_type::XS);
    const auto heavy = ligand_heavy_indices_by_serial(ligand);
    require(!heavy.empty(), "missing ligand heavy atoms");
    double nearest = wide_range;
    for (sz i = 0; i < all.size(); ++i) {
        require(!all[i].empty() && std::isfinite(all[i][0]) && (!i || all[i-1][0] <= all[i][0]),
                "native energies must be finite and ranked");
        require(xyz[i].size() == 3 * ligand.num_atoms(), "incomplete coordinate vector");
        for (sz j = 0; j < i; ++j) {
            double squared = 0;
            for (const auto& atom : heavy) {
                for (sz dim = 0; dim < 3; ++dim) {
                    const sz k = 3 * atom.second + dim;
                    const double delta = xyz[i][k] - xyz[j][k];
                    squared += delta * delta;
                }
            }
            const double distance = std::sqrt(squared / heavy.size());
            require(distance + 1e-10 >= spacing, "native outputs violate the requested heavy-atom RMSD");
            nearest = std::min(nearest, distance);
        }
    }
    require(all.front()[0] < all.back()[0], "fixture must exercise energy-window exclusion");
    for (double range : {0., 1., 3., 100., wide_range}) {
        for (int limit : {0, 1, 4, 9, 15}) {
            const auto scores = vina.get_poses_energies(limit, range);
            const auto coordinates = vina.get_poses_coordinates(limit, range);
            sz n = 0;
            while (n < all.size() && n < sz(limit) && all[n][0] <= all[0][0] + range) ++n;
            require(scores.size() == n && coordinates.size() == n, "energy/coordinate getter count mismatch");
            for (sz i = 0; i < n; ++i)
                require(scores[i] == all[i] && coordinates[i] == xyz[i], "getter rank or coordinates changed");

            std::istringstream lines(vina.get_poses(limit, range));
            std::string line;
            sz pose = 0;
            bool has_score = false;
            std::set<unsigned> seen;
            while (std::getline(lines, line)) {
                if (line.find("REMARK VINA RESULT:") == 0) {
                    require(pose < n && !has_score && std::abs(scores[pose][0] - std::stod(line.substr(19))) < 0.00051,
                            "API/PDBQT affinity mismatch");
                    has_score = true;
                } else if (line.find("ATOM  ") == 0 || line.find("HETATM") == 0) {
                    const auto serial = static_cast<unsigned>(std::stoul(line.substr(6, 5)));
                    const auto atom = heavy.find(serial);
                    if (atom == heavy.end()) continue;
                    require(pose < n && seen.insert(serial).second, "unexpected PDBQT model or repeated heavy atom");
                    for (sz dim = 0; dim < 3; ++dim)
                        require(std::abs(coordinates[pose][3 * atom->second + dim] -
                                         std::stod(line.substr(30 + 8 * dim, 8))) < 0.00051,
                                "API/PDBQT heavy-atom coordinate mismatch");
                } else if (line == "ENDMDL") {
                    require(has_score && seen.size() == heavy.size(), "incomplete PDBQT model");
                    ++pose;
                    has_score = false;
                    seen.clear();
                }
            }
            require(pose == n && !has_score && seen.empty(), "PDBQT model count mismatch");
        }
    }
    std::cout << "Rank-window cases and 25 API/coordinate/PDBQT combinations passed for "
              << all.size() << " poses; nearest pair " << nearest << " Angstrom.\n";
}
