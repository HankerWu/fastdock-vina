#include "vina.h"
#include "coords.h"

#include <cmath>
#include <iostream>
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
    vina.global_search(4, 9, 1.0, 3000);
    const auto all = vina.get_poses_energies(9, 100);
    const auto xyz = vina.get_poses_coordinates(9, 100);
    require(all.size()==8 && xyz.size()==8, "default native output should retain eight distinct poses");
    for (double range : {0., 1., 3., 100.}) {
        for (int limit : {0, 1, 4, 9, 15}) {
            const auto scores = vina.get_poses_energies(limit, range);
            const auto coordinates = vina.get_poses_coordinates(limit, range);
            std::istringstream lines(vina.get_poses(limit, range));
            std::string line;
            std::vector<double> written;
            while (std::getline(lines,line))
                if (line.find("REMARK VINA RESULT:")==0) written.push_back(std::stod(line.substr(19)));
            sz n=0;
            while (n<all.size() && n<sz(limit) && all[n][0]<=all[0][0]+range) ++n;
            require(scores.size()==n && coordinates.size()==n && written.size()==n, "getter count mismatch");
            for (sz i=0;i<n;++i) {
                require(scores[i]==all[i] && coordinates[i]==xyz[i], "getter rank or coordinates changed");
                require(std::abs(scores[i][0]-written[i])<0.00051, "API/PDBQT affinity mismatch");
            }
        }
    }
    std::cout << "Rank-window cases and 20 API/coordinate/PDBQT combinations passed.\n";
}
