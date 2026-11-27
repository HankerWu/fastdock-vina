// Offline single-trajectory observer. The shadow model owns all written coordinates.
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "model.h"

struct path_observer {
    model shadow;
    std::ofstream stream;
    std::vector<std::pair<int, sz>> serials;
    path_observer(const model& source) : shadow(source), stream(std::getenv("FASTDOCK_PATH_TRACE"), std::ios::app) {
        VINA_CHECK(stream.good() && source.num_ligands() == 1 && !source.num_flex());
        for (const auto& line : source.get_ligand(0).cont) {
            if (!line.second || source.get_atom(*line.second).is_hydrogen()) continue;
            serials.emplace_back(std::stoi(line.first.substr(6, 5)), *line.second);
        }
        std::sort(serials.begin(), serials.end());
        stream << std::setprecision(17) << "serial";
        for (const auto& item : serials) stream << ',' << item.first;
        stream << '\n';
    }
    void write(sz step, const char* stage, const output_type& pose, int evals, const std::string& detail) {
        stream << step << ',' << stage << ',' << evals << ',' << pose.e << ',' << detail;
        shadow.set(pose.c);
        for (const auto& item : serials) {
            const vec x = shadow.get_coords(item.second);
            for (sz k=0; k<3; ++k) stream << ',' << x[k];
        }
        stream << '\n';
    }
};

static std::string mutation_kind(const conf& before, const conf& after) {
    const auto& a = before.ligands.front();
    const auto& b = after.ligands.front();
    if (vec_distance_sqr(a.rigid.position,b.rigid.position) > 0) return "translation";
    if (a.rigid.orientation != b.rigid.orientation) return "rotation";
    for (sz i=0; i<a.torsions.size(); ++i)
        if (a.torsions[i] != b.torsions[i]) return "torsion_" + std::to_string(i);
    return "none";
}
