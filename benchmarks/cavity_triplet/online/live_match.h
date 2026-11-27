// POSIX research adapter. Exchange one real, post-mutation proposal over a
// socketpair with the invoking process. No RNG draw and no extra Vina eval.
#pragma once
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <unistd.h>

class live_match {
    int fd = -1;
    unsigned trajectory = 0, step = 0, current = 0;
public:
    live_match() {
        if (const char* path = std::getenv("FD_TRIPLET_SOCKET")) {
            fd = std::stoi(path);
            trajectory = std::stoul(std::getenv("FD_TRIPLET_TRAJECTORY"));
            step = std::stoul(std::getenv("FD_TRIPLET_STEP"));
        }
    }
    void begin() { ++current; }
    void apply(unsigned n, const model& original, conf& candidate) {
        if (fd < 0 || current != trajectory || n != step) return;
        VINA_CHECK(candidate.ligands.size() == 1 && candidate.flex.empty());
        model state(original);
        state.set(candidate);
        const int copy = dup(fd);
        VINA_CHECK(copy >= 0);
        FILE* channel = fdopen(copy, "r+");
        VINA_CHECK(channel);
        auto& rigid = candidate.ligands.front().rigid;
        const auto& q = rigid.orientation;
        std::fprintf(channel, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                     rigid.position[0], rigid.position[1], rigid.position[2],
                     q.R_component_1(), q.R_component_2(), q.R_component_3(), q.R_component_4());
        std::set<sz> heavy;
        for (const auto& a : state.chemical_heavy_atom_elements()) heavy.insert(a.first);
        std::fprintf(channel, "%zu\n", heavy.size());
        std::vector<std::pair<unsigned long, sz>> by_serial;
        for (const auto& line : state.get_ligand(0).cont)
            if (line.second && heavy.count(*line.second))
                by_serial.emplace_back(std::stoul(line.first.substr(6, 5)), *line.second);
        std::sort(by_serial.begin(), by_serial.end());
        for (sz i = 0; i < by_serial.size(); ++i) {
            VINA_CHECK(i == 0 || by_serial[i-1].first < by_serial[i].first);
            const vec& p = state.movable_coords(by_serial[i].second);
            std::fprintf(channel, "%.17g %.17g %.17g\n", p[0], p[1], p[2]);
        }
        VINA_CHECK(by_serial.size() == heavy.size() && std::fflush(channel) == 0);
        double value[7];
        for (double& v : value) VINA_CHECK(std::fscanf(channel, "%lf", &v) == 1 && std::isfinite(v));
        std::fclose(channel);
        for (sz i = 0; i < 3; ++i) rigid.position[i] = value[i];
        rigid.orientation = qt(value[3], value[4], value[5], value[6]);
        VINA_CHECK(quaternion_is_normalized(rigid.orientation));
    }
};
