// Research-only, one CPU and one dry ligand. No reference coordinates enter this observer.
#pragma once
#include "model.h"
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <set>

namespace basin_trace {
enum Event : uint32_t { start = 1, gradient = 2, adopted = 3, hunt_return = 4, mc_accept = 5,
    second_return = 6, saved = 7, pooled = 8, explicit_input = 9, explicit_output = 10 };
enum Phase : uint32_t { other = 0, hunt = 1, second = 2, refine = 3 };
struct Writer {
    std::ofstream out;
    uint32_t trajectory = 0, step = 0, index = 0;
    Phase phase = other;
    std::array<uint32_t, 5> tags{};
    std::vector<sz> atom_indices;
    std::vector<double> values;
    bool pending = false;
    Writer() {
        const char* path = std::getenv("FD_BASIN_TRACE");
        VINA_CHECK(path);
        out.open(path, std::ios::binary);
        VINA_CHECK(out.good());
    }
    void flush() {
        if (!pending) return;
        out.write(reinterpret_cast<const char*>(tags.data()), sizeof(tags));
        out.write(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(double));
        VINA_CHECK(out.good());
        pending = false;
    }
    ~Writer() { flush(); }
    void initialize(const model& m) {
        if (!atom_indices.empty()) return;
        VINA_CHECK(m.num_ligands() == 1 && m.num_flex() == 0);
        std::set<sz> heavy;
        for (const auto& entry : m.chemical_heavy_atom_elements()) heavy.insert(entry.first);
        std::vector<uint32_t> serials;
        for (const auto& line : m.get_ligand(0).cont) {
            if (line.second && heavy.count(*line.second)) {
                serials.push_back(std::stoul(line.first.substr(6, 5)));
                atom_indices.push_back(*line.second);
            }
        }
        const uint32_t count = serials.size();
        VINA_CHECK(count == heavy.size() && count > 0);
        out.write(reinterpret_cast<const char*>(&count), sizeof(count));
        out.write(reinterpret_cast<const char*>(serials.data()), sizeof(uint32_t) * count);
        values.resize(1 + 3 * count);
    }
    void begin(Phase value) { phase = value; index = 0; }
    void write(Event event, const model& m, fl energy) {
        flush();
        initialize(m);
        if (event == gradient) ++index;
        tags = {{trajectory, step, phase, event, index}};
        values[0] = energy;
        sz n = 1;
        for (sz i : atom_indices)
            for (sz d = 0; d < 3; ++d) values[n++] = m.movable_coords(i)[d];
        pending = true;
    }
    void adopt() {
        VINA_CHECK(pending && tags[3] == gradient && tags[4] > 1);
        tags[3] = adopted;
    }
};
inline Writer& writer() { static Writer w; return w; }
}
