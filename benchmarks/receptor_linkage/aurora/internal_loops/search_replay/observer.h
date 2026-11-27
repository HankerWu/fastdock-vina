// Offline observer for one rigid-receptor ligand, one CPU. Read-only with
// respect to the production model, RNG, counters and minimizer.
#ifndef FASTDOCK_SEARCH_REPLAY_OBSERVER_H
#define FASTDOCK_SEARCH_REPLAY_OBSERVER_H
#include "model.h"
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>

namespace replay_trace {
struct State {
    std::ofstream out;
    unsigned trajectory = 0, step = 0, eval = 0;
    const char* phase = nullptr;
    uint64_t fingerprint = 0;
    unsigned focus[2][2] = {};
    State() {
        const char* path = std::getenv("FD_REPLAY_TRACE");
        VINA_CHECK(path);
        out.open(path);
        VINA_CHECK(out.good());
        out << std::setprecision(17);
        const char* steps = std::getenv("FD_REPLAY_HUNTS");
        if (steps) VINA_CHECK(std::sscanf(steps, "%u,%u,%u,%u",
            &focus[0][0], &focus[0][1], &focus[1][0], &focus[1][1]) == 4);
    }
    bool focused() const {
        return trajectory >= 1 && trajectory <= 2 && (step == focus[trajectory-1][0] || step == focus[trajectory-1][1]);
    }
    void prefix(const char* event) {
        out << trajectory << ' ' << step << ' ' << event;
    }
    void conf_values(const conf& c) {
        VINA_CHECK(c.ligands.size() == 1 && c.flex.empty());
        const auto& ligand = c.ligands.front();
        const auto& pose = ligand.rigid;
        for (sz i = 0; i < 3; ++i) out << ' ' << pose.position[i];
        const auto& q = pose.orientation;
        out << ' ' << q.R_component_1() << ' ' << q.R_component_2()
            << ' ' << q.R_component_3() << ' ' << q.R_component_4();
        for (fl t : ligand.torsions) out << ' ' << t;
    }
    void coordinates(const model& m) {
        for (const vec& atom : m.get_chemical_heavy_atom_movable_coords())
            for (sz i = 0; i < 3; ++i) out << ' ' << atom[i];
    }
    void position(const char* event, const conf& c) {
        prefix(event);
        conf_values(c);
        out << '\n';
    }
    void begin(const char* stage) {
        phase = stage;
        fingerprint = 14695981039346656037ULL;
        eval = 0;
    }
    void mix(fl value) {
        uint64_t bits;
        static_assert(sizeof(value) == sizeof(bits), "This trace expects binary64");
        std::memcpy(&bits, &value, sizeof(bits));
        fingerprint = (fingerprint ^ bits) * 1099511628211ULL;
    }
    void local(const char* event, const output_type& candidate, int calls, const model& m) {
        prefix("path");
        out << ' ' << phase << ' ' << fingerprint << '\n';
        phase = nullptr;
        prefix(event);
        out << ' ' << calls << ' ' << candidate.e;
        conf_values(candidate.c);
        out << '\n';
        if (focused()) {
            prefix("end_xyz");
            out << ' ' << event;
            coordinates(m);
            out << '\n';
        }
    }
    void value(const conf& c, fl energy, const flv& grad, const model& m) {
        if (!phase) return;
        mix(energy);
        for (fl g : grad) mix(g);
        if (!focused()) return;
        prefix("eval");
        out << ' ' << phase << ' ' << ++eval << ' ' << energy;
        conf_values(c);
        out << " grad";
        for (fl g : grad) out << ' ' << g;
        out << " xyz";
        coordinates(m);
        out << '\n';
    }
};
inline State& current() {
    static State state;
    return state;
}
}
#endif
