// Single-worker frozen-bank experiment. Only the selected rigid transform is
// transferred back to search. The model used for type permutation is private.
#pragma once

#include "cache.h"
#include "inject.h"
#include <limits>

unsigned long long observed_scalar_grid_positions();

struct native_field_rank {
    static void select(rigid_injection& injection, unsigned trajectory, unsigned step,
                       const model& original, const igrid& field,
                       const conf& source, fl cap) {
        const char* path = std::getenv("FD_FIELD_BANK");
        if (!path || trajectory != injection.trajectory || step != injection.step) return;
        VINA_CHECK(source.ligands.size() == 1 && source.flex.empty());
        VINA_CHECK(original.atom_typing_used() == atom_type::XS && dynamic_cast<const cache*>(&field));
        const char* arm = std::getenv("FD_FIELD_ARM");
        VINA_CHECK(arm);
        const bool permute = std::string(arm) == "perm";
        VINA_CHECK(permute || std::string(arm) == "true");
        model probe(original);
        std::set<sz> heavy;
        for (const auto& atom : original.chemical_heavy_atom_elements()) heavy.insert(atom.first);
        std::vector<std::pair<unsigned long, sz>> serials;
        for (const auto& line : original.get_ligand(0).cont) {
            if (line.second && heavy.count(*line.second))
                serials.emplace_back(std::stoul(line.first.substr(6, 5)), *line.second);
        }
        std::sort(serials.begin(), serials.end());
        VINA_CHECK(serials.size() == heavy.size() && serials.size() >= 3);
        const sz n = serials.size();
        sz changed = 0;
        for (sz i = 0; i < n; ++i) {
            VINA_CHECK(!i || serials[i-1].first < serials[i].first);
            const sz dst = serials[i].second;
            const sz src = serials[permute ? (i + n/2) % n : i].second;
            const sz t = original.atoms[src].xs;
            VINA_CHECK(t < XS_TYPE_SIZE && cache::map_types({t}).count() == 1);
            changed += probe.atoms[dst].xs != t;
            probe.atoms[dst].xs = t;
            injection.log << "field_type " << serials[i].first << ' '
                          << original.atoms[dst].xs << ' ' << t << '\n';
        }
        std::ifstream bank(path);
        VINA_CHECK(bank.is_open());
        std::ofstream trace;
        if (const char* dest = std::getenv("FD_FIELD_TRACE")) {
            trace.open(dest);
            VINA_CHECK(trace.is_open());
            trace << "index,score\n" << std::setprecision(17);
        }
        conf placement(source);
        auto& rigid = placement.ligands[0].rigid;
        fl best = std::numeric_limits<fl>::infinity();
        sz selected = 0;
        const auto before = observed_scalar_grid_positions();
        for (sz i = 0; i < 1024; ++i) {
            std::array<fl, 7> v;
            for (fl& component : v) VINA_CHECK((bank >> component) && std::isfinite(component));
            rigid.position = vec(v[0], v[1], v[2]);
            rigid.orientation = qt(v[3], v[4], v[5], v[6]);
            VINA_CHECK(quaternion_is_normalized(rigid.orientation));
            probe.set(placement);
            const fl score = field.eval(probe, cap);
            VINA_CHECK(std::isfinite(score));
            if (trace.is_open()) trace << i << ',' << score << '\n';
            if (score < best) {
                best = score;
                selected = i;
                injection.value = v;
            }
        }
        bank >> std::ws;
        VINA_CHECK(bank.eof());
        const auto queries = observed_scalar_grid_positions() - before;
        VINA_CHECK(queries == 1024 * n);
        if (trace.is_open()) { trace.flush(); VINA_CHECK(trace.good()); }
        injection.log << "field_rank " << arm << " selected " << selected << " energy " << best
                      << " heavy " << n << " changed " << changed << " positions " << queries
                      << " cap " << cap << '\n';
    }
};
