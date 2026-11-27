#ifndef VINA_FIXED_RING_H
#define VINA_FIXED_RING_H

#include "model.h"

// Keep the prepared offsets on the union of the opened ring paths. Entity
// indices address MC moves (two rigid choices), whereas gradient indices
// address BFGS coordinates (six rigid components).
class fixed_ring {
    std::vector<std::pair<sz,sz>> angles;
    szv gradient_indices;
    szv moves;

    static unsigned mark(const branch& b, const interacting_pair& pair, sz& next,
                         std::vector<bool>& locked) {
        const sz i = next++;
        auto inside = [&](sz n){ return b.node.begin <= n && n < b.node.end; };
        unsigned hit = (inside(pair.a) ? 1u : 0u) | (inside(pair.b) ? 2u : 0u);
        for (const auto& child : b.children) hit |= mark(child, pair, next, locked);
        if (hit == 1 || hit == 2) locked[i] = true;
        return hit;
    }

public:
    explicit fixed_ring(const model& m) {
        const auto pairs = m.get_glue_pairs();
        if (pairs.empty()) return;
        const auto size = m.get_size();
        sz entity = 0, derivative = 0;
        for (sz i = 0; i < m.num_ligands(); ++i) {
            const auto lig = m.get_ligand(i);
            std::vector<bool> locked(size.ligands[i]);
            for (const auto& pair : pairs) {
                if (pair.a < lig.begin || pair.a >= lig.end || pair.b < lig.begin || pair.b >= lig.end)
                    continue;
                sz next = 0;
                for (const auto& child : lig.children) mark(child,pair,next,locked);
                VINA_CHECK(next == locked.size());
            }
            moves.push_back(entity++);
            moves.push_back(entity++);
            derivative += conf_size::rigid_degrees_of_freedom;
            for (sz j = 0; j < locked.size(); ++j, ++entity, ++derivative) {
                if (locked[j]) {
                    angles.emplace_back(i,j);
                    gradient_indices.push_back(derivative);
                } else moves.push_back(entity);
            }
        }
        for (sz n : size.flex)
            for (sz j = 0; j < n; ++j) moves.push_back(entity++);
    }
    bool active() const { return !angles.empty(); }
    const szv& allowed_moves() const { return moves; }
    void prepare(conf& c) const {
        for (const auto& a : angles) c.ligands[a.first].torsions[a.second] = 0;
    }
    void project(flv& g) const {
        for (sz i : gradient_indices) g[i] = 0;
    }
};

#endif
