/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

/* Modified for FastDock Vina: build lookup tables for the type pairs used
   by the receptor index and molecular interaction lists, and share the
   distance-to-sample scale across each set of tables.
   See docs/PROVENANCE.md for upstream and experimental source identities. */

#ifndef VINA_PRECALCULATE_H
#define VINA_PRECALCULATE_H

#include <bitset>
#include <cstddef>

#include "model.h"
#include "scoring_function.h"

inline flv sample_distances(fl factor, fl cutoff_sqr) {
    // Interpolation and centered derivatives need the two samples above the cutoff.
    flv distances(sz(factor * cutoff_sqr) + 3);
    VINA_FOR_IN(i, distances)
        distances[i] = std::sqrt(i / factor);
    return distances;
}

struct precalculate_element
{
public:
    precalculate_element() = default;

    template<typename Evaluate>
    precalculate_element(const flv& distances, fl cap, const Evaluate& evaluate)
        : smooth(distances.size()) {
        VINA_FOR_IN(i, distances)
            smooth[i][0] = (std::min)(cap, evaluate(distances[i]));
        // The centered radial derivative is zero at both ends of the table.
        for (sz i = 1; i + 1 < distances.size(); ++i) {
            const fl delta = distances[i + 1] - distances[i - 1];
            smooth[i][1] = (smooth[i + 1][0] - smooth[i - 1][0]) / (delta * distances[i]);
        }
    }

    bool initialized() const { return !smooth.empty(); }
    fl fast_sample(sz i) const {
        assert(i + 1 < smooth.size());
        return (smooth[i + 1][0] + smooth[i][0]) / 2;
    }
    pr eval_deriv(fl r2_factored) const {
        assert(r2_factored + 1 < smooth.size());
        // Valid sample offsets fit the signed difference type of the vector.
        const auto i = static_cast<std::ptrdiff_t>(r2_factored);
        const fl rem = r2_factored - i;
        assert(rem >= -epsilon_fl);
        assert(rem < 1 + epsilon_fl);
        const sample lower = smooth[i];
        const sample result = lower + rem * (smooth[i + 1] - lower);
        return pr(result[0], result[1]);
    };
private:
    typedef fl sample __attribute__((vector_size(16)));
    std::vector<sample> smooth; // energy and radial derivative interpolate independently
};

struct precalculate {
    precalculate() = default;

    explicit precalculate(fl cutoff) : m_cutoff_sqr(sqr(cutoff)) {}

    void populate(const ScoringFunction& sf, const std::bitset<XS_TYPE_SIZE>& receptor_types,
                  const std::bitset<XS_TYPE_SIZE>& ligand_types) {
        if (receptor_types.none() || ligand_types.none()) return;
        if (m_data.dim() == 0)
            m_data = triangular_matrix<precalculate_element>(XS_TYPE_SIZE, precalculate_element());
        flv distances;
        VINA_FOR(t2, XS_TYPE_SIZE) {
            VINA_RANGE(t1, 0, t2 + 1) {
                const bool used = (receptor_types[t1] && ligand_types[t2]) ||
                                  (receptor_types[t2] && ligand_types[t1]);
                if (!used || m_data(t1, t2).initialized()) continue;
                if (distances.empty()) distances = sample_distances(m_factor, m_cutoff_sqr);
                m_data(t1, t2) = precalculate_element(distances, max_fl,
                    [&](fl r) { return sf.eval(t1, t2, r); });
            }
        }
    }

    fl eval_fast(sz type_pair_index, fl r2) const {
        return element(type_pair_index).fast_sample(fast_index(r2));
    }
    sz fast_index(fl r2) const {
        assert(r2 <= m_cutoff_sqr);
        return sz(m_factor * r2);
    }
    // Isolate a linear component using the same midpoint samples as eval_fast.
    template<class Potential>
    fl eval_fast_component(fl r2, const Potential& potential) const {
        const sz i = fast_index(r2);
        return (potential(std::sqrt(i / m_factor)) + potential(std::sqrt((i + 1) / m_factor))) / 2;
    }
    const precalculate_element& element(sz type_pair_index) const {
        assert(m_data(type_pair_index).initialized());
        return m_data(type_pair_index);
    }
    pr eval_deriv(sz type_pair_index, fl r2) const {
        assert(r2 <= m_cutoff_sqr);
        return element(type_pair_index).eval_deriv(m_factor * r2);
    }
    fl cutoff_sqr() const { return m_cutoff_sqr; }

private:
    fl m_cutoff_sqr = 0;
    static constexpr fl m_factor = 32;
    triangular_matrix<precalculate_element> m_data;
};

struct precalculate_byatom {
    precalculate_byatom() = default;

    precalculate_byatom(const ScoringFunction& sf, const model& m, fl cap = max_fl, fl factor = 32)
        : m_cutoff_sqr(sqr(sf.get_cutoff())), m_max_cutoff_sqr(sqr(sf.get_max_cutoff())),
          m_factor(factor),
          m_data(sf.get_atom_typing() == atom_type::XS ? XS_TYPE_SIZE : m.num_atoms(), precalculate_element()) {
        VINA_CHECK(factor > epsilon_fl);
        const flv distances = sample_distances(factor, m_cutoff_sqr);

        const auto build_tables = [&](const interacting_pairs& pairs, const flv& samples) {
            for (const interacting_pair& pair : pairs) {
                auto& table = m_data(pair.table_index);
                if (!table.initialized())
                    table = precalculate_element(samples, cap,
                        [&](fl r) { return sf.eval(m.atoms[pair.a], m.atoms[pair.b], r); });
            }
        };

        // Build long-range tables first, so any shared pair has all needed samples.
        if (m_max_cutoff_sqr > m_cutoff_sqr && !m.glue_pairs.empty())
            build_tables(m.glue_pairs, sample_distances(factor, m_max_cutoff_sqr));
        else
            build_tables(m.glue_pairs, distances);
        for (const auto& ligand : m.ligands)
            build_tables(ligand.pairs, distances);
        build_tables(m.inter_pairs, distances);
        build_tables(m.other_pairs, distances);
    }

    fl eval_fast(const interacting_pair& pair, fl r2) const {
        return element(pair).fast_sample(sz(m_factor * r2));
    }
    pr eval_deriv(const interacting_pair& pair, fl r2) const {
        return element(pair).eval_deriv(m_factor * r2);
    }
    fl cutoff_sqr() const { return m_cutoff_sqr; }
    fl max_cutoff_sqr() const { return m_max_cutoff_sqr; }

private:
    const precalculate_element& element(const interacting_pair& pair) const {
        assert(m_data(pair.table_index).initialized());
        return m_data(pair.table_index);
    }

    fl m_cutoff_sqr = 0;
    fl m_max_cutoff_sqr = 0;
    fl m_factor = 0;
    triangular_matrix<precalculate_element> m_data;
};

#endif
