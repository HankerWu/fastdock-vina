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
   by the receptor index and molecular interaction lists.
   See docs/PROVENANCE.md for upstream and experimental source identities. */

#ifndef VINA_PRECALCULATE_H
#define VINA_PRECALCULATE_H

#include <bitset>

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
    precalculate_element(const flv& distances, fl factor_, fl cap, const Evaluate& evaluate)
        : smooth(distances.size()), fast(distances.size()), factor(factor_) {
        VINA_FOR_IN(i, distances)
            smooth[i].first = (std::min)(cap, evaluate(distances[i]));
        init_from_smooth_fst(distances);
    }

    bool initialized() const { return !smooth.empty(); }
    fl eval_fast(fl r2) const{
        assert(r2 * factor < fast.size());
        return fast[sz(factor * r2)];
    };
    const fl* fast_samples() const { return fast.data(); }
    pr eval_deriv(fl r2) const{
        fl r2_factored = factor * r2;
        assert(r2_factored + 1 < smooth.size());
        sz i1 = sz(r2_factored);
        sz i2 = i1 + 1; // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
        assert(i1 < smooth.size());
        assert(i2 < smooth.size());
        fl rem = r2_factored - i1;
        assert(rem >= -epsilon_fl);
        assert(rem < 1 + epsilon_fl);
        const pr &p1 = smooth[i1];
        const pr &p2 = smooth[i2];
        fl e = p1.first + rem * (p2.first - p1.first);
        fl dor = p1.second + rem * (p2.second - p1.second);
        return pr(e, dor);
    };
private:
    void init_from_smooth_fst(const flv &rs){
        sz n = smooth.size();
        VINA_CHECK(rs.size() == n);
        VINA_CHECK(fast.size() == n);
        VINA_FOR(i, n) {
            // calculate dor's
            fl &dor = smooth[i].second;
            if (i == 0 || i == n - 1)
                dor = 0;
            else {
                fl delta = rs[i + 1] - rs[i - 1];
                fl r = rs[i];
                dor = (smooth[i + 1].first - smooth[i - 1].first) / (delta * r);
            }
            // calculate fast's from smooth.first's
            fl f1 = smooth[i].first;
            fl f2 = (i + 1 >= n) ? 0 : smooth[i + 1].first;
            fast[i] = (f2 + f1) / 2;
        }
    };
    prv smooth; // [(e, dor)]
    flv fast;
    fl factor = 0;
};

struct precalculate {
    precalculate() = default;

    precalculate(const ScoringFunction& sf, const std::bitset<XS_TYPE_SIZE>& receptor_types,
                 fl cap = max_fl, fl factor = 32)
        : m_cutoff_sqr(sqr(sf.get_cutoff())), m_factor(factor),
          m_data(XS_TYPE_SIZE, precalculate_element()) {
        VINA_CHECK(factor > epsilon_fl);
        const flv distances = sample_distances(factor, m_cutoff_sqr);
        VINA_FOR(t2, XS_TYPE_SIZE) {
            VINA_RANGE(t1, 0, t2 + 1) {
                if (!receptor_types[t1] && !receptor_types[t2]) continue;
                m_data(t1, t2) = precalculate_element(distances, factor, cap,
                    [&](fl r) { return sf.eval(t1, t2, r); });
            }
        }
    }

    fl eval_fast(sz type_pair_index, fl r2) const {
        assert(r2 <= m_cutoff_sqr);
        return element(type_pair_index).eval_fast(r2);
    }
    sz fast_index(fl r2) const {
        assert(r2 <= m_cutoff_sqr);
        return sz(m_factor * r2);
    }
    const precalculate_element& element(sz type_pair_index) const {
        assert(m_data(type_pair_index).initialized());
        return m_data(type_pair_index);
    }
    pr eval_deriv(sz type_pair_index, fl r2) const {
        assert(r2 <= m_cutoff_sqr);
        return element(type_pair_index).eval_deriv(r2);
    }
    fl cutoff_sqr() const { return m_cutoff_sqr; }

private:
    fl m_cutoff_sqr = 0;
    fl m_factor = 0;
    triangular_matrix<precalculate_element> m_data;
};

struct precalculate_byatom {
    precalculate_byatom() = default;

    precalculate_byatom(const ScoringFunction& sf, const model& m, fl cap = max_fl, fl factor = 32)
        : m_cutoff_sqr(sqr(sf.get_cutoff())), m_max_cutoff_sqr(sqr(sf.get_max_cutoff())),
          m_use_type_tables(sf.get_atom_typing() == atom_type::XS),
          m_data(m_use_type_tables ? XS_TYPE_SIZE : m.num_atoms(), precalculate_element()) {
        VINA_CHECK(factor > epsilon_fl);
        const flv distances = sample_distances(factor, m_cutoff_sqr);

        const auto build_table = [&](precalculate_element& table, sz a, sz b, const flv& samples) {
            table = precalculate_element(samples, factor, cap,
                [&](fl r) { return sf.eval(m.atoms[a], m.atoms[b], r); });
        };
        if (!m_use_type_tables) {
            // AD4 depends on atom charges, so its lookup remains indexed by atom pair.
            VINA_FOR(j, m.num_atoms())
                VINA_RANGE(i, 0, j + 1)
                    build_table(m_data(i, j), i, j, distances);
            return;
        }

        const auto build_type_tables = [&](const interacting_pairs& pairs, const flv& samples) {
            for (const interacting_pair& pair : pairs) {
                auto& table = m_data(pair.type_pair_index);
                if (!table.initialized())
                    build_table(table, pair.a, pair.b, samples);
            }
        };

        // Build long-range tables first, so any shared pair has all needed samples.
        if (!m.glue_pairs.empty())
            build_type_tables(m.glue_pairs, sample_distances(factor, m_max_cutoff_sqr));
        for (const auto& ligand : m.ligands)
            build_type_tables(ligand.pairs, distances);
        build_type_tables(m.inter_pairs, distances);
        build_type_tables(m.other_pairs, distances);
    }

    fl eval_fast(const interacting_pair& pair, fl r2) const {
        return element(pair).eval_fast(r2);
    }
    pr eval_deriv(const interacting_pair& pair, fl r2) const {
        return element(pair).eval_deriv(r2);
    }
    fl cutoff_sqr() const { return m_cutoff_sqr; }
    fl max_cutoff_sqr() const { return m_max_cutoff_sqr; }

private:
    const precalculate_element& element(const interacting_pair& pair) const {
        return m_use_type_tables ? m_data(pair.type_pair_index)
                                 : m_data(pair.a, pair.b);
    }

    fl m_cutoff_sqr = 0;
    fl m_max_cutoff_sqr = 0;
    bool m_use_type_tables = false;
    triangular_matrix<precalculate_element> m_data;
};

#endif
