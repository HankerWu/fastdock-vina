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

#ifndef VINA_SCORING_FUNCTION_H
#define VINA_SCORING_FUNCTION_H

#include <stdlib.h>
#include <list>
#include <cmath>
#include <stdexcept>
#include <string>
#include "atom.h"
#include "conf_independent.h"
#include "potentials.h"
#include "common.h"


//Forward declaration
struct model;

enum scoring_function_choice {SF_VINA, SF_AD42, SF_VINARDO};

// FastDock: one ordinary interaction radius for all XS potential terms.
constexpr fl default_xs_interaction_cutoff = 8.0;
// Ordinary XS contacts and tabulated closure lookups are bounded; the closure restraint itself is not.
constexpr fl xs_table_range_limit = 20.0;

class ScoringFunction {
public:
    ScoringFunction() {
        m_num_potentials = 0;
        m_num_conf_independents = 0;
    }
    ScoringFunction(const scoring_function_choice sf_choice, const flv& weights,
                    fl xs_cutoff = default_xs_interaction_cutoff, bool smooth_cutoff = false)
        : m_smooth_cutoff(smooth_cutoff) {
        if (smooth_cutoff && (sf_choice != SF_VINA || xs_cutoff <= 1))
            throw std::invalid_argument("smooth_cutoff requires Vina and an interaction_cutoff greater than 1 Angstrom");
        if (sf_choice != SF_AD42 && (!std::isfinite(xs_cutoff) || xs_cutoff <= 0 || xs_cutoff > xs_table_range_limit))
            throw std::invalid_argument("interaction_cutoff must be positive and at most " +
                                        std::to_string(static_cast<int>(xs_table_range_limit)) + " Angstrom");
        switch (sf_choice)
        {
            case SF_VINA:
            {
                m_potentials.push_back(new vina_gaussian(0, 0.5, xs_cutoff));
                m_potentials.push_back(new vina_gaussian(3, 2.0, xs_cutoff));
                m_potentials.push_back(new vina_repulsion(0.0, xs_cutoff));
                m_potentials.push_back(new vina_hydrophobic(0.5, 1.5, xs_cutoff));
                m_vina_hbond_index = m_potentials.size();
                m_potentials.push_back(new vina_non_dir_h_bond(-0.7, 0, xs_cutoff));
                m_potentials.push_back(new linearattraction());
                m_conf_independents.push_back(new num_tors_div());
                m_atom_typing = atom_type::XS;
                m_cutoff = xs_cutoff;
                m_max_cutoff = xs_table_range_limit;
                break;
            }
            case SF_VINARDO:
            {
                m_potentials.push_back(new vinardo_gaussian(0, 0.8, xs_cutoff));
                m_potentials.push_back(new vinardo_repulsion(0, xs_cutoff));
                m_potentials.push_back(new vinardo_hydrophobic(0, 2.5, xs_cutoff));
                m_potentials.push_back(new vinardo_non_dir_h_bond(-0.6, 0, xs_cutoff));
                m_potentials.push_back(new linearattraction());
                m_conf_independents.push_back(new num_tors_div());
                m_atom_typing = atom_type::XS;
                m_cutoff = xs_cutoff;
                m_max_cutoff = xs_table_range_limit;
                break;
            }
            case SF_AD42:
            {
                constexpr fl ad4_short_cutoff = 8.0;
                constexpr fl ad4_long_cutoff = 20.48;
                m_potentials.push_back(new ad4_vdw(0.5, 100000, ad4_short_cutoff));
                m_potentials.push_back(new ad4_hb(0.5, 100000, ad4_short_cutoff));
                m_potentials.push_back(new ad4_electrostatic(100, ad4_long_cutoff));
                m_potentials.push_back(new ad4_solvation(3.6, 0.01097, true, ad4_long_cutoff));
                m_potentials.push_back(new linearattraction());
                m_conf_independents.push_back(new ad4_tors_add());
                m_atom_typing = atom_type::AD;
                m_cutoff = ad4_long_cutoff;
                m_max_cutoff = ad4_long_cutoff;
                break;
            }
            default:
            {
                std::cout << "INSIDE everything::everything()   sfchoice = " << sf_choice << "\n";
                VINA_CHECK(false);
                break;
            }
        }

        m_num_potentials = m_potentials.size();
        m_num_conf_independents = m_conf_independents.size();
        m_weights = weights;
    };
    void Destroy()
    {
        for (auto p : m_potentials)
        {
            delete p;
        }
        m_potentials.clear();
        m_num_potentials = 0;
        for (auto p : m_conf_independents)
        {
            delete p;
        }
        m_conf_independents.clear();
        m_num_conf_independents = 0;
    }
    ~ScoringFunction() {
        Destroy();
    }

    fl eval(const atom& a, const atom& b, fl r) const{ // intentionally not checking for cutoff
        return pair_energy(a, b, r);
    }
    fl eval(sz t1, sz t2, fl r) const {
        return pair_energy(t1, t2, r);
    }
    fl vina_hbond(sz t1, sz t2, fl r) const {
        const fl e = m_weights[m_vina_hbond_index] * m_potentials[m_vina_hbond_index]->eval(t1, t2, r);
        return m_smooth_cutoff ? e * switch_factor(r) : e;
    }
    bool smooth_cutoff() const { return m_smooth_cutoff; }
    fl conf_independent(const model& m, fl e) const{
        // Iterator for weights
        flv::const_iterator it = m_weights.begin() + m_num_potentials;
        conf_independent_inputs in(m); // FIXME quite inefficient, but I think speed is irrelevant here, right?
        VINA_FOR (i, m_num_conf_independents)
        {
            // We don't accumulate energy. Why? I don't know...
            e = m_conf_independents[i]->eval(in, e, it);
        }
        assert(it == m_weights.end());
        return e;
    };
    fl get_cutoff() const { return m_cutoff; }
    fl get_max_cutoff() const { return m_max_cutoff; }
    fl glue_weight() const { return m_weights[m_num_potentials - 1]; }
    atom_type::t get_atom_typing() const { return m_atom_typing; }
    szv get_atom_types() const{
        szv tmp;
        VINA_FOR(i, num_atom_types(m_atom_typing))
        {
          tmp.push_back(i);
        }
        return tmp;
    };
    sz get_num_atom_types() const { return num_atom_types(m_atom_typing); }
    flv get_weights() const { return m_weights; }

private:
    // The transition width is fixed at one Angstrom. Evaluate from the outer
    // endpoint to retain precision as the interaction vanishes.
    fl switch_factor(fl r) const {
        const fl t = m_cutoff - r; // (cutoff - r) / 1 Angstrom
        if (t >= 1) return 1;
        if (t <= 0) return 0;
        return t * t * t * (10 + t * (6 * t - 15));
    }

    template<typename AtomOrType>
    fl pair_energy(const AtomOrType& a, const AtomOrType& b, fl r) const {
        // For switched Vina, the last potential is the macrocycle closure.
        // Add it separately so its physical range and strength are preserved.
        const int count = m_num_potentials - int(m_smooth_cutoff);
        fl e = 0;
        VINA_FOR(i, count)
            e += m_weights[i] * m_potentials[i]->eval(a, b, r);
        if (m_smooth_cutoff)
            e = e * switch_factor(r) + m_weights[count] * m_potentials[count]->eval(a, b, r);
        return e;
    }

    bool m_smooth_cutoff = false;
    std::vector<Potential*> m_potentials;
    std::vector<ConfIndependent*> m_conf_independents;
    sz m_vina_hbond_index = 0;
    flv m_weights;
    fl m_cutoff;
    fl m_max_cutoff;
    int m_num_potentials;
    int m_num_conf_independents;
    atom_type::t m_atom_typing;
};

#endif
