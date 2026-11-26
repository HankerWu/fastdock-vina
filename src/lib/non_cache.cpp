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

// FastDock: the explicit receptor owns its index, tables and shared energy pass.

#include "non_cache.h"
#include "cache.h"
#include "curl.h"
#include "net_polar.h"

non_cache::non_cache(const model& receptor, const grid_dims& bounds,
                     std::shared_ptr<const ScoringFunction> scoring_, fl slope_)
	: slope(slope_), sgrid(receptor, szv_grid_dims(bounds), sqr(scoring_->get_cutoff())),
	  gd(bounds), p(scoring_->get_cutoff()), scoring(std::move(scoring_)) {}

void non_cache::prepare_for_model(const model& m) {
	if (scoring)
		prepare_for_model(m, *scoring, p);
}

void non_cache::prepare_for_model(const model& m, const ScoringFunction& sf, precalculate& tables) const {
	std::bitset<XS_TYPE_SIZE> types;
	for (sz type : m.get_movable_atom_types(atom_type::XS))
		if (!is_glue_type(type)) types.set(type);
	tables.populate(sf, sgrid.atom_types(), types);
}

void non_cache::populate_maps(cache& maps, const model& receptor, const std::bitset<XS_TYPE_SIZE>& types,
                              bool retain_tables) {
	if (retain_tables) {
		p.populate(*scoring, sgrid.atom_types(), types);
		maps.populate(receptor, p, sgrid, types);
	} else {
		populate_maps_from(maps, receptor, types, *scoring);
	}
}

void non_cache::populate_maps_from(cache& maps, const model& receptor, const std::bitset<XS_TYPE_SIZE>& types,
                                   const ScoringFunction& scoring_function) const {
	precalculate temporary(scoring_function.get_cutoff());
	temporary.populate(scoring_function, sgrid.atom_types(), types);
	maps.populate(receptor, temporary, sgrid, types);
}

fl non_cache::eval(const model& m, fl v) const {
	return eval_impl<false>(m, v, p, false, nullptr, nullptr, nullptr);
}

fl non_cache::eval(const model& m, fl v, const precalculate& scores) const {
	return eval_impl<false>(m, v, scores, false, nullptr, nullptr, nullptr);
}

fl non_cache::eval_intra(model& m, fl v) const {
	if (m.num_flex() == 0) return 0;
	return eval_impl<false>(m, v, p, true, nullptr, nullptr, nullptr);
}

fl non_cache::directional_eval(const model& m, fl v, const ScoringFunction& sf,
                               const net_polar& polar, net_polar_terms& terms) const {
	return eval_impl<true>(m, v, p, false, &sf, &polar, &terms);
}

template<bool Directional>
fl non_cache::eval_impl(const model& m, fl v, const precalculate& scores, bool flex_only,
                         const ScoringFunction* sf, const net_polar* polar, net_polar_terms* terms) const {
	fl e = 0;
	const fl cutoff_sqr = scores.cutoff_sqr();
	const auto& receptor = m.grid_atoms();

	sz n = num_atom_types(atom_type::XS);

	VINA_FOR(i, m.num_movable_atoms()) {
		if (flex_only && m.is_atom_in_ligand(i)) continue;
		fl this_e = 0;
		fl direction_delta = 0;
		fl out_of_bounds_penalty = 0;
		const atom& a = m.atoms[i];
		if (a.xs >= n || is_glue_type(a.xs)) continue;

		const vec& a_coords = m.coords[i];
		vec adjusted_a_coords; adjusted_a_coords = a_coords;
		VINA_FOR_IN(j, gd) {
			if(gd[j].n_voxels > 0) {
				if     (a_coords[j] < gd[j].begin) { adjusted_a_coords[j] = gd[j].begin; out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin); }
				else if(a_coords[j] > gd[j].end  ) { adjusted_a_coords[j] = gd[j]  .end; out_of_bounds_penalty += std::abs(a_coords[j] - gd[j]  .end); }
			}
		}
		out_of_bounds_penalty *= slope;

		const szv& possibilities = sgrid.possibilities(adjusted_a_coords);

		VINA_FOR_IN(possibilities_j, possibilities) {
			const sz j = possibilities[possibilities_j];
			const atom& b = receptor[j];
			sz t2 = b.get(atom_type::XS);
			if(t2 >= n) continue;
			vec r_ba; r_ba = adjusted_a_coords - b.coords; // FIXME why b-a and not a-b ?
			fl r2 = sqr(r_ba);
			if(r2 < cutoff_sqr) {
				sz type_pair_index = get_type_pair_index(atom_type::XS, a, b);
				this_e += scores.eval_fast(type_pair_index, r2);
				if (Directional && (a.el == EL_TYPE_N || a.el == EL_TYPE_O) &&
				    (b.el == EL_TYPE_N || b.el == EL_TYPE_O) && xs_h_bond_possible(a.xs, t2)) {
					const fl native_hbond = scores.eval_fast_component(r2, [&](fl distance) {
						return sf->vina_hbond(a.xs, t2, distance);
					});
					if (native_hbond != 0) {
						const fl directed_hbond = native_hbond * polar->direction(m, i, j, adjusted_a_coords, b.coords);
						direction_delta += directed_hbond - native_hbond;
						terms->native_hbond += native_hbond;
						terms->directed_hbond += directed_hbond;
					}
				}
			}
		}
		if (Directional) {
			fl replaced = this_e + direction_delta;
			curl(replaced, v);
			curl(this_e, v);
			terms->direction_delta += replaced - this_e;
		} else {
			curl(this_e, v);
		}
		e += this_e + out_of_bounds_penalty;
	}
	return e;
}

bool non_cache::within(const model& m, fl margin) const {
	VINA_FOR(i, m.num_movable_atoms()) {
		if(m.atoms[i].is_hydrogen()) continue;
		const vec& a_coords = m.coords[i];
		VINA_FOR_IN(j, gd)
			if(gd[j].n_voxels > 0)
				if(a_coords[j] < gd[j].begin - margin || a_coords[j] > gd[j].end + margin) 
					return false;
	}
	return true;
}

fl non_cache::eval_deriv(      model& m, fl v) const { // clean up
	fl e = 0;
	const fl cutoff_sqr = p.cutoff_sqr();
	const auto& receptor = m.grid_atoms();

	sz n = num_atom_types(atom_type::XS);

	VINA_FOR(i, m.num_movable_atoms()) {
		fl this_e = 0;
		vec deriv(0, 0, 0);
		vec out_of_bounds_deriv(0, 0, 0);
		fl out_of_bounds_penalty = 0;
		const atom& a = m.atoms[i];
		if (a.xs >= n || is_glue_type(a.xs)) {
			// The pair pass adds closure forces after each receptor evaluation.
			m.minus_forces[i].assign(0);
			continue;
		}
		const vec& a_coords = m.coords[i];
		vec adjusted_a_coords; adjusted_a_coords = a_coords;
		VINA_FOR_IN(j, gd) {
			if(gd[j].n_voxels > 0) {
				if     (a_coords[j] < gd[j].begin) { adjusted_a_coords[j] = gd[j].begin; out_of_bounds_deriv[j] = -1; out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].begin); }
				else if(a_coords[j] > gd[j].end  ) { adjusted_a_coords[j] = gd[j]  .end; out_of_bounds_deriv[j] =  1; out_of_bounds_penalty += std::abs(a_coords[j] - gd[j]  .end); }
			}
		}
		out_of_bounds_penalty *= slope;
		out_of_bounds_deriv *= slope;

		const szv& possibilities = sgrid.possibilities(adjusted_a_coords);

		VINA_FOR_IN(possibilities_j, possibilities) {
			const sz j = possibilities[possibilities_j];
			const atom& b = receptor[j];
			sz t2 = b.get(atom_type::XS);
			if(t2 >= n) continue;
			vec r_ba; r_ba = adjusted_a_coords - b.coords; // FIXME why b-a and not a-b ?
			fl r2 = sqr(r_ba);
			if(r2 < cutoff_sqr) {
				sz type_pair_index = get_type_pair_index(atom_type::XS, a, b);
				pr e_dor =  p.eval_deriv(type_pair_index, r2);
				this_e += e_dor.first;
				deriv += e_dor.second * r_ba;
			}
		}
		curl(this_e, deriv, v);
		m.minus_forces[i] = deriv + out_of_bounds_deriv;
		e += this_e + out_of_bounds_penalty;
	}
	return e;
}
