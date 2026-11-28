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

// FastDock: reuse the trial conformation across Monte Carlo steps.

#include "monte_carlo.h"
#include "coords.h"
#include "mutate.h"
#include "quasi_newton.h"
#include "trajectory_archive.h"

output_type monte_carlo::operator()(model& m, const precalculate_byatom& p, const igrid& ig, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator) const {
	output_container tmp;
	this->operator()(m, tmp, p, ig, corner1, corner2, increment_me, generator); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
	if(new_f < old_f) return true;
	const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
	return random_fl(0, 1, generator) < acceptance_probability;
}

// out is sorted
void monte_carlo::operator()(model& m, output_container& out, const precalculate_byatom& p, const igrid& ig, const vec& corner1, const vec& corner2, incrementable* increment_me, rng& generator, trajectory_archive* archive, search_trajectory_counters* counters,
                             const boost::optional<vec>& pocket_start, bool prepared_torsions) const {
    int evalcount = 0;
	vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
	const conf_size s = m.get_size();
	output_type tmp(s, 0);
	tmp.c.randomize(corner1, corner2, generator);
	// Zero offsets reproduce the prepared internal geometry. Draw all native
	// random values first to preserve the placement and subsequent RNG stream.
	if (prepared_torsions)
		for (auto& lig : tmp.c.ligands)
			std::fill(lig.torsions.begin(), lig.torsions.end(), 0.0);
	if (pocket_start || (counters && m.num_ligands() == 1 && m.num_flex() == 0)) {
		auto& ligand = tmp.c.ligands.front();
		vec center = m.ligand_heavy_centroid(ligand);
		if (pocket_start) {
			ligand.rigid.position += *pocket_start - center;
			if (counters)
				center = m.ligand_heavy_centroid(ligand);
		}
		if (counters) {
			counters->initial_centroid = center;
			counters->pocket_anchored = bool(pocket_start);
		}
	}
	output_type candidate(s, 0);
	fl best_e = max_fl;
	bool reference_from_save = false;
	quasi_newton quasi_newton_par(s, local_steps);
	VINA_U_FOR(step, global_steps) {
		if(increment_me)
			++(*increment_me);
		if((max_evals > 0) & (evalcount > max_evals))
			break;
		candidate.c = tmp.c;
		const bool reoriented = mutate_conf(candidate.c, m, mutation_amplitude, generator,
		                                    mixed_reorientation, separate_hunt_chain, centered_rotation, balanced_torsion);
		const int before_hunt = evalcount;
		quasi_newton_par(m, p, ig, candidate, hunt_cap, evalcount, nullptr, 0, true);
		if (counters) {
			++counters->hunt_calls;
			counters->hunt_evals += evalcount - before_hunt;
			counters->rigid_reorientations += reoriented;
			counters->mixed_accept_checks += reference_from_save;
		}
		if(step == 0 || metropolis_accept(tmp.e, candidate.e, temperature, generator)) {
			std::swap(tmp.c, candidate.c);
			tmp.e = candidate.e;
			reference_from_save = false;
			if (counters)
				++counters->accepted_steps;

			// Hunt minimization left m at the accepted conformation. The shadow
			// archive must observe it before the historical second minimization.
			if (archive)
				archive->capture(m, tmp, min_rmsd, num_saved_mins);

			// FIXME only for very promising ones
			if(tmp.e < best_e || out.size() < num_saved_mins) {
				if (archive)
					++archive->native_refinements;
				const int before_second = evalcount;
				// The separated chain keeps the accepted hunt state and its
				// energy. Reuse the trial buffer for the saved minimum.
				output_type& saved = separate_hunt_chain ? candidate : tmp;
				if (separate_hunt_chain)
					saved.c = tmp.c;
				quasi_newton_par(m, p, ig, saved, authentic_v, evalcount);
				if (counters) {
					++counters->second_calls;
					counters->second_evals += evalcount - before_second;
				}
				saved.coords = m.get_heavy_atom_movable_coords();
				add_to_output_container(out, saved, min_rmsd, num_saved_mins); // 20 - max size
				reference_from_save = !separate_hunt_chain;
				if(tmp.e < best_e)
					best_e = tmp.e;
			}
		}
	}
	if (archive)
		archive->native_evaluations = evalcount;
	VINA_CHECK(!out.empty());
	VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order
}
