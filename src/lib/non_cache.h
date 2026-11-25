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

#ifndef VINA_NON_CACHE_H
#define VINA_NON_CACHE_H

#include "igrid.h"
#include "szv_grid.h"
#include "precalculate.h"

class net_polar;
struct net_polar_terms;

struct non_cache : public igrid {
	non_cache() = default;
	non_cache(const model& receptor, const grid_dims& bounds,
	          std::shared_ptr<const ScoringFunction> scoring, fl slope);
	void prepare_for_model(const model& m);
	void prepare_for_model(const model& m, const ScoringFunction& sf, precalculate& tables) const;
	void populate_maps(cache& maps, const model& receptor, const std::bitset<XS_TYPE_SIZE>& types,
	                   bool retain_tables);
	virtual fl eval      (const model& m, fl v) const; // needs m.coords // clean up
	fl eval(const model& m, fl v, const precalculate& scores) const;
	fl directional_eval(const model& m, fl v, const ScoringFunction& sf,
	                    const net_polar& polar, net_polar_terms& terms) const;
	virtual fl eval_intra(      model& m, fl v) const;
	virtual fl eval_deriv(      model& m, fl v) const; // needs m.coords, sets m.minus_forces // clean up
	bool within(const model& m, fl margin = 0.0001) const;
	fl slope = 0;
private:
	template<bool Directional>
	fl eval_impl(const model& m, fl v, const precalculate& scores, bool flex_only,
	             const ScoringFunction* sf, const net_polar* polar, net_polar_terms* terms) const;
	szv_grid sgrid;
	grid_dims gd;
	precalculate p;
	std::shared_ptr<const ScoringFunction> scoring;
};

#endif

