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

#include "coords.h"

fl rmsd_upper_bound(const vecv& a, const vecv& b) {
	VINA_CHECK(a.size() == b.size());
	fl acc = 0;
	VINA_FOR_IN(i, a) 
		acc += vec_distance_sqr(a[i], b[i]);
	return (a.size() > 0) ? std::sqrt(acc / a.size()) : 0;
}

std::pair<sz, fl> find_closest(const vecv& a, const output_container& b) {
	std::pair<sz, fl> tmp(b.size(), max_fl);
	VINA_FOR_IN(i, b) {
		fl res = rmsd_upper_bound(a, b[i].coords);
		if(i == 0 || res < tmp.second)
			tmp = std::pair<sz, fl>(i, res);
	}
	return tmp;
}

void add_to_output_container(output_container& out, const output_type& t, fl min_rmsd, sz max_size) {
	VINA_CHECK(out.empty() || t.coords.size() == out.front().coords.size());
	if (max_size == 0 || (out.size() == max_size && !(t.e < out.back().e)))
		return;
	const auto overlaps = [&](const output_type& resident) {
		return rmsd_upper_bound(t.coords, resident.coords) < min_rmsd;
	};
	// Residents are already separated and sorted by energy. Earlier residents
	// win ties. A single better or tied neighbour blocks the newcomer.
	auto position = out.begin();
	while (position != out.end() && !(t.e < position->e)) {
		if (overlaps(*position)) return;
		++position;
	}
	// A new representative can bridge two distinct basins. Remove every
	// lower-priority neighbour to maintain pairwise separation after insertion.
	const auto inserted = out.insert(position, new output_type(t));
	out.erase_if(inserted + 1, out.end(), overlaps);
	if (out.size() > max_size) out.pop_back();
}

output_container select_distinct_poses(const output_container& poses, fl min_rmsd, sz max_size) {
	std::vector<const output_type*> by_energy;
	by_energy.reserve(poses.size());
	for(const auto& pose : poses)
		by_energy.push_back(&pose);
	std::stable_sort(by_energy.begin(), by_energy.end(), [](const output_type* a, const output_type* b) {
		return a->e < b->e;
	});

	output_container selected;
	selected.reserve(std::min(max_size, poses.size()));
	for(const auto* pose : by_energy) {
		if(selected.size() == max_size)
			break;
		if(find_closest(pose->coords, selected).second >= min_rmsd)
			selected.push_back(new output_type(*pose));
	}
	return selected;
}
