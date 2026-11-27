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

#ifndef VINA_COORDS_H
#define VINA_COORDS_H

#include <boost/container/small_vector.hpp>
#include "conf.h"
#include "atom.h" // for atomv

fl rmsd_upper_bound(const vecv& a, const vecv& b);
std::pair<sz, fl> find_closest(const vecv& a, const output_container& b);
// Update a sorted, separated online pool using the current residents and one arrival.
void add_to_output_container(output_container& out, const output_type& t, fl min_rmsd, sz max_size);
using pose_indices = boost::container::small_vector<sz, 9>;
// Apply the requested energy window before testing distances in existing rank order.
pose_indices ranked_output_indices(const output_container& poses, sz limit, fl min_rmsd, boost::optional<fl> max_energy);

// Keep the lowest-energy poses at least min_rmsd apart, breaking score ties by input order.
output_container select_distinct_poses(const output_container& poses, fl min_rmsd, sz max_size);


#endif
