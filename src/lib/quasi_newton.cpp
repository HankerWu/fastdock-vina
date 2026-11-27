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

// Modified for FastDock Vina: reuse coordinates from the final BFGS evaluation.

#include "quasi_newton.h"

namespace {
boost::optional<bfgs_rigid_metric> centroid_metric(model& m, const conf& base) {
    m.set(base);
    const vec centroid = m.ligand_heavy_centroid(0);
    const auto atoms = m.get_chemical_heavy_atom_movable_coords();
    fl squared_radius = 0;
    for (const auto& atom : atoms)
        squared_radius += vec_distance_sqr(atom, centroid);
    if (squared_radius > 0) {
        const fl scale = std::sqrt(2 * squared_radius / (3 * atoms.size()));
        if (std::isfinite(scale) && scale > epsilon_fl)
            return bfgs_rigid_metric{centroid - base.ligands.front().rigid.position, 1 / scale};
    }
    return boost::none;
}

} // namespace

void quasi_newton::operator()(model& m, const precalculate_byatom& p, const igrid& ig,
                              output_type& out, const vec& v, int& evalcount, pose_template* reference, int max_evals, search_backtracking* monitor) {
    const auto rigid = rigid_metric ? centroid_metric(m, out.c) : boost::none;
    auto evaluate = [&](const conf& c, flv& gradient) {
        m.set(c);
        const fl energy = m.eval_deriv(p, ig, v, gradient, reference);
        if (rigid)
            rigid->to_centroid_gradient(gradient);
        return energy;
    };
    VINA_CHECK(!geometry_metric || strict_armijo);
    const flv inverse = geometry_metric ? m.geometry_inverse_metric(out.c) : flv();
    const bfgs_result result = strict_armijo
        ? workspace.minimize<true>(evaluate, out.c, max_steps, evalcount,
                                   geometry_metric ? &inverse : nullptr, max_evals, rigid ? &*rigid : nullptr)
        : workspace.minimize(evaluate, out.c, max_steps, evalcount, nullptr, 0, nullptr, monitor);
    if (result.restore_coordinates)
        m.set(out.c);
    out.e = result.energy;
}
