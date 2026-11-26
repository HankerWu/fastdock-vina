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
// Diagonal of the pullback of the Cartesian RMS displacement metric, J^T J/N.
// The central coordinate probes do not evaluate the scoring function. This is
// an experimental preconditioner, not a force-field energy or entropy term.
flv geometry_inverse_metric(model& m, const conf& base) {
    const conf_size size = m.get_size();
    VINA_CHECK(size.ligands.size() == 1 && size.flex.empty());
    flv inverse(size.num_degrees_of_freedom(), 1);
    flv delta(inverse.size(), 0);
    conf trial(size);
    const fl eps = 1e-4;
    // Translational entries are analytically one. Angular entries can have
    // very different lever arms and are estimated in the actual conf chart.
    for (sz k = 3; k < inverse.size(); ++k) {
        delta[k] = 1;
        trial.set_incremented(base, delta, eps);
        m.set(trial);
        const vecv plus = m.get_chemical_heavy_atom_movable_coords();
        trial.set_incremented(base, delta, -eps);
        m.set(trial);
        const vecv minus = m.get_chemical_heavy_atom_movable_coords();
        VINA_CHECK(!plus.empty() && plus.size() == minus.size());
        fl sum = 0;
        VINA_FOR_IN(i, plus) sum += vec_distance_sqr(plus[i], minus[i]);
        const fl g = sum / (4 * eps * eps * plus.size());
        VINA_CHECK(std::isfinite(g) && g >= 0);
        // Bound the inverse to [1/64, 4] to regularize weak/large DOFs.
        inverse[k] = 1 / std::min(fl(64), std::max(fl(0.25), g));
        delta[k] = 0;
    }
    m.set(base);
    return inverse;
}
} // namespace

void quasi_newton::operator()(model& m, const precalculate_byatom& p, const igrid& ig,
                              output_type& out, const vec& v, int& evalcount, pose_template* reference) {
    auto evaluate = [&](const conf& c, flv& gradient) {
        m.set(c);
        return m.eval_deriv(p, ig, v, gradient, reference);
    };
    VINA_CHECK(!geometry_metric || strict_armijo);
    const flv inverse = geometry_metric ? geometry_inverse_metric(m, out.c) : flv();
    const bfgs_result result = strict_armijo
        ? workspace.minimize<true>(evaluate, out.c, max_steps, evalcount,
                                   geometry_metric ? &inverse : nullptr)
        : workspace.minimize(evaluate, out.c, max_steps, evalcount);
    if (result.restore_coordinates)
        m.set(out.c);
    out.e = result.energy;
}
