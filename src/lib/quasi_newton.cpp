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

void quasi_newton::operator()(model& m, const precalculate_byatom& p, const igrid& ig,
                              output_type& out, const vec& v, int& evalcount) {
    auto evaluate = [&](const conf& c, change& gradient) {
        m.set(c);
        return m.eval_deriv(p, ig, v, gradient);
    };
    const bfgs_result result = workspace.minimize(evaluate, out.c, max_steps, evalcount);
    if (result.reverted)
        m.set(out.c);
    out.e = result.energy;
}
