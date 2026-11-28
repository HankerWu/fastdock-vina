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

// Modified for FastDock Vina: reuse the BFGS workspace across local minimizations.

#ifndef VINA_QUASI_NEWTON_H
#define VINA_QUASI_NEWTON_H

#include "bfgs.h"
#include "model.h"

class quasi_newton {
    bfgs_workspace workspace;
    unsigned max_steps;

public:
    bool strict_armijo = false;
    bool geometry_metric = false;
    bool rigid_metric = false;
    quasi_newton(const conf_size& size, unsigned steps) : workspace(size), max_steps(steps) {}
    void operator()(model& m, const precalculate_byatom& p, const igrid& ig,
                    output_type& out, const vec& v, int& evalcount, pose_template* reference = nullptr,
                    int max_evals = 0, fl final_trial_spacing = 0);
};

#endif
