#include "bfgs.h"
#include <stdexcept>

static void check(fl actual, fl expected) {
    if (!(std::abs(actual - expected) < 1e-13))
        throw std::runtime_error("principal chart trial mismatch");
}

int main() {
    conf_size layout;
    layout.ligands = {2, 1};
    layout.flex = {0, 2};
    conf state(layout);
    flv p(layout.num_degrees_of_freedom(), 0);
    p[0] = 1000;
    p[9] = -400;
    check(principal_angle_step(state, p), 1);
    p[3] = 3;
    p[4] = -4;
    check(principal_angle_step(state, p), pi / 5);
    p[11] = 12;
    p[12] = 5;
    check(principal_angle_step(state, p), pi / 13);
    p[6] = -22;
    check(principal_angle_step(state, p), pi / 22);
    p[14] = 30;
    check(principal_angle_step(state, p), pi / 30);
    p[16] = -40;
    check(principal_angle_step(state, p), pi / 40);

    // Exercise the same acceptance calculation that grid search calls, including
    // one rejected full trial. The returned alpha must describe the actual step.
    state.set_to_null();
    conf trial(layout);
    flv trial_gradient(p.size(), 0);
    fl trial_energy = 0;
    int evaluations = 0;
    auto objective = [&](const conf&, flv&) { return evaluations == 0 ? fl(1) : fl(-1); };
    fl alpha = line_search(objective, state, 0, p, -1, trial, trial_gradient,
                           trial_energy, evaluations, 0, true);
    check(alpha, pi / 80);
    check(trial.flex.back().torsions[1], -pi / 2);
    check(trial.ligands[0].rigid.position[0], alpha * p[0]);
    if (evaluations != 2) throw std::runtime_error("wrong number of evaluated trials");
    evaluations = 0;
    auto reject = [](const conf&, flv&) { return fl(1); };
    check(line_search(reject, state, 0, p, -1, trial, trial_gradient,
                      trial_energy, evaluations, 0, true), 0);
    check(trial.flex.back().torsions[1], -pi / 512);
    if (evaluations != 10) throw std::runtime_error("wrong exhausted trial count");
}
