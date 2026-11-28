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
}
