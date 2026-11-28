// Use the public candidate header to test the tangent and accepted-step contract.
#include "bfgs.h"
#include <iostream>
#include <stdexcept>

void check(fl value, fl expected) {
    if (!(std::abs(value - expected) < 1e-12))
        throw std::runtime_error("hunt tangent differs from analytical result");
}

int main() {
    const conf_size layout{{2, 1}, {0, 2}};
    conf state(layout);
    state.set_to_null();
    flv p(layout.num_degrees_of_freedom(), 0), g(p.size(), 0);
    p[0] = 17;
    p[9] = -4;
    g[0] = -1;
    g[9] = 1;
    auto q = p;
    limit_hunt_direction(state, g, q);
    if (q != p) throw std::runtime_error("pure translation changed");
    p[3] = 0.1;
    p[16] = -0.5;
    q = p;
    limit_hunt_direction(state, g, q);
    if (q != p) throw std::runtime_error("short angular step changed");

    // Second rigid body's rotation and terminal flexible-residue torsion set the bound.
    p[11] = 12;
    p[12] = 5;
    g[11] = -6;
    q = p;
    limit_hunt_direction(state, g, q);
    check(q[0], 17); check(q[9], -4);
    check(std::hypot(q[11], q[12]), pi); check(q[16], -0.5 * pi / 13);
    p[16] = -40;
    q = p;
    limit_hunt_direction(state, g, q);
    check(q[16], -pi); check(q[11], 12 * pi / 40);
    check(scalar_product(q, g), (pi / 40) * (-72) - 21);

    // Positive combined translational work with downhill total slope uses full scaling.
    g[0] = 1;
    q = p;
    limit_hunt_direction(state, g, q);
    for (sz i = 0; i < q.size(); ++i) check(q[i], p[i] * (pi / 40));
    check(scalar_product(q, g), (-72 + 17 - 4) * (pi / 40));

    // Run through BFGS. First accept the real half step; reject all later trials.
    const conf_size single{{1}, {}};
    conf origin(single);
    origin.set_to_null();
    int evaluations = 0;
    bfgs_workspace workspace(single);
    auto landscape = [&](const conf& c, flv& grad) {
        const fl shift = c.ligands.front().rigid.position[0];
        grad[0] = 8;
        grad[6] = 20;
        if (evaluations == 0) return fl(0);
        if (evaluations == 2) {
            check(shift, -4); check(c.ligands[0].torsions[0], -pi / 2);
            return fl(-0.5);
        }
        return fl(1);
    };
    auto result = workspace.minimize(landscape, origin, 3, evaluations, nullptr, 0, nullptr, true);
    check(result.energy, -0.5); check(origin.ligands[0].rigid.position[0], -4);
    check(origin.ligands[0].torsions[0], -pi / 2);
    if (evaluations != 13 || !result.restore_coordinates)
        throw std::runtime_error("hunt did not restore accepted half-step after ten failures");
    std::cout << "pure translation, short arc, multiple ligands/residues, mixed descent, fallback, half-step and failure restoration passed\n";
}
