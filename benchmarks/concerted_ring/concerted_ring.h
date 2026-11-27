// A single cut bond is encoded by two opposite real/dummy endpoint matches.
// Each search worker owns its model. Branch axes below refer to that worker
// and remain valid while model::set updates the existing atom frames.
#ifndef VINA_CONCERTED_RING_H
#define VINA_CONCERTED_RING_H
#include "model.h"
#include <array>

namespace concerted_ring {
using six = std::array<fl, 6>;
using ends = std::array<sz, 4>;

struct axis {
    const segment* frame;
    sz torsion;
    unsigned mask;
};

inline unsigned collect(const branch& b, const ends& atoms, sz& next,
                        std::vector<axis>& axes) {
    const sz torsion = next++;
    unsigned mask = 0;
    for (sz i = 0; i < atoms.size(); ++i)
        if (b.node.begin <= atoms[i] && atoms[i] < b.node.end)
            mask |= 1u << i;
    for (const auto& child : b.children) mask |= collect(child, atoms, next, axes);
    const unsigned p = mask & 3, q = mask >> 2;
    if (p == 1 || p == 2 || q == 1 || q == 2)
        axes.push_back({&b.node, torsion, mask});
    return mask;
}

inline bool topology(const model& m, sz index, ends& atoms, std::vector<axis>& axes) {
    const auto& lig = m.ligand_tree(index);
    sz count = 0;
    for (const auto& pair : m.get_glue_pairs()) {
        if (pair.a < lig.begin || pair.a >= lig.end || pair.b < lig.begin || pair.b >= lig.end) continue;
        if (count == 4) return false;
        atoms[count++] = pair.a;
        atoms[count++] = pair.b;
    }
    if (count != 4) return false;
    sz next = 0;
    for (const auto& child : lig.children) collect(child, atoms, next, axes);
    return true;
}

inline six geometry(const model& m, const ends& atoms, const std::vector<axis>& axes,
                    std::vector<six>& jac) {
    std::array<vec, 4> xyz;
    for (sz i = 0; i < 4; ++i) xyz[i] = m.get_coords(atoms[i]);
    six residual;
    for (sz p = 0; p < 2; ++p) {
        const vec delta = xyz[2*p] - xyz[2*p+1];
        for (sz d = 0; d < 3; ++d) residual[3*p+d] = delta[d];
    }
    for (sz j = 0; j < axes.size(); ++j) {
        const axis& a = axes[j];
        const vec& origin = a.frame->get_origin();
        for (sz p = 0; p < 2; ++p) {
            vec lever(0, 0, 0);
            if (a.mask & (1u << (2*p))) lever += xyz[2*p] - origin;
            if (a.mask & (2u << (2*p))) lever -= xyz[2*p+1] - origin;
            const vec velocity = cross_product(a.frame->get_axis(), lever);
            for (sz d = 0; d < 3; ++d) jac[j][3*p+d] = velocity[d];
        }
    }
    return residual;
}

// Solve (J J^T + lambda I)y=b. The cut bond's fixed length makes the six
// Cartesian constraints rank deficient. Relative damping also handles
// near-collinear axes. J is in angstrom/radian with dimensionless angles.
inline bool solve(const std::vector<six>& jac, six& rhs) {
    six lower[6]{};
    for (const auto& col : jac)
        for (sz i = 0; i < 6; ++i)
            for (sz j = 0; j <= i; ++j) lower[i][j] += col[i]*col[j];
    fl trace = 0;
    for (sz i = 0; i < 6; ++i) trace += lower[i][i];
    if (!(trace > 0)) return false;
    for (sz i = 0; i < 6; ++i) {
        lower[i][i] += trace * std::sqrt(epsilon_fl) / 6;
        for (sz j = 0; j <= i; ++j) {
            for (sz k = 0; k < j; ++k) lower[i][j] -= lower[i][k]*lower[j][k];
            if (i == j) {
                if (!(lower[i][i] > 0)) return false;
                lower[i][i] = std::sqrt(lower[i][i]);
            } else lower[i][j] /= lower[j][j];
        }
        for (sz j = 0; j < i; ++j) rhs[i] -= lower[i][j]*rhs[j];
        rhs[i] /= lower[i][i];
    }
    for (int i = 5; i >= 0; --i) {
        for (int j = i+1; j < 6; ++j) rhs[i] -= lower[j][i]*rhs[j];
        rhs[i] /= lower[i][i];
    }
    return true;
}

inline fl dot(const six& a, const six& b) {
    fl value = 0;
    for (sz i = 0; i < 6; ++i) value += a[i]*b[i];
    return value;
}

inline bool propose(conf& c, model& m, sz index, sz driver, fl draw, fl amplitude) {
    ends atoms;
    std::vector<axis> axes;
    if (!topology(m, index, atoms, axes)) return false;
    const auto chosen = std::find_if(axes.begin(), axes.end(),
                                     [driver](const axis& a){ return a.torsion == driver; });
    if (chosen == axes.end()) return false;
    const sz selected = chosen - axes.begin();
    m.set(c);
    std::vector<six> jac(axes.size());
    const six initial = geometry(m, atoms, axes, jac);
    // Only move from an already closed ring. Unsuccessful proposals are nulls.
    constexpr fl tolerance_sq = 1e-12;
    if (dot(initial, initial) > tolerance_sq) return false;
    const fl radius = m.gyration_radius(index);
    if (!(radius > epsilon_fl)) return false;
    const fl angle = draw * std::min(fl(1), amplitude/(pi*radius));
    six rhs = jac[selected];
    if (!solve(jac, rhs)) return false;
    auto& angles = c.ligands[index].torsions;
    const flv saved = angles;
    for (sz j = 0; j < axes.size(); ++j)
        angles[axes[j].torsion] += angle*((j == selected ? 1 : 0) - dot(jac[j],rhs));

    fl previous = max_fl;
    for (unsigned iteration = 0; iteration <= 6; ++iteration) {
        m.set(c);
        rhs = geometry(m, atoms, axes, jac);
        const fl residual_sq = dot(rhs, rhs);
        if (residual_sq <= tolerance_sq) return true;
        if (!(residual_sq < previous) || iteration == 6 || !solve(jac, rhs)) break;
        previous = residual_sq;
        for (sz j = 0; j < axes.size(); ++j) angles[axes[j].torsion] -= dot(jac[j], rhs);
    }
    angles = saved;
    m.set(c);
    return false;
}
}
#endif
