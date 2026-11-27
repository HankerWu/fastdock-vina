// Compare the production table's force both to its returned energy and to the
// underlying Vina potential. The two references differ at finite resolution.
#include "precalculate.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {

struct errors {
    size_t n = 0;
    double sum_table = 0, max_table = 0, sum_direct = 0, max_direct = 0;

    void add(double force, double table, double direct) {
        ++n;
        const double dt = std::abs(force - table), dd = std::abs(force - direct);
        sum_table += dt;
        max_table = std::max(max_table, dt);
        sum_direct += dd;
        max_direct = std::max(max_direct, dd);
    }
};

void require_close(double observed, double expected, double tolerance) {
    if (!std::isfinite(observed) || std::abs(observed - expected) > tolerance)
        throw std::runtime_error("Pair lookup disagrees with reference");
}

void audit(double factor, sz a, sz b, const ScoringFunction& sf) {
    const double cutoff = sf.get_cutoff();
    const auto nodes = sample_distances(factor, cutoff * cutoff);
    const auto direct = [&](double r) { return sf.eval(a, b, r); };
    const precalculate_element lookup(nodes, max_fl, direct);
    const auto sample = [&](double q) {
#ifdef PAIR_AUDIT_NATIVE
        return lookup.eval_deriv(q);
#else
        return lookup.eval_deriv(q, factor);
#endif
    };
    const auto evaluate = [&](double r) { return sample(factor * r * r); };
    const double h = 1e-6;
    std::array<errors, 4> groups;
    const char* labels[] = {"overlap_lt2A", "contact_2_to_6A", "tail", "last_two_bins"};
    const size_t bins = static_cast<size_t>(factor * cutoff * cutoff);
    for (size_t i = 0; i < bins; ++i) {
        const double e0 = direct(nodes[i]), e1 = direct(nodes[i + 1]);
        require_close(lookup.fast_sample(i), (e0 + e1) / 2, 0);
        for (double fraction : {0.125, 0.5, 0.875}) {
            const double r = std::sqrt((i + fraction) / factor);
            const auto pair = evaluate(r);
#ifdef PAIR_AUDIT_MONOTONE
            if (pair.first < std::min(e0, e1) - 2e-12 || pair.first > std::max(e0, e1) + 2e-12 ||
                (e1 - e0) * pair.second < -2e-12)
                throw std::runtime_error("Hermite interpolation overshot or reversed");
#else
            require_close(pair.first, e0 + fraction * (e1 - e0), 5e-12);
#endif
            const double numerical = (evaluate(r + h).first - evaluate(r - h).first) / (2*h);
            const double from_potential = (direct(r + h) - direct(r - h)) / (2*h);
            const double returned = r * pair.second;
            const size_t group = r < 2 ? 0 : r < 6 ? 1 : i + 2 >= bins ? 3 : 2;
            groups[group].add(returned, numerical, from_potential);
#if !defined(PAIR_AUDIT_NATIVE)
            // Finite difference of a cubic has a truncation error growing
            // with the third derivative near the hard cutoff and r=0.
            require_close(returned, numerical, 3e-5);
#endif
        }
    }
    for (size_t i = 0; i < groups.size(); ++i) {
        const auto& g = groups[i];
        std::cout << factor << ',' << a << ',' << b << ',' << labels[i] << ',' << g.n
                  << ',' << g.sum_table/g.n << ',' << g.max_table << ','
                  << g.sum_direct/g.n << ',' << g.max_direct << '\n';
    }
#ifdef PAIR_AUDIT_MONOTONE
    for (size_t i = 1; i <= bins; ++i) {
        const auto at = sample(i);
        require_close(at.first, direct(nodes[i]), 0);
        const auto before = sample(std::nextafter(double(i), 0.0));
        require_close(at.first, before.first, 3e-12);
        require_close(at.second, before.second, 3e-9);
    }
#elif !defined(PAIR_AUDIT_NATIVE)
    // At nodes the energy can have a kink. Query the integer bin index
    // directly to avoid rounding a representable right boundary to its left.
    const size_t node = 100;
    const auto at_node = sample(node);
    require_close(at_node.first, direct(nodes[node]), 0);
    require_close(at_node.second, 2 * factor * (direct(nodes[node + 1]) - direct(nodes[node])), 0);
    require_close(0 * evaluate(0).second, 0, 0);
#endif
}
} // namespace

int main() {
    const flv weights{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439,
                      50., 5 * 0.05846 / 0.1 - 1};
    ScoringFunction sf(SF_VINA, weights);
    std::cout << std::setprecision(12)
              << "factor,xs1,xs2,region,n,mean_abs_dE_error,max_abs_dE_error,"
                 "mean_abs_direct_error,max_abs_direct_error\n";
    for (double factor : {16., 32., 64.}) {
        audit(factor, XS_TYPE_C_H, XS_TYPE_C_H, sf);
        audit(factor, XS_TYPE_N_D, XS_TYPE_O_A, sf);
        audit(factor, XS_TYPE_O_A, XS_TYPE_O_A, sf);
    }
}
