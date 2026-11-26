#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "precalculate.h"

namespace {
void require_close(double observed, double expected, double tolerance, const char* label) {
    if (!std::isfinite(observed) || std::abs(observed - expected) > tolerance) {
        std::cerr << label << ": observed " << observed << ", expected " << expected << '\n';
        throw std::runtime_error(label);
    }
}

// Independent polynomial in distance from the onset, not from the cutoff.
double envelope(double distance, double cutoff) {
    if (distance <= cutoff - 1) return 1;
    if (distance >= cutoff) return 0;
    const double t = distance - (cutoff - 1);
    return 1 - t*t*t*(10 + t*(-15 + 6*t));
}

void check_energy(const ScoringFunction& original, const ScoringFunction& switched,
                  double cutoff, sz first, sz second) {
    atom a, b;
    a.xs = first;
    b.xs = second;
    for (double offset : {-1.5, -1.0, -0.95, -0.75, -0.5, -0.1, -0.001, 0., 0.3}) {
        const double r = cutoff + offset;
        const double expected = original.eval(first, second, r) * envelope(r, cutoff);
        require_close(switched.eval(first, second, r), expected, 1e-12, "XS type pair");
        require_close(switched.eval(a, b, r), expected, 1e-12, "atom pair");
    }
}

void check_table(const ScoringFunction& sf, bool switched) {
    const auto distances = sample_distances(32, 64);
    const precalculate_element table(distances, max_fl, [&](double r) {
        return sf.eval(XS_TYPE_C_H, XS_TYPE_C_H, r);
    });
    const auto evaluate = [&](double r) { return table.eval_deriv(32*r*r); };
    for (double r : {7.997070, 7.999023}) {
        const double h = 1.e-7;
        const double slope = (evaluate(r+h).first - evaluate(r-h).first) / (2*h);
        const double returned = r * evaluate(r).second;
        const double residual = std::abs(returned - slope);
        std::cout << (switched ? "switched" : "native") << " r=" << r
                  << " dE/dr=" << slope << " returned=" << returned
                  << " residual=" << residual << '\n';
        if (switched) {
            require_close(residual, 0., 1.e-6, "switched terminal-bin force residual");
            require_close(returned, 0., 2.e-6, "switched terminal-bin radial force");
        } else if (!(residual > 0.4)) {
            throw std::runtime_error("Native control no longer reproduces cutoff artifact");
        }
    }
}
} // namespace

int main() {
    const flv weights{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439,
                      50., 5 * 0.05846 / 0.1 - 1};
    std::cout << std::setprecision(9);
    for (double cutoff : {8.0, 6.234, 3.5}) {
        ScoringFunction original(SF_VINA, weights, cutoff);
        ScoringFunction switched(SF_VINA, weights, cutoff, true);
        check_energy(original, switched, cutoff, XS_TYPE_C_H, XS_TYPE_C_H);
        check_energy(original, switched, cutoff, XS_TYPE_N_D, XS_TYPE_O_A);
        if (cutoff == 8.) {
            check_table(original, false);
            check_table(switched, true);
        }
    }

    // The long-range closure restraint serves a separate purpose from ordinary atom contacts.
    ScoringFunction glue(SF_VINA, flv{0., 0., 0., 0., 0., 50., 0.}, 8., true);
    atom anchor, pseudoatom;
    anchor.xs = XS_TYPE_C_H_CG0;
    pseudoatom.xs = XS_TYPE_G0;
    for (double r : {7.5, 8., 12., 19.9, 20.}) {
        const double expected = r < 20. ? 50.*r : 0.;
        require_close(glue.eval(anchor, pseudoatom, r), expected, 1.e-12, "atom-pair closure");
        require_close(glue.eval(anchor.xs, pseudoatom.xs, r), expected, 1.e-12, "type-pair closure");
    }
}
