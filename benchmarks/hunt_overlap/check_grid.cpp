// Verify the research interpolant itself, including the box-wall chain rule.
#include "grid.h"
#include <iostream>
#include <stdexcept>

void close(fl observed, fl expected, fl tol) {
    if (!std::isfinite(observed) || std::abs(observed - expected) > tol)
        throw std::runtime_error("grid check: " + std::to_string(observed) + " vs " + std::to_string(expected));
}

int main() {
    grid_dims d;
    for (auto& axis : d) { axis.begin = 0; axis.end = 2; axis.n_voxels = 2; }
    grid field(d);
    for (sz z = 0; z < 3; ++z)
        for (sz y = 0; y < 3; ++y)
            for (sz x = 0; x < 3; ++x)
                field.m_data(x,y,z) = 20.0 * x + 2.0 * y - 4.0 * z;
    // Negative, weak positive, exactly at the join, high, and outside a face.
    for (const auto& at : {vec(0.01,0.5,0.5), vec(0.1,0.5,0.5), vec(0.55,0.5,0.5),
                           vec(1.4,0.7,0.6), vec(2.3,0.7,0.6)}) {
        const fl e = 20 * std::min(at[0],2.0) + 2 * at[1] - 4 * at[2];
        const fl cap = 10, wall = 17;
        const fl expected = e <= 0 ? e : (e <= cap ? cap * e / (cap + e) : (e + cap) / 4);
        const fl slope = e <= 0 ? 1 : (e <= cap ? sqr(cap/(cap+e)) : 0.25);
        vec gradient;
        const fl value = field.evaluate_hunt_tangent(at,wall,cap,&gradient);
        close(value, expected + wall * std::max(at[0]-2.0,0.0), 1e-12);
        close(value, field.evaluate_hunt_tangent(at,wall,cap,nullptr), 1e-12);
        close(gradient[0], at[0] > 2 ? wall : 20*slope, 1e-12);
        close(gradient[1], 2*slope, 1e-12);
        close(gradient[2], -4*slope, 1e-12);
        for (sz k = 0; k < 3; ++k) {
            vec a=at, b=at;
            const fl h = 1e-6;
            a[k] += h; b[k] -= h;
            const fl fd = (field.evaluate_hunt_tangent(a,wall,cap,nullptr)
                         - field.evaluate_hunt_tangent(b,wall,cap,nullptr)) / (2*h);
            close(gradient[k], fd, 4e-6);
        }
    }
    std::cout << "passed: five energy/gradient cases and 15 central differences\n";
}
