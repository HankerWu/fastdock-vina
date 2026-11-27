// Verify the explicit rigid-receptor table's one-sided endpoint, unchanged
// midpoint samples, and fallback under custom cutoffs, switching and Vinardo.
#include "precalculate.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
void close(fl a, fl b, fl tol) {
    if (!std::isfinite(a) || std::abs(a-b) > tol)
        throw std::runtime_error("shifted receptor table mismatch");
}

void check(const ScoringFunction& sf, sz a, sz b, bool shifted) {
    std::bitset<XS_TYPE_SIZE> lhs, rhs;
    lhs.set(a);
    rhs.set(b);
    precalculate p(sf.get_cutoff());
    p.populate(sf, lhs, rhs);
    const fl cutoff = sf.get_cutoff(), rc2 = cutoff*cutoff;
    const sz index = triangular_matrix_index_permissive(XS_TYPE_SIZE, a, b);
    const fl limit = shifted ? sf.eval(a,b,std::nextafter(cutoff, 0.0)) : 0;
    const auto potential = [&](fl r){ return sf.eval(a,b,r); };
    const auto gridpoint = [](sz i){ return std::sqrt(i/32.); };
    for (sz i = 0; i < sz(32*rc2); ++i) {
        const fl left = potential(gridpoint(i)), right = potential(gridpoint(i+1));
        const bool final = shifted && (i+1 == sz(32*rc2));
        const fl delta = (final ? limit : right) - left;
        const fl r2 = (i + 0.375) / 32;
        const auto eg = p.eval_deriv(index, r2);
        close(p.eval_fast(index,r2), (left+right)/2, 0);
        close(eg.first, left-limit+0.375*delta, 2e-13);
        close(eg.second, 64*delta, 2e-13);
        const fl r = std::sqrt(r2), h = 1e-6;
        const fl fd = (p.eval_deriv(index,sqr(r+h)).first-p.eval_deriv(index,sqr(r-h)).first)/(2*h);
        close(r * eg.second, fd, 4e-8);
    }
    if (shifted) {
        const fl r = 7.999335847431335, h = 1e-6;
        const auto eg = p.eval_deriv(index,r*r);
        const fl direct = (potential(r+h)-potential(r-h))/(2*h);
        close(p.eval_deriv(index,rc2).first,0,0);
        close(p.eval_deriv(index,rc2).second,0,0);
        close(eg.second*r, direct, 1e-5);
        std::cout << a << ',' << b << ',' << r << ',' << eg.first << ',' << limit
                  << ',' << r*eg.second << ',' << direct << '\n';
    }
}
} // namespace
int main() {
    const flv w{-0.035579,-0.005156,0.840245,-0.035069,-0.587439,50.,5*0.05846/0.1-1};
    const ScoringFunction vina(SF_VINA,w), custom(SF_VINA,w,6.234), switched(SF_VINA,w,8.,true);
    const ScoringFunction vinardo(SF_VINARDO, flv{-0.045,0.8,-0.035,-0.6,50.,5*0.05846/0.1-1});
    std::cout << std::setprecision(16) << "xs1,xs2,r,shifted_energy,left_limit,radial_derivative,direct_derivative\n";
    for (auto ab : {std::make_pair(XS_TYPE_C_P,XS_TYPE_C_P), std::make_pair(XS_TYPE_N_D,XS_TYPE_O_A),
                    std::make_pair(XS_TYPE_O_A,XS_TYPE_O_A)}) {
        check(vina, ab.first, ab.second, true);
        check(custom, ab.first, ab.second, false);
        check(switched, ab.first, ab.second, false);
        check(vinardo, ab.first, ab.second, false);
    }
}
