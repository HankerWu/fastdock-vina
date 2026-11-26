// Build with: c++ -std=c++14 -O2 -I src/lib benchmarks/rigid_window/analytic.cpp -o /tmp/rigid-analytic
#include "printable_box.h"
#include "rigid_window.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
void close_to(fl actual, fl expected, fl tolerance, const char *label) {
	if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
		std::cerr << label << " observed=" << actual << " expected=" << expected << '\n';
		throw std::runtime_error(label);
	}
}

void integration(unsigned count) {
	const auto nodes = rigid_window::make_samples(count);
	rigid_window::integral flat, harmonic, large_positive;
	fl sum = 0, second_r = 0, first_theta = 0, cross_moment = 0;
	const long double a = rigid_window::angle;
	const long double expected_theta = (a * a / 2 - a * std::sin(a) + 1 - std::cos(a)) / (a - std::sin(a));
	// A uniform 3D ball has E[(r/R)^2]=3/5, independent of Haar angle.
	for (const auto &q : nodes) {
		const fl r2 = sqr(q.shift);
		const fl theta = std::atan2(q.sine, q.cosine);
		close_to(q.axis.norm(), 1, 4e-15, "axis unit norm");
		if (r2 > 1. + 1e-12 || theta > a + 1e-12 || q.weight <= 0 || q.weight > 3.0000000001)
			throw std::runtime_error("outside proposal support");
		sum += q.weight;
		second_r += q.weight * r2;
		first_theta += q.weight * theta;
		cross_moment += q.weight * r2 * theta;
		flat.add(0, q.weight);
		large_positive.add(5000, q.weight);
		harmonic.add(2 * rigid_window::thermal * r2, q.weight);
	}
	const fl c = 2.;
	const fl harmonic_integral =
		3 * (std::sqrt(pi) * std::erf(std::sqrt(c)) / (4 * std::pow(c, 1.5)) - std::exp(-c) / (2 * c));
	const fl expected_harmonic = -rigid_window::thermal * std::log(harmonic_integral);
	close_to(flat.correction(), 0, 2e-12, "flat potential");
	close_to(large_positive.correction(), 5000, 2e-9, "constant energy offset");
	close_to(sum / nodes.size(), 1., 0.003, "absolute proposal normalization");
	close_to(second_r / sum, 0.6, 0.0015, "mean squared radius");
	close_to(first_theta / sum, expected_theta, 0.0006, "mean Haar angle");
	close_to(cross_moment / sum, 0.6 * expected_theta, 0.0004, "translation-rotation independence");
	close_to(harmonic.correction(), expected_harmonic, 0.0015, "analytic radial free energy");
	std::cout << "nodes=" << nodes.size() << " mass/N=" << sum / nodes.size() << " mean_r2=" << second_r / sum
			  << " mean_theta=" << first_theta / sum << " analytic_theta=" << expected_theta
			  << " harmonic=" << harmonic.correction() << " analytic_harmonic=" << expected_harmonic
			  << " flat=" << flat.correction() << " flat_ess=" << flat.effective_samples() << '\n';
}

void printable_coordinates() {
	// A pose can be wider than the naive "shrink each side by half a milli-Angstrom" box.
	const vec lo(0.0004, -2.0004, -1), hi(10.0004, 3.0004, 1);
	const vecv points = {vec(0.00039, -2.00041, 0), vec(9.99990, 2.9999, 0.1)};
	const auto shift = printable_box_shift(points, lo, hi);
	if (!(shift[0] > 0 && shift[0] < 0.00012 && shift[1] > 0 && shift[2] == 0))
		throw std::runtime_error("feasible narrow margin was not repaired");
	for (const auto &p : points)
		for (sz d = 0; d < 3; ++d) {
			const fl actual = p[d] + shift[d];
			std::ostringstream text;
			text << std::fixed << std::setprecision(3) << actual;
			const fl printed = std::stod(text.str());
			if (actual < lo[d] - 1e-14 || actual > hi[d] + 1e-14 || printed < lo[d] || printed > hi[d])
				throw std::runtime_error("raw or printed coordinate remained outside the wall");
		}
	close_to(sqr(printable_box_shift({vec(2, 0, 0), vec(4, 2, 0)}, lo, hi)), 0, 0, "interior pose moved");
	close_to(sqr(printable_box_shift({vec(-1, 0, 0)}, lo, hi)), 0, 0, "genuine exterior pose moved");
	close_to(sqr(printable_box_shift({vec(0.00039, 0, 0), vec(10.0006, 0, 0)}, lo, hi)), 0, 0,
			 "impossible rigid pose moved");
	close_to(sqr(printable_box_shift({}, lo, hi)), 0, 0, "empty pose moved");
}

void frames() {
	auto centroid = [](const vecv &v) {
		vec c(0, 0, 0);
		for (const auto &p : v)
			c += p;
		return (1. / v.size()) * c;
	};
	auto rotate = [](const vec &p) { return vec(-p[2], p[0], -p[1]); };
	const vec offset(8, -11, 4);
	const std::vector<vecv> shapes = {{vec(0, 0, 0)},
									  {vec(0, 0, 0), vec(0, 0, 0)},
									  {vec(1, 2, 3), vec(3, 6, 9), vec(-1, -2, -3)},
									  {vec(0, 0, 0), vec(2, 0.2, 0), vec(-1, 1, 0), vec(1, -0.1, 3)}};
	for (const auto &shape : shapes) {
		const auto F = rigid_window::frame(shape, centroid(shape));
		close_to(F.x.norm(), 1, 2e-14, "frame X");
		close_to(F.y.norm(), 1, 2e-14, "frame Y");
		close_to(F.z.norm(), 1, 2e-14, "frame Z");
		close_to(F.x * F.y, 0, 2e-14, "frame orthogonality");
		close_to(F.z * cross_product(F.x, F.y), 1, 2e-14, "frame handedness");
		if (shape.size() != 4)
			continue;
		vecv other;
		for (const auto &p : shape)
			other.push_back(rotate(p) + offset);
		const auto other_F = rigid_window::frame(other, centroid(other));
		const vec axis(0.2, 0.7, -0.4);
		close_to((rotate(F.apply(axis)) - other_F.apply(axis)).norm(), 0, 2e-14, "rigid covariance");
	}
}

} // namespace
int main() {
	try {
		std::cout << std::setprecision(10);
		for (const unsigned n : {2048U, 8192U})
			integration(n);
		frames();
		printable_coordinates();
		for (const unsigned n : {0U, 8U, 17U, 131072U}) {
			try {
				rigid_window::make_samples(n);
			} catch (const std::invalid_argument &) {
				continue;
			}
			throw std::runtime_error("invalid node count accepted");
		}
		std::cout << "Analytic moments, frame degeneracies, printable walls and invalid counts passed\n";
	} catch (const std::exception &error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}
