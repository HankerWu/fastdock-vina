// Numerical regression for strict Armijo acceptance and the accepted BFGS state.
#include "bfgs.h"

#include <iostream>
#include <stdexcept>

struct gradient {
	fl value;
	sz num_floats() const { return 1; }
	fl operator()(sz) const { return value; }
	fl& operator()(sz) { return value; }
};

struct point {
	fl x;
	void increment(const gradient& direction, fl alpha) { x += alpha * direction.value; }
};

void require(bool condition, const char* message) {
	if (!condition) throw std::runtime_error(message);
}

int main() {
	const point origin{0};
	const gradient positive{1};
	const gradient down{-1};
	const gradient up{1};
	point trial{0};
	gradient trial_g{0};
	fl trial_f = 0;
	int evaluations = 0;
	auto uphill = [](const point&, gradient& g) {
		g.value = 1;
		return fl(1);
	};
	const fl legacy_alpha = line_search(uphill, 1, origin, positive, 0, down, trial, trial_g, trial_f, evaluations);
	require(evaluations == 10 && trial.x == -1.0 / 512 && legacy_alpha == 1.0 / 1024,
	        "Legacy line search did not reproduce the untested half-step");
	evaluations = 0;
	require(line_search<true>(uphill, 1, origin, positive, 0, down, trial, trial_g, trial_f, evaluations) == 0
	        && evaluations == 10, "Ten Armijo failures must return no accepted step");
	evaluations = 0;
	require(line_search<true>(uphill, 1, origin, positive, 0, up, trial, trial_g, trial_f, evaluations) == 0
	        && evaluations == 0, "An ascent direction must not be evaluated");

	auto step_then_uphill = [](const point& p, gradient& g) {
		g.value = 1;
		return p.x == 0 ? fl(0) : p.x == -1 ? fl(-1) : fl(-0.5);
	};
	point current = origin;
	gradient current_g{0};
	evaluations = 0;
	const fl legacy_f = bfgs(step_then_uphill, current, current_g, 2, 0, 10, evaluations);
	require(legacy_f == -0.5 && current.x < -1, "Legacy BFGS did not reproduce an uphill commit");
	current = origin;
	evaluations = 0;
	require(bfgs<true>(step_then_uphill, current, current_g, 2, 0, 10, evaluations) == -1
	        && current.x == -1 && current_g.value == 1 && evaluations == 12,
	        "After Armijo failure BFGS must retain the last accepted point and gradient");

	auto steep_quadratic = [](const point& p, gradient& g) {
		g.value = 3 * p.x;
		return 1.5 * p.x * p.x;
	};
	evaluations = 0;
	const fl accepted_alpha = line_search<true>(steep_quadratic, 1, point{1}, gradient{3}, 1.5,
	                                            gradient{-3}, trial, trial_g, trial_f, evaluations);
	require(accepted_alpha == 0.5 && trial.x == -0.5 && trial_f == 0.375
	        && trial_g.value == -1.5 && evaluations == 2,
	        "The accepted half-step, coordinates, value and gradient must agree");

	auto mixed_nonfinite = [](const point& p, gradient& g) {
		g.value = p.x == -0.5 ? std::numeric_limits<fl>::quiet_NaN() : 1;
		return p.x == -1 ? std::numeric_limits<fl>::infinity() : p.x;
	};
	evaluations = 0;
	require(line_search<true>(mixed_nonfinite, 1, origin, positive, 0, down, trial, trial_g, trial_f, evaluations) == 0.25
	        && trial.x == -0.25 && trial_f == -0.25 && trial_g.value == 1 && evaluations == 3,
	        "Armijo must reject nonfinite values and nonfinite gradients");

	auto stationary_quadratic = [](const point& p, gradient& g) {
		g.value = p.x;
		return 0.5 * p.x * p.x;
	};
	current = origin;
	evaluations = 0;
	require(bfgs<true>(stationary_quadratic, current, current_g, 5, 0, 10, evaluations) == 0
	        && evaluations == 1 && current.x == 0 && current_g.value == 0,
	        "A stationary starting point must stop before trials");
	current = point{1};
	evaluations = 0;
	require(bfgs<true>(stationary_quadratic, current, current_g, 5, 0, 10, evaluations) == 0
	        && evaluations == 2 && current.x == 0 && current_g.value == 0,
	        "A stationary accepted point must stop using the new gradient");

	auto negative_curvature = [](const point& p, gradient& g) {
		g.value = 1 - p.x;
		return p.x - 0.5 * p.x * p.x;
	};
	current = origin;
	evaluations = 0;
	require(bfgs<true>(negative_curvature, current, current_g, 2, 0, 10, evaluations) == -7.5
	        && current.x == -3 && current_g.value == 4 && evaluations == 3,
	        "Negative curvature must not give the initial inverse Hessian a negative scale");

	std::cout << "Armijo numerical regressions passed.\n";
}
