// Numerical regression for strict Armijo acceptance and the accepted BFGS state.
#include "bfgs.h"

#include <iostream>
#include <stdexcept>

const conf_size size{{0}, {}};

conf point(fl x) {
	conf result(size);
	result.set_to_null();
	result.ligands.front().rigid.position[0] = x;
	return result;
}

fl coordinate(const conf& position) {
	return position.ligands.front().rigid.position[0];
}

flv direction(fl x) {
	flv result(size.num_degrees_of_freedom(), 0);
	result.front() = x;
	return result;
}

void require(bool condition, const char* message) {
	if (!condition) throw std::runtime_error(message);
}

int main() {
	const conf origin = point(0);
	const flv down = direction(-1);
	const flv up = direction(1);
	conf trial(size);
	flv trial_g(size.num_degrees_of_freedom(), 0);
	fl trial_f = 0;
	int evaluations = 0;
	auto uphill = [](const conf&, flv& g) {
		g.front() = 1;
		return fl(1);
	};
	const fl legacy_alpha = line_search(uphill, origin, 0, down, -1, trial, trial_g, trial_f, evaluations);
	require(evaluations == 10 && coordinate(trial) == -1.0 / 512 && legacy_alpha == 1.0 / 1024,
	        "Legacy line search did not reproduce the untested half-step");
	evaluations = 0;
	require(line_search<true>(uphill, origin, 0, down, -1, trial, trial_g, trial_f, evaluations) == 0
	        && evaluations == 10, "Ten Armijo failures must return no accepted step");
	evaluations = 0;
	require(line_search<true>(uphill, origin, 0, up, 1, trial, trial_g, trial_f, evaluations) == 0
	        && evaluations == 0, "An ascent direction must not be evaluated");

	auto step_then_uphill = [](const conf& p, flv& g) {
		g.front() = 1;
		return coordinate(p) == 0 ? fl(0) : coordinate(p) == -1 ? fl(-1) : fl(-0.5);
	};
	conf current = origin;
	bfgs_workspace workspace(size);
	evaluations = 0;
	const fl legacy_f = workspace.minimize(step_then_uphill, current, 2, evaluations).energy;
	require(legacy_f == -0.5 && coordinate(current) < -1, "Legacy BFGS did not reproduce an uphill commit");
	current = origin;
	evaluations = 0;
	require(workspace.minimize<true>(step_then_uphill, current, 2, evaluations).energy == -1
	        && coordinate(current) == -1 && evaluations == 12,
	        "After Armijo failure BFGS must retain the last accepted point");

	auto steep_quadratic = [](const conf& p, flv& g) {
		g.front() = 3 * coordinate(p);
		return 1.5 * coordinate(p) * coordinate(p);
	};
	evaluations = 0;
	const fl accepted_alpha = line_search<true>(steep_quadratic, point(1), 1.5, direction(-3), -9, trial, trial_g, trial_f, evaluations);
	require(accepted_alpha == 0.5 && coordinate(trial) == -0.5 && trial_f == 0.375
	        && trial_g.front() == -1.5 && evaluations == 2,
	        "The accepted half-step, coordinates, value and gradient must agree");

	auto mixed_nonfinite = [](const conf& p, flv& g) {
		g.front() = coordinate(p) == -0.5 ? std::numeric_limits<fl>::quiet_NaN() : 1;
		return coordinate(p) == -1 ? std::numeric_limits<fl>::infinity() : coordinate(p);
	};
	evaluations = 0;
	require(line_search<true>(mixed_nonfinite, origin, 0, down, -1, trial, trial_g, trial_f, evaluations) == 0.25
	        && coordinate(trial) == -0.25 && trial_f == -0.25 && trial_g.front() == 1 && evaluations == 3,
	        "Armijo must reject nonfinite values and nonfinite gradients");

	auto stationary_quadratic = [](const conf& p, flv& g) {
		g.front() = coordinate(p);
		return 0.5 * coordinate(p) * coordinate(p);
	};
	current = origin;
	evaluations = 0;
	require(workspace.minimize<true>(stationary_quadratic, current, 5, evaluations).energy == 0
	        && evaluations == 1 && coordinate(current) == 0,
	        "A stationary starting point must stop before trials");
	current = point(1);
	evaluations = 0;
	require(workspace.minimize<true>(stationary_quadratic, current, 5, evaluations).energy == 0
	        && evaluations == 2 && coordinate(current) == 0,
	        "A stationary accepted point must stop using the new gradient");

	auto negative_curvature = [](const conf& p, flv& g) {
		g.front() = 1 - coordinate(p);
		return coordinate(p) - 0.5 * coordinate(p) * coordinate(p);
	};
	current = origin;
	evaluations = 0;
	require(workspace.minimize<true>(negative_curvature, current, 2, evaluations).energy == -7.5
	        && coordinate(current) == -3 && evaluations == 3,
	        "Negative curvature must not give the initial inverse Hessian a negative scale");

	// An anisotropic quadratic models coordinate scaling with unequal lever arms.
	auto anisotropic = [](const conf& p, flv& g) {
		std::fill(g.begin(), g.end(), 0);
		const auto& x = p.ligands.front().rigid.position;
		g[0] = x[0];
		g[1] = 16 * x[1];
		return 0.5 * (x[0] * x[0] + 16 * x[1] * x[1]);
	};
	current = point(1);
	current.ligands.front().rigid.position[1] = 1;
	conf unscaled = current;
	flv inverse(size.num_degrees_of_freedom(), 1);
	inverse[1] = 1.0 / 16;
	evaluations = 0;
	require(workspace.minimize<true>(anisotropic, current, 1, evaluations, &inverse).energy == 0
	        && evaluations == 2, "Positive diagonal preconditioner must solve separable scaled quadratic");
	evaluations = 0;
	require(workspace.minimize<true>(anisotropic, unscaled, 1, evaluations).energy > 0
	        && evaluations > 2, "Identity control must not solve scaled quadratic in one step");

	std::cout << "Armijo and geometry-preconditioned numerical regressions passed.\\n";
}
