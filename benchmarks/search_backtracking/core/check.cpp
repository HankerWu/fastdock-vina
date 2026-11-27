// Focused scalar tests for the actual MC line search and BFGS exhaustion path.
#include "bfgs.h"
#include <iostream>
#include <stdexcept>

static void require(bool passed, const char* message) {
	if (!passed) throw std::runtime_error(message);
}

int main() {
	const fl infinity = std::numeric_limits<fl>::infinity();
	require(quadratic_contraction(0.5, -2, 100) == 0.1, "lower guard");
	require(quadratic_contraction(0.5, -2, -0.5) == 0.5, "upper guard");
	require(quadratic_contraction(0.5, -2, infinity) == 0.5, "infinite trial");
	require(quadratic_contraction(0.5, 1, 2) == 0.5, "non-descent fallback");

	const conf_size layout{{0}, {}};
	conf start(layout), trial(layout);
	start.set_to_null();
	start.ligands[0].rigid.position[0] = 1;
	const auto x = [](const conf& c) { return c.ligands[0].rigid.position[0]; };
	const auto energy = [](fl position) { return 4 * position * position; };
	auto quadratic = [&](const conf& c, flv& gradient) {
		gradient[0] = 8 * x(c);
		return energy(x(c));
	};
	const fl h = 1e-6, p = -8, pg = -64;
	const fl finite_difference = (energy(1 + h * p) - energy(1 - h * p)) / (2 * h);
	require(std::abs(pg - finite_difference) < 1e-7, "directional slope");
	flv direction(layout.num_degrees_of_freedom(), 0), gradient(direction);
	direction[0] = p;
	for (auto mode : {backtracking_mode::native, backtracking_mode::halving, backtracking_mode::quadratic}) {
		search_backtracking monitor{mode};
		int evaluations = 0;
		fl value = 0;
		const fl step = line_search(quadratic, start, energy(1), direction, pg, trial,
		                            gradient, value, evaluations, 0, &monitor);
		const bool interpolated = mode == backtracking_mode::quadratic;
		require(step == 0.125 && x(trial) == 0 && value == 0 && gradient[0] == 0,
		        "accepted coordinates, value and derivative must agree");
		require(evaluations == (interpolated ? 2 : 4) && monitor.searches == 1 &&
		        monitor.rejected == unsigned(evaluations - 1) && monitor.exhausted == 0, "trial accounting");

		// BFGS accepts the first unit step and encounters a nondecreasing
		// local model on every later trial; corrected arms restore that point.
		auto reject_after_first = [&](const conf& c, flv& g) {
			g[0] = 1;
			return x(c) == 1 ? fl(0) : x(c) == 0 ? fl(-1) : fl(-0.5);
		};
		conf current = start;
		bfgs_workspace workspace(layout);
		monitor = search_backtracking{mode};
		evaluations = 0;
		const auto result = workspace.minimize(reject_after_first, current, 2, evaluations,
		                                       nullptr, 0, nullptr, &monitor);
		require(monitor.searches == 2 && monitor.rejected == 10 && monitor.exhausted == 1 &&
		        evaluations == 12, "ten exhausted trials");
		if (mode == backtracking_mode::native)
			require(result.energy == -0.5 && x(current) == -1.0 / 512, "native control");
		else
			require(result.energy == -1 && x(current) == 0 && result.restore_coordinates,
			        "restore last accepted point and model coordinates");
	}
	std::cout << "Backtracking numeric checks passed\n";
}
