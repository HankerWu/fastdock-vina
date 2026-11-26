#ifndef VINA_PRINTABLE_BOX_H
#define VINA_PRINTABLE_BOX_H

#include "common.h"

// Saving PDBQT rounds Cartesian coordinates to 0.001 A. A converged
// rigid ligand just outside a wall can otherwise be scored with a slope of
// 1e6 and be impossible to re-score from its published coordinates.
inline vec printable_box_shift(const vecv &xyz, const vec &lo, const vec &hi) {
	constexpr fl half_unit = 0.0005;
	vec shift(0, 0, 0);
	if (xyz.empty())
		return shift;
	vec minimum = xyz.front(), maximum = xyz.front();
	for (const auto &p : xyz)
		for (sz d = 0; d < 3; ++d) {
			minimum[d] = std::min(minimum[d], p[d]);
			maximum[d] = std::max(maximum[d], p[d]);
		}
	for (sz d = 0; d < 3; ++d) {
		// A real out-of-box solution is a different failure. Do not move it.
		if (minimum[d] < lo[d] - half_unit || maximum[d] > hi[d] + half_unit)
			return zero_vec;
		// Intersect exact in-box coordinates with the rounding cells of
		// representable in-box decimals. Keep a small guard at rounding ties.
		const fl guard = 1e-9;
		const fl first = std::ceil(lo[d] * 1000) / 1000;
		const fl last = std::floor(hi[d] * 1000) / 1000;
		if (first > last)
			return zero_vec;
		const fl lower = std::max(lo[d], first - half_unit + guard) - minimum[d];
		const fl upper = std::min(hi[d], last + half_unit - guard) - maximum[d];
		if (lower > upper)
			return zero_vec;
		if (lower > 0)
			shift[d] = std::min(upper, lower + guard);
		else if (upper < 0)
			shift[d] = std::max(lower, upper - guard);
	}
	return shift;
}

#endif
