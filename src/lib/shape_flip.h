// Research proposal for a single ligand in a rigid receptor.
#ifndef VINA_SHAPE_FLIP_H
#define VINA_SHAPE_FLIP_H

#include "model.h"

namespace shape_flip {

// A half-turn about an eigenvector preserves all second moments of the
// heavy-atom envelope. Eigenvector signs cancel in a half-turn.
inline bool principal_axis(const vecv& xyz, const vec& center, const vec& draw, vec& axis) {
	mat cov(0,0,0, 0,0,0, 0,0,0), axes(1,0,0, 0,1,0, 0,0,1);
	for (const auto& point : xyz) {
		const vec d = point - center;
		VINA_FOR(i, 3) VINA_FOR(j, 3) cov(i,j) += d[i] * d[j];
	}
	const fl trace = cov(0,0) + cov(1,1) + cov(2,2);
	if (!(trace > 0) || !std::isfinite(trace) || sqr(draw) == 0) return false;
	const fl tol = 64 * std::numeric_limits<fl>::epsilon() * trace;

	// Cyclic Jacobi diagonalization; decline unresolved or repeated roots.
	for (unsigned sweep = 0; sweep < 12; ++sweep) {
		bool diagonal = true;
		for (sz p = 0; p < 2; ++p) for (sz q = p + 1; q < 3; ++q) {
			const fl off = cov(p,q);
			if (std::abs(off) <= tol) continue;
			diagonal = false;
			const fl tau = (cov(q,q) - cov(p,p)) / (2 * off);
			const fl t = std::copysign(1. / (std::abs(tau) + std::hypot(1., tau)), tau);
			const fl c = 1. / std::hypot(1., t), s = t * c;
			cov(p,p) -= t * off;
			cov(q,q) += t * off;
			cov(p,q) = cov(q,p) = 0;
			for (sz k = 0; k < 3; ++k) {
				if (k != p && k != q) {
					const fl a = cov(k,p), b = cov(k,q);
					cov(p,k) = cov(k,p) = c * a - s * b;
					cov(q,k) = cov(k,q) = s * a + c * b;
				}
				const fl x = axes(k,p), y = axes(k,q);
				axes(k,p) = c * x - s * y;
				axes(k,q) = s * x + c * y;
			}
		}
		if (diagonal) break;
	}
	fl best = -1;
	for (sz i = 0; i < 3; ++i) {
		for (sz j = 0; j < i; ++j)
			if (std::abs(cov(i,j)) > tol || std::abs(cov(i,i) - cov(j,j)) <= tol) return false;
		const vec v(axes(0,i), axes(1,i), axes(2,i));
		const fl projection = std::abs(draw * v);
		if (projection > best) {
			best = projection;
			axis = v;
		}
	}
	return true;
}

inline bool apply(conf& c, model& m, const vec& draw) {
	m.set(c); // the last trial may have been rejected
	const vec center = m.ligand_heavy_centroid(0);
	vec axis;
	if (!principal_axis(m.get_chemical_heavy_atom_movable_coords(), center, draw, axis)) return false;
	qt q(0, axis[0], axis[1], axis[2]);
	quaternion_normalize(q);
	auto& rigid = c.ligands.front().rigid;
	rigid.position = center + quaternion_to_r3(q) * (rigid.position - center);
	rigid.orientation = q * rigid.orientation;
	quaternion_normalize_approx(rigid.orientation);
	return true;
}
}
#endif
