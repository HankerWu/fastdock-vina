/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

/* Modified for FastDock Vina (selected D): column-major symmetric BFGS storage
   and two-double CPU SIMD, preserving each output accumulation order.
   See docs/PROVENANCE.md for upstream and experimental source identities. */

#ifndef VINA_BFGS_H
#define VINA_BFGS_H

#include "conf.h"
#include <cstring>

// Only the local BFGS inverse Hessian uses this full symmetric storage.
class flmat {
	flv m_data;
	sz m_dim;
public:
	explicit flmat(sz n) : m_data(n * n), m_dim(n) {}
	void clear() { std::fill(m_data.begin(), m_data.end(), 0); }
	sz dim() const { return m_dim; }
	sz index(sz i, sz j) const { return i + j * m_dim; }
	const fl& operator()(sz i) const { return m_data[i]; }
	fl& operator()(sz i) { return m_data[i]; }
	const fl& operator()(sz i, sz j) const { return m_data[index(i, j)]; }
	fl& operator()(sz i, sz j) { return m_data[index(i, j)]; }
};

typedef fl bfgs_lane2 __attribute__((vector_size(16)));

inline void minus_mat_vec_product(const flmat& m, const flv& values, flv& out) {
	const sz n = m.dim();
	sz i = 0;
	for(; i + 1 < n; i += 2) {
		bfgs_lane2 sums = {0, 0};
		VINA_FOR(j, n) {
			bfgs_lane2 coefficients;
			std::memcpy(&coefficients, &m(i + j * n), sizeof(coefficients));
			const bfgs_lane2 input = {values[j], values[j]};
			sums += coefficients * input;
		}
		fl lanes[2];
		std::memcpy(lanes, &sums, sizeof(sums));
		out[i] = -lanes[0];
		out[i + 1] = -lanes[1];
	}
	for(; i < n; ++i) {
		fl sum = 0;
		VINA_FOR(j, n)
			sum += m(i + j * n) * values[j];
		out[i] = -sum;
	}
}

inline fl scalar_product(const flv& a, const flv& b) {
	fl sum = 0;
	VINA_FOR_IN(i, a)
		sum += a[i] * b[i];
	return sum;
}

inline void bfgs_update(flmat& h, const flv& direction, const flv& y, fl alpha, fl yp, flv& minus_hy) {
	if(alpha * yp < epsilon_fl) return;
	minus_mat_vec_product(h, y, minus_hy);
	const fl yhy = -scalar_product(y, minus_hy);
	const fl r = 1 / (alpha * yp);
	VINA_FOR(j, h.dim())
		for(sz i = 0; i <= j; ++i) {
			h(i, j) += alpha * r * (minus_hy[i] * direction[j] + minus_hy[j] * direction[i])
			         + alpha * alpha * (r*r * yhy + r) * direction[i] * direction[j];
			h(j, i) = h(i, j);
		}
}

template<bool StrictArmijo = false, typename F>
fl line_search(F& f, const conf& x, fl f0, const flv& p, fl pg, conf& x_new, flv& g_new, fl& f1, int& evalcount,
               int max_evals = 0) {
	if (StrictArmijo && (!std::isfinite(pg) || pg >= 0))
		return 0;
	const fl c0 = 0.0001;
	const unsigned max_trials = 10;
	const fl multiplier = 0.5;
	fl alpha = 1;

	VINA_U_FOR(trial, max_trials) {
		if (StrictArmijo && max_evals > 0 && evalcount >= max_evals)
			return 0;
		x_new.set_incremented(x, p, alpha);
		f1 = f(x_new, g_new);
		++evalcount;
		if (StrictArmijo) {
			if (std::isfinite(f1) && f1 < f0 && f1 - f0 <= c0 * alpha * pg &&
			    std::all_of(g_new.begin(), g_new.end(), [](fl x) { return std::isfinite(x); }))
				return alpha;
		} else if (f1 - f0 < c0 * alpha * pg) {
			break;
		}
		alpha *= multiplier;
	}
	return StrictArmijo ? 0 : alpha;
}

inline void set_diagonal(flmat& m, fl x) {
	VINA_FOR(i, m.dim())
		m(i, i) = x;
}

struct bfgs_result {
	fl energy;
	bool restore_coordinates;
};

// Fixed centroid and length scale for rigid-body steps in a single ligand.
struct bfgs_rigid_metric {
	vec centroid_offset;
	fl inverse_rotation_scale;

	void to_conformation_step(flv& step) const {
		const vec angle = inverse_rotation_scale * vec(step[3], step[4], step[5]);
		const vec shift = vec(step[0], step[1], step[2]) - cross_product(angle, centroid_offset);
		VINA_FOR(i, 3) {
			step[i] = shift[i];
			step[i + 3] = angle[i];
		}
	}

	void to_centroid_gradient(flv& gradient) const {
		const vec force(gradient[0], gradient[1], gradient[2]);
		const vec torque(gradient[3], gradient[4], gradient[5]);
		const vec scaled = inverse_rotation_scale * (torque - cross_product(centroid_offset, force));
		VINA_FOR(i, 3)
			gradient[i + 3] = scaled[i];
	}
};

// Each search or refinement loop owns a workspace for its fixed conformation
// layout. Every minimization starts with a fresh positive diagonal inverse Hessian.
class bfgs_workspace {
	conf trial, original;
	flv gradient, trial_gradient, direction, minus_hy;
	flmat h;

public:
	explicit bfgs_workspace(const conf_size& size)
		: trial(size), original(size), gradient(size.num_degrees_of_freedom()),
		  trial_gradient(gradient.size()), direction(gradient.size()),
		  minus_hy(gradient.size()), h(gradient.size()) {}

	template<bool StrictArmijo = false, typename F>
	bfgs_result minimize(F& f, conf& x, unsigned max_steps, int& evalcount,
	                     const flv* inverse_metric = nullptr, int max_evals = 0,
	                     const bfgs_rigid_metric* rigid_metric = nullptr) {
		h.clear();
		if (inverse_metric) {
			VINA_CHECK(StrictArmijo && inverse_metric->size() == h.dim());
			VINA_FOR(i, h.dim()) {
				VINA_CHECK(std::isfinite((*inverse_metric)[i]) && (*inverse_metric)[i] > 0);
				h(i, i) = (*inverse_metric)[i];
			}
		} else set_diagonal(h, 1);
		original = x;
		fl energy = f(x, gradient);
		++evalcount;
		const fl initial_energy = energy;
		if (StrictArmijo && (!std::isfinite(energy) ||
		    !std::all_of(gradient.begin(), gradient.end(), [](fl x) { return std::isfinite(x); }) ||
		    std::sqrt(scalar_product(gradient, gradient)) < 1e-5))
			return {energy, false};

		VINA_U_FOR(step, max_steps) {
			minus_mat_vec_product(h, gradient, direction);
			const flv* conformation_step = &direction;
			if (StrictArmijo && rigid_metric) {
				minus_hy = direction;
				rigid_metric->to_conformation_step(minus_hy);
				conformation_step = &minus_hy;
			}
			fl trial_energy = 0;
			const fl alpha = line_search<StrictArmijo>(f, x, energy, *conformation_step, scalar_product(direction, gradient),
			                                           trial, trial_gradient, trial_energy, evalcount, max_evals);
			if (StrictArmijo && alpha == 0)
				return {energy, true}; // The model may reflect the last rejected trial.
			energy = trial_energy;
			std::swap(x, trial);
			// The inverse Hessian is only used by subsequent steps.
			if (step + 1 == max_steps) break;
			const auto& convergence_gradient = StrictArmijo ? trial_gradient : gradient;
			if (!(std::sqrt(scalar_product(convergence_gradient, convergence_gradient)) >= 1e-5)) break;
			gradient.swap(trial_gradient);
			// The previous gradient is now scratch until the next trial.
			flv& y = trial_gradient;
			VINA_FOR_IN(i, y)
				y[i] = gradient[i] - y[i];

			const fl yp = scalar_product(y, direction);
			const fl sy = alpha * yp;
			if (StrictArmijo && (!std::isfinite(sy) || sy < epsilon_fl))
				continue;
			if (step == 0) {
				fl yy = scalar_product(y, y);
				if (inverse_metric) {
					yy = 0;
					VINA_FOR_IN(i, y) yy += y[i] * y[i] * (*inverse_metric)[i];
				}
				if (StrictArmijo) {
					if (std::isfinite(yy) && yy > epsilon_fl) {
						const fl scale = sy / yy;
						if (std::isfinite(scale) && scale > 0) {
							if (inverse_metric) {
								// Preserve the relative geometric scales during the first
								// secant rescaling, rather than replacing D by identity.
								bool valid = true;
								VINA_FOR(i, h.dim()) {
									const fl value = scale * (*inverse_metric)[i];
									valid = valid && std::isfinite(value) && value > 0;
								}
								if (valid) VINA_FOR(i, h.dim()) h(i, i) = scale * (*inverse_metric)[i];
							} else set_diagonal(h, scale);
						}
					}
				} else if (std::abs(yy) > epsilon_fl) {
					set_diagonal(h, sy / yy);
				}
			}
			bfgs_update(h, direction, y, alpha, yp, minus_hy);
		}

		if (!(energy <= initial_energy)) { // also reverts when the comparison is unordered
			std::swap(x, original);
			return {initial_energy, true};
		}
		return {energy, false};
	}
};

#endif
