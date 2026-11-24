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

#include "matrix.h"
#include <cstring>

// Only the local BFGS inverse Hessian uses this full symmetric storage.
class flmat {
	flv m_data;
	sz m_dim;
public:
	flmat() : m_dim(0) {}
	flmat(sz n, fl filler) : m_data(n * n, filler), m_dim(n) {}
	sz dim() const { return m_dim; }
	sz index(sz i, sz j) const { return i + j * m_dim; }
	sz index_permissive(sz i, sz j) const { return index(i, j); }
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

template<typename F, typename Conf, typename Change>
fl line_search(F& f, const Conf& x, const fl f0, const Change& p, fl pg, Conf& x_new, Change& g_new, fl& f1, int& evalcount) { // returns alpha
	const fl c0 = 0.0001;
	const unsigned max_trials = 10;
	const fl multiplier = 0.5;
	fl alpha = 1;

	VINA_U_FOR(trial, max_trials) {
		x_new = x; x_new.increment(p, alpha);
		f1 = f(x_new, g_new);
		evalcount++;
		if(f1 - f0 < c0 * alpha * pg) // FIXME check - div by norm(p) ? no?
			break;
		alpha *= multiplier;
	}
	return alpha;
}

inline void set_diagonal(flmat& m, fl x) {
	VINA_FOR(i, m.dim())
		m(i, i) = x;
}

template<typename F, typename Conf, typename Change>
fl bfgs(F& f, Conf& x, Change& g, const unsigned max_steps, int& evalcount) { // x is I/O, final value is returned
	sz n = g.num_floats();
	flmat h(n, 0);
	set_diagonal(h, 1);

	Change g_new(g);
	Conf x_new(x);
	fl f0 = f(x, g);
	evalcount++;

	fl f_orig = f0;
	Change g_orig(g);
	Conf x_orig(x);

	Change p(g);

	// Keep the optimizer's linear algebra flat. Only the model and the
	// conformational increment need the ligand/residue layout.
	flv gradient(n), direction(n), y(n), minus_hy(n);
	VINA_FOR(i, n) gradient[i] = g(i);

	VINA_U_FOR(step, max_steps) {
		minus_mat_vec_product(h, gradient, direction);
		VINA_FOR(i, n) p(i) = direction[i];
		fl f1 = 0;
		const fl alpha = line_search(f, x, f0, p, scalar_product(direction, gradient),
		                             x_new, g_new, f1, evalcount);
		f0 = f1;
		x = x_new;
		if(!(std::sqrt(scalar_product(gradient, gradient)) >= 1e-5)) break;
		g = g_new;
		VINA_FOR(i, n) {
			const fl next = g(i);
			y[i] = next - gradient[i];
			gradient[i] = next;
		}

		const fl yp = scalar_product(y, direction);
		if(step == 0) {
			const fl yy = scalar_product(y, y);
			if(std::abs(yy) > epsilon_fl)
				set_diagonal(h, alpha * yp / yy);
		}
		bfgs_update(h, direction, y, alpha, yp, minus_hy);
	}

	if(!(f0 <= f_orig)) { // succeeds for nans too
		f0 = f_orig;
		x = x_orig;
		g = g_orig;
	}
	return f0;
}

#endif
