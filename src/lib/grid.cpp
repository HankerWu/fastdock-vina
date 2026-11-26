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

#include "grid.h"

#include <cstddef>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

void grid::init(const grid_dims& gd) {
	m_data.resize(gd[0].n_voxels + 1, gd[1].n_voxels + 1, gd[2].n_voxels + 1); // number of sample points == n_voxels + 1
	m_init = vec(gd[0].begin, gd[1].begin, gd[2].begin);
	m_range = vec(gd[0].span(), gd[1].span(), gd[2].span());
	assert(m_range[0] > 0);
	assert(m_range[1] > 0);
	assert(m_range[2] > 0);
	m_dim_fl_minus_1 = vec(m_data.dim0() - 1, 
	                       m_data.dim1() - 1,
			               m_data.dim2() - 1);
	VINA_FOR(i, 3) {
		m_factor[i] = m_dim_fl_minus_1[i] / m_range[i];
		m_factor_inv[i] = 1 / m_factor[i];
	}
}

grid_dims grid::read_map(const path& filename) {
	ifile input(filename);
	const auto invalid = [&](const char* reason) {
		return std::runtime_error("Invalid affinity map " + filename.string() + ": " + reason);
	};
	std::string line;
	for (int i = 0; i < 3; ++i) {
		if (!std::getline(input, line))
			throw invalid("incomplete header");
	}

	std::string label;
	fl spacing;
	if (!std::getline(input, line) ||
	    !(std::istringstream(line) >> label >> spacing) ||
	    label != "SPACING" || !(spacing > 0) || !std::isfinite(spacing))
		throw invalid("expected positive SPACING");

	if (!std::getline(input, line))
		throw invalid("missing NELEMENTS");
	std::istringstream counts(line);
	if (!(counts >> label) || label != "NELEMENTS")
		throw invalid("missing NELEMENTS");
	grid_dims dimensions;
	for (auto& axis : dimensions) {
		std::ptrdiff_t voxels;
		if (!(counts >> voxels) || voxels <= 0 || voxels % 2 != 0)
			throw invalid("NELEMENTS must contain three positive even counts");
		axis.n_voxels = voxels;
	}

	if (!std::getline(input, line))
		throw invalid("missing CENTER");
	std::istringstream centers(line);
	if (!(centers >> label) || label != "CENTER")
		throw invalid("missing CENTER");
	for (auto& axis : dimensions) {
		fl center;
		if (!(centers >> center))
			throw invalid("CENTER must contain three coordinates");
		const fl halfspan = axis.n_voxels * spacing / 2.0;
		axis.begin = center - halfspan;
		axis.end = center + halfspan;
		if (!std::isfinite(axis.begin) || !std::isfinite(axis.end) || !(axis.begin < axis.end))
			throw invalid("invalid grid bounds");
	}

	grid loaded(dimensions);
	for (fl& sample : loaded.m_data) {
		if (!std::getline(input, line))
			throw invalid("too few samples");
		const char* begin = line.c_str();
		char* end;
		sample = std::strtod(begin, &end);
		if (end == begin || !std::isfinite(sample) ||
		    line.find_first_not_of(" \t\r\f\v", end - begin) != std::string::npos)
			throw invalid("expected one finite number per sample");
	}
	while (std::getline(input, line))
		if (line.find_first_not_of(" \t\r\f\v") != std::string::npos)
			throw invalid("too many samples");
	if (input.bad())
		throw invalid("read failed");

	*this = std::move(loaded);
	return dimensions;
}

// Catmull-Rom cardinal weights with centered node slopes. A linearly
// extrapolated ghost node gives a one-sided secant at either box face.
namespace {
struct cubic_axis {
	sz index[4];
	fl weight[4], derivative[4];

	cubic_axis(sz low, sz size, fl t) {
		const fl t2 = t * t, t3 = t2 * t;
		index[0] = low ? low - 1 : low;
		index[1] = low;
		index[2] = low + 1;
		index[3] = std::min(low + 2, size - 1);
		weight[0] = (-t + 2*t2 - t3) / 2;
		weight[1] = 1 + (-5*t2 + 3*t3) / 2;
		weight[2] = (t + 4*t2 - 3*t3) / 2;
		weight[3] = (-t2 + t3) / 2;
		derivative[0] = (-1 + 4*t - 3*t2) / 2;
		derivative[1] = (-10*t + 9*t2) / 2;
		derivative[2] = (1 + 8*t - 9*t2) / 2;
		derivative[3] = (-2*t + 3*t2) / 2;
		if (low == 0) fold_ghost(0, 1, 2);
		if (low + 2 == size) fold_ghost(3, 2, 1);
	}

	void fold_ghost(int ghost, int edge, int next) {
		weight[edge] += 2 * weight[ghost];
		weight[next] -= weight[ghost];
		derivative[edge] += 2 * derivative[ghost];
		derivative[next] -= derivative[ghost];
		weight[ghost] = derivative[ghost] = 0;
	}
};

fl cubic_sample(const array3d<fl>& data, const boost::array<sz,3>& cell, const vec& t, vec* grad) {
	const cubic_axis x(cell[0], data.dim0(), t[0]), y(cell[1], data.dim1(), t[1]), z(cell[2], data.dim2(), t[2]);
	fl energy = 0;
	if (grad) grad->assign(0);
	for (int k = 0; k < 4; ++k) {
		fl plane = 0, dx_plane = 0, dy_plane = 0;
		for (int j = 0; j < 4; ++j) {
			const fl* row = &data(0, y.index[j], z.index[k]);
			fl value = 0, dx = 0;
			for (int i = 0; i < 4; ++i) {
				const fl sample = row[x.index[i]];
				value += x.weight[i] * sample;
				if (grad) dx += x.derivative[i] * sample;
			}
			plane += y.weight[j] * value;
			if (grad) {
				dx_plane += y.weight[j] * dx;
				dy_plane += y.derivative[j] * value;
			}
		}
		energy += z.weight[k] * plane;
		if (grad) {
			(*grad)[0] += z.weight[k] * dx_plane;
			(*grad)[1] += z.weight[k] * dy_plane;
			(*grad)[2] += z.derivative[k] * plane;
		}
	}
	return energy;
}
} // namespace

template<bool Cubic>
fl grid::evaluate_aux(const vec& location, fl slope, fl v, vec* deriv) const { // sets *deriv if not NULL
	vec s  = elementwise_product(location - m_init, m_factor); 

	vec miss(0, 0, 0);
	boost::array<int, 3> region;
	boost::array<sz, 3> a;

	VINA_FOR(i, 3) {
		if(s[i] < 0) {
			miss[i] = -s[i];
			region[i] = -1;
			a[i] = 0; 
			s[i] = 0;
		}       
		else if(s[i] >= m_dim_fl_minus_1[i]) {
			miss[i] = s[i] - m_dim_fl_minus_1[i];
			region[i] = 1;
			assert(m_data.dim(i) >= 2);
			a[i] = m_data.dim(i) -  2; 
			s[i] = 1;
		}
		else {
			region[i] = 0; // now that region is boost::array, it's not initialized
			// A valid grid index fits the signed difference type of the storage.
			const auto index = static_cast<std::ptrdiff_t>(s[i]);
			a[i] = index;
			s[i] -= index;
		}
		assert(s[i] >= 0);
		assert(s[i] <= 1);
		assert(a[i] >= 0);
		assert(a[i]+1 < m_data.dim(i));
	}
	const fl penalty = slope * (miss * m_factor_inv); // FIXME check that inv_factor is correctly initialized and serialized
	assert(penalty > -epsilon_fl);

	if (Cubic) {
		fl value = cubic_sample(m_data, a, s, deriv);
		if (deriv) {
			curl(value, *deriv, v);
			VINA_FOR(i, 3)
				(*deriv)[i] = region[i] ? slope * region[i] : m_factor[i] * (*deriv)[i];
		} else {
			curl(value, v);
		}
		return value + penalty;
	}

	// Samples are stored x-fastest. Locate the first corner once and
	// use the row and plane strides for its seven neighbors.
	const sz row = m_data.dim0();
	const fl* lower = &m_data(a[0], a[1], a[2]);
	const fl* upper = lower + row * m_data.dim1();
	const fl f000 = lower[0];
	const fl f100 = lower[1];
	const fl f010 = lower[row];
	const fl f110 = lower[row + 1];
	const fl f001 = upper[0];
	const fl f101 = upper[1];
	const fl f011 = upper[row];
	const fl f111 = upper[row + 1];

	const fl x = s[0];
	const fl y = s[1];
	const fl z = s[2];

	const fl mx = 1-x;
	const fl my = 1-y;
	const fl mz = 1-z;

	fl f = 
		f000 *  mx * my * mz  +
		f100 *   x * my * mz  +
		f010 *  mx *  y * mz  + 
		f110 *   x *  y * mz  +
		f001 *  mx * my *  z  +
		f101 *   x * my *  z  +
		f011 *  mx *  y *  z  +
		f111 *   x *  y *  z  ;

	if(deriv) { // valid pointer
		const fl x_g = 
			f000 * (-1)* my * mz  +
			f100 *   1 * my * mz  +
			f010 * (-1)*  y * mz  + 
			f110 *   1 *  y * mz  +
			f001 * (-1)* my *  z  +
			f101 *   1 * my *  z  +
			f011 * (-1)*  y *  z  +
			f111 *   1 *  y *  z  ;


		const fl y_g = 
			f000 *  mx *(-1)* mz  +
			f100 *   x *(-1)* mz  +
			f010 *  mx *  1 * mz  + 
			f110 *   x *  1 * mz  +
			f001 *  mx *(-1)*  z  +
			f101 *   x *(-1)*  z  +
			f011 *  mx *  1 *  z  +
			f111 *   x *  1 *  z  ;


		const fl z_g =  
			f000 *  mx * my *(-1) +
			f100 *   x * my *(-1) +
			f010 *  mx *  y *(-1) + 
			f110 *   x *  y *(-1) +
			f001 *  mx * my *  1  +
			f101 *   x * my *  1  +
			f011 *  mx *  y *  1  +
			f111 *   x *  y *  1  ;

		vec gradient(x_g, y_g, z_g);
		curl(f, gradient, v);
		vec gradient_everywhere;

		VINA_FOR(i, 3) {
			gradient_everywhere[i] = ((region[i] == 0) ? gradient[i] : 0);
			(*deriv)[i] = m_factor[i] * gradient_everywhere[i] + slope * region[i];
		}

		return f + penalty;
	}
	else {
		curl(f, v);
		return f + penalty;
	}
} 

template fl grid::evaluate_aux<true>(const vec&, fl, fl, vec*) const;
template fl grid::evaluate_aux<false>(const vec&, fl, fl, vec*) const;
