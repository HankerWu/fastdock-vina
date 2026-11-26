// Local translation/rotation integral on a frozen explicit-receptor pose.
#ifndef VINA_RIGID_WINDOW_H
#define VINA_RIGID_WINDOW_H

#include <boost/random/sobol.hpp>
#include "precalculate.h"
#include "curl.h"

namespace rigid_window {
constexpr fl radius = 1.0;
const fl angle = 15.0 * pi / 180.0;
constexpr fl thermal = 0.00198720425864 * 298.15;

inline vec unit_sphere(fl u, fl v) {
	const fl z = 2 * u - 1, xy = std::sqrt(std::max(0., 1 - z * z)), az = 2 * pi * v;
	return vec(xy * std::cos(az), xy * std::sin(az), z);
}

// Haar measure in the SO(3) geodesic ball is proportional to (1 - cos(theta)) dtheta.
inline fl angular_volume(fl x) {
	if (x < 0.01)
		return x * x * x * (1. / 6 - x * x / 120 + x * x * x * x / 5040);
	return x - std::sin(x);
}

struct sample {
	vec shift, axis;
	fl sine, cosine, weight;
};

inline std::vector<sample> make_samples(unsigned count) {
	if (count < 16 || count > 65536 || (count & (count - 1)))
		throw std::invalid_argument("rigid_window_samples requires a power of two between 16 and 65536");
	boost::random::sobol engine(6);
	const fl scales[] = {1., 0.5, 0.25};
	std::vector<sample> nodes;
	nodes.reserve(3 * count);
	for (unsigned i = 0; i < count; ++i) {
		fl u[6];
		for (fl &value : u)
			value = std::ldexp(static_cast<fl>(engine()),
							   -std::numeric_limits<boost::random::sobol::result_type>::digits);
		for (fl scale : scales) {
			const fl r = radius * scale * std::cbrt(u[0]);
			fl lower = 0, upper = angle * scale;
			const fl target = u[3] * angular_volume(upper);
			for (unsigned k = 0; k < 38; ++k) {
				const fl mid = 0.5 * (lower + upper);
				if (angular_volume(mid) < target)
					lower = mid;
				else
					upper = mid;
			}
			const fl theta = 0.5 * (lower + upper);
			fl density = 0;
			for (fl s : scales) {
				if (r <= radius * s && theta <= angle * s)
					density += angular_volume(angle) / (3 * s * s * s * angular_volume(angle * s));
			}
			nodes.push_back({r * unit_sphere(u[1], u[2]), unit_sphere(u[4], u[5]), std::sin(theta),
							 std::cos(theta), 1 / density});
		}
	}
	return nodes;
}

struct local_frame {
	vec x, y, z;
	vec apply(const vec &p) const { return p[0] * x + p[1] * y + p[2] * z; }
};

inline local_frame frame(const vecv &xyz, const vec &center) {
	vec x(0, 0, 0), y(0, 0, 0);
	for (const vec &p : xyz)
		if (sqr(p - center) > sqr(x))
			x = p - center;
	if (x.norm() < 1e-8)
		x = vec(1, 0, 0);
	else
		x = (1 / x.norm()) * x;
	for (const vec &p : xyz) {
		const vec delta = p - center, transverse = delta - (x * delta) * x;
		if (sqr(transverse) > sqr(y))
			y = transverse;
	}
	if (y.norm() < 1e-8) {
		const vec axis = std::abs(x[0]) < 0.8 ? vec(1, 0, 0) : vec(0, 1, 0);
		y = axis - (axis * x) * x;
	}
	y = (1 / y.norm()) * y;
	return {x, y, cross_product(x, y)};
}

struct integral {
	fl maximum = -std::numeric_limits<fl>::infinity(), total = 0, squared = 0, mass = 0;
	void add(fl delta, fl weight) {
		const fl log_term = std::log(weight) - delta / thermal;
		if (log_term > maximum) {
			const fl scale = std::exp(maximum - log_term);
			total *= scale;
			squared *= scale * scale;
			maximum = log_term;
		}
		const fl term = std::exp(log_term - maximum);
		total += term;
		squared += term * term;
		mass += weight;
	}
	fl correction() const { return -thermal * (maximum + std::log(total / mass)); }
	fl effective_samples() const { return total * total / squared; }
};

struct neighbor {
	vec coords;
	const precalculate_element *table;
};
struct probe_atom {
	vec displacement;
	std::vector<neighbor> neighbors;
};
struct result {
	fl inter, correction, half, uniform, effective;
};

class evaluator {
public:
	evaluator(const model &m, const ScoringFunction &sf, unsigned count)
		: nodes(make_samples(count)), tables(sf.get_cutoff()), cutoff(sf.get_cutoff()) {
		std::bitset<XS_TYPE_SIZE> rec, lig;
		for (const atom &a : m.grid_atoms())
			if (a.xs < XS_TYPE_SIZE)
				rec.set(a.xs);
		for (sz t : m.get_movable_atom_types(atom_type::XS))
			if (!is_glue_type(t))
				lig.set(t);
		tables.populate(sf, rec, lig);
	}

	result evaluate(const model &m, fl slope) const {
		const vecv xyz = m.get_chemical_heavy_atom_movable_coords();
		if (xyz.empty())
			throw std::invalid_argument("rigid_window_pool requires a ligand heavy atom");
		vec center(0, 0, 0);
		for (const vec &p : xyz)
			center += p;
		center *= 1. / xyz.size();
		const auto basis = frame(xyz, center);
		std::vector<probe_atom> atoms;
		VINA_FOR(i, m.num_movable_atoms()) {
			const atom_base &a = m.movable_atom(i);
			if (a.el == EL_TYPE_Dummy || is_glue_type(a.xs) || (a.ad >= AD_TYPE_CG0 && a.ad <= AD_TYPE_CG3))
				throw std::invalid_argument("rigid_window_pool does not support closure or dummy atoms");
			if (a.xs >= XS_TYPE_SIZE)
				continue;
			const vec &p = m.movable_coords(i);
			probe_atom probe{p - center, {}};
			const fl reach = cutoff + radius + 2 * std::sin(angle / 2) * probe.displacement.norm();
			for (const atom &b : m.grid_atoms()) {
				if (b.xs < XS_TYPE_SIZE && vec_distance_sqr(p, b.coords) <= reach * reach)
					probe.neighbors.push_back(
						{b.coords, &tables.element(get_type_pair_index(atom_type::XS, a, b))});
			}
			atoms.push_back(std::move(probe));
		}
		const auto energy = [&](const vec &shift, const vec &axis, fl sn, fl cs) {
			fl energy = 0;
			for (const auto &a : atoms) {
				const vec &v = a.displacement;
				const vec pos =
					center + shift + cs * v + sn * cross_product(axis, v) + ((1 - cs) * (axis * v)) * axis;
				fl e = 0;
				for (const neighbor &b : a.neighbors) {
					const fl d2 = vec_distance_sqr(pos, b.coords);
					if (d2 < tables.cutoff_sqr())
						e += b.table->fast_sample(tables.fast_index(d2));
				}
				curl(e, 1000);
				energy += e;
			}
			return energy;
		};
		const fl initial = energy(vec(0, 0, 0), vec(0, 0, 0), 0, 1);
		integral full, half, uniform;
		VINA_FOR_IN(i, nodes) {
			const auto &p = nodes[i];
			const fl delta =
				slope * (energy(basis.apply(p.shift), basis.apply(p.axis), p.sine, p.cosine) - initial);
			full.add(delta, p.weight);
			if (i < nodes.size() / 2)
				half.add(delta, p.weight);
			if (i % 3 == 0)
				uniform.add(delta, 1);
		}
		return {initial, full.correction(), half.correction(), uniform.correction(), full.effective_samples()};
	}

private:
	std::vector<sample> nodes;
	precalculate tables;
	fl cutoff;
};
} // namespace rigid_window
#endif
