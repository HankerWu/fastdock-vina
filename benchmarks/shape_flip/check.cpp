#include "shape_flip.h"
#include "parse_pdbqt.h"
#include "mutate.h"
#include <iostream>
#include <iomanip>

static vec mean(const vecv& xyz) {
	vec sum(0,0,0);
	for (const auto& p : xyz) sum += p;
	return (1./xyz.size()) * sum;
}
static mat covariance(const vecv& xyz) {
	mat c(0,0,0, 0,0,0, 0,0,0);
	const vec center = mean(xyz);
	for (const auto& p : xyz) {
		const vec d = p - center;
		VINA_FOR(i,3) VINA_FOR(j,3) c(i,j) += d[i]*d[j]/xyz.size();
	}
	return c;
}

int main(int argc, char** argv) {
	VINA_CHECK(argc == 2);
	vec rejected;
	const vec draw(0.7,-0.1,0.2), zero(0,0,0);
	VINA_CHECK(!shape_flip::principal_axis({zero}, zero, draw, rejected));
	VINA_CHECK(!shape_flip::principal_axis({vec(-1,0,0),vec(1,0,0)}, zero, draw, rejected));
	model m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
	rng states(173), native(129), flipped(129);
	size_t opportunities=0, paired=0;
	fl com=0, cov=0, back=0, separation=0;
	for (unsigned k=0; k<1200; ++k) {
		conf c=m.get_initial_conf(), stale=c;
		c.randomize(vec(-5,-5,-5),vec(5,5,5),states);
		stale.randomize(vec(-5,-5,-5),vec(5,5,5),states);
		m.set(c);
		const auto before = m.get_chemical_heavy_atom_movable_coords();
		conf a=c, b=c;
		m.set(stale); mutate_conf(a,m,2,native,false);
		m.set(stale); const bool turned=mutate_conf(b,m,2,flipped,false,false,false,false,true);
		VINA_CHECK(native==flipped);
		m.set(b);
		const auto after=m.get_chemical_heavy_atom_movable_coords();
		if (!turned) {
			m.set(a);
			const auto parity=m.get_chemical_heavy_atom_movable_coords();
			VINA_FOR_IN(i,after) VINA_CHECK(vec_distance_sqr(parity[i],after[i])==0);
			++paired;
			continue;
		}
		++opportunities;
		com=std::max(com,(mean(before)-mean(after)).norm());
		const auto old=covariance(before), current=covariance(after);
		VINA_FOR(i,9) cov=std::max(cov,std::abs(old.data[i]-current.data[i]));
		VINA_FOR_IN(i,after) VINA_FOR(j,i) separation=std::max(separation,
		    std::abs((before[i]-before[j]).norm()-(after[i]-after[j]).norm()));
		// A common laboratory-frame draw selects the same axis after a half-turn.
		conf repeated=c;
		const vec direction=random_inside_sphere(states);
		VINA_CHECK(shape_flip::apply(repeated,m,direction));
		VINA_CHECK(shape_flip::apply(repeated,m,direction));
		m.set(repeated);
		const auto twice=m.get_chemical_heavy_atom_movable_coords();
		VINA_FOR_IN(i,twice) back=std::max(back,(before[i]-twice[i]).norm());
	}
	VINA_CHECK(opportunities>0 && paired>0 && com<1e-9 && cov<1e-9 && back<1e-8 && separation<1e-9);
	std::cout<<argv[1]<<","<<opportunities<<","<<paired<<","<<std::setprecision(12)
	         <<com<<","<<cov<<","<<back<<","<<separation<<"\n";
}
