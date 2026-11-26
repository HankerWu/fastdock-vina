// Validate finite-angle proposals against a direct Rodrigues sweep. Read axes
// from PDBQT BRANCH serials, independent of the engine's torsion-tree lookup.
#include "parse_pdbqt.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

static void require(bool ok, const char* message) {
	if (!ok) throw std::runtime_error(message);
}

static vec turn(const vec& x, const vec& o, const vec& u, fl t) {
	const vec r = x - o;
	return o + std::cos(t)*r + std::sin(t)*cross_product(u,r) + ((1-std::cos(t))*(u*r))*u;
}

int main(int argc, char** argv) {
	require(argc == 2, "pass a prepared flexible ligand PDBQT");
	const auto filename = argv[1];
	model m = parse_ligand_pdbqt_from_file(filename, atom_type::XS);
	conf start = m.get_initial_conf();
	std::map<unsigned, sz> atom_indices;
	for (const auto& l : m.get_ligand(0).cont)
		if (l.second) atom_indices.emplace(std::stoul(l.first.substr(6,5)), *l.second);
	std::vector<std::pair<sz, sz>> bonds;
	std::ifstream f(filename);
	for (std::string line; std::getline(f,line);) {
		if (line.compare(0,7,"BRANCH ") != 0) continue;
		std::istringstream fields(line.substr(7));
		unsigned a,b;
		fields >> a >> b;
		bonds.emplace_back(atom_indices.at(a), atom_indices.at(b));
	}
	auto& l = start.ligands.at(0);
	require(bonds.size() == l.torsions.size() && !bonds.empty(),"fixture axes not in one-to-one correspondence");
	for (sz i=0; i<bonds.size(); ++i) l.torsions[i] = 0.17*i-0.4;
	quaternion_increment(l.rigid.orientation, vec(0.57,-0.41,0.92));
	l.rigid.position += vec(-3,7,12);
	unsigned tested=0;
	fl worst_reverse=0, worst_shape=0, greatest_drop=0, worst_scan=0;
	for (sz i=0; i<bonds.size(); ++i) {
		for (fl delta : {-pi, -2.6, -0.18, 1e-7, 0.71, 1.9, pi}) {
			m.set(start);
			const auto before=m.get_chemical_heavy_atom_movable_coords();
			const vec o=m.get_coords(bonds[i].first);
			vec u=m.get_coords(bonds[i].second)-o;
			u=(1/u.norm())*u;
			auto naive=start;
			naive.ligands[0].torsions[i]+=delta;
			m.set(naive);
			const auto raw=m.get_chemical_heavy_atom_movable_coords();
			auto proposed=start; // Deliberately leave the worker at the naive trial.
			m.set_balanced_torsion(proposed.ligands[0],0,i,l.torsions[i]+delta);
			m.set(proposed);
			const auto actual=m.get_chemical_heavy_atom_movable_coords();
			fl sse=0, base=0, mesh=max_fl;
			for (sz j=0; j<actual.size(); ++j) {
				sse+=vec_distance_sqr(actual[j],before[j]);
				base+=vec_distance_sqr(raw[j],before[j]);
				for (sz k=0; k<j; ++k)
					worst_shape=std::max(worst_shape,std::abs((actual[j]-actual[k]).norm()-(raw[j]-raw[k]).norm()));
			}
			for (int angle=-180; angle<=180; ++angle) {
				fl loss=0;
				for(sz j=0; j<raw.size(); ++j)
					loss+=vec_distance_sqr(turn(raw[j],o,u,angle*pi/180),before[j]);
				mesh=std::min(mesh,loss);
			}
			worst_scan=std::max(worst_scan,sse-mesh);
			greatest_drop=std::max(greatest_drop,std::sqrt(base/before.size())-std::sqrt(sse/before.size()));
			m.set_balanced_torsion(proposed.ligands[0],0,i,l.torsions[i]);
			m.set(proposed);
			const auto reversed=m.get_chemical_heavy_atom_movable_coords();
			for (sz j=0;j<before.size();++j)
				worst_reverse=std::max(worst_reverse,(before[j]-reversed[j]).norm());
			++tested;
		}
	}
	require(worst_scan < 1e-7 && worst_reverse < 1e-8 && worst_shape < 1e-8,"finite-angle proposal failed independent numerical check");
	std::cout<<std::setprecision(9)<<filename<<" proposals="<<tested<<" excess_scan_SSE="<<worst_scan
	         <<" max_reverse_A="<<worst_reverse<<" max_shape_A="<<worst_shape<<" max_RMSD_reduction_A="<<greatest_drop<<"\n";
}
