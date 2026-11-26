// Independent scalar reference: ordinary non_cache grid on each rounded PDBQT, no swept neighbor list.
#include "non_cache.h"
#include "parse_pdbqt.h"

int main(int argc, char **argv) {
	try {
		if (argc != 4)
			throw std::runtime_error("reference receptor.pdbqt prepared.pdbqt transformed.pdbqt");
		const model rec = parse_receptor_pdbqt(argv[1]);
		model m = rec;
		m.append(parse_ligand_pdbqt_from_file(argv[2], atom_type::XS));
		const auto sf = std::make_shared<ScoringFunction>(
			SF_VINA, flv{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50, 5 * 0.05846 / 0.1 - 1});
		grid_dims gd;
		for (sz j = 0; j < 3; ++j) {
			gd[j].begin = max_fl;
			gd[j].end = -max_fl;
			for (const atom &a : rec.grid_atoms()) {
				gd[j].begin = std::min(gd[j].begin, a.coords[j] - 9);
				gd[j].end = std::max(gd[j].end, a.coords[j] + 9);
			}
			gd[j].n_voxels = static_cast<sz>(std::ceil((gd[j].end - gd[j].begin) / 0.5));
		}
		non_cache grid(rec, gd, sf, 0);
		grid.prepare_for_model(m);
		std::ifstream in(argv[3]);
		std::string block, line;
		std::cout << std::setprecision(14);
		while (std::getline(in, line)) {
			block += line + "\n";
			if (line == "ENDMDL") {
				replace_fixed_pose_coordinates(m, block);
				std::cout << grid.eval(m, 1000) << "\n";
				block.clear();
			}
		}
		if (in.bad() || !block.empty())
			throw std::runtime_error("Bad reference input");
	} catch (const std::exception &e) {
		std::cerr << e.what() << '\n';
		return 1;
	}
}
