#include "parse_pdbqt.h"
#include "polar_geometry.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {

std::string atom_record(sz serial, const char* type, fl x, fl y = 0, fl z = 0) {
	std::ostringstream line;
	line << "ATOM  " << std::setw(5) << serial << "  X   UNL     1    "
	     << std::fixed << std::setprecision(3) << std::setw(8) << x << std::setw(8) << y << std::setw(8) << z
	     << "  1.00  0.00     0.000 " << type << '\n';
	return line.str();
}

std::string donor(sz serial, fl x, fl y = 0) {
	return atom_record(serial, "N", x, y) + atom_record(serial + 1, "HD", x, y, 1.02);
}

std::string acceptor(sz serial, fl x, fl bonded_carbon) {
	return atom_record(serial, "OA", x) + atom_record(serial + 1, "C", bonded_carbon);
}

polar_geometry::terms measure(const std::string& receptor_path, const std::string& ligand,
                             const std::string& receptor) {
	{
		std::ofstream out(receptor_path);
		out << receptor;
		if (!out)
			throw std::runtime_error("Could not write the temporary receptor fixture");
	}
	model combined = parse_receptor_pdbqt(receptor_path);
	combined.append(parse_ligand_pdbqt_from_string("ROOT\n" + ligand + "ENDROOT\nTORSDOF 0\n", atom_type::XS));
	return polar_geometry(combined).evaluate(combined);
}

void require(polar_geometry::terms actual, fl radial, fl angular, const char* label) {
	if (std::abs(actual.radial - radial) > 1e-12 || std::abs(actual.angular - angular) > 1e-12) {
		std::ostringstream problem;
		problem << label << " expected (" << radial << ", " << angular << ") but got ("
		        << actual.radial << ", " << actual.angular << ')';
		throw std::runtime_error(problem.str());
	}
}

fl radial_weight(fl distance) {
	return std::exp(-sqr((distance - 2.9) / 0.45) / 2);
}

} // namespace

int main(int argc, char** argv) {
	if (argc != 2)
		throw std::runtime_error("Pass an unused receptor scratch PDBQT path");
	const std::string path = argv[1];
	const auto ligand_acceptor = acceptor(1, 0, 1.30);
	require(measure(path, ligand_acceptor, donor(1, -2.90)), 1, 1, "ligand acceptor");
	require(measure(path, ligand_acceptor, donor(1, 2.90)), 1, 0, "angular lower clamp");
	require(measure(path, ligand_acceptor, donor(1, -4.0)), radial_weight(4), radial_weight(4), "4-Angstrom edge");
	require(measure(path, ligand_acceptor, donor(1, -4.001)), 0, 0, "outside edge");
	require(measure(path, ligand_acceptor, donor(1, 0)), 0, 0, "coincident donor and acceptor");

	const fl distance = std::hypot(-0.253, 2.889);
	const fl gate = (std::cos(80 * pi / 180) + 0.253 / distance) /
	                (std::cos(80 * pi / 180) - std::cos(110 * pi / 180));
	require(measure(path, ligand_acceptor, donor(1, -0.253, 2.889)),
	        radial_weight(distance), radial_weight(distance) * gate, "interior cosine gate");
	const auto two_neighbors = ligand_acceptor + atom_record(3, "C", -0.65, 1.126);
	require(measure(path, two_neighbors, donor(1, -2.90)), 1, 0, "smallest of two bonded angles");

	require(measure(path, donor(1, 0), acceptor(1, -2.90, -4.20)), 1, 1, "receptor acceptor");
	const auto amphoteric = ligand_acceptor + atom_record(3, "HD", 0, 0, 0.98);
	const auto partners = donor(1, 2.90) + acceptor(3, -3.35, -4.65);
	require(measure(path, amphoteric, partners), 1, radial_weight(3.35), "separate maxima across both roles");

	for (const char* kind : {"HD", "W", "G0", "Zn"})
		require(measure(path, atom_record(1, "OA", 0) + atom_record(2, kind, 1.30), donor(1, -2.90)),
		        0, 0, "hydrogen, dummy or metal alone is not a heavy covalent anchor");
	require(measure(path, ligand_acceptor, atom_record(1, "Zn", -2.90)), 0, 0, "metal is not a donor");
	std::cout << "Polar geometry: distances, cosine gates, smallest angle, both-role saturation, and atom restrictions passed.\n";
}
