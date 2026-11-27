#include "parse_pdbqt.h"
#include "vina.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

std::string convert_XS_to_string(sz t);

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "types") {
        const model m = parse_ligand_pdbqt_from_file(argv[2], atom_type::XS);
        for (const auto& record : m.get_ligand(0).cont) {
            if (!record.second) continue;
            const atom& a = m.get_atom(*record.second);
            if (a.is_hydrogen()) continue;
            if (a.xs >= XS_TYPE_SIZE) throw std::runtime_error("unsupported atom type");
            std::cout << record.first.substr(6, 5) << " " << convert_XS_to_string(a.xs) << "\n";
        }
        return 0;
    }
    if (argc < 8 || (std::string(argv[1]) != "refine" && std::string(argv[1]) != "fixed"))
        throw std::runtime_error("types prepared.pdbqt | refine receptor ligand cx cy cz start.pdbqt ...");
    Vina v("vina", 1, 1, 0);
    v.set_receptor(argv[2]);
    v.set_ligand_from_file(argv[3]);
    v.compute_vina_maps(std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]),
                        24, 24, 24, 0.375, false, true);
    std::cout << std::setprecision(16);
    for (int i=7; i<argc; ++i) {
        if (std::string(argv[1]) == "fixed") {
            setenv("FD_SCAN_SHOW_UNCAPPED","1",1);
            v.set_pose_coordinates_for_scoring(argv[i]);
            const auto score = v.score(0);
            std::cout << "fixed " << argv[i] << " " << score[1] << " " << score[5] << "\n";
            continue;
        }
        v.set_ligand_from_file(argv[i]);
        const auto before = v.score(0);
        const auto result = v.optimize();
        v.write_pose(std::string(argv[i])+".refined");
        std::cout << "result " << argv[i] << " " << before[1]+before[5] << " "
                  << result[1]+result[5] << "\n";
    }
}
