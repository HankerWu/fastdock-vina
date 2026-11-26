// Stream real benchmark coordinates through the production overlap implementation.
#include <iomanip>
#include <iostream>
#include "pose_template.h"

int main() {
    std::size_t donor_count, query_count;
    std::cout << std::setprecision(16);
    while (std::cin >> donor_count >> query_count) {
        std::string line;
        std::getline(std::cin, line);
        std::stringstream donor;
        for (std::size_t i = 0; i < donor_count; ++i) {
            if (!std::getline(std::cin, line)) return 2;
            donor << line << '\n';
        }
        pose_template::ligand_atoms elements;
        vecv coords;
        for (std::size_t i = 0; i < query_count; ++i) {
            std::string element;
            fl x, y, z;
            if (!(std::cin >> element >> x >> y >> z)) return 2;
            elements.emplace_back(i, element);
            coords.emplace_back(x, y, z);
        }
        pose_template model(donor, elements);
        std::cout << model.overlap(coords) << '\n';
    }
    return 0;
}
