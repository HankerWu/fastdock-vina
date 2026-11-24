#include <iostream>
#include <stdexcept>

#include "grid.h"

int main(int argc, char** argv) {
    if (argc != 2)
        return 2;
    const path filename(argv[1]);
    grid affinity;
    const auto write_fixture = [&](int samples, bool invalid_number = false) {
        ofile output(filename);
        output << "GRID_PARAMETER_FILE none\nGRID_DATA_FILE none\nMACROMOLECULE none\n"
               << "SPACING 0.5\nNELEMENTS 2 2 2\nCENTER 0 0 0\n";
        for (int i = 0; i < samples; ++i) {
            if (invalid_number && i == 7)
                output << "nan\n";
            else
                output << i << '\n';
        }
    };
    const auto has_expected_samples = [&] {
        return affinity.evaluate(vec(0, -0.5, 0.5), 1e6, max_fl) == 19
            && affinity.evaluate(vec(-0.5, 0.5, 0), 1e6, max_fl) == 15;
    };
    write_fixture(27);
    affinity.read_map(filename);
    if (!has_expected_samples())
        return 1;

    const auto rejects = [&](int samples, bool invalid_number, const std::string& reason) {
        write_fixture(samples, invalid_number);
        try {
            affinity.read_map(filename);
        } catch (const std::runtime_error& error) {
            return std::string(error.what()).find(reason) != std::string::npos
                && has_expected_samples();
        }
        return false;
    };
    if (!rejects(26, false, "too few samples") ||
        !rejects(28, false, "too many samples") ||
        !rejects(27, true, "expected one finite number"))
        return 1;

    std::cout << "Affinity map sample order and malformed map rejection passed\n";
}
