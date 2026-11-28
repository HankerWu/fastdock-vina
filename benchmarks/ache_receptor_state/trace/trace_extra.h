#pragma once
#include "observer.h"
#include <map>
inline std::vector<double> coord_key(const output_type& pose) {
    std::vector<double> key;
    for (const auto& point : pose.coords)
        for (sz d = 0; d < 3; ++d) key.push_back(point[d]);
    key.push_back(pose.e);
    return key;
}
