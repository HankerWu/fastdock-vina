// Verify successive updates against stable greedy selection on the current
// residents plus one arrival. Discarded history is intentionally inaccessible.
#include "coords.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <stdexcept>

struct sample { double x, energy; sz id; };

static output_type make_pose(sample value) {
    output_type pose(conf{}, value.energy);
    pose.lb = pose.ub = pose.intra = pose.inter = pose.conf_independent = pose.unbound = pose.total = 0;
    pose.native_rank = value.id;
    pose.coords = {vec(value.x, 0, 0), vec(value.x, 2, 1)};
    return pose;
}

static void update_reference(std::vector<sample>& pool, sample arrival, fl radius, sz capacity) {
    auto candidates = pool;
    candidates.push_back(arrival);
    std::stable_sort(candidates.begin(), candidates.end(), [](const sample& a, const sample& b) {
        return a.energy < b.energy;
    });
    pool.clear();
    for (const auto& candidate : candidates) {
        if (pool.size() == capacity) break;
        if (std::none_of(pool.begin(), pool.end(), [&](const sample& resident) {
                return std::abs(resident.x - candidate.x) < radius;
            })) pool.push_back(candidate);
    }
}

int main() {
    output_container bridge;
    add_to_output_container(bridge, make_pose({0, -3, 1}), 1, 3);
    add_to_output_container(bridge, make_pose({1.625, -1, 2}), 1, 3);
    // The less favourable neighbour at 1.625 is closer, but the better one
    // at 0 must block the arrival without evicting any existing resident.
    add_to_output_container(bridge, make_pose({0.875, -2, 3}), 1, 3);
    if (bridge.size() != 2 || bridge[0].native_rank != 1 || bridge[1].native_rank != 2)
        throw std::runtime_error("better non-nearest neighbour failed to block insertion");
    add_to_output_container(bridge, make_pose({0.875, -4, 4}), 1, 3);
    if (bridge.size() != 1 || bridge[0].native_rank != 4)
        throw std::runtime_error("best bridge failed to remove both weaker neighbours");

    // Includes identical coordinates, an exact-radius pair, a bridge, a far
    // pose and equal energies. Two translated atoms give |dx| as RMSD.
    const std::array<double, 5> xs{{0, 0.75, 1.625, 1.625, 2.625}};
    std::array<int, 5> order{{0, 1, 2, 3, 4}};
    sz updates = 0;
    do {
        for (unsigned pattern = 0; pattern < 243; ++pattern) {
            unsigned code = pattern;
            std::array<sample, 5> samples;
            for (sz i = 0; i < samples.size(); ++i) {
                samples[i] = {xs[i], double(code % 3), i};
                code /= 3;
            }
            for (sz capacity : {0u, 1u, 2u, 3u, 5u}) {
                for (fl radius : {0., 1.}) {
                    output_container actual;
                    std::vector<sample> expected;
                    for (int id : order) {
                        const auto& arrival = samples[id];
                        add_to_output_container(actual, make_pose(arrival), radius, capacity);
                        update_reference(expected, arrival, radius, capacity);
                        if (actual.size() != expected.size()) throw std::runtime_error("incorrect pool size");
                        for (sz i = 0; i < expected.size(); ++i)
                            if (actual[i].native_rank != expected[i].id)
                                throw std::runtime_error("incorrect representative or tie ordering");
                        ++updates;
                    }
                }
            }
        }
    } while (std::next_permutation(order.begin(), order.end()));
    std::cout << updates << " online pose archive updates matched the batch reference.\n";
}
