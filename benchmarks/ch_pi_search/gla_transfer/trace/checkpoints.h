// Single-worker diagnostic. The production binary does not include this file.
#pragma once
#include "observer.h"

inline sz& traced_trajectory() {
    static thread_local sz current = 0;
    return current;
}

inline void checkpoint(const char* stage, sz trajectory, sz index, fl energy, const vecv& xyz) {
    auto& stream = archive_trace();
    stream << stage << ' ' << trajectory << ' ' << index << ' ' << energy;
    for (const vec& atom : xyz)
        stream << ' ' << atom[0] << ' ' << atom[1] << ' ' << atom[2];
    stream << '\n';
}
