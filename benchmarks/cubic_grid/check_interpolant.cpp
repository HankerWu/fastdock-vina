#include "grid.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void close(fl actual, fl expected, fl tol) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tol)
        throw std::runtime_error("Interpolation or gradient check failed");
}

void check_field(sz nx, sz ny, sz nz) {
    grid_dims box;
    const sz n[] = {nx, ny, nz};
    const vec spacing(0.25, 0.5, 0.75);
    for (sz j=0; j<3; ++j) {
        box[j].begin = -1;
        box[j].n_voxels = n[j];
        box[j].end = -1 + n[j] * spacing[j];
    }
    grid g(box);
    const vec linear(2, -3, 0.25);
    for (sz z=0; z<=nz; ++z)
        for (sz y=0; y<=ny; ++y)
            for (sz x=0; x<=nx; ++x)
                g.m_data(x,y,z) = 4 + linear * g.index_to_argument(x,y,z);
    const fl wall = 13;
    for (fl f : {0.0, 0.031, 0.25, 0.51, 0.993}) {
        vec p;
        for (sz j=0; j<3; ++j) p[j] = box[j].begin + f * box[j].span();
        vec gradient;
        close(g.evaluate_cubic(p, wall, max_fl, &gradient), 4+linear*p, 1e-12);
        for (sz j=0; j<3; ++j) close(gradient[j], linear[j], 1e-11);
    }

    // Degree-two tensor polynomial in the interior; exact node and centered slopes.
    if (nx >= 4 && ny >= 4 && nz >= 4) {
        for (sz z=0; z<=nz; ++z)
            for (sz y=0; y<=ny; ++y)
                for (sz x=0; x<=nx; ++x) {
                    const vec p = g.index_to_argument(x,y,z);
                    g.m_data(x,y,z) = 2 + p[0]*p[0]*p[1] - p[2]*p[1]*p[1];
                }
        const vec p = g.index_to_argument(2,2,2) + vec(.03,.07,.05);
        vec grad;
        close(g.evaluate_cubic(p,wall,max_fl,&grad), 2+p[0]*p[0]*p[1]-p[2]*p[1]*p[1],1e-12);
        close(grad[0],2*p[0]*p[1],1e-11);
        close(grad[1],p[0]*p[0]-2*p[2]*p[1],1e-11);
        close(grad[2],-p[1]*p[1],1e-11);
    }

    // Arbitrary nodal data, including points straddling each interior face.
    sz counter=0;
    for (fl& v : g.m_data) v=3+std::sin(++counter*1.52);
    for (const fl x : {-1.12, -0.915, -0.75, box[0].end+0.11}) {
        const vec p(x,-0.72,-0.51);
        vec grad;
        close(g.evaluate_cubic(p,wall,10,&grad),g.evaluate_cubic(p,wall,10,nullptr),1e-12);
        for (sz j=0; j<3; ++j) {
            if (p[j] == box[j].end) continue; // linear box wall has a kink at its face
            vec delta(0,0,0); delta[j]=1e-6;
            const fl fd=(g.evaluate_cubic(p+delta,wall,10,nullptr)-g.evaluate_cubic(p-delta,wall,10,nullptr))/(2e-6);
            close(grad[j],fd,1e-4);
        }
    }
}
}
int main() {
    check_field(1,1,1);
    check_field(2,3,1);
    check_field(6,5,4);
    std::cout << "Tricubic nodal field, anisotropic derivatives, curl and box-wall checks passed\n";
}
