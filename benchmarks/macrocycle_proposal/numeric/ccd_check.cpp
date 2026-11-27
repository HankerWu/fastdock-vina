// A closure-only forward-kinematics oracle. It recomputes every endpoint for
// each angle at 0, pi/2 and pi. It shares neither the projected lever-arm
// calculation nor subtree endpoint membership with the production sweep.
#include "mutate.h"
#include "parse_pdbqt.h"

#include <iomanip>
#include <iostream>

namespace {
fl objective(model& m, const conf& c, const interacting_pairs& pairs) {
    m.set(c);
    fl sum = 0;
    for (const auto& pair : pairs) sum += vec_distance_sqr(m.get_coords(pair.a), m.get_coords(pair.b));
    return sum;
}

void postorder(const branch& b, sz& next, szv& order) {
    const sz index = next++;
    for (const auto& child : b.children) postorder(child, next, order);
    order.push_back(index);
}

sz reference_sweep(model& m, conf& c, const interacting_pairs& pairs) {
    sz used = 0;
    for (sz i = 0; i < c.ligands.size(); ++i) {
        sz next = 0;
        szv order;
        const auto ligand = m.get_ligand(i);
        for (const auto& child : ligand.children) postorder(child, next, order);
        for (const sz j : order) {
            fl& angle = c.ligands[i].torsions[j];
            const fl saved = angle;
            const fl f0 = objective(m,c,pairs);
            angle = saved + pi/2;
            const fl f90 = objective(m,c,pairs);
            angle = saved + pi;
            const fl f180 = objective(m,c,pairs);
            const fl cosine = (f180 - f0)/2, sine = (f0 + f180)/2 - f90;
            angle = saved;
            // Summing whole-system squared distances makes flat-axis roundoff
            // proportional to the total objective, including remote ligands.
            if (std::hypot(cosine,sine) <= 512 * std::numeric_limits<fl>::epsilon() *
                                            std::max<fl>({1, f0, f90, f180})) continue;
            angle += std::atan2(sine, cosine);
            normalize_angle(angle);
            ++used;
        }
    }
    return used;
}

vecv positions(model& m, const conf& c) {
    m.set(c);
    vecv xyz;
    for (sz i = 0; i < m.num_movable_atoms(); ++i) xyz.push_back(m.get_coords(i));
    return xyz;
}

fl delta_angles(const conf& first, const conf& second) {
    fl max_error = 0;
    for (sz i = 0; i < first.ligands.size(); ++i)
        for (sz j = 0; j < first.ligands[i].torsions.size(); ++j)
            max_error = std::max(max_error, std::abs(std::remainder(
                first.ligands[i].torsions[j] - second.ligands[i].torsions[j], 2*pi)));
    return max_error;
}
}

int main(int argc, char** argv) {
    if (argc < 3) return 2;
    model input = parse_ligand_pdbqt_from_file(argv[2], atom_type::XS);
    for (int i = 3; i < argc; ++i)
        input.append(parse_ligand_pdbqt_from_file(argv[i], atom_type::XS));
    const auto pairs = input.get_glue_pairs();
    const conf prepared = input.get_initial_conf();
    for (unsigned seed : {191009u,458921u,740317u,384529u,834721u}) {
        rng stream(seed);
        rng generator(static_cast<rng::result_type>(random_int(0, 1000000, stream)));
        model m = input;
        conf first = prepared;
        first.randomize(vec(-12,-12,-12), vec(12,12,12), generator);
        mutate_conf(first, m, 2, generator, false);
        const fl initial = objective(m, first, pairs);
        conf projected = first, reference = first, transformed = first;
        m.project_first_ring_proposal(projected);
        const sz used = reference_sweep(m, reference, pairs);
        sz n_torsions = 0;
        for (const auto& ligand : first.ligands) n_torsions += ligand.torsions.size();
        const auto actual = positions(m, projected), target = positions(m, reference);
        const fl final = objective(m, projected, pairs), ref_final = objective(m, reference, pairs);
        const qt q = angle_to_quaternion(vec(0.387,-0.681,0.214));
        const mat rotation = quaternion_to_r3(q);
        const vec translation(18,-23,14);
        for (auto& ligand : transformed.ligands) {
            ligand.rigid.position = translation + rotation * ligand.rigid.position;
            ligand.rigid.orientation = q * ligand.rigid.orientation;
        }
        m.project_first_ring_proposal(transformed);
        const auto shifted = positions(m, transformed);
        fl atom_error = 0, equivariance_error = 0;
        for (sz i = 0; i < actual.size(); ++i) {
            atom_error = std::max(atom_error, std::sqrt(vec_distance_sqr(actual[i],target[i])));
            equivariance_error = std::max(equivariance_error,std::sqrt(
                vec_distance_sqr(translation + rotation * actual[i],shifted[i])));
        }
        const fl angle_error = delta_angles(projected, reference);
        std::cout << std::setprecision(12) << argv[1] << ',' << seed << ',' << pairs.size() << ','
                  << used << ',' << n_torsions - used << ',' << initial << ',' << final << ',' << ref_final << ',' << angle_error
                  << ',' << atom_error << ',' << equivariance_error << '\n';
        if (angle_error > 1e-8 || atom_error > 1e-8 || equivariance_error > 1e-8 ||
            final > initial + 1e-8 || std::abs(final - ref_final) > 1e-8) return 1;
    }
}
