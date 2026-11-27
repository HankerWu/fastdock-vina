"""After prediction, examine all original probes for the two motivating failures."""

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np

from run import paired, parse_pose, rms, save
from scan import read_ligand, receptor_maps, translate_search, write_start


def main():
    p = argparse.ArgumentParser()
    for f in ("dev", "work", "probe", "vina"):
        p.add_argument("--" + f, type=Path, required=True)
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    meta = json.loads((a.dev / "metadata.json").read_text())
    rows = []
    oriented = []
    for query, rec in [("5EPS", "5C7N"), ("3WKE", "4OD0")]:
        center = np.array(meta[rec]["center"])
        inp = a.dev / query / "ligand.pdbqt"
        fields, origin, spacing = receptor_maps(
            a.vina, a.dev / rec / "receptor.pdbqt", center, a.work / "maps" / rec
        )
        ligand = read_ligand(a.probe, inp)
        rid = {serial: k for k, (_, serial, _) in enumerate(ligand[1])}
        indices = [rid[i] for i in meta[query]["serial"]]
        rot, shift, *_ = paired.alignment(meta[query], meta[rec])
        refs = np.array(meta[query]["mapped"]) @ rot + shift
        for seed in (43103, 43117):
            for arm in ("lattice", "random"):
                probe_xyz = []
                scores = []
                rmsds = []
                rots = []

                # Orthogonal decomposition into heavy-atom centroid translation and centered shape.
                def collect(i, rotated, pos, energy):
                    ref_centers = refs.mean(axis=1)
                    observed = rotated[indices]
                    centroid = observed.mean(0)
                    shape = (
                        np.square(
                            refs - ref_centers[:, None, :] - (observed - centroid)
                        )
                        .sum(axis=-1)
                        .mean(axis=1)
                    )
                    squared = shape[:, None] + np.square(
                        ref_centers[:, None, :] - (pos + centroid)
                    ).sum(axis=-1)
                    distances = np.sqrt(squared.min(axis=0))
                    assert np.isclose(
                        distances[0], rms(observed + pos[0], refs), atol=1e-9
                    )
                    scores.append(energy)
                    rmsds.append(distances)
                    probe_xyz.append(pos)
                    rots.append(rotated)
                    low = int(energy.argmin())
                    oriented.append(
                        {
                            "query": query,
                            "receiver": rec,
                            "seed": seed,
                            "arm": arm,
                            "orientation": i,
                            "probes": len(pos),
                            "coarse_near": int((distances <= 2).sum()),
                            "minimum_rmsd": distances.min(),
                            "lowest_energy": energy[low],
                            "rmsd_at_lowest": distances[low],
                            "best_energy_if_near": energy[distances <= 2].min()
                            if np.any(distances <= 2)
                            else "",
                        }
                    )

                _winners, selected, _n = translate_search(
                    fields, origin, spacing, center, ligand, seed, arm, collect
                )
                en = np.concatenate(scores)
                ds = np.concatenate(rmsds)
                mask = ds <= 2
                mins = np.cumsum([0] + [len(x) for x in scores])
                nearbest = (
                    int(np.argmin(np.where(mask, en, np.inf))) if np.any(mask) else None
                )
                nearidx = int(ds.argmin())
                near_or_winner = sum(
                    bool(np.greater_equal(2, distances[values.argmin()]))
                    for distances, values in zip(rmsds, scores)
                )
                near_after_filter = sum(
                    rms(xyz[indices], refs) <= 2 for _, _, xyz in selected
                )
                dest = a.work / f"{query}_{rec}_{seed}_{arm}"
                dest.mkdir(exist_ok=True)
                paths = []
                for title, j in [
                    ("closest", nearidx),
                    ("best_energy_within2", nearbest),
                ]:
                    if j is None:
                        continue
                    i = int(np.searchsorted(mins[1:], j, side="right"))
                    local_idx = j - mins[i]
                    path = dest / (title + ".pdbqt")
                    before = rots[i] + probe_xyz[i][local_idx]
                    write_start(path, *ligand[:2], before)
                    paths.append((title, path, j, i))
                files = [str(item[1]) for item in paths]
                cmd = [
                    str(a.probe),
                    "refine",
                    str(a.dev / rec / "receptor.pdbqt"),
                    str(inp),
                    *map(str, center),
                    *files,
                ]
                out = subprocess.run(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
                )
                (dest / "oracle.log").write_text(out.stdout)
                if out.returncode:
                    raise RuntimeError(out.stdout)
                refined = [
                    l.split()
                    for l in out.stdout.splitlines()
                    if l.startswith("result ")
                ]
                for (title, path, j, i), line in zip(paths, refined):
                    xyz = parse_pose(
                        Path(str(path) + ".refined"), meta[query]["serial"]
                    )
                    rows.append(
                        {
                            "query": query,
                            "receiver": rec,
                            "seed": seed,
                            "arm": arm,
                            "diagnostic": title,
                            "all_probes": len(ds),
                            "probes_within2": int(mask.sum()),
                            "orientations_with_near": sum(
                                np.any(d <= 2) for d in rmsds
                            ),
                            "near_among_orientation_winners": near_or_winner,
                            "near_among_final_starts": near_after_filter,
                            "proj_rmsd": ds[j],
                            "coarse_score": en[j],
                            "global_rank": 1 + int((en < en[j]).sum()),
                            "same_orientation_rank": 1 + int((scores[i] < en[j]).sum()),
                            "explicit_before": float(line[2]),
                            "explicit_after": float(line[3]),
                            "after_rmsd": rms(xyz, refs),
                        }
                    )
    save(a.work / "trace.csv", rows)
    save(a.work / "orientation.csv", oriented)


if __name__ == "__main__":
    main()
