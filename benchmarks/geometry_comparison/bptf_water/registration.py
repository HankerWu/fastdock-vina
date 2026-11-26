"""Confirm or relocate independently recurrent sites using receptor water coordinates."""

from pathlib import Path
import numpy as np


def match_observed(sites, receptor_pdb, dry_pdbqt, chain="A"):
    """Return one-to-one matches within 1.25 Å, closest first with deterministic ties."""
    heavy = np.array(
        [
            [float(line[k : k + 8]) for k in (30, 38, 46)]
            for line in Path(dry_pdbqt).read_text().splitlines()
            if line.startswith(("ATOM  ", "HETATM"))
            and line[77:].strip() not in ("H", "HD")
        ]
    )
    observed = []
    for line in Path(receptor_pdb).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if (
            line.startswith("HETATM")
            and line[21] == chain
            and line[17:20] in ("HOH", "DOD")
            and line[76:78].strip() == "O"
            and float(line[54:60]) >= 0.8
        ):
            xyz = np.array([float(line[k : k + 8]) for k in (30, 38, 46)])
            if np.linalg.norm(heavy - xyz, axis=1).min() >= 2.5:
                observed.append((line[21:27], xyz))
    candidates = sorted(
        (float(np.linalg.norm(xyz - s["xyz"])), i, j)
        for i, s in enumerate(sites)
        for j, (_, xyz) in enumerate(observed)
    )
    used_sites, used_waters, matched = set(), set(), []
    for distance, i, j in candidates:
        if distance > 1.25:
            break
        water_id, xyz = observed[j]
        if i not in used_sites and water_id not in used_waters:
            used_sites.add(i)
            used_waters.add(water_id)
            matched.append(
                dict(
                    sites[i],
                    observed_xyz=xyz.tolist(),
                    observed_id=water_id,
                    registration_distance=distance,
                )
            )
    return sorted(matched, key=lambda s: s["id"])
