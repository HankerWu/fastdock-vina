"""Generate rigid starts using receptor-only Vina fields; independent of benchmark labels."""

import subprocess

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial.transform import Rotation


def receptor_maps(vina, receptor, center, output):
    output.mkdir(parents=True, exist_ok=True)
    empty = output / "empty"
    empty.mkdir(exist_ok=True)
    if any(empty.iterdir()):
        raise ValueError("Map generation requires the empty ligand list")
    cmd = [
        str(vina),
        "--receptor",
        str(receptor),
        "--spacing",
        "0.75",
        "--write_maps",
        str(output / "grid"),
        "--cpu",
        "1",
        "--batch",
        str(empty),
        "--dir",
        str(empty),
    ]
    for axis, val in zip("xyz", center):
        cmd.extend(["--center_" + axis, str(val), "--size_" + axis, "24"])
    with (output / "build.log").open("w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    fields = {}
    for path in sorted(output.glob("grid.*.map")):
        with path.open() as f:
            header = [next(f) for _ in range(6)]
            spacing = float(header[3].split()[1])
            n = np.array(header[4].split()[1:], dtype=int) + 1
            origin = (
                np.array(header[5].split()[1:], dtype=float) - spacing * (n - 1) / 2
            )
            array = np.fromstring(f.read(), sep=" ").reshape(n[::-1]).transpose(2, 1, 0)
            fields[path.name.split(".")[1]] = array
    return fields, origin, spacing


def read_ligand(probe, path):
    text = path.read_text().splitlines()
    rows = [
        (i, int(l[6:11]), [float(l[a : a + 8]) for a in (30, 38, 46)])
        for i, l in enumerate(text)
        if l.startswith(("ATOM", "HETATM"))
    ]
    types = subprocess.check_output([str(probe), "types", str(path)], text=True)
    types = {int(l.split()[0]): l.split()[1] for l in types.splitlines()}
    heavy = [i for i, (_, serial, _) in enumerate(rows) if serial in types]
    xyz = np.array([pos for _, _, pos in rows])
    xyz -= xyz[heavy].mean(0)
    return text, rows, xyz, np.array(heavy), [types[rows[i][1]] for i in heavy]


def write_start(path, text, rows, coords):
    out = text.copy()
    for (index, _, _), p in zip(rows, coords):
        line = out[index]
        out[index] = line[:30] + "".join(f"{v:8.3f}" for v in p) + line[54:]
    path.write_text("\n".join(out) + "\n")


def translate_search(fields, origin, spacing, center, ligand, seed, arm, observe=None):
    _, _, base, heavy, types = ligand
    phase = Rotation.random(random_state=seed)
    rotations = (Rotation.create_group("O") * phase).as_matrix()
    rng = np.random.default_rng(seed)
    proposals = []
    probes = 0
    for index, r in enumerate(rotations):
        rotated = base @ r.T
        offsets = rotated[heavy]
        # Keep rounding and map-header serialization away from either physical edge.
        low = center - 12.0 - offsets.min(0) + 0.002
        high = center + 12.0 - offsets.max(0) - 0.002
        if np.any(low >= high):
            continue
        n = np.ceil((high - low) / 1.5).astype(int)
        count = int(n.prod())
        if arm == "lattice":
            axes = [
                a + (np.arange(k) + 0.5) * (b - a) / k for a, b, k in zip(low, high, n)
            ]
            positions = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(
                -1, 3
            )
        elif arm == "random":
            positions = rng.uniform(low, high, (count, 3))
        else:
            raise ValueError(arm)
        energy = np.zeros(count)
        for off, typ in zip(offsets, types):
            indices = ((positions + off - origin) / spacing).T
            energy += map_coordinates(
                fields[typ],
                indices,
                order=1,
                prefilter=False,
                mode="constant",
                cval=np.inf,
            )
        if observe is not None:
            observe(index, rotated, positions, energy)
        j = int(np.argmin(energy))
        probes += count
        proposals.append((float(energy[j]), index, rotated + positions[j]))
    selected = []
    for candidate in sorted(proposals, key=lambda x: (x[0], x[1])):
        if all(
            np.square(candidate[2][heavy] - old[2][heavy]).sum(axis=1).mean() >= 1
            for old in selected
        ):
            selected.append(candidate)
        if len(selected) == 9:
            break
    return proposals, selected, probes
