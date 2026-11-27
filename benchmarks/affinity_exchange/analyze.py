"""Replay frozen first poses and test the effect of arbitrary additive references."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


BENCH = Path(__file__).resolve().parents[1]
ARCHIVE = BENCH / "affinity_reference_ensemble" / "results"
RT = 0.00198720425864083 * 298.15
ARMS = ["native", "rigid", "hinge", "separate", "joint", "size"]
PANELS = {
    "galectin": ("Gal1", "Gal3"),
    "galectin_5_16": ("Gal1", "Gal3"),
    "kinase": ("AURKA", "CHEK1", "S6K1"),
}


def labels():
    gal = pd.read_csv(BENCH / "galectin_affinity/table1.csv")
    parts = []
    for target in PANELS["galectin"]:
        parts.append(
            pd.DataFrame(
                {
                    "target": target,
                    "ligand": gal.compound,
                    "kd_nM": gal[f"{target}_Kd_uM"] * 1000,
                    "sem_nM": gal[f"{target}_error_uM"] * 1000,
                }
            )
        )
    kin = pd.read_csv(BENCH / "torsion_kinase_affinity/panel.csv")
    kin = kin[kin.target.isin(["CHEK1", "AURKA"])].rename(
        columns={"kd_nm": "kd_nM", "cid": "ligand"}
    )
    s6k = pd.read_csv(ARCHIVE / "s6k_panel.csv")
    s6k = s6k.loc[s6k.decision.eq("included")].rename(columns={"cid": "ligand"})
    parts += [kin, s6k.assign(target="S6K1")]
    return pd.concat(parts, ignore_index=True)[["target", "ligand", "kd_nM", "sem_nM"]]


def scores(nodes, paths, published):
    cells = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
    cells = cells[cells.nodes.eq(nodes)].copy()
    assert (
        len(cells) == 500
        and not cells.duplicated(["target", "ligand", "state", "seed"]).any()
    )
    cells["separate"] = cells.rigid + cells.hinge
    for arm in ["rigid", "hinge", "separate", "joint"]:
        cells[arm] += cells.native
    cells["size"] = -cells.n_heavy
    data = cells.groupby(["target", "ligand", "state"])[ARMS].median()
    data = data.groupby(["target", "ligand"]).mean()
    archived = pd.read_csv(published).set_index(["target", "ligand"])[ARMS].sort_index()
    assert data.index.equals(archived.index)
    np.testing.assert_allclose(data, archived, rtol=0, atol=6e-8)
    obs = cells.groupby(["target", "ligand"]).kd_nM.agg(["min", "max"])
    data = data.reset_index().merge(
        labels(), on=["target", "ligand"], how="left", validate="one_to_one"
    )
    np.testing.assert_array_equal(data.kd_nM, obs["min"])
    np.testing.assert_array_equal(data.kd_nM, obs["max"])
    return data


def projection(frame, arm):
    # One dependent indicator in a connected panel is harmless in least squares.
    design = pd.concat(
        [pd.get_dummies(frame[col]).astype(float) for col in ["target", "ligand"]],
        axis=1,
    ).to_numpy()
    energy = RT * np.log(frame.kd_nM.to_numpy() * 1e-9)
    error = energy - frame[arm].to_numpy()
    offsets, _, rank, _ = np.linalg.lstsq(design, error, rcond=None)
    residual = error - design @ offsets
    np.testing.assert_allclose(design.T @ residual, 0.0, atol=1e-9)
    return dict(
        entries=len(frame),
        ligands=frame.ligand.nunique(),
        single_target_ligands=int(frame.ligand.value_counts().eq(1).sum()),
        residual_df=len(frame) - rank,
        best_additive_cell_rmse=float(np.sqrt(np.mean(residual**2))),
        max_abs_residual=float(np.max(np.abs(residual))),
    )


def rectangles(left, right, arms=ARMS):
    common = left.index.intersection(right.index).sort_values()
    for i, j in itertools.combinations(common, 2):
        # Signs are fixed by ascending ligand ID; positive = ligand i loses
        # relative to j on the left target compared to the right target.
        a, b, c, d = left.loc[i], left.loc[j], right.loc[i], right.loc[j]
        x, y, z, w = [row.kd_nM for row in (a, b, c, d)]
        numerator, denominator = x * w, y * z
        measured = RT * np.log(numerator / denominator)
        reversal = (x >= 3 * y and w >= 3 * z) or (y >= 3 * x and z >= 3 * w)
        row = dict(
            i=i,
            j=j,
            k_left_i=x,
            k_left_j=y,
            k_right_i=z,
            k_right_j=w,
            experimental=measured,
            exchange_ge3=max(numerator, denominator) >= 3 * min(numerator, denominator),
            reversal_ge3=reversal,
            experimental_lo=np.nan,
            experimental_hi=np.nan,
        )
        if np.isfinite(a.sem_nM):
            lo = (x - a.sem_nM) * (w - d.sem_nM) / ((y + b.sem_nM) * (z + c.sem_nM))
            hi = (x + a.sem_nM) * (w + d.sem_nM) / ((y - b.sem_nM) * (z - c.sem_nM))
            row.update(experimental_lo=RT * np.log(lo), experimental_hi=RT * np.log(hi))
        for arm in arms:
            prediction = (a[arm] - b[arm]) - (c[arm] - d[arm])
            if abs(prediction) < 1e-12:
                prediction = 0.0
            yield dict(**row, arm=arm, predicted=prediction)


def diagnostics(frame):
    err = frame.predicted - frame.experimental
    signs = np.sign(frame.predicted) * np.sign(frame.experimental)
    eligible = frame.exchange_ge3
    return dict(
        n=len(frame),
        mae=np.abs(err).mean(),
        rmse=np.sqrt((err**2).mean()),
        n_exchange_ge3=int(eligible.sum()),
        credit_ge3=((1 + signs[eligible]) / 2).sum(),
        reversals=int(frame.reversal_ge3.sum()),
        impossible=int((frame.reversal_ge3 & (signs <= 0)).sum()),
    )


def run(dest):
    dest.mkdir(parents=True, exist_ok=True)
    rows, floors, summary = [], [], []
    inputs = [
        BENCH / "galectin_affinity/table1.csv",
        BENCH / "torsion_kinase_affinity/panel.csv",
        ARCHIVE / "s6k_panel.csv",
        Path(__file__).with_name("PROTOCOL.md"),
    ]
    for count, label, s6k_label in [
        (1024, "primary", "s6k_cells"),
        (8192, "sensitivity", "s6k_sensitivity"),
    ]:
        files = [
            ARCHIVE / f"{label}_historical_cells.csv.xz",
            ARCHIVE / f"{s6k_label}.csv.xz",
        ]
        published = ARCHIVE / f"{label}_ligands.csv"
        inputs.extend([*files, published])
        data = scores(count, files, published)
        for panel, targets in PANELS.items():
            block = data[data.target.isin(targets)]
            if panel == "galectin_5_16":
                block = block[block.ligand.between(5, 16)]
            for arm in ARMS:
                floors.append(
                    dict(nodes=count, panel=panel, arm=arm, **projection(block, arm))
                )
            for u, v in itertools.combinations(targets, 2):
                pair = pd.DataFrame(
                    rectangles(
                        *(block[block.target.eq(t)].set_index("ligand") for t in (u, v))
                    )
                )
                ids = sorted(set(pair.i) | set(pair.j))
                pair["target_left"], pair["target_right"] = u, v
                rows.append(pair.assign(nodes=count, panel=panel))
                for arm, f in pair.groupby("arm", sort=False):
                    stats = diagnostics(f)
                    leave_one = [
                        diagnostics(f[(f.i != i) & (f.j != i)])["mae"] for i in ids
                    ]
                    summary.append(
                        dict(
                            nodes=count,
                            panel=panel,
                            left=u,
                            right=v,
                            arm=arm,
                            ligands=len(ids),
                            **stats,
                            loo_mae_min=min(leave_one),
                            loo_mae_max=max(leave_one),
                        )
                    )
    for name, frame in [
        ("contrasts", pd.concat(rows)),
        ("summary", pd.DataFrame(summary)),
        ("reference_floor", pd.DataFrame(floors)),
    ]:
        suffix = ".csv.xz" if name == "contrasts" else ".csv"
        frame.to_csv(dest / (name + suffix), index=False, float_format="%.10g")
    # Post hoc: the historical Vina unbound reference depends on the pool.
    # Replacing that with a single prepared-conformer reference is distinguishable
    # from adding a common per-ligand offset to every already-reported score.
    readback_path = BENCH / "galectin_affinity/results/readback.csv"
    inputs.append(readback_path)
    readback = (
        pd.read_csv(readback_path)
        .query("arm == 'native'")
        .rename(columns={"compound": "ligand"})
    )
    assert len(readback) == 208
    assert (
        (readback.groupby("ligand")[["common_reference", "scaling"]].nunique() == 1)
        .all()
        .all()
    )
    np.testing.assert_allclose(
        readback.common_score - readback.replayed,
        readback.scaling * (readback.pool_reference - readback.common_reference),
        rtol=0,
        atol=1e-7,
    )
    readback["target"] = readback.state.map(
        {"1GZW": "Gal1", "4Q26": "Gal1", "1KJL": "Gal3", "3T1L": "Gal3"}
    )
    assert readback.target.notna().all()
    modes = ["saved", "replayed", "common_score"]
    med = readback.groupby(["target", "ligand", "state"])[modes].median()
    block = (
        med.groupby(["target", "ligand"])
        .mean()
        .reset_index()
        .merge(labels(), on=["target", "ligand"], validate="one_to_one")
    )
    reference_rows = pd.DataFrame(
        rectangles(
            *(
                block[block.target.eq(t)].set_index("ligand")
                for t in PANELS["galectin"]
            ),
            arms=modes,
        )
    )
    reference_rows.to_csv(
        dest / "galectin_reference.csv", index=False, float_format="%.10g"
    )
    hashes = {
        str(p.relative_to(BENCH)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(inputs)
    }
    (dest / "inputs.json").write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dest", type=Path)
    run(parser.parse_args().dest)
