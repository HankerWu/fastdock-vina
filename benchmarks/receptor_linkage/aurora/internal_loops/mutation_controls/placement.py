"""Crystal-labelled rigid-placement lower bounds on logged proposals and local returns."""

import argparse
import hashlib
import struct
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import check_predictions as audit
import replay


def analyze(aurora, streams, results, dest):
    dest.mkdir(exist_ok=False)
    cases = audit.rows(results / "runs.csv")
    rows, best = [], []
    for row in cases:
        drug, graph, seed = row["drug"], row["graph"], row["seed"]
        stream = streams / f"{drug}_{graph}_{seed}.bin"
        with stream.open("rb") as f:
            assert hashlib.file_digest(f, "sha256").hexdigest() == row["raw_sha256"]
        with stream.open("rb") as f:
            (n,) = struct.unpack("<I", f.read(4))
            serials = list(struct.unpack(f"<{n}I", f.read(4 * n)))
        dtype = np.dtype(
            [
                ("tr", "<u4"),
                ("step", "<u4"),
                ("phase", "<u4"),
                ("event", "<u4"),
                ("index", "<u4"),
                ("energy", "<f8"),
                ("xyz", "<f8", (n, 3)),
            ]
        )
        raw = np.memmap(stream, dtype=dtype, mode="r", offset=4 * (n + 1))
        query, by_query = audit.frame.read_query(aurora, drug)
        ref, _ = audit.refmol(aurora, drug)
        maps = ref.GetSubstructMatches(query, uniquify=False, maxMatches=10000)
        assert maps and len(maps) < 10000
        ix = [by_query.index(s) for s in serials]
        references = ref.GetConformer().GetPositions()[np.array(maps)[:, ix]]
        ref_center = references[0].mean(0)
        for stage, mask in (
            (
                "proposal",
                (raw["phase"] == 1) & (raw["event"] == 2) & (raw["index"] == 1),
            ),
            ("hunt_return", raw["event"] == 4),
        ):
            indices = np.flatnonzero(mask)
            x = raw["xyz"][indices]
            center = x.mean(1)
            dx = x - center[:, None]
            c2 = np.sum((center - ref_center) ** 2, axis=1)
            t2 = np.full(len(x), np.inf)
            s2 = t2.copy()
            for y in references:
                dy = y - y.mean(0)
                t2 = np.minimum(t2, np.mean(np.sum((dx - dy) ** 2, axis=2), axis=1))
                u, _, vt = np.linalg.svd(np.swapaxes(dx, 1, 2) @ dy)
                sign = np.linalg.det(u @ vt)
                u[:, :, -1] *= sign[:, None]
                s2 = np.minimum(
                    s2, np.mean(np.sum((dx @ (u @ vt) - dy) ** 2, axis=2), axis=1)
                )
            d2 = c2 + t2
            assert np.all(t2 - s2 > -1e-8)
            values = dict(
                as_observed=np.sqrt(d2),
                translate_only=np.sqrt(t2),
                rotate_only=np.sqrt(c2 + s2),
                rigid_refit=np.sqrt(s2),
                centroid=np.sqrt(c2),
                orientation_residual=np.sqrt(np.maximum(t2 - s2, 0)),
            )
            expected = min(
                float(r["min_A"])
                for r in audit.rows(results / "stages.csv")
                if r["drug"] == drug
                and r["seed"] == seed
                and r["graph"] == graph
                and r["stage"] == ("hunt_start" if stage == "proposal" else stage)
            )
            assert abs(values["as_observed"].min() - expected) < 1e-8
            key = dict(drug=drug, seed=seed, graph=graph, stage=stage, n=len(x))
            rows.append(
                dict(
                    **key,
                    **{f"min_{k}_A": v.min() for k, v in values.items()},
                    **{
                        f"le2_{k}": int((v <= 2).sum())
                        for k, v in values.items()
                        if k != "orientation_residual"
                    },
                )
            )
            for k in (
                "as_observed",
                "translate_only",
                "rotate_only",
                "rigid_refit",
                "centroid",
            ):
                i = values[k].argmin()
                event = raw[indices[i]]
                best.append(
                    dict(
                        **key,
                        criterion=k,
                        row=indices[i],
                        trajectory=event["tr"],
                        step=event["step"],
                        **{f"{t}_A": v[i] for t, v in values.items()},
                    )
                )
    replay.geom.write_csv(dest, "placement.csv", rows)
    replay.geom.write_csv(dest, "placement_witnesses.csv", best)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("aurora", "streams", "results", "dest"):
        p.add_argument(name, type=Path)
    analyze(*(v.resolve() for v in vars(p.parse_args()).values()))
