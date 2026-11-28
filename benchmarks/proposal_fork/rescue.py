"""Count successes/losses of one paired tenth trial on stable paths and solve each original proposal anew."""

import argparse
import csv
import json
import lzma
import os
from pathlib import Path
import sys
import tarfile
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from benchmarks.proposal_fork.replay import (
    ARMS,
    HERE,
    TRACE,
    argv,
    csvfile,
    digest,
    fold,
    native,
    records,
    restore,
    run_panel,
    structures,
)


def fixtures(work, historical):
    if historical is None:
        historical = work / "history"
        restore(HERE.parent / "bfgs_source_replay/raw", historical)
    ext = work / "external"
    with tarfile.open(
        HERE.parent / "stage_local_structure/external/inputs.tar.xz"
    ) as archive:
        archive.extractall(ext, filter="data")
    e, skipped = structures(ext)
    assert not skipped
    fpath = historical / "fixtures/reciprocal"
    f = [
        json.loads((fpath / c / "preparation.json").read_text())
        for c in ("2NNQ", "3FR2")
    ]
    for inputs, q, r, seed, drivers in [
        (fpath, *f, 85237, ("e7031e9", "2326430")),
        (ext, *e, 85231, ARMS),
        (ext, *e[::-1], 85231, ARMS),
    ]:
        case = f"{q['pdb']}-{r['pdb']}-{seed}"
        for driver in drivers:
            if inputs == ext:
                path = work / f"{case}-{driver}.native"
                path.write_bytes(
                    lzma.decompress(
                        (TRACE / f"reciprocal-{case}-{driver}.xz").read_bytes()
                    )
                )
            else:
                path = historical / case / f"{driver}.native"
            yield case, driver, inputs, q, r, seed, fold(native(path)[0])


def main(work, out, binary, historical):
    live = binary is not None
    baseline = HERE / "raw"
    if live:
        out.mkdir(parents=True, exist_ok=False)
    previous = json.loads((baseline / "costs.json").read_text())
    forks = list(csv.DictReader((baseline / "forks.csv").open()))
    sums, receipts = [], []
    for case, driver, inputs, q, r, seed, proposals in fixtures(work, historical):
        key = f"{case}-{driver}"
        selected = set(range(0, len(proposals), 97))
        selected.update(
            int(row[driver + "_index"]) for row in forks if row["case"] == case
        )
        orig = next(item for item in previous if item["key"] == key + "-stable")
        states = work / f"{key}.inputs"
        indices = work / f"{key}.indices"
        states.write_text("".join(s["P"][0][3] + "\n" for s in proposals))
        indices.write_text("".join(f"{i}\n" for i in sorted(selected)))
        assert (
            run_panel.digest(states) == orig["inputs"]["states"]
            and run_panel.digest(indices) == orig["inputs"]["indices"]
        )
        stable = records(baseline / f"{key}-stable.shadow.xz")
        archived_trials = lzma.decompress(
            (baseline / f"{key}-stable.trials.xz").read_bytes()
        )
        shadow = work / f"{key}.shadow"
        trials = work / f"{key}.trials"
        rescues = work / f"{key}.rescues"
        if live:
            for mode in ("isolated", "production"):
                os.environ.update(
                    FD_BFGS_SHADOW=str(states),
                    FD_BFGS_TRACE=str(shadow),
                    FD_FORK_TRIALS=str(trials),
                    FD_FORK_INDICES=str(indices),
                )
                if mode == "isolated":
                    os.environ.update(
                        FD_FORK_COUNTERFACTUAL="1", FD_FORK_RESCUES=str(rescues)
                    )
                else:
                    os.environ.pop("FD_FORK_COUNTERFACTUAL", None)
                cmd = argv(binary, inputs, q, r, seed, shadow.with_suffix(".unused"))
                stats = run_panel.execute(cmd, out / f"{key}-{mode}.log", 180)
                assert stats["status"] == "ok", (key, mode, stats)
                hashes = {}
                if mode == "isolated":
                    assert trials.read_bytes() == archived_trials, key
                    assert shadow.read_bytes() == lzma.decompress(
                        (baseline / f"{key}-stable.shadow.xz").read_bytes()
                    ), key
                    files = [rescues]
                    hashes["reproduced_stable_trace"] = digest(shadow.read_bytes())
                    hashes["reproduced_stable_trials"] = digest(trials.read_bytes())
                else:
                    files = [shadow, trials]
                for source in files:
                    name = key + source.suffix + ".xz"
                    dest = out / name
                    dest.write_bytes(lzma.compress(source.read_bytes(), preset=3))
                    hashes[name] = run_panel.digest(dest)
                receipts.append(
                    dict(
                        key=key,
                        mode=mode,
                        argv=cmd,
                        binary=run_panel.digest(binary),
                        inputs=orig["inputs"],
                        outputs=hashes,
                        **stats,
                    )
                )
        else:
            for row in json.loads((out / "costs.json").read_text()):
                if row["key"] == key:
                    assert row["inputs"] == orig["inputs"]
                    for name, sha in row["outputs"].items():
                        if name.startswith("reproduced_"):
                            data = (
                                archived_trials
                                if name.endswith("trials")
                                else lzma.decompress(
                                    (baseline / f"{key}-stable.shadow.xz").read_bytes()
                                )
                            )
                            assert digest(data) == sha
                        else:
                            assert run_panel.digest(out / name) == sha
        events = [
            line.split(",")
            for line in lzma.decompress((out / f"{key}.rescues.xz").read_bytes())
            .decode()
            .splitlines()
        ]
        losses, rescued, old_failed, old_pass, changed = 0, 0, 0, 0, 0
        for event in events:
            i, pg, f0, old_a, new_a, old_e, new_e, old, new = map(float, event)
            assert i == int(i) and old_a == 2**-9 and 0 < new_a <= old_a
            assert old == int(old_e - f0 < 1e-4 * old_a * pg) and new == int(
                new_e - f0 < 1e-4 * new_a * pg
            )
            changed += new_a != old_a
            old_failed += not old
            old_pass += bool(old)
            rescued += bool(new) and not old
            losses += bool(old) and not new
        assert old_failed == sum(int(row[5]) for row in stable)
        result = records(out / f"{key}.shadow.xz")
        assert len(result) == len(stable) and all(
            a[9:11] == b[9:11] for a, b in zip(result, stable)
        )
        ref = run_panel.replay.reciprocal.reference(q, r)
        serials = [
            int(l[6:11])
            for l in (inputs / q["pdb"] / "ligand.pdbqt").read_text().splitlines()
            if l.startswith(("ATOM  ", "HETATM"))
            and l[77:].strip() not in {"H", "HD", "HS"}
        ]
        order = [serials.index(i) for i in q["serial"]]
        errors = []
        for row in result:
            xyz = np.fromstring(row[8], dtype=int, sep=" ").reshape(-1, 3)[order] / 1000
            errors.append(
                float(np.sqrt(np.square(xyz[None] - ref).sum(-1).mean(-1).min()))
            )
        de = np.array([float(b[2]) - float(a[2]) for a, b in zip(stable, result)])
        no_tenth = {i for i in range(len(proposals))} - {int(row[0]) for row in events}
        assert all(result[i] == stable[i] for i in no_tenth)
        sums.append(
            dict(
                key=key,
                proposals=len(result),
                isolated_tenths=len(events),
                alpha_changed=changed,
                old_exhaustions=old_failed,
                rescued=rescued,
                old_passed_tenths=old_pass,
                old_pass_lost=losses,
                no_tenth_both_identical=len(no_tenth),
                stable_calls=sum(int(r[4]) for r in stable),
                new_calls=sum(int(r[4]) for r in result),
                new_exhaustions=sum(int(r[5]) for r in result),
                returned_E_lower=int((de < -1e-7).sum()),
                returned_E_higher=int((de > 1e-7).sum()),
                new_near=sum(e <= 2 for e in errors),
                new_best_A=min(errors),
            )
        )
        print(sums[-1], flush=True)
    csvfile(out / "comparison.csv", sums, live)
    if live:
        (out / "costs.json").write_text(json.dumps(receipts, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--binary", type=Path)
    p.add_argument("--history", type=Path)
    a = p.parse_args()
    with tempfile.TemporaryDirectory() as work:
        main(Path(work), a.out.resolve(), a.binary, a.history)
