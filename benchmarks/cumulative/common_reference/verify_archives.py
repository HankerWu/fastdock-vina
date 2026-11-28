"""Reconstruct cross-seed choices from raw released PDBQT and scorer logs."""

import csv
from pathlib import Path
import json
import tarfile

HERE = Path(__file__).resolve().parent
MAIN = HERE.parents[1] / "geometry_comparison/upstream/four_versions/results"
EXT = HERE.parent / "external/results"


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def verify():
    result = HERE / "results"
    pairs = rows(result / "pairs.csv")
    count = chosen = 0
    with (
        tarfile.open(MAIN / "raw.tar.xz") as docking,
        tarfile.open(EXT / "raw.tar.xz") as external,
        tarfile.open(result / "logs.tar.xz") as scored,
    ):
        for p in pairs:
            cohort, arm, query, receiver = [
                p[k] for k in ("cohort", "arm", "query", "receiver")
            ]
            archive = external if cohort == "external" else docking
            prefix = "-".join((cohort, arm, query, receiver))
            refs = {}
            for seed in (85231, 85237):
                key = f"{cohort}-{query}-{receiver}-{seed}-{arm}.pdbqt"
                if key not in archive.getnames():
                    continue
                data = archive.extractfile(key).read().decode()
                for model in data.split("ENDMDL"):
                    tokens = model.splitlines()
                    score_lines = [
                        line
                        for line in tokens
                        if line.startswith("REMARK VINA RESULT:")
                    ]
                    if not score_lines:
                        continue
                    rank = int(
                        next(
                            line for line in tokens if line.startswith("MODEL")
                        ).split()[1]
                    )
                    original = float(score_lines[0].split()[3])
                    raw = float(
                        next(
                            line
                            for line in tokens
                            if line.startswith("REMARK INTER + INTRA")
                        ).split()[-1]
                    )
                    refs[(seed, rank)] = dict(original=original, archived=raw)
            if p["status"] == "unavailable_search":
                continue
            ids = sorted(refs)
            lines = (
                scored.extractfile(prefix + "/fixed.log").read().decode().splitlines()
            )
            for l in lines:
                if l.startswith("POSE "):
                    _, index, value, inter, intra = l.split()
                    refs[ids[int(index)]]["fixed"] = float(value)
            for s in (("batch", "minus3") if cohort == "external" else ("batch",)):
                log = scored.extractfile(f"{prefix}/{s}.log").read().decode()
                for l in log.splitlines():
                    if l.startswith("Scoring ligand:"):
                        identity = tuple(
                            map(int, Path(l.split(": ", 1)[1]).stem.split("-"))
                        )
                    elif l.startswith("Estimated Free Energy of Binding"):
                        refs[identity][s] = float(l.split(":")[1].split()[0])
            assert int(p["modes"]) == len(ids)
            count += len(ids)
            for strategy in ("original", "archived", "fixed", "batch", "minus3"):
                winner = p.get(strategy + "_winner")
                if all(strategy in refs[i] for i in ids):
                    best = min(ids, key=lambda i: (refs[i][strategy], *i))
                    assert winner == f"{best[0]}-{best[1]}.pdbqt"
                    chosen += 1
                else:
                    assert not winner
    return dict(
        compared_seed_pairs=len(pairs),
        raw_modes=count,
        independently_reconstructed_winners=chosen,
    )


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
