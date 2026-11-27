"""Verify archived end-to-end predictions, native counters, summaries and code provenance."""

import hashlib
import re
import subprocess
import tarfile
import tempfile
from pathlib import Path

import summarize

HERE = Path(__file__).resolve().parent
COUNTERS = (
    "hunt_calls",
    "hunt_evals",
    "second_calls",
    "second_evals",
    "explicit_calls",
    "explicit_evals",
)


def verify_runs(tarpath, tables):
    wanted = {}
    for prefix, table in tables:
        for row in summarize.read(table):
            key = prefix + "/raw/" + row["key"]
            if key in wanted:
                raise ValueError("Duplicate prediction " + key)
            wanted[key] = row
    counts = {"pose": 0, "log": 0}
    with tarfile.open(tarpath, "r|xz") as archive:
        for member in archive:
            key, ext = member.name.rsplit(".", 1) if "." in member.name else ("", "")
            if key not in wanted:
                continue
            row = wanted[key]
            if ext == "pdbqt":
                data = archive.extractfile(member).read()
                digest = hashlib.sha256(data).hexdigest()
                if digest != row.get("sha256", row.get("pool_sha256")):
                    raise ValueError("Pool hash differs: " + key)
                if data.count(b"ENDMDL") != int(row["modes"]):
                    raise ValueError("Mode count differs: " + key)
                counts["pose"] += 1
            elif ext == ("counts" if row["arm"] == "upstream" else "log"):
                data = archive.extractfile(member).read().decode()
                for name in COUNTERS:
                    if int(re.search(r"\b" + name + r"=(\d+)", data)[1]) != int(
                        row[name]
                    ):
                        raise ValueError("Counter differs: " + key + " " + name)
                if int(row["joint"]) != sum(
                    int(row[n]) for n in COUNTERS if n.endswith("_evals")
                ):
                    raise ValueError("Total differs: " + key)
                counts["log"] += 1
    if counts != {"pose": len(wanted), "log": len(wanted)}:
        raise ValueError(f"{tarpath}: expected {len(wanted)} pairs, found {counts}")
    print(
        f"{tarpath.parent.name}: {len(wanted)} prediction files and counter logs verified"
    )
    return len(wanted)


def main():
    h, e, k, t = (
        HERE / n for n in ("historical", "external", "fkbp_affinity", "holdout")
    )
    checked = verify_runs(
        h / "poses.tar.xz",
        [("first", h / "runs.csv"), ("persistent", h / "persistent/runs.csv")],
    )
    checked += verify_runs(
        e / "outputs.tar.xz",
        [
            (part, e / part / "runs.csv")
            for part in ("paired", "repeat", "lower", "persistent", "persistent_repeat")
        ],
    )
    checked += verify_runs(
        k / "evidence.tar.xz", [("observed", k / "results/runs.csv")]
    )
    checked += verify_runs(
        t / "observations.tar.xz", [("results", t / "results/runs.csv")]
    )
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        summarize.main(out)
        for name in ("target_summary.csv", "valid_first.csv"):
            if (out / name).read_bytes() != (HERE / name).read_bytes():
                raise ValueError("Aggregation differs: " + name)
    null = HERE / "null_control"
    wall = e / "wall_match"
    verify_runs(
        null / "raw.tar.xz",
        [(p, null / "results" / p / "runs.csv") for p in ("5QCT", "5K0X")],
    )
    verify_runs(wall / "raw.tar.xz", [("validation", wall / "validation/runs.csv")])
    for root, file in (
        (e, "SHA256SUMS"),
        (k, "evidence.sha256"),
        (null, "raw.sha256"),
        (wall, "SHA256SUMS"),
    ):
        subprocess.run(["sha256sum", "--quiet", "-c", file], cwd=root, check=True)
    print(
        f"{checked} primary and lower-budget records plus null and wall-time controls verified. Target and frozen-pool tables reproduced."
    )


if __name__ == "__main__":
    main()
