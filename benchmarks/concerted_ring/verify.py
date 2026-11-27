"""Offline integrity and outcome checks; reconstruct summaries from untouched CSVs."""

import hashlib
import tarfile
from pathlib import Path
from summarize import read, summarize

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
summary = (RESULTS / "summary.csv").read_bytes()
summarize(RESULTS)
assert summary == (RESULTS / "summary.csv").read_bytes()
models = 0
with tarfile.open(HERE / "coordinates.tar.xz", "r:xz") as archive:
    for group in ("panel", "transfer", "catd", "fresh", "lower"):
        runs = read(RESULTS / group / "runs.csv")
        metrics = read(RESULTS / group / "models.csv")
        for row in runs:
            key = row["key"]
            block = [m for m in metrics if m["key"] == key]
            contents = archive.extractfile(group + "/raw/" + key + ".pdbqt").read()
            assert hashlib.sha256(contents).hexdigest() == row["sha256"]
            assert len(block) == int(row["modes"]) == contents.count(b"ENDMDL")
            assert (
                max(
                    abs(float(m["rmsd_A"]) - float(m["independent_sdf_rmsd_A"]))
                    for m in block
                )
                < 0.0015
            )
            models += len(block)
frozen = read(RESULTS / "frozen.csv")
assert len(frozen) == 710 and all(int(r["applied"]) for r in frozen)
assert max(float(r["max_endpoint_A"]) for r in frozen) < 1e-6
assert max(float(r["max_jacobian_error_A_per_rad"]) for r in frozen) < 1e-6
assert len(read(RESULTS / "stages2" / "first_difference.csv")) == 9
assert len(read(RESULTS / "held-stages" / "first_difference.csv")) == 6
clean = read(RESULTS / "clean_replay.csv")
assert len(clean) == 32 and clean[-2]["sha256"] == clean[-1]["sha256"]
print(
    f"{models} output poses, {len(frozen)} common ring draws and 15 observed first divergences verified."
)
