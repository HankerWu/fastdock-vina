"""Verify a real round-start or restored production execution against every fixed assay pool."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

FIELDS = ("modes", "seed", "cap", "receptor", "compound", "assay", "top_score", "output_sha256",
          "hunt_evals", "second_evals", "explicit_evals", "late_prescore")
EXPECTED_BINARY = "f6ff9e0d9d6f0bc6e77e65409f9c7f88fccd826eca064a49e9c91a9a73496f9f"


def verify(work, reference):
    def inputs(d):
        return json.loads((d / "inputs.json").read_text())["files"]
    assert inputs(work) == inputs(reference)
    sub = Path("outputs/stable/nominal")
    names = {p.name for p in (reference / sub).glob("*.json")}
    assert len(names) == 140 and names == {p.name for p in (work / sub).glob("*.json")}
    total = {k: 0 for k in ("modes", "hunt_evals", "second_evals", "explicit_evals", "late_prescore")}
    cpu = wall = 0.0
    peak = 0
    for name in sorted(names):
        actual, old = (base / sub / name for base in (work, reference))
        now, before = (json.loads(f.read_text()) for f in (actual, old))
        assert now["valid"] and before["valid"]
        assert now["binary_sha256"] == before["binary_sha256"] == EXPECTED_BINARY
        assert all(now[k] == before[k] for k in FIELDS), name
        body = actual.with_suffix(".pdbqt").read_bytes()
        assert body == old.with_suffix(".pdbqt").read_bytes()
        assert hashlib.sha256(body).hexdigest() == now["output_sha256"]
        for k in total:
            total[k] += now[k]
        cpu += now["user_s"] + now["sys_s"]
        wall += now["wall_s"]
        peak = max(peak, now["maxrss_kib"])
    return dict(binary_sha256=EXPECTED_BINARY, identical_pools=len(names), **total,
                cpu_s=cpu, sum_job_wall_s=wall, maxrss_kib=peak)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("reference", type=Path, help="study.py setup workspace")
    p.add_argument("out", type=Path)
    p.add_argument("--archive", action="store_true")
    a = p.parse_args()
    evidence = verify(a.work, a.reference)
    a.out.mkdir(parents=True, exist_ok=True)
    body = json.dumps(evidence, indent=2) + "\n"
    if a.archive:
        (a.out/"verification.json").write_text(body)
        with tarfile.open(a.out/"raw.tar.xz", "w:xz") as tar:
            for subpath in ("inputs.json", "outputs/stable/nominal"):
                tar.add(a.work/subpath, arcname=subpath)
    else:
        assert (a.out/"verification.json").read_text() == body
    print(body)
