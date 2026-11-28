"""Verify the prior FKBP refinement claim directly from the native logs and coordinate hashes."""

import csv
import hashlib
import json
import re
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np

SOURCE = Path(__file__).resolve().parents[2] / "affinity_conditioning/results"


def verify():
    archive = SOURCE / "refine_raw.tar.xz"
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    assert digest == json.loads((SOURCE / "refine_receipt.json").read_text())["sha256_raw"]
    native, count = {}, 0
    with (SOURCE / "refine_jobs.csv").open() as file, tarfile.open(archive) as tar:
        for job in csv.DictReader(file):

            def read(name, identity=job["job"]):
                return tar.extractfile(f"evidence/{identity}/{name}").read()

            assert read("explicit.graph") == read("inferred.graph")
            for name, field in (
                ("input.pdbqt", "input_sha256"),
                ("output.pdbqt", "output_sha256"),
                ("inferred.graph", "graph_sha256"),
            ):
                assert hashlib.sha256(read(name)).hexdigest() == job[field]
            observed = re.findall(rb"REFINE_OBSERVE\s+\S+\s+\S+\s+(\d+)", read("observer.log"))
            calls = sum(int(n) for n in observed)
            assert calls == int(job["observed_gradient_calls"])
            count += calls
            values = re.findall(rb"Estimated Free Energy of Binding\s+: (\S+)", read("native.log"))
            assert len(values) == 2
            native[job["job"]] = np.array(list(map(float, values)))
    groups = defaultdict(dict)
    with (SOURCE / "refine_legs.csv").open() as file:
        for row in csv.DictReader(file):
            key = tuple(row[k] for k in ("panel", "arm", "seed", "receiver"))
            groups[key][row["leg"]] = native[row["job"]]
    result = {
        "verified_raw_sha256": digest,
        "unique_jobs_with_verified_input_output_and_graph_hashes": len(native),
        "observed_gradient_calls": count,
        "cells": len(groups),
    }
    for label, leg in (("net", "native5"), ("substitution", "feasible5")):
        a = np.array([v["native6"] - v[leg] for v in groups.values()])
        result["negative_" + label + "_before_after"] = (a < 0).sum(0).tolist()
        result[label + "_direction_flips" if label == "net" else "substitution_flips"] = int(
            (np.sign(a[:, 0]) != np.sign(a[:, 1])).sum()
        )
    return result


if __name__ == "__main__":
    data = verify()
    path = Path(__file__).resolve().parent / "results/previous_replays.json"
    assert data == json.loads(path.read_text())["previous_fkbp12_refinement_log_check"]
    print(json.dumps(data, indent=2))
