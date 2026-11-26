"""Replay every control with the previous release and compare exact model bytes."""

import csv
import hashlib
import json
import sys
import tempfile
from pathlib import Path
from run import run

work, binary = map(lambda p: Path(p).resolve(), sys.argv[1:])
with tempfile.TemporaryDirectory() as folder:
    root = Path(folder)
    for name in ("panel.csv", "prepared"):
        (root / name).symlink_to(work / name)
    run(root, binary, 3, ["native"], ["SYK", "FAK"])
    count = 0
    for p in (root / "pools").glob("*.pdbqt"):
        assert p.read_bytes() == (work / "pools" / p.name).read_bytes(), p.name
        count += 1
    assert count == 52
(work / "native_parity.json").write_text(
    json.dumps(
        dict(
            cells=count,
            byte_identical=count,
            baseline="ce9a89fcb41c7424a30340546be50a972570b876",
            binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        ),
        indent=2,
    )
    + "\n"
)
