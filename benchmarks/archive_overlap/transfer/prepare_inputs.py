"""Prepare the registered human complexes and preserve failures and source hashes."""

import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "pair_gradient/pose"))
import prepare

SYSTEMS = {
    "3P6D": ("FABP4", "ZGB", "A"),
    "3P6E": ("FABP4", "ZGC", "A"),
    "4GQ3": ("MEN1", "0RO", "A"),
    "4GQ4": ("MEN1", "0RT", "A"),
    "4HG7": ("MDM2", "NUT", "A"),
    "1UY6": ("HSP90", "PU3", "A"),
    "3H82": ("EPAS1", "020", "A"),
    "3F1O": ("EPAS1", "2XY", "A"),
}


def main(work, systems=SYSTEMS):
    import os

    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    prepare.CASES.update(systems)
    sources = []
    for code, (target, component, chain) in systems.items():
        start = time.monotonic()
        path = work / code / "preparation.json"
        try:
            prepare.prepare(work, code)
            meta = json.loads(path.read_text())
            status = meta["receptor_status"]
        except (ValueError, RuntimeError, OSError) as exc:
            status = f"preparation failed: {exc}"
        raw = work / "raw"
        files = [raw / f"{code}.pdb", raw / f"{component}.cif"]
        # Chain A is the target protein in both prespecified panels. Record its
        # deposited source so an incorrect organism guess cannot become metadata.
        record = (
            re.search(
                r"^SOURCE.+?ORGANISM_SCIENTIFIC: ([^;]+);",
                files[0].read_text(),
                re.MULTILINE,
            )
            if files[0].exists()
            else None
        )
        sources.append(
            dict(
                pdb=code,
                target=target,
                deposited_organism=record.group(1) if record else None,
                component=component,
                chain=chain,
                status=status,
                preparation_wall_s=time.monotonic() - start,
                source=f"https://www.rcsb.org/structure/{code}",
                ccd=f"https://files.rcsb.org/ligands/download/{component}.cif",
                sha256={
                    p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in files
                    if p.exists()
                },
            )
        )
        (work / "sources.json").write_text(json.dumps(sources, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    main(parser.parse_args().work.resolve())
