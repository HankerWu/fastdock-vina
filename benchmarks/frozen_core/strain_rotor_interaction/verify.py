"""Independently check contact counts with the compiled CLI and the unscaled joint-score formula."""
import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

from replay import BENCH, HERE, index, load

sys.path.insert(0, str(BENCH / "terminal_rotor_score"))
from assay import assay_inputs


def check(binary, directory):
    keys = ["target", "compound", "receptor", "seed"]
    row = index(load(directory / "transfer.csv"), keys)
    combined = index([r for r in load(directory / "first_poses.csv.xz")
                      if r["panel"] == "transfer_archived"], keys)
    max_error, n, unsupported = 0., 0, 0
    with TemporaryDirectory() as tmp:
        work = Path(tmp)
        for target, ligand, value, receptor, seed, protein, ligand_input, pool in assay_inputs():
            if target not in {"S6K1", "Gal3"}:
                continue
            key = target, ligand, receptor, str(seed)
            first = pool.decode().split("ENDMDL", 1)[0]
            if hashlib.sha256(pool).hexdigest() != row[key]["pool_sha256"]:
                raise ValueError("Pose pool identity mismatch")
            match = lambda field: float(re.search(r"^REMARK " + field + r":\s+(-?[\d.]+)", first, re.M)[1])
            inter, unbound = match("INTER"), match("UNBOUND")
            (work / "rec").write_bytes(protein)
            (work / "lig").write_bytes(ligand_input)
            (work / "pose").write_text("\n".join(l for l in first.splitlines() if not l.startswith("MODEL")) + "\n")
            cmd = [str(binary), "--score_only", "--receptor", "rec", "--ligand", "lig", "--pose_coordinates", "pose",
                   "--contact_rotor_score", "--unbound_energy", str(unbound), "--cpu", "1", "--seed", "1"]
            result = subprocess.run(cmd, text=True, capture_output=True, cwd=work)
            if row[key]["rotor_status"] != "ok":
                if result.returncode == 0 or "hydration or macrocycle pseudo-atoms" not in result.stderr or key in combined:
                    raise ValueError(f"Expected C++ domain rejection for {key}")
                unsupported += 1
                continue
            result.check_returncode()
            out = result.stdout
            contact = re.search(r"^REMARK VINA CONTACT ROTOR SCORE:\s+\S+\s+(\d+)\s+(\d+)", out, re.M)
            if contact is None or int(contact[1]) != int(row[key]["constrained"]):
                raise ValueError(f"C++ contact count mismatch at {key}: {contact}")
            strain = float(row[key]["fixed_mmff"]) - float(row[key]["free_mmff"])
            direct = inter + strain + int(contact[1]) * .00198720425864 * 298.15 * math.log(3)
            max_error = max(max_error, abs(direct - float(combined[key]["joint"])))
            n += 1
    if n + unsupported != len(row) or n != len(combined) or max_error > .00151:
        raise ValueError(f"Incomplete or inconsistent first-pose reconstruction: {n}, {max_error}")
    print(json.dumps(dict(pools_attempted=len(row), pools_verified=n, expected_cpp_abstentions=unsupported, cpp_contact_mismatches=0, max_abs_direct_joint_error=max_error,
                          binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest()), indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("data", nargs="?", type=Path, default=HERE)
    a = p.parse_args()
    check(a.binary.resolve(), a.data)
