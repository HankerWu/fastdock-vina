"""Run a real single-case regression using existing example inputs; not efficacy evidence."""
import argparse
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

def run(binary, args, expect_ok=True):
    result = subprocess.run([str(binary), *args], capture_output=True, text=True)
    if expect_ok and result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    if not expect_ok and not result.returncode:
        raise AssertionError("expected rejection: " + str(args))
    return result.stdout + result.stderr

def counters(log):
    return {int(m.group(1)): m.group(0) for m in re.finditer(
        r"Search trajectory: index=(\d+)[^\n]*", log)}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    binary = parser.parse_args().binary.resolve()
    rec = ROOT / "examples/basic/1iep_receptor.pdbqt"
    lig = ROOT / "examples/basic/1iep_ligand.pdbqt"
    box = ROOT / "examples/basic/1iep_receptor.box.txt"
    with tempfile.TemporaryDirectory() as temp:
        temp = Path(temp)
        common = ["--receptor", str(rec), "--ligand", str(lig), "--config", str(box),
                  "--pose_selection", "late", "--exhaustiveness", "4", "--max_evals", "3000",
                  "--num_modes", "9", "--refine_candidates", "18", "--energy_range", "100",
                  "--seed", "20261126", "--search_counters"]
        def dock(label, flags):
            out = temp / (label + ".pdbqt")
            log = run(binary, common + flags + ["--out", str(out)])
            data = out.read_text()
            assert "REMARK VINA RESULT:" in data and "MODEL" in data
            return log, data
        baseline, _ = dock("baseline", ["--cpu", "1"])
        mixed, poses = dock("mixed", ["--cpu", "1", "--mixed_vinardo_search"])
        repeated, again = dock("repeat", ["--cpu", "2", "--mixed_vinardo_search"])
        assert "Mixed Vinardo search: vina=2 vinardo=2" in mixed
        a, b, c = map(counters, (baseline, mixed, repeated))
        assert set(a) == set(b) == set(c) == {1, 2, 3, 4}
        assert all(a[i] == b[i] for i in (1, 3)), (a, b)
        assert b == c and poses == again, "thread-count reproducibility"
        scores = [float(x) for x in re.findall(r"REMARK VINA RESULT:\s+([-+0-9.eE]+)", poses)]
        assert scores == sorted(scores)
        assert len(scores) > 0
        for flags in (["--mixed_cutoff_search"], ["--scoring", "vinardo"], ["--smooth_cutoff"]):
            log = run(binary, common + ["--cpu", "1", "--mixed_vinardo_search", *flags,
                                         "--out", str(temp / "invalid.pdbqt")], False)
            assert "mixed_vinardo_search" in log, log
        print("Mixed Vinardo regression passed: native odd trajectories, thread reproducibility, score order, incompatible modes.")

if __name__ == "__main__":
    main()
