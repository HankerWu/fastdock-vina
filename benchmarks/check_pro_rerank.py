#!/usr/bin/env python3
"""Regression smoke for fixed-pool Vina Pro on the bundled 1IEP ligand."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import math
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory


EXAMPLE = Path(__file__).resolve().parents[1] / "examples/basic"
LIGAND = EXAMPLE / "1iep_ligand.pdbqt"
BOX = EXAMPLE / "1iep_receptor.box.txt"
COMMON = ["--receptor", str(EXAMPLE / "1iep_receptor.pdbqt"), "--cpu", "1",
          "--seed", "20261022", "--verbosity", "0"]


def models(path):
    return [model.splitlines() for model in path.read_text().strip().split("ENDMDL") if model.strip()]


def remark(model, key):
    prefix = f"REMARK {key}:"
    return next(line[len(prefix):].split() for line in model if line.startswith(prefix))


def atoms(model):
    return [line for line in model if line.startswith(("ATOM  ", "HETATM"))]


def heavy_coords(model):
    return [tuple(float(line[start:start + 8]) for start in (30, 38, 46))
            for line in atoms(model) if line[77:].strip() not in {"H", "HD"}]


def run(binary, args):
    result = subprocess.run([str(binary), *COMMON, *args], capture_output=True, text=True)
    return result


def check(binary, baseline, work):
    def dock(executable, label, *args, window=100, config=BOX, ligand=LIGAND, modes=9):
        output = work / f"{label}.pdbqt"
        result = run(executable, ["--ligand", str(ligand), "--config", str(config),
                                  "--exhaustiveness", "1", "--num_modes", str(modes), "--energy_range",
                                  str(window), "--out", str(output), *args])
        result.check_returncode()
        return output

    config = work / "pro.conf"
    config.write_text(BOX.read_text() + "\nrerank_gauss2_scale = 0.5\n")

    # Keep the MODEL body exactly as vina_split would, with prior Pro and AD4
    # score remarks added to exercise every reserved diagnostic label.
    split_ligand = work / "split_1iep.pdbqt"
    old_model = models(EXAMPLE / "1iep_docked.pdbqt")[0][1:]
    annotations = ["REMARK VINA PRO RESULT: -4.200 0.0 0.0", "REMARK VINA NATIVE: -13.264 1",
                   "REMARK VINA PRO GAUSS2 SCALE: 0.5", "REMARK CONF_INDEPENDENT: 1.5",
                   "REMARK VINA provenance: split 1IEP", "REMARK INTERACTION NOTE: keep this comment"]
    split_ligand.write_text("\n".join([*annotations, *old_model]) + "\n")
    quick_redock = {"ligand": split_ligand, "modes": 1}
    with ThreadPoolExecutor(max_workers=4) as pool:
        original = pool.submit(dock, baseline, "original")
        native = pool.submit(dock, binary, "native")
        identity = pool.submit(dock, binary, "identity", "--rerank_gauss2_scale", "1")
        zero = pool.submit(dock, binary, "zero", "--rerank_gauss2_scale", "0")
        half = pool.submit(dock, binary, "half", config=config)
        narrow = pool.submit(dock, binary, "narrow", "--rerank_gauss2_scale", "0", window=0.1, config=config)
        old_redock = pool.submit(dock, baseline, "old_redock", "--max_evals", "1000", **quick_redock)
        native_redock = pool.submit(dock, binary, "native_redock", "--max_evals", "1000", **quick_redock)
        identity_redock = pool.submit(dock, binary, "identity_redock", "--max_evals", "1000",
                                      "--rerank_gauss2_scale", "1", **quick_redock)
        pro_redock = pool.submit(dock, binary, "pro_redock", "--max_evals", "1000",
                                 "--rerank_gauss2_scale", "0", **quick_redock)

    assert original.result().read_bytes() == native.result().read_bytes() == identity.result().read_bytes()
    originals = models(native.result())
    assert len(originals) > 1

    for scale, path in ((0.0, zero.result()), (0.5, half.result())):
        ranked = models(path)
        pro_scores = []
        native_ranks = []
        best_coords = heavy_coords(ranked[0])
        unbound = remark(ranked[0], "INTRA")
        for index, model in enumerate(ranked, start=1):
            native_score, rank_text = remark(model, "VINA NATIVE")
            rank = int(rank_text)
            pro_score, lower, upper = map(float, remark(model, "VINA PRO RESULT"))
            assert f"MODEL {index}" in model
            assert float(remark(model, "VINA PRO GAUSS2 SCALE")[0]) == scale
            assert native_score == remark(originals[rank - 1], "VINA RESULT")[0]
            assert atoms(model) == atoms(originals[rank - 1])
            assert remark(model, "UNBOUND") == unbound
            coords = heavy_coords(model)
            rmsd = math.sqrt(sum(sum((a - b) ** 2 for a, b in zip(pt, ref))
                                 for pt, ref in zip(coords, best_coords)) / len(coords))
            assert 0 <= lower <= upper and abs(upper - rmsd) < 0.002
            pro_scores.append(pro_score)
            native_ranks.append(rank)
        assert pro_scores == sorted(pro_scores)
        assert sorted(native_ranks) == list(range(1, len(originals) + 1))

    zero_models = models(zero.result())
    assert int(remark(zero_models[0], "VINA NATIVE")[1]) != 1
    best = float(remark(zero_models[0], "VINA PRO RESULT")[0])
    expected_window = [model for model in zero_models
                       if float(remark(model, "VINA PRO RESULT")[0]) <= best + 0.1]
    assert models(narrow.result()) == expected_window

    assert old_redock.result().read_bytes() == native_redock.result().read_bytes() == identity_redock.result().read_bytes()
    pro_marks = ("REMARK VINA PRO RESULT:", "REMARK VINA NATIVE:", "REMARK VINA PRO GAUSS2 SCALE:",
                 "REMARK INTER + INTRA:", "REMARK INTER:", "REMARK INTRA:", "REMARK UNBOUND:")
    diagnostics = (*pro_marks, "REMARK VINA RESULT:", "REMARK CONF_INDEPENDENT:")
    redocked, = models(pro_redock.result())
    marks_seen = tuple(line.partition(":")[0] + ":" for line in redocked if line.startswith(diagnostics))
    assert marks_seen == pro_marks
    native_redocked, = models(native_redock.result())
    assert [line for line in redocked if not line.startswith(diagnostics)] == [
        line for line in native_redocked if not line.startswith(diagnostics)
    ]
    assert all(line in redocked for line in annotations[-2:])
    assert any(line.startswith("REMARK SMILES IDX ") for line in redocked)

    for flag in (["--no_refine"], ["--score_only"], ["--local_only"], ["--randomize_only"],
                 ["--scoring", "vinardo"], ["--flex", str(LIGAND)], ["--maps", "unused"],
                 ["--ligand", str(LIGAND), str(LIGAND)]):
        ligand = [] if flag[0] == "--ligand" else ["--ligand", str(LIGAND)]
        result = run(binary, ["--config", str(BOX), "--rerank_gauss2_scale", "0", *ligand, *flag])
        assert result.returncode and "rerank_gauss2_scale" in result.stderr
    for value in ("-0.1", "1.1", "nan", "inf"):
        result = run(binary, ["--config", str(BOX), "--ligand", str(LIGAND), "--rerank_gauss2_scale", value])
        assert result.returncode and "rerank_gauss2_scale" in result.stderr


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vina", type=Path, help="Vina Pro build")
    parser.add_argument("baseline", type=Path, help="pre-Pro Vina build")
    args = parser.parse_args()
    with TemporaryDirectory(prefix="fastdock-pro-smoke-") as work:
        check(args.vina.resolve(), args.baseline.resolve(), Path(work))
    print("Vina Pro smoke passed: native identity, fixed pool, rank, RMSD, window, redocking, and unsupported options")


if __name__ == "__main__":
    main()
