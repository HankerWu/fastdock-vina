"""Fixed four-compound PARP1 SPR set, two preselected receptors, two seeds."""

import argparse
import hashlib
import json
import shutil
import subprocess
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import gemmi
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
LIGANDS = {
    "talazoparib": ("4UND", "2YQ"),
    "olaparib": ("7KK4", "09L"),
    "niraparib": ("7KK5", "3JD"),
    "veliparib": ("7KK6", "78P"),
}
RECEIVERS = ("4UND", "7KK4")
SEEDS = (98207, 98213)


def prep(parent, work):
    chem = []
    for name, (pdb, ccd) in LIGANDS.items():
        dest = work / (name + ".pdbqt")
        if pdb in RECEIVERS:
            shutil.copyfile(parent / pdb / "ligand.pdbqt", dest)
            smi = json.loads((parent / pdb / "metadata.json").read_text())[
                "canonical_smiles"
            ]
            source = parent / pdb / (ccd + ".cif")
        else:
            source = work / (ccd + ".cif")
            if not source.exists():
                with urllib.request.urlopen(
                    f"https://files.rcsb.org/ligands/download/{ccd}.cif", timeout=60
                ) as res:
                    source.write_bytes(res.read())
            d = (
                gemmi.cif.read_file(str(source))
                .sole_block()
                .get_mmcif_category("_pdbx_chem_comp_descriptor.")
            )
            smi = next(
                s for t, s in zip(d["type"], d["descriptor"]) if t == "SMILES_CANONICAL"
            )
            graph = Chem.MolFromSmiles(smi)
            if Chem.GetFormalCharge(graph) == 0:
                ring_size = 6 if name == "niraparib" else 5
                ids = graph.GetSubstructMatches(
                    Chem.MolFromSmarts(f"[N;H1;+0;r{ring_size};X3;!a]")
                )
                assert len(ids) == 1
                n = graph.GetAtomWithIdx(ids[0][0])
                n.SetFormalCharge(1)
                n.SetNumExplicitHs(2)
                n.SetNoImplicit(True)
                Chem.SanitizeMol(graph)
            assert Chem.GetFormalCharge(graph) == 1
            smi = Chem.MolToSmiles(graph)
            mol = Chem.AddHs(Chem.MolFromSmiles(smi))
            params = AllChem.ETKDGv3()
            params.randomSeed = 43103
            assert (
                AllChem.EmbedMolecule(mol, params) == 0
                and AllChem.MMFFOptimizeMolecule(mol, maxIters=1000) == 0
            )
            (setup,) = MoleculePreparation().prepare(mol)
            block, ok, reason = PDBQTWriterLegacy.write_string(
                setup, add_index_map=True
            )
            assert ok, reason
            dest.write_text(block)
        chem.append(
            dict(
                compound=name,
                pdb=pdb,
                ccd=ccd,
                smiles=smi,
                charge=Chem.GetFormalCharge(Chem.MolFromSmiles(smi)),
                ccd_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                pdbqt_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),
            )
        )
    (work / "chemistry.json").write_text(json.dumps(chem, indent=2) + "\n")


def dock(parent, work, binary):
    out = work / "pools"
    out.mkdir(exist_ok=True)

    def one(job):
        name, rec, seed = job
        center = json.loads((parent / rec / "metadata.json").read_text())["center"]
        prefix = out / f"{name}-{rec}-{seed}"
        cmd = [
            str(binary.resolve()),
            "--receptor",
            str(parent / rec / "receptor.pdbqt"),
            "--ligand",
            str(work / (name + ".pdbqt")),
            "--exhaustiveness",
            "4",
            "--max_evals",
            "60000",
            "--energy_range",
            "100",
            "--num_modes",
            "9",
            "--seed",
            str(seed),
            "--cpu",
            "1",
            "--contact_rotor_score",
            "--out",
            str(prefix.with_suffix(".pdbqt")),
        ]
        for a, v in zip("xyz", center):
            cmd += ["--center_" + a, str(v), "--size_" + a, "24"]
        with prefix.with_suffix(".log").open("w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
        assert prefix.with_suffix(".pdbqt").read_bytes().endswith(b"ENDMDL\n")
        return prefix.name

    with ThreadPoolExecutor(4) as p:
        for name in p.map(
            one, [(l, r, s) for l in LIGANDS for r in RECEIVERS for s in SEEDS]
        ):
            print(name, flush=True)
    (work / "binary.sha256").write_text(
        hashlib.sha256(binary.read_bytes()).hexdigest() + "  vina\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("parent", type=Path)
    p.add_argument("--binary", type=Path)
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    prep(a.parent, a.work)
    dock(a.parent, a.work, a.binary)
