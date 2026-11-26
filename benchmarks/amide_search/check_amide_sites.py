"""Regression checks for graph-based amide site generation and pose ranking."""

import sys
import unittest
from pathlib import Path

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from amide_rotamers import active_amide_sites, amides, search_amide_records
from rank_amide_pool import rank


def prepare(smiles, flexible=True):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=360) == 0
    AllChem.UFFOptimizeMolecule(mol)
    setup = MoleculePreparation(flexible_amides=flexible)(mol)[0]
    text, ok, error = PDBQTWriterLegacy.write_string(setup)
    assert ok, error
    return text


def pool_model(ligand, number, coords):
    records = [f"MODEL {number}", "REMARK VINA RESULT: -5.0 0.0 0.0"]
    for line in ligand.splitlines():
        if line.startswith(("ATOM ", "HETATM ")):
            serial = int(line[6:11])
            if serial in coords:
                line = (
                    line[:30] + "".join(f"{x:8.3f}" for x in coords[serial]) + line[54:]
                )
        records.append(line)
    return "\n".join(records + ["ENDMDL"]) + "\n"


class AmideSites(unittest.TestCase):
    def test_recognition(self):
        # Arity of the measured nitrogen side of each flexible site.
        for smiles, arms in (
            ("CC(=O)N(C)C", [2]),  # equivalent methyls
            ("CC(=O)N(C)CC", [2]),
            ("CC(=O)NC", [1]),
            ("CC(=O)N1CCCCC1", [2]),  # exocyclic acyl C–N bond
            ("CC(=O)N(C)CCNC(C)=O", [1, 2]),
            ("CC(=O)N", []),  # primary
            ("CC(=O)[NH+](C)C", []),  # charged
            ("O=C1NCCC1", []),  # lactam
            ("COC(=O)N(C)C", []),  # carbamate
            ("CNC(=O)NC", []),  # urea
            ("CC(=O)N(C)C(=O)C", []),  # N-substituted imide
            ("CCC(=O)NC(=O)C", []),  # NH imide
        ):
            with self.subTest(smiles=smiles):
                sites = active_amide_sites(prepare(smiles), "all-arms")
                self.assertEqual(sorted(len(site) - 3 for site in sites), arms)
                self.assertTrue(
                    all(tuple(sorted(site[3:])) == site[3:] for site in sites)
                )

        fixed = prepare("CC(=O)NC", flexible=False)
        self.assertTrue(amides(fixed))
        self.assertEqual(active_amide_sites(fixed, "all-arms"), [])
        with self.assertRaisesRegex(ValueError, "No eligible"):
            search_amide_records(fixed, 5.0)

    def test_symmetric_site_energy_and_serials(self):
        ligand = prepare("CC(=O)N(C)C")
        [site] = active_amide_sites(ligand, "all-arms")
        o, c, n, x1, x2 = site
        self.assertEqual(active_amide_sites(ligand, "distinguished"), [])
        self.assertEqual(
            search_amide_records(ligand, 8.0).split(), [*map(str, site), "8.0"]
        )
        # One arm planar and the other perpendicular. Swapping the labels has no effect.
        coords = {
            c: (0, 0, 0),
            n: (1, 0, 0),
            o: (0, 1, 0),
            x1: (1, 1, 0),
            x2: (1, 0, 1),
        }
        pool = pool_model(ligand, 1, coords)
        coords[x1], coords[x2] = coords[x2], coords[x1]
        pool += pool_model(ligand, 2, coords)
        self.assertEqual([r["score"] for r in rank(ligand, pool, 8.0)], [-5.0, -5.0])
        for row in rank(ligand, pool, 8.0, "all-arms"):
            self.assertAlmostEqual(row["amide_penalty"], 4.0)
            self.assertAlmostEqual(row["score"], -1.0)
        self.assertEqual(
            [r["score"] for r in rank(ligand, pool, 0.0, "all-arms")], [-5.0, -5.0]
        )
        for bad in (-1.0, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                search_amide_records(ligand, bad)


if __name__ == "__main__":
    unittest.main()
