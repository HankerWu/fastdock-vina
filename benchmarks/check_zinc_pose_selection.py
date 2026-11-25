#!/usr/bin/env python3
"""Exercise zinc pose selection through its CLI using an ideal tetrahedral site."""

import json
import math
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


TOOL = Path(__file__).resolve().parents[1] / "tools/select_zinc_pose.py"


def atom(serial, atom_type, xyz, residue="LIG"):
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {atom_type:>4} {residue:>3} A   1    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00    0.000 {atom_type:>3}\n")


def receptor(plane_height):
    rows = [atom(10, "Zn", (0, 0, 0), "ZN")]
    for index in range(3):
        phi = 2 * math.pi * index / 3
        radius = math.sqrt(2 ** 2 - plane_height ** 2)
        rows.append(atom(11 + index, "NA", (radius * math.cos(phi), radius * math.sin(phi), plane_height), "HIS"))
    return "".join(rows)


def pose(number, xyz, score, extra="", pro=None, donor_type="OA"):
    if pro is None:
        header = f"REMARK VINA RESULT: {score} 0.000 0.000\n"
    else:
        pro_score, rank = pro
        header = (f"REMARK VINA PRO RESULT: {pro_score} 0.000 0.000\n"
                  f"REMARK VINA NATIVE: {score} {rank}\nREMARK VINA PRO GAUSS2 SCALE: 0.5\n")
    return f"MODEL {number}\n{header}{extra}{atom(5, donor_type, xyz)}ENDMDL\n"


def check(work):
    receptor_path = work / "receptor.pdbqt"
    poses_path = work / "poses.pdbqt"
    tetrahedral = receptor(-2 / 3)
    receptor_path.write_text(tetrahedral)

    def select(poses, mode=None):
        poses_path.write_text(poses)
        command = [sys.executable, str(TOOL), "--receptor", str(receptor_path),
                   "--poses", str(poses_path), "--zinc-serial", "10", "--donor-serial", "5"]
        if mode is not None:
            command.extend(["--mode", mode])
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return json.loads(result.stdout)

    bad_direction = pose(1, (2.10, 0, 0), -12)
    perfect = (0, 0, 2.10)
    result = select(bad_direction + pose(2, perfect, -8, extra="REMARK VINA RESULT: 50 0 0\n") +
                    pose(3, (0, 0, 2.25), -7) + pose(4, (0, 0, 3.00), -13) +
                    pose(5, perfect, -10, extra=atom(6, "NA", (1.95, 0, 0))))
    assert result["status"] == "selected" and result["selected_model"] == 2
    assert result["site"]["applicable"] and result["site"]["vacant_direction"] == [0, 0, 1]
    assert [row["reason"] for row in result["poses"]] == [
        "donor_outside_open_direction", None, None, "donor_distance_outside_window", "other_ligand_coordinators_in_shell"
    ]
    selected = result["poses"][1]
    assert selected["native_score"] == -8 and selected["ignored_trailing_score_headers"] == 1
    assert selected["zinc_distance_angstrom"] == 2.10 and selected["direction_error_deg"] == 0
    assert selected["virtual_point_error_angstrom"] == 0

    stale = "REMARK VINA PRO RESULT: 99.00 0 0\nREMARK VINA NATIVE: 99.00 99\nREMARK VINA RESULT: 99.00 0 0\n"
    pro_runs = pose(1, perfect, -8, pro=(-20, 2)) + pose(2, perfect, -9, extra=stale, pro=(-4, 1))
    result = select(pro_runs)
    assert result["status"] == "selected" and result["selected_model"] == 2
    assert result["poses"][1]["pro_score"] == -4 and result["poses"][1]["ignored_trailing_score_headers"] == 3
    assert result == select(pro_runs, "hard")

    result = select(bad_direction)
    assert result["status"] == "abstain" and result["reason"] == "no_geometry_pass" and result["selected_model"] is None
    assert select("")["reason"] == "no_candidates"

    receptor_path.write_text(receptor(0))
    result = select(pose(1, perfect, -8))
    assert result["status"] == "abstain" and result["reason"] == "no_applicable_site"
    assert result["site"]["reason"] == "zinc_too_close_to_coordinator_plane"

    receptor_path.write_text(tetrahedral + atom(15, "OA", perfect, "HOH"))
    result = select(pose(1, perfect, -8))
    assert result["status"] == "abstain" and len(result["site"]["coordinators"]) == 4

    receptor_path.write_text(tetrahedral + atom(16, "Cl", (0, 0, 2.28), "CL"))
    result = select(pose(1, perfect, -8))
    assert result["status"] == "abstain" and result["site"]["reason"] == "nearby_receptor_halogen"
    assert result["site"]["nearby_halogen_serials"] == [16]

    receptor_path.write_text(tetrahedral)
    result = select(pose(1, perfect, -10, donor_type="Na"))
    assert result["status"] == "abstain"
    assert result["poses"][0]["reason"] == "donor_missing_or_not_n_o_s"

    result = select(pose(1, perfect, -8, extra=atom(6, "Na", (1.95, 0, 0)), donor_type="NA"))
    assert result["status"] == "selected" and result["selected_model"] == 1
    assert result["poses"][0]["donor_element"] == "N"
    assert result["poses"][0]["other_nearby_coordinator_serials"] == []

    receptor_path.write_text(tetrahedral + atom(17, "Na", (0, 2.20, 0), "NA"))
    result = select(pose(1, perfect, -8))
    assert result["status"] == "selected"
    assert [row["serial"] for row in result["site"]["coordinators"]] == [11, 12, 13]

    receptor_path.write_text(tetrahedral)
    nitrogen_poses = (
        pose(1, (2.00, 0, 0), -9.4, donor_type="NA")
        + pose(2, (0, 0, 2.564), -9.0, donor_type="N")
        + pose(3, (1.20, 0, 1.60), -8.6, donor_type="NA")
        + pose(4, (0, 0, 9.00), -9.2, donor_type="NA")
        + pose(5, (0, 0, 1.649), -30, donor_type="NA")
        + pose(6, (0, 0, 2.564), -30, donor_type="NA", extra=atom(6, "OA", (2.50, 0, 0)))
        + pose(7, (0, 0, 2.00), -30, donor_type="OA")
        + pose(8, (0, 0, 2.00), -30, donor_type="Na")
        + pose(9, (0, 0, 1.65), -8.5, donor_type="NA")
    )
    continuous = select(nitrogen_poses, "continuous")
    radial = select(nitrogen_poses, "radial")
    assert continuous["status"] == "selected" and continuous["selected_model"] == 2
    assert radial["status"] == "selected" and radial["selected_model"] == 1
    assert continuous["poses"][0]["reason"] == "donor_outside_open_direction"
    assert continuous["poses"][1]["reason"] == "donor_distance_outside_window"
    assert continuous["poses"][5]["reason"] == "donor_distance_outside_window"
    assert select(nitrogen_poses)["selected_model"] == 7
    for result, mode in [(continuous, "continuous"), (radial, "radial")]:
        assert result["mode"] == mode and result["reason"] is None
        assert [row["model"] for row in result["poses"]] == list(range(1, 10))
        assert [row["scoring_eligible"] for row in result["poses"]] == [True] * 4 + [False] * 4 + [True]
        assert [row["scoring_reason"] for row in result["poses"][4:8]] == [
            "donor_distance_below_minimum", "other_ligand_coordinators_in_shell",
            "donor_missing_or_not_nitrogen", "donor_missing_or_not_nitrogen",
        ]
        for row in result["poses"][4:8]:
            assert row["scoring_delta_angstrom"] is row["geometry_reward"] is row["ranking_score"] is None
        far_pose = result["poses"][3]
        assert 0 < far_pose["geometry_reward"] < 1e-9
        assert far_pose["ranking_score"] < far_pose["native_score"]

    for index, directional_delta, radial_delta in [(0, math.sqrt(8), 0), (1, 0.564, 0.564),
                                                    (2, math.sqrt(1.6), 0)]:
        for result, delta in [(continuous, directional_delta), (radial, radial_delta)]:
            row = result["poses"][index]
            reward = 3 * math.exp(-(delta ** 2) / 2)
            assert math.isclose(row["scoring_delta_angstrom"], delta, abs_tol=1e-12)
            assert math.isclose(row["geometry_reward"], reward, rel_tol=1e-12)
            assert math.isclose(row["ranking_score"], row["native_score"] - reward, rel_tol=1e-12)

    tied_pro_poses = (pose(1, (0, 0, 2.00), -8, pro=(-21, 2), donor_type="NA")
                     + pose(2, (0, 0, 2.00), -8, pro=(-4, 1), donor_type="N"))
    for mode in ("continuous", "radial"):
        assert select(tied_pro_poses, mode)["selected_model"] == 2
        no_nitrogen = select(pose(1, perfect, -8), mode)
        assert no_nitrogen["status"] == "abstain" and no_nitrogen["selected_model"] is None
        assert no_nitrogen["reason"] == "no_scoring_eligible_pose"
        assert select("", mode)["reason"] == "no_candidates"

    receptor_path.write_text(receptor(0))
    invalid_site = select(pose(1, perfect, -8, donor_type="NA"), "continuous")
    assert invalid_site["reason"] == "no_applicable_site"
    assert invalid_site["poses"][0]["scoring_reason"] == "site_inapplicable"


if __name__ == "__main__":
    with TemporaryDirectory(prefix="fastdock-zinc-smoke-") as directory:
        check(Path(directory))
    print("Zinc pose selection smoke passed: hard geometry, continuous and radial scores, native ties, and abstention")
