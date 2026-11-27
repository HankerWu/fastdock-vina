"""Audit frozen candidate ranks, AD4 trajectories, map-derived public scores and Table II directions."""
import argparse
from pathlib import Path
import re
import sys

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from frozen_assay import ARMS, CASES, ROOT, old, table
from check_assay import rms, rows

TRACE_ARMS = ("match", *ARMS)
ALL_ARMS = ("official", "base", "match", *ARMS, "high", "official_high")


def parse_xyz(p, expected):
    lines = p.read_text().splitlines()
    assert [int(v) for v in lines[0][1:].split()] == list(expected)
    for line in lines[1:]:
        tag, scope, pos, text = line.split(" ", 3)
        yield tag, int(scope), int(pos), np.fromstring(text, sep=" ").reshape(-1,3)


def rank_array(x):
    idx = np.argsort(x, kind="stable")
    rank = np.empty(len(idx), dtype=int)
    rank[idx] = np.arange(1,len(x)+1)
    return rank


def trace(work, out):
    refs = {r:old.reference(work,r,ccd) for r,ccd in (("1FKG","SB3"),("1FKH","SBX"))}
    pool_stats, details, stage_rows, candidate_rows, contacts = [], [], [], [], []
    for rec, n in CASES:
        stem = f"{rec}_{n}"
        lig = (work / "affinity" / f"{n}.pdbqt").read_text()
        mol, chemical_order, _ = old.molecule(lig)
        ccd_source = {8:"1FKG", 9:"1FKH"}.get(n)
        native = old.crystal_positions(refs[ccd_source], refs[rec], mol) if ccd_source else None
        recxyz = np.array([[float(l[i:i+8]) for i in (30,38,46)]
            for l in (work / rec / "receptor.pdbqt").read_text().splitlines()
            if l.startswith(("ATOM", "HETATM")) and l[77:].strip() not in ("H","HD","HS")])
        tree = cKDTree(recxyz)
        with np.load(work/"plans"/(stem+".npz")) as bank, np.load(work/"ranks"/(stem+".npz")) as scores:
            serials = bank["serials"].tolist()
            order = [serials.index(i) for i in chemical_order]
            xyz = bank["xyz"]
            poses = (xyz-xyz.mean(0)) @ bank["matched_rotation"].transpose(0,2,1) + bank["center"][:,None]
            shape = bank["matched_score"]
            dist = rms(poses[:, order], native) if ccd_source else None
            ranks = {arm:rank_array(np.round(shape,10) if arm=="match" else scores[arm]) for arm in TRACE_ARMS}
            winners = {arm:int(np.argmin(v)) for arm,v in ranks.items()}
            assert winners["match"] == int(bank["matched_best"])
            near_contact = tree.query(poses)[0].min(axis=1)
            for i in range(len(poses)):
                candidate_rows.append(dict(receiver=rec, compound=n, index=i,
                    rmsd_A=dist[i] if ccd_source else "", min_heavy_contact_A=near_contact[i],
                    shape_score=shape[i], typed_kcal_mol=scores["typed"][i], perm_kcal_mol=scores["perm"][i],
                    **{arm+"_rank":rank[i] for arm,rank in ranks.items()}))
            for arm in TRACE_ARMS:
                i = winners[arm]
                p = work/"pools"/arm/stem
                snapshots = list(parse_xyz(p.with_suffix(".xyz"), serials))
                merged = list(parse_xyz(work / "merge" / f"{stem}_{arm}.xyz", serials))
                stages = {tag:x for tag,scope,pos,x in snapshots if tag in ("selected","hunt","second")}
                np.testing.assert_allclose(stages["selected"], poses[i],rtol=0,atol=1e-10)
                history = {int(l.split()[0]):[int(v) for v in l.split()[4:]]
                           for l in p.with_suffix(".history").read_text().splitlines()}
                released = np.array([old.xyz(b,serials) for b in p.with_suffix(".pdbqt").read_text().split("ENDMDL") if b.strip()])
                def eq(x): return np.array_equal(x,stages["second"])
                end_at = [pos for tag,tr,pos,x in snapshots if tag=="end_trajectory" and tr==1 and eq(x)]
                merged_at = [pos for _,_,pos,x in merged if eq(x)]
                final_at = (1 + np.flatnonzero(np.max(abs(released-stages["second"]), axis=(1,2))<=.00051)).tolist()
                assert len(end_at)<=1 and len(merged_at)<=1 and len(final_at)<=1 and bool(1 in history[max(history)])==bool(end_at)
                assert not merged_at or end_at
                assert not final_at or merged_at
                d = dist if ccd_source else np.array([])
                accepted = [(pos,rms(x[order], native)) for tag,tr,pos,x in snapshots if ccd_source and tag=="offered" and tr==1]
                def first(mask):
                    return int(min(ranks[arm][mask])) if np.any(mask) else ""
                info = dict(receiver=rec, compound=n, arm=arm, selected_index=i,
                    changed_types=int(np.sum(scores["true_types"]!=scores["perm_types"])),
                    nheavy=len(order), min_heavy_contact_A=near_contact[i],
                    generated_le2=sum(d<=2) if ccd_source else "", generated_le3=sum(d<=3) if ccd_source else "",
                    first_le2_rank=first(d<=2), first_le3_rank=first(d<=3),
                    best_generated_A=min(d) if ccd_source else "", selected_A=d[i] if ccd_source else "",
                    hunt_A=float(rms(stages["hunt"][order],native)) if ccd_source else "",
                    second_A=float(rms(stages["second"][order],native)) if ccd_source else "",
                    entered_local_id1=int(1 in history[1]), first_evict_step=next((step for step,h in history.items() if 1 not in h), ""),
                    end_trajectory_rank=end_at[0] if end_at else "", merged_rank=merged_at[0] if merged_at else "",
                    published_rank=final_at[0] if final_at else "",
                    later_first_offered_le2_step=next((step for step,err in accepted if step>1 and err<=2),""),
                    later_first_offered_le3_step=next((step for step,err in accepted if step>1 and err<=3),""))
                details.append(info)
                for stage,scope,rank,x in [*snapshots,*merged,*[("published",0,k+1,x) for k,x in enumerate(released)]]:
                    if stage=="published":
                        contacts.append(dict(receiver=rec, compound=n, arm=arm, rank=rank,
                                             min_heavy_distance_A=tree.query(x)[0].min()))
                    if ccd_source:
                        stage_rows.append(dict(receiver=rec, compound=n, arm=arm, stage=stage,
                            scope=scope, step_or_rank=rank, rmsd_A=float(rms(x[order],native))))
        for arm in ARMS:
            lines = (work/"pools"/arm/(stem+".log")).read_text()
            counts = {(int(c),int(d)):int(x) for c,d,x in re.findall(
                        r"FIELD_COUNT combined=(\d+) derivative=(\d+) positions=(\d+)", lines)}
            assert len(counts)==4 and counts[0,0]==counts[0,1]==0
            pool_stats.append(dict(receiver=rec, compound=n, arm=arm,
                ad4_scalar_positions=counts[1,0], ad4_gradient_positions=counts[1,1],
                downstream_channel_interpolations=sum(counts.values())*3))
    table(out/"stages.csv",details)
    table(out/"downstream_field.csv",pool_stats)
    table(out/"min_contacts.csv",contacts)
    table(work/"candidate_audit.csv",candidate_rows)
    table(work/"stage_audit.csv",stage_rows)


def check(work, out):
    old.inspect(work,out, arms=ALL_ARMS, comparison="typed")
    trace(work,out)
    scores = {(r["arm"],int(r["compound"])):float(r["top1_score"]) for r in rows(out/"scores.csv") if r["receiver"]=="mean"}
    desc = rows(out / "descriptors.csv")
    comparisons = []
    for arm, control in (("typed","base"),("typed","perm"),("typed","match"),("perm","match")):
        dx = [scores[arm,int(d["compound"])]-scores[control,int(d["compound"])] for d in desc]
        comparisons.append(dict(arm=arm,control=control,n=len(desc),
          rho_size=spearmanr(dx,[int(d["heavy_atoms"]) for d in desc]).statistic,
          rho_torsdof=spearmanr(dx,[int(d["torsdof"]) for d in desc]).statistic,
          n_distinct_charges=len(set(d["formal_charge"] for d in desc))))
    table(out / "bias_contrasts.csv",comparisons)
    cells = rows(out / "pose_cells.csv")
    by_mol = []
    for arm in ALL_ARMS:
        for n in (8,9):
            sub = [r for r in cells if r["arm"]==arm and int(r["compound"])==n]
            a = np.array([float(r["top1"]) for r in sub])
            q = np.quantile(a,[.25,.5,.75,.9,.95])
            by_mol.append(dict(arm=arm,compound=n,n_receivers=len(a),top1_mean=a.mean(),
                **dict(zip(("q25","median","q75","p90","p95"),q,strict=True)),
                **{k: np.mean([float(r[k]) for r in sub]) for k in ("best5","oracle","regret")},
                gt5=sum(a>5),gt10=sum(a>10)))
    table(out/"per_molecule.csv",by_mol)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    check(a.work.resolve(), a.out.resolve())
