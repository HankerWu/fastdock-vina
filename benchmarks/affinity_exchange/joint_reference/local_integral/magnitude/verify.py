"""Recompute witnesses and a full continuous-domain cover from raw seed scores."""
import argparse
import csv
import json
import sys
from collections import Counter
from decimal import Decimal as D, getcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
# Load the previous independent verifier under an unambiguous name.
from importlib.util import module_from_spec, spec_from_file_location
spec = spec_from_file_location("parent_verifier", HERE.parent / "verify.py")
raw = module_from_spec(spec)
spec.loader.exec_module(raw)

getcontext().prec = 50
RT = raw.RT


def verify(r, source, sem):
    energies, kds = source
    ts = r["structures"]
    arm = r["arm"]
    axes = r["axes"]
    assert set(axes) == {t for t in ts if len(ts[t]) == 2}
    keys = {k for k in kds if k[0] in ts}
    assert len(keys) == r["cells"]
    for t in ts:
        assert sorted({s for _, target, _, s in energies if target == t}) == ts[t]

    def energies_at(t, ligand, p):
        e = [D(0) if arm == "zero" else energies[arm, t, ligand, s] for s in ts[t]]
        if len(e) == 1:
            return e[0]
        f = p[axes.index(t)]
        return -RT * (((1-f)*(-e[0]/RT).exp()) + f*(-e[1]/RT).exp()).ln()

    def observed(key):
        kd = kds[key]
        err = D(0) if r["observation"] == "mean" else sem[key]
        return tuple(RT*(x*D("1e-9")).ln() for x in (kd-err, kd+err))

    upper = D(str(r["upper"]))
    for t, i in keys:
        y1, y2 = observed((t, i))
        model = (energies_at(t, i, r["population"]) + r["target_offsets"][t]
                 + r["ligand_offsets"][str(i)])
        assert max(y1-model, model-y2, D(0)) <= upper

    paths = set()
    cover_area = D(0)
    for b in [r["full_domain"], *r["boxes"]]:
        path, terms = b["path"], b["cycle"]
        if b is not r["full_domain"]:
            assert path not in paths
            paths.add(path)
        else:
            assert path == ""
        left, right = [D(0)]*len(axes), [D(1)]*len(axes)
        for n, bit in enumerate(path):
            col = n % len(axes)
            assert bit in "01"
            (left if bit == "1" else right)[col] = (left[col]+right[col])/2
        if b is not r["full_domain"]:
            cover_area += D(1) / (2**len(path))
        assert len(terms) == len({t for t, _, _ in terms})
        assert Counter(i for _, i, _ in terms) == Counter(j for _, _, j in terms)
        bottom = top = y_lower = y_upper = model_lower = model_upper = D(0)
        for t, i, j in terms:
            assert i != j
            ilo, ihi = observed((t, i))
            jlo, jhi = observed((t, j))
            d1, d2 = (energies_at(t, i, p) - energies_at(t, j, p) for p in (left, right))
            y_lower += ilo-jhi
            y_upper += ihi-jlo
            model_lower += min(d1,d2)
            model_upper += max(d1,d2)
            bottom += ilo-jhi-max(d1,d2)
            top += ihi-jlo-min(d1,d2)
        assert max(bottom, -top, D(0)) / (2*len(terms)) >= b["lower"]
        if "positive_scale_lower" in b:
            assert not path
            assert y_upper < 0 < model_lower or model_upper < 0 < y_lower
            assert min(abs(y_lower), abs(y_upper)) / (2*len(terms)) >= b["positive_scale_lower"]
    assert cover_area == 1
    assert all(path[:i] not in paths for path in paths for i in range(len(path)))
    assert r["lower"] == min(b["lower"] for b in r["boxes"])
    assert D(0) <= upper-r["lower"] < D(".001")


def main(dest):
    reports = json.loads((dest / "witnesses.json").read_text(), parse_float=D)
    raw_data = {str(n): raw.read_raw(str(n)) for n in (1024,8192)}
    raw.verify_native_sources(raw_data)
    sem = {}
    for row in raw.read(HERE.parents[3] / "galectin_affinity/table1.csv"):
        for t in ("Gal1","Gal3"):
            sem[t, int(row["compound"])] = D(1000) * D(row[t+"_error_uM"])
    with (dest / "summary.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(reports) == len(rows) == 12
    for r, row in zip(reports, rows, strict=True):
        for key, value in row.items():
            assert str(r[key]) == value or (key in ("lower","upper") and D(value) == r[key])
        verify(r, raw_data[str(r["nodes"])], sem)
    print(f'{len(reports)} primal witnesses and {sum(len(r["boxes"]) for r in reports)} covering boxes verified.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    main(parser.parse_args().results)
