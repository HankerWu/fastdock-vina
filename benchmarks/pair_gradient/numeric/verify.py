"""Check archived coordinates/energies, torque decomposition, and terminal reference."""
import csv
import gzip
import math
from pathlib import Path
import sys


def require_close(x, y, tol):
    if not math.isclose(x, y, rel_tol=0, abs_tol=tol):
        raise ValueError(f"{x} differs from {y} by more than {tol}")


def verify(directory):
    with gzip.open(directory / "model.csv.gz", "rt") as f:
        rows = list(csv.DictReader(f))
    arms = {arm: {} for arm in ("native", "linear", "hermite", "shifted")}
    for row in rows:
        key = (row["ligand"], row["receptor"], row["stage"], row["dof"], row["h"])
        values = tuple(float(row[k]) for k in ("energy", "central_difference", "returned"))
        if key in arms[row["arm"]] or not all(map(math.isfinite, values)):
            raise ValueError("Duplicate or non-finite full-model result")
        arms[row["arm"]][key] = values
    if any(set(v) != set(arms["native"]) or len(v) != 1128 for v in arms.values()):
        raise ValueError("Incomplete or mismatched model probes")
    for key, native in arms["native"].items():
        linear = arms["linear"][key]
        if native[:2] != linear[:2] or (key[2] == "hunt" and arms["shifted"][key] != linear):
            raise ValueError("Control energy or search evaluation changed")
        if float(key[-1]) == 1e-7:
            for name, tol in (("linear", 5e-4), ("hermite", 1e-5), ("shifted", 5e-4)):
                _, fd, returned = arms[name][key]
                require_close(fd, returned, tol)
    key = next(k for k in arms["native"] if k[:4] == ("2YWP","2YWP","explicit","4") and float(k[4]) == 1e-7)
    with gzip.open(directory / "2ywp_torque.csv.gz", "rt") as f:
        trace = list(csv.DictReader(f))
    for arm in ("native", "linear"):
        summed = math.fsum(float(row[arm+"_torque"]) for row in trace)
        require_close(summed, arms[arm][key][2], 1e-10)
    largest = max(trace, key=lambda r: abs(float(r["native_torque"])-float(r["linear_torque"])))
    if (largest["xs1"], largest["xs2"]) != ("1", "1"):
        raise ValueError("Expected a C_P–C_P pair")
    xyz = lambda k: tuple(float(largest[k+"_"+axis]) for axis in "xyz")
    r = math.dist(xyz("ligand"), xyz("rec"))
    require_close(r, 7.999335847431335, 1e-12)
    # C_P–C_P: Vina's interior repulsion has already vanished at this r.
    # Both Gaussian terms survive. Neither polar carbon is hydrophobic.
    def gaussian_tail(s):
        d = s-3.8
        return -0.035579*math.exp(-(d/0.5)**2) - 0.005156*math.exp(-((d-3)/2)**2)
    h = 1e-6
    derivative = (gaussian_tail(r+h)-gaussian_tail(r-h))/(2*h)
    with (directory / "shift.csv").open() as f:
        entry = next(row for row in csv.DictReader(f) if (row["xs1"],row["xs2"]) == ("1","1"))
    require_close(derivative, float(entry["direct_derivative"]), 2e-12)
    require_close(derivative, float(entry["radial_derivative"]), 1e-5)


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name("results"))
