"""Build three offline path observers using the parent's exact compile and link flags."""

import argparse
import concurrent.futures
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent


def insert(s, anchor, text):
    if s.count(anchor) != 1:
        raise ValueError(anchor)
    return s.replace(anchor, text + anchor)


def instrument(src, stem, arm):
    if stem == "vina":
        s = insert(src, '#include "vina.h"', '#include "observer.h"\n')
        s = insert(
            s,
            "const sz refinement_inputs = poses.size();",
            "path_observer state_trace(m_model);\n\t\t\t\t",
        )
        s = insert(
            s,
            "VINA_FOR(p, 5){",
            'state_trace.write(i, "explicit_in", poses[i], evalcount, "-");\n\t\t\t\t\t',
        )
        return insert(
            s,
            "++refinement_calls;",
            'state_trace.write(i, "explicit_out", poses[i], evalcount, std::to_string(p));\n\t\t\t\t\t\t',
        )
    s = insert(src, '#include "monte_carlo.h"', '#include "observer.h"\n')
    s = insert(s, "VINA_U_FOR(step, global_steps)", "path_observer state_trace(m);\n\t")
    s = insert(
        s,
        "candidate.c = tmp.c;",
        'state_trace.write(step, "source", tmp, evalcount, reference_from_save ? "save" : "hunt");\n\t\t',
    )
    if arm != "parent":
        s = insert(
            s,
            "const bool eligible = SILENT_TORSION_ALWAYS",
            'state_trace.write(step, "native_proposal", candidate, evalcount, mutation_kind(tmp.c,candidate.c));\n        ',
        )
    s = insert(
        s,
        "const int before_hunt = evalcount;",
        'state_trace.write(step, "proposal", candidate, evalcount, mutation_kind(tmp.c,candidate.c));\n\t\t',
    )
    s = insert(
        s,
        "if(step == 0 || metropolis_accept",
        'state_trace.write(step, "hunt", candidate, evalcount, "-");\n\t\t',
    )
    s = insert(
        s,
        "std::swap(tmp.c, candidate.c);",
        'state_trace.write(step, "accepted", candidate, evalcount, "-");\n\t\t\t',
    )
    s = insert(
        s,
        "saved.coords = m.get_heavy_atom_movable_coords();",
        'state_trace.write(step, "second", saved, evalcount, "-");\n\t\t\t\t',
    )
    return insert(
        s,
        "reference_from_save = !separate_hunt_chain;",
        'for (sz slot=0;slot<out.size();++slot) state_trace.write(step, "pool", out[slot], evalcount, std::to_string(slot));\n\t\t\t\t',
    )


def main(a):
    source = a.checkout.resolve()
    for stem in ("monte_carlo", "vina"):
        name = "src/lib/" + stem + ".cpp"
        expected = subprocess.check_output(
            ["git", "-C", source, "show", "3d9dfb9:" + name]
        )
        if (source / name).read_bytes() != expected:
            raise ValueError("use the frozen parent 3d9dfb9 for " + name)
    commands = [
        shlex.split(line)
        for line in a.build_log.read_text().splitlines()
        if line.strip()
    ]
    link = next(
        c for c in commands if "-o" in c and Path(c[c.index("-o") + 1]).name == "vina"
    )
    out = a.output.resolve()
    out.mkdir(exist_ok=True, parents=True)
    with tempfile.TemporaryDirectory() as tmp:
        w = Path(tmp)

        def obj(stem, arm):
            d = w / (stem + arm) / "src/lib"
            d.mkdir(parents=True)
            cpp = d / (stem + ".cpp")
            shutil.copyfile(source / "src/lib" / cpp.name, cpp)
            if stem == "monte_carlo" and arm != "parent":
                subprocess.run(
                    ["patch", "-p1", "-i", str(a.patch.resolve())],
                    cwd=d.parents[1],
                    check=True,
                )
            cpp.write_text(instrument(cpp.read_text(), stem, arm))
            o = d / (stem + ".o")
            c = next(
                c
                for c in commands
                if "-c" in c and Path(c[c.index("-c") + 1]).name == cpp.name
            ).copy()
            c[c.index("-c") + 1] = str(cpp)
            c[c.index("-o") + 1] = str(o)
            c = [
                ("-I" + str(source / "src/lib"))
                if x.startswith("-I") and x.endswith("/src/lib")
                else x
                for x in c
            ]
            c.extend(
                [
                    "-I" + str(HERE),
                    "-I" + str(source / "src/lib"),
                    "-DSILENT_TORSION_ALWAYS=" + str(int(arm == "always")),
                ]
            )
            subprocess.run(c, check=True)
            return o

        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            v = pool.submit(obj, "vina", "common")
            jobs = {
                arm: pool.submit(obj, "monte_carlo", arm)
                for arm in ("parent", "always", "adaptive")
            }
            for arm, job in jobs.items():
                c = link.copy()
                for stem, o in [("vina", v.result()), ("monte_carlo", job.result())]:
                    old = next(x for x in c if Path(x).name == stem + ".o")
                    c[c.index(old)] = str(o)
                c[c.index("-o") + 1] = str(out / arm)
                subprocess.run(c, check=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("checkout", "build_log", "patch", "output"):
        p.add_argument(name, type=Path)
    main(p.parse_args())
