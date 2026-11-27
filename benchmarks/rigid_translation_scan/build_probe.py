"""Build an observer with the baseline optimizer; link to tools/build.py objects."""

import argparse
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("objects", type=Path, help="obj directory built using tools/build.py")
p.add_argument("output", type=Path)
p.add_argument("--boost-prefix", type=Path)
p.add_argument(
    "--source-root", type=Path, help="same source checkout as the supplied object files"
)
a = p.parse_args()
here = Path(__file__).resolve().parent
lib = (a.source_root or here.parents[1]) / "src/lib"
out = a.output.resolve()
out.parent.mkdir(parents=True, exist_ok=True)
flags = ["c++", "-O3", "-DNDEBUG", "-std=c++14", "-pthread", "-I" + str(lib)]
link = []
if a.boost_prefix:
    flags += ["-I" + str(a.boost_prefix / "include")]
    link += [
        "-L" + str(a.boost_prefix / "lib"),
        "-Wl,--disable-new-dtags,-rpath," + str(a.boost_prefix / "lib"),
    ]
src = (lib / "vina.cpp").read_text()
needle = "energies_after_opt = score();"
assert src.count(needle) == 1
observer = out.parent / "translation_vina.cpp"
src = src.replace(
    needle, 'std::cout << "gradient_evals " << evalcount << "\\n";\n\t' + needle
)
needle = "lig_intra = m_model.evali(m_precalculated_byatom, authentic_v);"
idx = src.index(needle, src.index("std::vector<double> Vina::score(double"))
pos = idx + len(needle)
# Only fixed-coordinate diagnostics request the uncapped energy; optimization is untouched.
src = (
    src[:pos]
    + """
        if (std::getenv("FD_SCAN_SHOW_UNCAPPED")) {
            std::cout << "uncapped " << receptor.eval(m_model,max_fl) << " "
                      << m_model.evali(m_precalculated_byatom,max_vec) << std::endl;
        }"""
    + src[pos:]
)
observer.write_text(src)
obj = out.parent / "translation_vina.o"
subprocess.run(flags + ["-c", str(observer), "-o", str(obj)], check=True)
others = [
    str(o)
    for o in sorted(a.objects.glob("*.o"))
    if o.name not in ("main.o", "split.o", "vina.o")
]
subprocess.run(
    flags
    + [str(here / "probe.cpp"), str(obj)]
    + others
    + link
    + [
        "-lboost_thread",
        "-lboost_serialization",
        "-lboost_filesystem",
        "-lboost_program_options",
        "-o",
        str(out),
    ],
    check=True,
)
