"""Build an observation-only native MC object against a completed production build."""
import argparse
from pathlib import Path
import subprocess

pre = '''			std::cerr << std::setprecision(17) << "POSE_EVENT " << step << ' ' << saved.e;
			for (sz k = 0; k < m.num_movable_atoms(); ++k) {
				const auto& a = m.movable_atom(k);
				if (a.el == EL_TYPE_H) continue;
				const vec& xyz = m.movable_coords(k);
				std::cerr << ' ' << int(ad_is_closure_dummy(a.ad));
				for (sz d = 0; d < 3; ++d) std::cerr << ' ' << xyz[d];
			}
			std::cerr << '\\n';
'''
post = '''			std::cerr << "RESIDENTS";
			for (const auto& resident : out) std::cerr << ' ' << resident.e;
			std::cerr << '\\n';
'''

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("source", "build", "boost", "output"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    source = a.source.resolve() / "src/lib"
    anchor = 'add_to_output_container(out, saved, min_rmsd, num_saved_mins); // 20 - max size'
    code = (source / "monte_carlo.cpp").read_text()
    assert code.count(anchor) == 1
    a.output.mkdir(parents=True,exist_ok=True)
    mc = a.output / "monte_carlo.cpp"
    mc.write_text('#include <iomanip>\n' + code.replace(anchor, pre + anchor + '\n' + post))
    flags = ["c++", "-O3", "-DNDEBUG", "-std=c++14", "-pthread", "-fvisibility=hidden",
             '-DVERSION="v1.2.7-27-g3c65c0b"', "-I"+str(source), "-I"+str(a.boost / "include")]
    obj = a.output / "observer.o"
    subprocess.run([*flags, "-c", mc, "-o", obj],check=True)
    others = sorted(p for p in (a.build / "obj").glob("*.o") if p.stem not in ("monte_carlo","split"))
    lib = str(a.boost.resolve() / "lib")
    subprocess.run([*flags, obj, *others, "-L"+lib, "-Wl,-rpath,"+lib,
                    "-lboost_system", "-lboost_thread", "-lboost_serialization",
                    "-lboost_filesystem", "-lboost_program_options", "-o", a.output / "vina"],check=True)
