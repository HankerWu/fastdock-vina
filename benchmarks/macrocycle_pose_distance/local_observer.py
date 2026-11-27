"""Link the production executable with one print of native optimizer work at termination."""
import argparse
from pathlib import Path
import subprocess

if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("source","build","boost","out"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    src=a.source.resolve()/"src/lib"
    a.out.mkdir(exist_ok=True,parents=True)
    anchor="\treturn energies_after_opt;"
    code=(src/"vina.cpp").read_text()
    assert code.count(anchor)==1
    changed=a.out/"vina.cpp"
    changed.write_text(code.replace(anchor,'\tstd::cerr << "NATIVE_LOCAL_EVALUATIONS " << evalcount << "\\n";\n'+anchor))
    flags=["c++","-O3","-DNDEBUG","-std=c++14","-pthread","-fvisibility=hidden",
           '-DVERSION="v1.2.7-27-g3c65c0b"','-I'+str(src),'-I'+str(a.boost/"include")]
    obj=a.out/"vina.o"
    subprocess.run([*flags,"-c",changed,"-o",obj],check=True)
    others=sorted(q for q in (a.build/"obj").glob("*.o") if q.stem not in ("vina","split"))
    lib=str(a.boost.resolve()/"lib")
    subprocess.run([*flags,obj,*others,"-L"+lib,"-Wl,-rpath,"+lib,"-lboost_system","-lboost_thread",
                    "-lboost_serialization","-lboost_filesystem","-lboost_program_options","-o",a.out/"vina"],check=True)
