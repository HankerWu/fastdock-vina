"""Compare the actual AD4 ligand models against SMILES-mapped edges used for FKBP anchors."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import time

BENCH = Path(__file__).resolve().parents[2]
TEMPLATE = BENCH/"receptor_linkage/aurora/internal_loops/rigid_graph_control/read_graph.cpp"


def sha(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def build_probe(objects, source, dest, boost, work):
    objects = sorted(p for p in objects.glob("*.o") if p.stem not in ("main", "split"))
    if not objects:
        raise ValueError("supply an existing native Vina object directory")
    code = TEMPLATE.read_text()
    start = code.index("    model m = ")
    end = code.index("    std::ifstream pdbqt", start)
    code = code[:start] + "    model m = parse_ligand_pdbqt_from_file(argv[1], atom_type::AD);\n" + code[end:]
    code = code.replace("atom_xs_type", "atom_ad_type,partial_charge")
    code = code.replace("<< a.xs", "<< a.ad << ',' << std::setprecision(17) << a.charge")
    cpp = dest.with_suffix(".cpp")
    cpp.write_text(code)
    subprocess.run(["c++","-std=c++14","-O3","-DNDEBUG","-pthread","-I"+str(source),"-I"+str(boost/"include"),
                    cpp,*objects,"-L"+str(boost/"lib"),"-Wl,-rpath,"+str(boost/"lib"),"-lboost_system",
                    "-lboost_thread","-lboost_serialization","-lboost_filesystem","-lboost_program_options","-o",dest],
                   env=os.environ|{"TMPDIR":str(work)},check=True)


def main(a):
    a.work.mkdir(parents=True,exist_ok=True)
    a.out.mkdir(parents=True,exist_ok=True)
    data = {}
    start, cpu = time.perf_counter(), os.times()
    for label, source, objects in (("official",a.up_src,a.up_obj),("pro",a.pro_src,a.pro_obj)):
        program = a.work/label
        build_probe(objects,source,program,a.boost,a.work)
        data[label] = dict(parser_sha256=sha(source/"parse_pdbqt.cpp"),model_sha256=sha(source/"model.cpp"),
                           dumper_sha256=sha(program.with_suffix(".cpp")),binary_sha256=sha(program))
    results, allrows = [], []
    for n in range(4,12):
        path = a.work_inputs/"affinity"/f"{n}.pdbqt"
        heavy = {int(l[6:11]) for l in path.read_text().splitlines()
                 if l.startswith(("ATOM  ","HETATM")) and l.split()[-1] not in ("H","HD","HS")}
        authority = {tuple(sorted(map(int,l.split()))) for l in path.with_suffix(".bonds").read_text().splitlines()}
        models = {}
        for label in ("official","pro"):
            raw = subprocess.check_output([a.work/label,path],text=True)
            models[label] = list(csv.DictReader(io.StringIO(raw)))
            allrows.extend(dict(compound=n,model=label,**r) for r in models[label])
        edges = {tuple(sorted((int(r["atom"]),int(r["neighbor"])))) for r in models["pro"]
                 if int(r["atom"]) in heavy and int(r["neighbor"]) in heavy}
        equal = models["official"] == models["pro"]
        results.append(dict(compound=n,ligand_sha256=sha(path),sidecar_sha256=sha(path.with_suffix(".bonds")),
                            heavy=len(heavy),mapped_edges=len(authority),inferred_edges=len(edges),
                            missing=sorted(authority-edges),extra=sorted(edges-authority),
                            official_pro_parsed_rows_equal=equal))
    with gzip.GzipFile(str(a.out/"ad4_parsed.csv.gz"),"wb",mtime=0) as compressed:
        with io.TextIOWrapper(compressed,encoding="utf-8",newline="") as stream:
            w = csv.DictWriter(stream,allrows[0].keys())
            w.writeheader()
            w.writerows(allrows)
    end = os.times()
    data.update(cases=results,wall_s=time.perf_counter()-start,
                cpu_s=sum(end[:4])-sum(cpu[:4]))
    (a.out/"ad4_graphs.json").write_text(json.dumps(data,indent=2)+"\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("up_src","up_obj","pro_src","pro_obj","boost","work_inputs","work","out"):
        p.add_argument(name,type=lambda s:Path(s).resolve())
    main(p.parse_args())
