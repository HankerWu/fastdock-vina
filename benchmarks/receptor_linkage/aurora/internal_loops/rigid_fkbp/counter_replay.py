"""Count upstream work on the same AD4 calls and compare every output byte with the release."""

import argparse
from pathlib import Path
import re
import subprocess
import time
from run_panel import cpu_seconds, RECEIVERS, sha, table


def replay(work, binary):
    dest = work / "official_counter"
    dest.mkdir(exist_ok=True)
    rows = []
    for r in RECEIVERS:
        for n in range(4, 12):
            stem = f"{r}_{n}"
            release = work / "pools/official" / stem
            old_output = release.with_suffix(".pdbqt").relative_to(work)
            new_output = Path("official_counter") / old_output.name
            conf = dest / f"{stem}.conf"
            conf.write_text(
                release.with_suffix(".conf")
                .read_text()
                .replace(str(old_output), str(new_output))
            )
            t, c = time.perf_counter(), cpu_seconds()
            logfile = dest / f"{stem}.log"
            with logfile.open("w") as log:
                subprocess.run(
                    [str(binary), "--config", str(conf.relative_to(work))],
                    cwd=work,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            elapsed, cpu = time.perf_counter() - t, cpu_seconds() - c
            text = logfile.read_text()
            counts = {
                k: list(map(int, re.findall(r"\b" + k + r"=(\d+)", text)))
                for k in ("hunt_calls", "hunt_evals", "second_evals", "explicit_evals")
            }
            assert all(
                len(counts[k]) == 3
                for k in ("hunt_calls", "hunt_evals", "second_evals")
            )
            assert not any(
                counts["explicit_evals"]
            )  # AD4 uses maps, no final explicit refinement.
            eq = (work / old_output).read_bytes() == (work / new_output).read_bytes()
            rows.append(
                dict(
                    receiver=r,
                    compound=n,
                    byte_equal=int(eq),
                    wall_s=elapsed,
                    cpu_s=cpu,
                    steps=sum(counts["hunt_calls"]),
                    hunt_evals=sum(counts["hunt_evals"]),
                    second_evals=sum(counts["second_evals"]),
                    sha256=sha(work / new_output),
                    executable_sha256=sha(binary),
                )
            )
    table(work / "counter_replay.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("counter", type=Path)
    a = p.parse_args()
    replay(a.work.resolve(), a.counter.resolve())
