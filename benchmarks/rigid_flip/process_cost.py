"""Measure one Linux docking process without inheriting the analysis interpreter's high RSS."""
import json
import os
from pathlib import Path
import sys
import time

if __name__ == "__main__":
    log, stats, *command = sys.argv[1:]
    with open(log, "w") as stream:
        start = time.perf_counter()
        pid = os.posix_spawn(command[0], command, os.environ, file_actions=[
            (os.POSIX_SPAWN_DUP2, stream.fileno(), 1), (os.POSIX_SPAWN_DUP2, stream.fileno(), 2)
        ])
        _, status, usage = os.wait4(pid, 0)
        seconds = time.perf_counter() - start
    if os.waitstatus_to_exitcode(status):
        sys.exit(os.waitstatus_to_exitcode(status))
    Path(stats).write_text(json.dumps(dict(wall_s=seconds, cpu_s=usage.ru_utime+usage.ru_stime,
                                           peak_rss_KiB=usage.ru_maxrss)) + "\n")
