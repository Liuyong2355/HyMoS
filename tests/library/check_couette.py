"""Run library and legacy CLI checks without modifying baseline files."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

memory, cli, baseline, input_file = map(Path, sys.argv[1:])
source_run = next((baseline / "runs").glob("Couette_*"))
with tempfile.TemporaryDirectory(prefix="hymos-couette-") as temp:
    work = Path(temp)
    (work / ".control.1").write_text("999 test sentinel\n")
    subprocess.run([str(memory), str(source_run / "Sol1th.dat"),
                    str(source_run / "res_step"), str(source_run / "DisSol1th.dat"),
                    str(input_file)], cwd=work, check=True, timeout=120)
    runs = work / "runs"
    runs.mkdir()
    env = dict(os.environ, HYMOS_RUNS_ROOT=str(runs))
    proc = subprocess.run([str(cli), "run", "couette", str(baseline / "input.txt")],
                          cwd=work, env=env, capture_output=True, text=True, timeout=120)
    if proc.returncode:
        raise RuntimeError(proc.stdout + proc.stderr)
    actual = next(runs.glob("Couette_*"))
    for name in ("Sol1th.dat", "DisSol1th.dat", "DATA1th.dat"):
        if (actual / name).read_bytes() != (source_run / name).read_bytes():
            raise AssertionError("CLI baseline mismatch: " + name)
    def residual(path):
        return [line.split()[:2] for line in path.read_text().splitlines()
                if line and line[0].isdigit()]
    if residual(actual / "res_step") != residual(source_run / "res_step"):
        raise AssertionError("CLI residual mismatch")
    print("PASS CLI byte-identical field/distribution/restart and residual history")
