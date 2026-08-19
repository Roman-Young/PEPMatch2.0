"""Shared subprocess runner for the external aligner wrappers.

Every tool used to be invoked with os.system(), which does not raise and whose exit
code was discarded. A tool that died -- missing binary (exit 127), OOM kill (SIGKILL),
bad flag -- produced no output file, and the wrapper then failed later on a confusing
FileNotFoundError from pandas. On an unattended cluster run that is close to
undiagnosable after the fact.

run_tool() checks the exit code, tees stdout/stderr to a per-tool log next to the
results, and raises a RuntimeError naming the tool, exit code and stderr tail.
"""
import os
import subprocess
from pathlib import Path

LOG_DIR = Path(os.environ.get('PEPMATCH_BENCH_LOGDIR', '.'))

# How many external binaries this process has exec'd. The memory measurement needs to
# know: a tool that exec's a standalone binary (BLAST/DIAMOND/MMseqs2) has its real
# footprint in a child that carries NO Python interpreter, whereas an in-process or
# fork-parallel method (PEPMatch, brute force) shares this process's interpreter
# baseline. Those two cases must be reduced to a comparable number differently.
EXTERNAL_RUNS = 0


def run_tool(command, label):
  """Run `command` (a shell string), failing loudly on a non-zero exit.

  Returns the CompletedProcess. Raises RuntimeError if the tool failed, with enough
  context to diagnose it from the results table alone.
  """
  global EXTERNAL_RUNS
  EXTERNAL_RUNS += 1

  LOG_DIR.mkdir(parents=True, exist_ok=True)
  log_path = LOG_DIR / f'{label}.log'

  proc = subprocess.run(command, shell=True, capture_output=True, text=True)

  with open(log_path, 'a') as log:
    log.write(f'$ {command}\n')
    log.write(f'-- exit {proc.returncode}\n')
    if proc.stdout:
      log.write(proc.stdout)
    if proc.stderr:
      log.write('--- stderr ---\n')
      log.write(proc.stderr)
    log.write('\n')

  if proc.returncode != 0:
    stderr_tail = (proc.stderr or '').strip().splitlines()[-5:]
    hint = ''
    if proc.returncode == 127:
      hint = ' (binary not found on PATH)'
    elif proc.returncode < 0 or proc.returncode > 128:
      hint = ' (killed by signal -- often the OOM killer)'
    raise RuntimeError(
      f'{label} exited {proc.returncode}{hint}. See {log_path}. '
      f'stderr tail: {" | ".join(stderr_tail) if stderr_tail else "(empty)"}'
    )

  return proc
