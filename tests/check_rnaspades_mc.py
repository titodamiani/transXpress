#!/usr/bin/env python3
"""Check run_rnaspades of workflow/scripts/tools.sh with a fake rnaspades.py.

Usage: check_rnaspades_mc.py   (run from the clone)
The fake rnaspades.py acts as rnaSPAdes 4.3.0: it records its arguments and writes K21, K33 and K127 with
configs/config.info (correct_mismatches true only in the last K), the work files .bin_reads, configs and tmp/spades_x,
a spades.log, and a run_spades.sh that names the K configs, writes transcripts.fasta and prints the spades-core log to
stdout (not into spades.log). The check calls run_rnaspades for each value of rnaspades_mismatch_correction and for
the cases that must stop. It also checks that the rnaspades rule passes the key.
Exit code 0: all as expected. 1: a case differs.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

TOOLS = Path("workflow/scripts/tools.sh").resolve()
K_DIRS = ["K21", "K33", "K127"]
LAST_K = "K127"
ARGS = ["--dataset", "samples.yaml", "-t", "4", "-m", "10", "-o", "out", "--ss", "rf", "--only-assembler"]

# FAKE_NO_KEY=1: the last K config has no correct_mismatches line. FAKE_MC_STAGE=1: run_spades.sh prints the stage.
# FAKE_STALE=1: a K201 folder of an earlier run, which run_spades.sh does not name.
FAKE = r"""#!/bin/bash
set -euo pipefail
echo "$*" >> "$FAKE_CALLS"
out=""
prev=""
for a in "$@"; do
  [ "$prev" = "-o" ] && out=$a
  prev=$a
done
mkdir -p "$out/.bin_reads" "$out/configs" "$out/tmp/spades_x"
touch "$out/.bin_reads/lib_0" "$out/configs/config.info" "$out/tmp/spades_x/file"
echo "spades.py log" > "$out/spades.log"
if [ "${FAKE_STALE:-0}" = 1 ]; then
  mkdir -p "$out/K201/configs"
  printf 'correct_mismatches true\n' > "$out/K201/configs/config.info"
fi
script='#!/bin/bash\nset -e\n'
for k in K21 K33 K127; do
  mkdir -p "$out/$k/configs"
  value=false
  [ $k = K127 ] && value=true
  if [ $k = K127 ] && [ "${FAKE_NO_KEY:-0}" = 1 ]; then
    printf 'output_dir ./\n' > "$out/$k/configs/config.info"
  else
    printf 'output_dir ./\ncorrect_mismatches %s\nmax_threads 4\n' $value > "$out/$k/configs/config.info"
  fi
  script="$script""echo spades-core $out/$k/configs/config.info $out/$k/configs/rna_mode.info\n"
  script="$script""echo '== STAGE == Assembling ($k)'\n"
done
if [ "${FAKE_MC_STAGE:-0}" = 1 ]; then
  script="$script""echo '== STAGE == Mismatch Correction (K127)'\n"
fi
script="$script""echo run_spades >> $FAKE_CALLS\necho '>t1' > $out/transcripts.fasta\n"
printf "$script" > "$out/run_spades.sh"
case " $* " in
  *" --only-generate-config "*) ;;
  *) bash "$out/run_spades.sh" >> "$out/spades.log" ;;
esac
"""


def run(value, work, env_extra=None):
    """Run run_rnaspades in the folder work. Returns (exit code, calls, configs by K, output text)."""
    work = Path(work)
    bindir = work / "bin"
    bindir.mkdir()
    fake = bindir / "rnaspades.py"
    fake.write_text(FAKE)
    fake.chmod(0o755)
    env = dict(os.environ, PATH="%s:%s" % (bindir, os.environ["PATH"]), FAKE_CALLS=str(work / "calls.txt"))
    env.update(env_extra or {})
    script = 'set -euo pipefail; source "%s"; run_rnaspades %s out %s' % (TOOLS, value, " ".join(ARGS))
    proc = subprocess.run(["bash", "-c", script], cwd=work, env=env, capture_output=True, text=True)
    calls_file = work / "calls.txt"
    calls = calls_file.read_text().splitlines() if calls_file.exists() else []
    configs = {}
    for k in K_DIRS + ["K201"]:
        path = work / "out" / k / "configs" / "config.info"
        configs[k] = path.read_text() if path.exists() else None
    return proc.returncode, calls, configs, proc.stdout + proc.stderr


def main():
    status = 0

    def fail(text):
        nonlocal status
        print("FAILED: " + text)
        status = 1

    # true: one call with the arguments as given, no config edit
    with tempfile.TemporaryDirectory() as work:
        code, calls, configs, out = run("true", work)
        if code != 0:
            fail("true: exit code %d\n%s" % (code, out))
        if calls != [" ".join(ARGS), "run_spades"]:
            fail("true: calls %s" % calls)
        if "correct_mismatches true" not in (configs[LAST_K] or ""):
            fail("true: the last K config was edited")

    # false: config only, edit of the last K (sorted as numbers), run_spades.sh, its output in spades.log, cleanup
    with tempfile.TemporaryDirectory() as work:
        code, calls, configs, out = run("false", work)
        if code != 0:
            fail("false: exit code %d\n%s" % (code, out))
        if calls != [" ".join(ARGS) + " --only-generate-config", "run_spades"]:
            fail("false: calls %s" % calls)
        text = configs[LAST_K] or ""
        if "correct_mismatches false" not in text or "correct_mismatches true" in text:
            fail("false: %s not edited:\n%s" % (LAST_K, text))
        if "max_threads 4" not in text or "output_dir ./" not in text:
            fail("false: %s lost other lines:\n%s" % (LAST_K, text))
        log = (Path(work) / "out" / "spades.log").read_text()
        if "spades.py log" not in log or "STAGE == Assembling (K127)" not in log:
            fail("false: spades.log lacks the spades.py or the run_spades.sh part:\n%s" % log)
        for left in [".bin_reads", "configs", "tmp/spades_x"]:
            if (Path(work) / "out" / left).exists():
                fail("false: out/%s is still there" % left)
        if not (Path(work) / "out" / "tmp").is_dir():
            fail("false: out/tmp was deleted (a single rnaspades.py call keeps it)")

    # false: a K folder of an earlier run that run_spades.sh does not name is not edited
    with tempfile.TemporaryDirectory() as work:
        code, _, configs, out = run("false", work, {"FAKE_STALE": "1"})
        if code != 0 or "correct_mismatches false" not in (configs[LAST_K] or ""):
            fail("false with a stale K201: exit code %d, %s:\n%s" % (code, LAST_K, configs[LAST_K]))
        if "correct_mismatches true" not in (configs["K201"] or ""):
            fail("false with a stale K201: K201 was edited")

    # false must stop when the key is missing, before run_spades.sh, with its own message
    with tempfile.TemporaryDirectory() as work:
        code, calls, _, out = run("false", work, {"FAKE_NO_KEY": "1"})
        if code == 0 or "run_spades" in calls or "no line 'correct_mismatches true'" not in out:
            fail("false without the key: exit code %d, calls %s\n%s" % (code, calls, out))

    # false must stop when the run_spades.sh output has the stage
    with tempfile.TemporaryDirectory() as work:
        code, _, _, out = run("false", work, {"FAKE_MC_STAGE": "1"})
        if code == 0 or "has a Mismatch Correction stage" not in out:
            fail("false with a Mismatch Correction stage: exit code %d\n%s" % (code, out))

    # another value must stop before rnaspades.py runs
    with tempfile.TemporaryDirectory() as work:
        code, calls, _, _ = run("False", work)
        if code == 0 or calls:
            fail("value 'False': exit code %d, calls %s" % (code, calls))

    rule = Path("workflow/rules/assembly.smk").read_text().split("rule rnaspades:")[1].split("\nrule ")[0]
    if 'run_rnaspades {params.mismatch_correction} "$SPADES_OUT"' not in rule:
        fail("rule rnaspades does not call run_rnaspades with params.mismatch_correction")
    if 'mismatch_correction="true" if config["rnaspades_mismatch_correction"] else "false"' not in rule:
        fail("rule rnaspades does not set params.mismatch_correction from the key")

    print("OK: 7 cases" if status == 0 else "FAILED")
    return status


if __name__ == "__main__":
    sys.exit(main())
