#!/usr/bin/env python3
"""Check run_rnaspades of workflow/scripts/tools.sh with a fake rnaspades.py.

Usage: check_rnaspades_mc.py   (run from the clone)
The fake rnaspades.py records its arguments and writes the folders K21, K33 and K127, each with configs/config.info
(correct_mismatches true), and a run_spades.sh that writes spades.log. The check calls run_rnaspades for each value of
rnaspades_mismatch_correction and for the cases that must stop. It also checks that the rnaspades rule passes the key.
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

# FAKE_NO_KEY=1: the configs have no correct_mismatches line. FAKE_MC_STAGE=1: run_spades.sh logs the stage.
FAKE = r"""#!/bin/bash
set -euo pipefail
echo "$*" >> "$FAKE_CALLS"
out=""
prev=""
for a in "$@"; do
  [ "$prev" = "-o" ] && out=$a
  prev=$a
done
for k in K21 K33 K127; do
  mkdir -p "$out/$k/configs"
  if [ "${FAKE_NO_KEY:-0}" = 1 ]; then
    printf 'output_dir ./\n' > "$out/$k/configs/config.info"
  else
    printf 'output_dir ./\ncorrect_mismatches true\nmax_threads 4\n' > "$out/$k/configs/config.info"
  fi
done
stage='echo "== STAGE == Assembling" >> "$(dirname "$0")/spades.log"'
if [ "${FAKE_MC_STAGE:-0}" = 1 ]; then
  stage="$stage; echo '== STAGE == Mismatch Correction (K127)' >> \"\$(dirname \"\$0\")/spades.log\""
fi
printf '#!/bin/bash\necho run_spades >> "%s"\n%s\n' "$FAKE_CALLS" "$stage" > "$out/run_spades.sh"
case " $* " in
  *" --only-generate-config "*) ;;
  *) bash "$out/run_spades.sh" ;;
esac
"""


def run(value, env_extra=None):
    """Run run_rnaspades in a new folder. Returns (exit code, calls, configs by K, folder)."""
    work = Path(tempfile.mkdtemp())
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
    for k in K_DIRS:
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
    code, calls, configs, out = run("true")
    if code != 0:
        fail("true: exit code %d\n%s" % (code, out))
    if calls != [" ".join(ARGS), "run_spades"]:
        fail("true: calls %s" % calls)
    if any("correct_mismatches false" in (c or "") for c in configs.values()):
        fail("true: a config was edited")

    # false: config only, edit of the last K only (sorted as numbers), then run_spades.sh
    code, calls, configs, out = run("false")
    if code != 0:
        fail("false: exit code %d\n%s" % (code, out))
    if calls != [" ".join(ARGS) + " --only-generate-config", "run_spades"]:
        fail("false: calls %s" % calls)
    for k in K_DIRS:
        text = configs[k] or ""
        if k == LAST_K:
            if "correct_mismatches false" not in text or "correct_mismatches true" in text:
                fail("false: %s not edited:\n%s" % (k, text))
            if "max_threads 4" not in text or "output_dir ./" not in text:
                fail("false: %s lost other lines:\n%s" % (k, text))
        elif "correct_mismatches true" not in text:
            fail("false: %s was edited:\n%s" % (k, text))

    # false must stop when the key is missing, before run_spades.sh
    code, calls, _, _ = run("false", {"FAKE_NO_KEY": "1"})
    if code == 0 or "run_spades" in calls:
        fail("false without the key: exit code %d, calls %s" % (code, calls))

    # false must stop when the log still has the stage
    code, _, _, _ = run("false", {"FAKE_MC_STAGE": "1"})
    if code == 0:
        fail("false with a Mismatch Correction stage in the log: exit code 0")

    # another value must stop before rnaspades.py runs
    code, calls, _, _ = run("False")
    if code == 0 or calls:
        fail("value 'False': exit code %d, calls %s" % (code, calls))

    rule = Path("workflow/rules/assembly.smk").read_text().split("rule rnaspades:")[1].split("\nrule ")[0]
    if 'run_rnaspades {params.mismatch_correction} "$SPADES_OUT"' not in rule:
        fail("rule rnaspades does not call run_rnaspades with params.mismatch_correction")
    if 'mismatch_correction="true" if config["rnaspades_mismatch_correction"] else "false"' not in rule:
        fail("rule rnaspades does not set params.mismatch_correction from the key")

    print("OK: 6 cases" if status == 0 else "FAILED")
    return status


if __name__ == "__main__":
    sys.exit(main())
