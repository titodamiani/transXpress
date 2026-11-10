# Shell functions that several rules share. A rule sources this file in its shell block.

# The computer a job ran on: host name, processor and CPU count.
record_host() {
  hostname
  grep -m1 "model name" /proc/cpuinfo 2>/dev/null || uname -m
  nproc 2>/dev/null || true
}

# The folder of the Trinity install. The Trinity command in the env is a link into it.
trinity_home() {
  python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))'
}

# publish_part FILE LOG: FILE.part gets the name FILE when its last line is "# [ok]", the end of a finished
# hmmscan or cmscan. Otherwise the job fails.
publish_part() {
  tail -n 1 "$1.part" | grep -qx '# \[ok\]' || { echo "$1.part does not end with # [ok]" >> "$2"; exit 1; }
  mv "$1.part" "$1"
}

# run_rnaspades <true|false> <out dir> <rnaspades.py arguments, with -o <out dir>>
# true: one rnaspades.py call, with the mismatch correction of the last K.
# false: rnaspades.py writes its configs and run_spades.sh only (--only-generate-config); then correct_mismatches is
# set to false in the config of the last K and run_spades.sh runs. rnaSPAdes has no command-line option for this.
# The step makes most of its changes against the reads and depends on the thread count (README, "rnaSPAdes").
run_rnaspades() {
  local mc=$1 out=$2
  shift 2
  case "$mc" in
    true)
      rnaspades.py "$@"
      return
      ;;
    false) ;;
    *)
      echo "run_rnaspades: the first argument must be true or false, not '$mc'" >&2
      return 1
      ;;
  esac
  rnaspades.py "$@" --only-generate-config || return 1
  # the last K config that run_spades.sh reads; a K folder of an earlier run in the same folder does not count
  local configs config c
  configs=$(grep -oE 'K[0-9]+/configs/config\.info' "$out/run_spades.sh" | sort -u -t K -k 2 -n)
  config=$out/$(echo "$configs" | tail -n 1)
  if ! grep -qE '^correct_mismatches[[:space:]]+true' "$config"; then
    echo "run_rnaspades: no line 'correct_mismatches true' in $config" >&2
    return 1
  fi
  sed -i.bak -E 's/^(correct_mismatches[[:space:]]+)true/\1false/' "$config"
  rm -f "$config.bak"
  for c in $configs; do
    if grep -qE '^correct_mismatches[[:space:]]+true' "$out/$c"; then
      echo "run_rnaspades: $out/$c still has correct_mismatches true" >&2
      return 1
    fi
  done
  echo "run_rnaspades: correct_mismatches false in $config"
  # spades-core writes its log to stdout only; spades.py would copy it into spades.log
  bash "$out/run_spades.sh" 2>&1 | tee -a "$out/spades.log"
  [ "${PIPESTATUS[0]}" = 0 ] || return 1
  [ -s "$out/transcripts.fasta" ] || { echo "run_rnaspades: no $out/transcripts.fasta" >&2; return 1; }
  if grep -q 'STAGE == Mismatch Correction' "$out/spades.log"; then
    echo "run_rnaspades: $out/spades.log has a Mismatch Correction stage" >&2
    return 1
  fi
  # spades.py deletes these at the end of its stages; run_spades.sh does not
  rm -rf "$out/.bin_reads" "$out/configs" "$out"/tmp/*
}
