# Node-local work folder for the heavy rules; a rule sources this file in its shell block.
#
#   scratch_start VAR LOG stop|allow [FAILURE_FUNCTION]
#
# VAR names the environment variable of the scratch key ("" works in the run folder). With a folder, this makes
# $WORK in $LOCAL_ROOT, writes its name to LOG and deletes it when the shell ends ($LOCAL_ROOT is empty without one).
# stop: exit when the folder is in RAM (tmpfs), because the work files would count against the job's memory.
# FAILURE_FUNCTION runs before the delete after an error, with the exit code in $rc. A job that Slurm kills
# at the time limit or with scancel runs nothing.

scratch_start() {
  local var=$1 log=$2 tmpfs=$3
  SCRATCH_FAILURE=${4:-}
  LOCAL_ROOT=""
  if [ -n "$var" ]; then LOCAL_ROOT="${!var:-}"; fi
  if [ -z "$LOCAL_ROOT" ]; then return 0; fi
  if [ "$tmpfs" = "stop" ] && [ "$(stat -f -c %T "$LOCAL_ROOT" 2>/dev/null)" = "tmpfs" ]; then
    echo "ERROR: the scratch folder $LOCAL_ROOT is in RAM (tmpfs) on $(hostname). Run this rule on a node with a local disk." | tee -a "$log" >&2
    exit 1
  fi
  WORK=$(mktemp -d "$LOCAL_ROOT/transxpress.XXXXXX")
  echo "work folder: $WORK" >> "$log" 2>&1
  trap scratch_on_exit EXIT
  trap "exit 143" TERM
}

scratch_on_exit() {
  rc=$?
  set +e
  trap - EXIT TERM
  if [ $rc -ne 0 ] && [ -n "$SCRATCH_FAILURE" ]; then "$SCRATCH_FAILURE"; fi
  rm -rf "$WORK"
  exit $rc
}


# The work folder of align_reads and kallisto (align_and_estimate_abundance.pl).
#
#   abundance_start RULE VAR LOG OUT_DIR TRANSCRIPTOME GENE_MAP SAMPLES   ...   abundance_finish
#
# The tool writes the indexes and one folder per replicate into the folder it runs in: OUT_DIR, or a node-local
# folder. That folder holds links to the transcriptome and the gene map and a samples file with full read paths.
# Paths are relative to the run folder. abundance_start sets RUN, LOG, OUT and WORK and goes to WORK.
# abundance_finish removes the links and the samples file and, with a node-local folder, copies the rest to
# OUT_DIR. After an error the files go to OUT_DIR too.

abundance_start() {
  ABUNDANCE_RULE=$1
  local var=$2
  RUN=$(pwd)
  LOG="$RUN/$3"
  ABUNDANCE_OUT_DIR=$4
  : > "$LOG"
  OUT="$RUN/$4"
  mkdir -p "$OUT"
  scratch_start "$var" "$LOG" stop abundance_on_failure
  if [ -z "$LOCAL_ROOT" ]; then WORK="$OUT"; fi
  cd "$WORK"
  ln -sfn "$RUN/$5" transcriptome.fasta
  ln -sfn "$RUN/$6" transcriptome.gene_trans_map
  awk -v run="$RUN" 'BEGIN { OFS = "\t" } { for (i = 3; i <= NF; i++) $i = run "/" $i; print }' "$RUN/$7" > samples_trimmed.txt
}

abundance_copy_back() {
  ( cd "$WORK" && find . -mindepth 1 -maxdepth 1 ! -type l ! -name samples_trimmed.txt -exec cp -rp -t "$OUT" {} + )
}

abundance_on_failure() {
  cd "$RUN"
  echo "$ABUNDANCE_RULE stopped with code $rc, so the files in $WORK go to $ABUNDANCE_OUT_DIR" &>> "$LOG"
  abundance_copy_back &>> "$LOG"
}

abundance_finish() {
  rm -f transcriptome.fasta transcriptome.gene_trans_map samples_trimmed.txt
  if [ -n "$LOCAL_ROOT" ]; then
    local size start
    size=$(du -sh "$WORK" | cut -f1)
    start=$(date +%s)
    abundance_copy_back &>> "$LOG"
    echo "copied $size from $WORK to $ABUNDANCE_OUT_DIR in $(( $(date +%s) - start )) s" &>> "$LOG"
  fi
}
