# Node-local work folder for the heavy rules. A rule sources this file in its shell block.
#
#   scratch_start VAR LOG stop|allow [FAILURE_FUNCTION]
#
# VAR is the name of the environment variable that holds the node-local folder, from the scratch key
# of config.yaml. It is empty when the key is "", and then the rule works in the run folder.
# With a folder, this makes the work folder $WORK inside $LOCAL_ROOT, writes its name to LOG, and
# deletes it when the shell ends. $LOCAL_ROOT is empty without the node-local folder.
# stop: the rule stops with a message when the folder is in RAM (tmpfs), because the work files
#       count against the job's memory. allow: no check.
# FAILURE_FUNCTION runs before the folder is deleted when the shell ends with an error, with the exit
# code in $rc. A job that Slurm kills at its time limit or with scancel runs nothing.

scratch_start() {
  local var=$1 log=$2 tmpfs=$3
  SCRATCH_FAILURE=${4:-}
  LOCAL_ROOT=""
  if [ -n "$var" ]; then LOCAL_ROOT="${!var:-}"; fi
  if [ -z "$LOCAL_ROOT" ]; then return 0; fi
  # A scratch folder in RAM (tmpfs) counts against the memory of the job
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
