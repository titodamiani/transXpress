#! /bin/bash

# Tells Snakemake whether a Slurm job succeeded, failed or still runs.
# Used by the cluster-status key in config.yaml.
#
# Snakemake calls this script with one argument, the Slurm job id, and expects
# exactly one line on stdout: success, failed or running. Any other output, or a
# non-zero exit, stops the whole workflow. So every path here prints one of the
# three words and exits 0, and every message goes to stderr.
#
# Without this script Snakemake waits for the job to write a .jobfinished or
# .jobfailed marker file. A job killed by the node or by the kernel writes
# neither, so Snakemake waits for it for ever.
#
# Slurm job state names: https://slurm.schedmd.com/sacct.html (JOB STATE CODES)

JOBID="$1"

# sbatch --parsable prints "jobid;clustername" on a federated cluster, and
# Snakemake passes that string through unchanged.
JOBID="${JOBID%%;*}"

if [ -z "$JOBID" ]; then
  echo "slurm-status.sh: called without a job id" >&2
  echo running
  exit 0
fi

# The answers come from two listings of this user's jobs, each built at most once
# a minute: squeue for the jobs in the queue, which reads the controller's memory,
# and sacct for the jobs of the last hour, which reads the Slurm database. So a run
# of any size costs Slurm at most two calls a minute. The listings sit in the
# temporary folder of the node, not on the shared file system.
CACHE_DIR="${TMPDIR:-/tmp}"
QUEUE_CACHE="$CACHE_DIR/transxpress-slurm-status-$USER.squeue"
STATE_CACHE="$CACHE_DIR/transxpress-slurm-status-$USER.sacct"
CACHE_MAX_AGE=60   # seconds a listing is trusted

# First line of a listing is the time it was built, so no stat call is needed.
NOW=$(date +%s)

cache_age() {
  if [ ! -s "$1" ]; then
    echo 999999
    return
  fi
  BUILT=$(head -1 "$1" 2>/dev/null)
  case "$BUILT" in
    ''|*[!0-9]*) echo 999999 ;;
    *) echo $((NOW - BUILT)) ;;
  esac
}

# $1 is the listing file, the rest is the Slurm command. A failed command keeps
# the old listing.
build_cache() {
  CACHE="$1"
  shift
  LISTING=$("$@" 2>/dev/null)
  if [ $? -ne 0 ]; then
    return 1
  fi
  # Write then move, so a reader never sees half a file.
  { echo "$NOW"; echo "$LISTING"; } > "$CACHE.$$" 2>/dev/null && mv -f "$CACHE.$$" "$CACHE" 2>/dev/null
  return 0
}

# -u "$USER" rather than --me, which needs a newer Slurm.
if [ "$(cache_age "$QUEUE_CACHE")" -gt "$CACHE_MAX_AGE" ]; then
  build_cache "$QUEUE_CACHE" squeue -h -u "$USER" -o %i
fi

# An array task appears as "123_4", so match the plain id and that form.
if tail -n +2 "$QUEUE_CACHE" 2>/dev/null | grep -qx -e "$JOBID" -e "${JOBID}_[0-9]*"; then
  echo running
  exit 0
fi

# The job is not in the queue listing: it has left the queue, or it is newer than
# the listing. The sacct listing tells the two apart. -X gives the job itself and
# not its steps.
if [ "$(cache_age "$STATE_CACHE")" -gt "$CACHE_MAX_AGE" ]; then
  build_cache "$STATE_CACHE" sacct -X -n -P -u "$USER" -S now-1hours -o JobID,State
fi
STATE=$(tail -n +2 "$STATE_CACHE" 2>/dev/null | awk -F'|' -v id="$JOBID" '$1 == id { print $2; exit }')

# "CANCELLED by 12345" and "CANCELLED+" both mean cancelled.
STATE="${STATE%% *}"
STATE="${STATE%+}"

if [ -z "$STATE" ]; then
  # The job is in neither listing, so it is newer than both, or Slurm did not
  # answer. Treat it as "wait", because a cluster outage is not a pipeline
  # failure, and a later poll asks again.
  echo running
  exit 0
fi

case "$STATE" in
  COMPLETED)
    echo success
    ;;
  PENDING|RUNNING|REQUEUED|RESIZING|SUSPENDED|COMPLETING|CONFIGURING|SIGNALING|STAGE_OUT|REQUEUE_HOLD|REQUEUE_FED|RESV_DEL_HOLD|STOPPED)
    echo running
    ;;
  BOOT_FAIL|CANCELLED|DEADLINE|FAILED|NODE_FAIL|OUT_OF_MEMORY|PREEMPTED|REVOKED|TIMEOUT|SPECIAL_EXIT)
    echo failed
    ;;
  *)
    # An unknown state counts as failed, so that a state this script has not met
    # stops the job instead of hanging the workflow. restart-times in config.yaml
    # gives the job another try.
    echo "slurm-status.sh: job $JOBID has unknown state '$STATE', treating it as failed" >&2
    echo failed
    ;;
esac

exit 0
