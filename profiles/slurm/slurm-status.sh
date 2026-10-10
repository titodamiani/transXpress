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

# One squeue call lists every job of this user, so the answer for most jobs comes
# from this cache instead of a Slurm query. squeue reads the controller's memory;
# sacct reads the Slurm database and is the expensive one.
CACHE=".snakemake/slurm-status-cache"
CACHE_MAX_AGE=60   # seconds a cache file is trusted
CACHE_MIN_AGE=5    # a miss rebuilds the cache only if it is older than this

# First line of the cache is the time it was built, so no stat call is needed.
NOW=$(date +%s)

cache_age() {
  if [ ! -s "$CACHE" ]; then
    echo 999999
    return
  fi
  BUILT=$(head -1 "$CACHE" 2>/dev/null)
  case "$BUILT" in
    ''|*[!0-9]*) echo 999999 ;;
    *) echo $((NOW - BUILT)) ;;
  esac
}

build_cache() {
  # -u "$USER" rather than --me, which needs a newer Slurm.
  QUEUE=$(squeue -h -u "$USER" -o %i 2>/dev/null)
  if [ $? -ne 0 ]; then
    return 1
  fi
  mkdir -p "$(dirname "$CACHE")" 2>/dev/null
  # Write then move, so a reader never sees half a file.
  { echo "$NOW"; echo "$QUEUE"; } > "$CACHE.$$" 2>/dev/null && mv -f "$CACHE.$$" "$CACHE" 2>/dev/null
  return 0
}

in_cache() {
  # An array task appears as "123_4", so match the plain id and that form.
  tail -n +2 "$CACHE" 2>/dev/null | grep -qx -e "$JOBID" -e "${JOBID}_[0-9]*"
}

if [ "$(cache_age)" -gt "$CACHE_MAX_AGE" ]; then
  build_cache
fi

if in_cache; then
  echo running
  exit 0
fi

# A job submitted since the cache was built is not in it yet. Rebuild once and
# look again, so a wave of new jobs costs one squeue call and not one sacct call
# each.
if [ "$(cache_age)" -gt "$CACHE_MIN_AGE" ]; then
  build_cache
  if in_cache; then
    echo running
    exit 0
  fi
fi

# The job has left the queue, so ask the database for its final state. -X gives
# the job itself and not its steps.
STATE=$(sacct -X -n -P -j "$JOBID" -o State 2>/dev/null | head -1)

# "CANCELLED by 12345" and "CANCELLED+" both mean cancelled.
STATE="${STATE%% *}"
STATE="${STATE%+}"

if [ -z "$STATE" ]; then
  # Slurm did not answer. Treat it as "wait", because a cluster outage is not a
  # pipeline failure, and the next poll asks again.
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
