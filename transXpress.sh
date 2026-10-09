#! /bin/bash

echo "Running the transXpress pipeline using snakemake"

CLUSTER="NONE"

if [ ! -z `which sbatch` ]; then
  CLUSTER="SLURM"
fi

if [ ! -z `which bsub` ]; then
  CLUSTER="LSF"
fi

if [ ! -z `which qsub` ]; then
  CLUSTER="PBS"
fi

CONDA_PREFIX_ARGS=()
if [ -n "$TRANSXPRESS_CONDA_PREFIX" ]; then
  CONDA_PREFIX_ARGS=(--conda-prefix "$TRANSXPRESS_CONDA_PREFIX")
else
  echo "Warning: TRANSXPRESS_CONDA_PREFIX is not set. Snakemake builds the conda environments in .snakemake/conda in this folder."
fi

case "$CLUSTER" in
"LSF")
  echo "Submitting snakemake jobs to LSF cluster"
  snakemake "${CONDA_PREFIX_ARGS[@]}" --conda-frontend conda --use-conda --latency-wait 60 --restart-times 1 --jobs 10000 --cluster "bsub -oo {log}.bsub -n {threads} -R rusage[mem={params.memory}000] -R span[hosts=1]" "$@"
  ;;
"SLURM")
  echo "Submitting snakemake jobs to SLURM cluster"
  snakemake --profile "$(dirname "$0")/profiles/slurm" "${CONDA_PREFIX_ARGS[@]}" "$@"
  ;;
"PBS")
  echo "Submitting snakemake jobs to PBS/Torque cluster"
  snakemake "${CONDA_PREFIX_ARGS[@]}" --conda-frontend conda --use-conda --latency-wait 60 --restart-times 1 --jobs 10000 --cluster "qsub -j oe -o {log}.pbs -l select=1:ncpus={threads}:mem={params.memory}gb" "$@"
  ;;
*)
  snakemake "${CONDA_PREFIX_ARGS[@]}" --conda-frontend conda --use-conda --cores all "$@"
  ;;
esac


