#!/bin/bash
# Makes a T1 run folder: the test data, a samples file and a config.yaml.
# Usage: tests/make_run.sh <folder> <assembler> <db_dir>
#   folder     the new run folder, directly inside the repository folder. Run the pipeline there:
#              snakemake --snakefile ../workflow/Snakefile
#   assembler  rnaspades or trinity
#   db_dir     a database folder made by the download_databases target (see README.md)
# SignalP and TargetP are off (signalp "" and targetp ""), and Trinity runs Inchworm on one thread,
# so that two runs of the same folder give the same transcripts.
set -euo pipefail
folder=$1
assembler=$2
db_dir=$3
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$folder"
cp -r "$here/test_data" "$here/test_samples.txt" "$folder/"
grep -v -E '^(assembler|signalp|targetp|trinity_parameters|kmers|lineage|db_dir):' "$here/config.yaml" > "$folder/config.yaml"
cat >> "$folder/config.yaml" <<CONFIG

assembler: "$assembler"
trinity_parameters: "--seqType fq --inchworm_cpu 1"
kmers: ""
signalp: ""
targetp: ""
lineage: "viridiplantae_odb12"
db_dir: "$db_dir"
CONFIG
