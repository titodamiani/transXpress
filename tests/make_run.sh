#!/bin/bash
# Makes a T1 run folder directly inside the clone: tests/make_run.sh <folder> <assembler> <db_dir>
# SignalP and TargetP are off, and Trinity runs Inchworm on one thread, so two runs give the same transcripts.
# README.md "Tests"; CONTRIBUTING.md "Tests".
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
