#! /bin/bash

SNAKEFILE="../Snakefile"

echo "Running the transXpress-trinity pipeline using snakemake"
snakemake --snakefile $SNAKEFILE --cores 8 "$@"

if command -v dot > /dev/null; then
  echo "Making DAG file describing pipeline execution"
  snakemake --snakefile $SNAKEFILE --dag | dot -Tsvg > dag.svg
else
  echo "Graphviz (dot) is not installed, skipping the DAG file"
fi

