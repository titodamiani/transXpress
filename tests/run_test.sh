#! /bin/bash

PROFILE="../profiles/local"

echo "Running the transXpress-trinity pipeline using snakemake"
snakemake --profile $PROFILE --cores 8 "$@"

if command -v dot > /dev/null; then
  echo "Making DAG file describing pipeline execution"
  snakemake --profile $PROFILE --dag | dot -Tsvg > dag.svg
else
  echo "Graphviz (dot) is not installed, skipping the DAG file"
fi

