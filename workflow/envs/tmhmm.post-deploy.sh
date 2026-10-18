#!/bin/sh
# Snakemake runs this script in the tmhmm env after it builds the env from tmhmm.linux-64.pin.txt, which
# holds the conda packages only. It adds tmhmm.py 1.2.1, the source file from PyPI with its sha256. pip
# builds it against the numpy and Cython of the env, and its only dependency, numpy, is in the pin file.
set -e
REQUIREMENTS=$(mktemp)
echo "tmhmm.py==1.2.1 --hash=sha256:1264b497ecfa295b2e3046f462a345c1e545e81d06e133e37ca0c692d7e7e4c6" > "$REQUIREMENTS"
pip install --no-deps --no-build-isolation --no-cache-dir --require-hashes --requirement "$REQUIREMENTS"
rm "$REQUIREMENTS"
