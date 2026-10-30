# Changelog

One entry for each merge into `master` that changes the code. Each entry says if a file under `results/` changes.
The tool versions of a run are in its `results/summary/run_info.txt`; cite those. Before this file, the merges are in `git log --merges`.

## Unreleased
- `results/summary/run_info.txt` has the lines `pfam_md5`, `rfam_md5` and `sprot_md5`: the md5 of the database files the annotation read (merge `c3f100e`). Results change: only `run_info.txt`, by these three lines. Two runs on one database folder give equal lines.
