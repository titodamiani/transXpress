# Changelog

One entry for each merge into `master` that changes the code. Each entry says if a file under `results/` changes.
The tool versions of a run are in its `results/summary/run_info.txt`; cite those. Before this file, the merges are in `git log --merges`.

## Unreleased
- CI job `t1-image` runs on `ubuntu-24.04`. On `ubuntu-latest` (Ubuntu 26.04) the Apptainer 1.5.2 package cannot install, because it needs `libfuse3-3`. Results change: none (CI only).
- `tests/allowed_diffs.tsv` has the line `results/expression/bowtie2_rsem/*` (mode `skip`), so `tests/compare_runs.py` accepts a run that holds the folder of rules `align_reads` and `IGV` and a run that does not. Results change: none (test tool only). The folder is not compared by the script (merge `5b8c8b1`).
- Rule `prepare_samples_for_IGV` asks for 16 GB (`mem_mb=16384`). `samtools sort --threads 4 -m 2G` holds 8 GB of buffers in its 4 threads, and 10 GB was not enough on piper09 (9 libraries; the `sacct` records and the job logs show out-of-memory kills of the sort jobs). Results change: none. This is a resource limit; the shell block is unchanged and no output depends on it (merge `27bdaaa`).
- `results/summary/run_info.txt` has the lines `pfam_md5`, `rfam_md5` and `sprot_md5`: the md5 of the database files the annotation read (merge `c3f100e`). Results change: only `run_info.txt`, by these three lines. Two runs on one database folder give equal lines.
- Rule `kallisto` passes the strand flag (`--rf-stranded` or `--fr-stranded`) to kallisto when `strand_specific` is not empty (merge `025f358`). Before, `align_and_estimate_abundance.pl` (Trinity 2.13.2 in the image; the same lines are in Trinity 2.15.2 and on its master) dropped the flag: it appends `kallisto_add_opts` to itself, not the flag. Results change: for a run with non-empty `strand_specific`, the outputs of rule `kallisto` (`results/expression/kallisto/` and the six `results/summary/expr_*` tables). A run with `strand_specific: ""` is unchanged.
