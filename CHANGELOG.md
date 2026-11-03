# Changelog

One entry for each merge into `master` that changes the code. Each entry says if a file under `results/` changes.
The tool versions of a run are in its `results/summary/run_info.txt`; cite those. Before this file, the merges are in `git log --merges`.

## Unreleased
- `results/summary/run_info.txt` has the lines `pfam_md5`, `rfam_md5` and `sprot_md5`: the md5 of the database files the annotation read (merge `c3f100e`). Results change: only `run_info.txt`, by these three lines. Two runs on one database folder give equal lines.
- Rule `kallisto` passes the strand flag (`--rf-stranded` or `--fr-stranded`) to kallisto when `strand_specific` is not empty. Before, `align_and_estimate_abundance.pl` (Trinity 2.15.2, also on Trinity master) dropped the flag: it appends `kallisto_add_opts` to itself, not the flag. Results change: for a run with non-empty `strand_specific`, the outputs of rule `kallisto` (`results/expression/kallisto/` and the six `results/summary/expr_*` tables). A run with `strand_specific: ""` is unchanged.
