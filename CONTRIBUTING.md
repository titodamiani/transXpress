# Contributing

Read `AGENTS.md` first (layout, rules, traps). Run commands and config keys are in `README.md`.

## Workflow
- One branch per fix. Run a small test for each commit.
- Before a merge, run a fresh paired test of the fix: two runs of the same commit (T1 folders from `tests/make_run.sh`, for `rnaspades` and `trinity`), compared with `tests/compare_runs.py`. Compare a T1 run with a reference run of the previous release too.
- Merge with a merge commit.

## Tests
- `tests/make_run.sh <folder> <assembler> <db_dir>`: a T1 run folder directly inside the clone (not under `runs/`) with the test data, `test_samples.txt` and a config from `tests/config.yaml`. SignalP and TargetP are off and Trinity runs with `--inchworm_cpu 1`. `db_dir` is a folder made by the `download_databases` target.
- `tests/compare_runs.py A B` compares the `results/` of two runs. Every file must exist in both and match byte for byte; a gzip file matches when its decompressed bytes match. Before the comparison it replaces the run folder path with `@RUN@`, the node-local work folder with `@SCRATCH@` and a conda env folder with `@ENV@`. A file that differs passes only when a line of `tests/allowed_diffs.tsv` names it. Fields, separated by tabs: path pattern (fnmatch, `*` matches `/`), mode, reason. Modes: `skip` (any difference, or in one run only), `strip:REGEX` (compare without the lines REGEX matches), `sorted:REGEX` (the same, then runs of spaces count as one and the lines are sorted). Exit code 0 when every file matches or is allowed to differ. `--sequences FASTA` compares only the sorted sequences of one fasta file.
- Differences the test does not cover: TargetP probabilities depend on the CPU count of the job (6th decimal; `targetp_parallel` asks for 1 thread and the T1 runs have TargetP off). Trinity records and `trinity.Trinity.fasta.gene_trans_map` come in another order in each run (`allowed_diffs.tsv` skips them; `transcriptome.fasta` is sorted). rnaSPAdes writes the glibc version into `results/assembly/rnaspades/params.txt`, so an image run and a conda run on hosts with another glibc differ in that line.
- CI (`.github/workflows/test.yml`), on each push to master, pull request and manual start:

| Job | What it does |
|---|---|
| `keys` | `tests/check_keys.py`: `config.yaml`, `config.schema.yaml` and the README table name the same keys |
| `databases` | fills a database folder with `download_databases` at the releases of `test.yml` and caches it |
| `t1` | for `rnaspades` and `trinity`: two T1 runs in conda mode, then `compare_runs.py` |
| `image` | calls `image.yml`: pin files, env paths of the `Dockerfile`, the main image on ghcr.io |
| `t1-image` | one T1 run with `rnaspades` in conda mode and one in image mode (Apptainer 1.5.2), then `compare_runs.py` |

## The main image
`MAIN_IMAGE` in `workflow/Snakefile` names the image (`docker://ghcr.io/titodamiani/transxpress:<tag>`). It holds the 10 rule envs of the `Dockerfile`, each built from its pin file, at `/conda-envs/<hash>`. The hash covers the env file, its pin file and its post-deploy script.

`image.yml` has two jobs:
- `pins` makes the pin file of each env file that has none (`conda list --explicit --md5` in a clean linux-64 conda), uploads the pin files as the artifact `pins` and fails, so that you commit them. A manual run makes the pin files of all env files again and builds no image. An env file with a `name:` line (`signalp.yaml`, `targetp.yaml`) is a named env and has no pin file.
- `image` checks the env paths of the `Dockerfile` against `snakemake --containerize`, then builds the image and pushes it under the tag of `MAIN_IMAGE`. It builds nothing when ghcr.io holds that tag from the same `Dockerfile`, and stops when ghcr.io holds it from another `Dockerfile` (a `Dockerfile` change needs another tag).

## Adding a tool
A tool needs a rule, an env file and the pin file of that env. A tool with a setting also needs a config key.
1. A rule in `workflow/rules/<step>.smk`. A rule file for a step of its own needs an `include:` line in `workflow/Snakefile`. Give the rule `conda: "../envs/<tool>.yaml"`, a `log:`, `threads:` and the resources `mem_mb` and `runtime`. If the output is a final result, add it to `rule all`.
2. `workflow/envs/<tool>.yaml` with the channels `conda-forge` and `bioconda`, a fixed tool version and no `name:` line. A tool from PyPI goes in a post-deploy script, as in `workflow/envs/tmhmm.post-deploy.sh`.
3. The pin file `workflow/envs/<tool>.linux-64.pin.txt`: push the branch, the `pins` job fails and uploads it, commit that file.
4. A setting is a key in `config.yaml` (one comment line), in `workflow/schemas/config.schema.yaml` (type, default) and in the README table. The rule reads `config["<key>"]`.
5. For image mode: run `snakemake --containerize` in a T1 run folder, copy the block of the env into the `Dockerfile` (a `COPY` of the pin file, `conda create --file` of the pin file, the post-deploy script if any), and give `MAIN_IMAGE` another tag.
