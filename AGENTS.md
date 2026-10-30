# transXpress: notes for agents and contributors

Snakemake 9 pipeline: RNA-Seq reads to an annotated transcriptome. Results go into papers, so a run must give the same `results/` for the same input. Run commands and config keys: `README.md`. Adding a tool, the image workflow and the tests: `CONTRIBUTING.md`.

## Layout
- `workflow/Snakefile` includes `workflow/rules/*.smk`; `rule all` is the default target. `common.smk` holds config checks, helper functions, `run_info`, `clean`, `delete_intermediates`.
- `workflow/scripts/`: Python for `script:` rules; `scratch.sh` (node-local work folder) and `tools.sh` (shared shell functions) are sourced by rule shell blocks.
- `workflow/envs/<tool>.yaml` with `<tool>.linux-64.pin.txt`: one conda env per tool group. `Dockerfile` builds the main image from them.
- `config.yaml` (comments), `workflow/schemas/config.schema.yaml` (types, defaults). A key you add goes in both and in the README table.
- A run is `runs/<name>/` in a clone, started with `snakemake --profile profiles/<site> -d runs/<name>`.

## Rules
- A change must not alter any file under `results/`, unless the change is the point and the results stay reproducible. Prove it with two runs of the same commit compared by `tests/compare_runs.py`, and a comparison with a reference run.
- One branch per fix, one small test per commit. Merge with a merge commit.
- Run a Slurm job only when the task says so.

## Traps
- `workflow/envs/*.yaml`, the pin files and `tmhmm.post-deploy.sh` are hashed into the rule env folder (`/conda-envs/<hash>` in the image). A change, even a comment, changes the `Dockerfile` paths and needs another tag in `MAIN_IMAGE` (`workflow/Snakefile`). Any change to `Dockerfile`, even a comment, also needs another tag: `image.yml` stops when ghcr.io holds the tag from another `Dockerfile`.
- A change inside a `shell:` or `run:` block (a comment too) changes the rule code. Snakemake then lists the finished jobs of an old run folder for a rerun; pass `--rerun-triggers mtime input params software-env` to continue such a folder.
- `PERL_HASH_SEED=0` in the TransDecoder rules, `kallisto` and `trinity_stats` fixes the output order of Perl scripts. Keep it.
- Chunk jobs of pfam and rfam write `<out>.part` and take the name `<out>` only when the last line is `# [ok]` (`publish_part` in `tools.sh`). Blast, SignalP, TargetP and TMHMM rename after a clean exit.
- `retries: 0` on `rnaspades`, `trinity`, `kallisto`, `align_reads`: a retry asks for the same resources and fails the same way.
- `pfam_parallel` uses 1 thread: more threads change the name column of `--domtblout` between runs.
- In a profile or `--config`, `scratch` is written without `$` (`SCRATCH`): the sbatch line passes through a shell.
- A rule with a `script:` runs in its own env but loads the `snakemake` package of the main env, so both use Python 3.11.
- `trimmomatic_split` makes the folder that `trimmomatic_parallel` writes the trimmed reads into.
