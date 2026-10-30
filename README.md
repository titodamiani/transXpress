![transXpress](logo/Transxpress_Logo_RGB.png)

transXpress: a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for de novo transcriptome assembly and annotation.

It trims RNA-Seq reads, assembles them with Trinity or rnaSPAdes, predicts ORFs (open reading frames) with TransDecoder, annotates transcripts and proteins (Pfam, Rfam, Swiss-Prot, TMHMM, SignalP, TargetP), estimates expression with kallisto and runs BUSCO, FastQC and MultiQC.

## Repository structure

```
transXpress/
├── workflow/        # Snakefile, rules/, scripts/, envs/ (one conda env per tool group), schemas/
├── profiles/        # Snakemake profiles: local, slurm, aurum
├── config.yaml      # Config of a run, one comment line per key. Copy it into each run folder
├── runs/<name>/     # One run folder per run: config.yaml, samples.txt, results/, logs/ (not tracked)
├── tests/           # Test data (T1), make_run.sh, compare_runs.py, allowed_diffs.tsv
├── containers/      # Recipes of the private SignalP and TargetP images
└── Dockerfile       # The main image with every rule env
```

## Installation

1. Clone the repository and install [Miniforge3](https://github.com/conda-forge/miniforge).
   ```bash
   git clone https://github.com/transXpress/transXpress.git
   ```
2. Make the main env. It holds Snakemake 9.27.0 and the executor plugins.
   ```bash
   unset R_LIBS_SITE    # the rules must use the R packages of their own envs
   conda config --set channel_priority strict
   conda env create --name transxpress --file workflow/envs/default.yaml
   conda activate transxpress
   ```
3. Get the rule envs, in one of two modes (set in the profile with `software-deployment-method`):
   - Conda mode (`profiles/local`, `profiles/slurm`): Snakemake builds the envs from the pin files in `workflow/envs/` during the first run. To build them before, run `snakemake --profile profiles/local --conda-prefix /path/to/envs --conda-create-envs-only`. Give the same `--conda-prefix` (or `conda-prefix:` in the profile) to every run.
   - Image mode (`profiles/aurum`): Apptainer runs each rule in the main image. Pull it once (about 5 GB, and `mksquashfs` needs a lot of memory, so pull in a cluster job):
     ```bash
     IMAGES=/path/to/images
     URL=docker://ghcr.io/titodamiani/transxpress:1.0.1    # MAIN_IMAGE in workflow/Snakefile
     apptainer pull "$IMAGES/$(printf %s "$URL" | md5sum | cut -d ' ' -f 1).simg" "$URL"
     ```
     In the profile: `software-deployment-method: [conda, apptainer]`, `apptainer-prefix: /path/to/images` and `apptainer-args: "--bind <clone>,<databases>,<reads>"` (a job sees only the run folder and the bound folders). `profiles/aurum/config.yaml` is an example.
4. Install SignalP and TargetP, or turn them off with `signalp: ""` and `targetp: ""`. See [SignalP and TargetP](#signalp-and-targetp).

## SignalP and TargetP

SignalP 6.0 fast and TargetP 2.0 have a DTU licence that forbids sharing them. Download both from DTU Health Tech and unpack them. Each tool runs either in an Apptainer image (needs `apptainer` in `software-deployment-method`) or in a conda env named `signalp` or `targetp` (steps in `containers/README.md`).

```bash
apptainer build --fakeroot --build-arg SIGNALP_PACKAGE=/path/to/signalp6_fast/signalp-6-package signalp.sif containers/signalp.def
apptainer build --fakeroot --build-arg TARGETP_DIR=/path/to/targetp-2.0 targetp.sif containers/targetp.def
apptainer exec signalp.sif signalp6 --fastafile <in.pep> --organism eukarya --output_dir <dir> --format none --mode fast    # test run
apptainer exec targetp.sif targetp -fasta <in.pep> -format short -org <pl|non-pl> -prefix <name>                           # test run
```
Set the image paths in `config.yaml`: `signalp_container: "/path/to/signalp.sif"` and `targetp_container: "/path/to/targetp.sif"`. With `""`, the rule uses the named conda env. Keep the images private.

## Running transXpress

Each run has a folder `runs/<name>/` in the clone, with a `config.yaml` (copy of the repository one) and a `samples.txt`: tab-separated, one library per line, two or more lines, no spaces in file names. Columns: condition, replicate name (unique), left reads, right reads (leave out for single-end reads). `tests/test_samples.txt` is an example.

```bash
conda activate transxpress
cd /path/to/transXpress
snakemake --profile profiles/local -d runs/<name>        # or profiles/slurm, profiles/aurum
```

| Profile | Jobs run | Rule envs |
|---|---|---|
| `profiles/local` | on this computer | conda mode |
| `profiles/slurm` | as Slurm jobs (50 at once). Start the main Snakemake process in its own job, with the `transxpress` env active | conda mode |
| `profiles/aurum` | as Slurm jobs on the Aurum cluster of IOCB Prague | image mode |

For another scheduler, copy `profiles/slurm` and use the matching [Snakemake executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/). Only Slurm is tested.

Targets after the profile: a rule or a file (`results/qc/multiqc/before_trim/multiqc_fastqc.txt` runs only the read QC); `align_reads` and `IGV` align the reads and sort the BAM files for IGV (not run by default); `delete_intermediates` deletes the trimmed reads and the kallisto, bowtie2 and RSEM indexes; `clean` deletes `results/` and `logs/`.

## Configuration

`config.yaml` has one comment line per key. `workflow/schemas/config.schema.yaml` gives the type and default of each key, and a wrong type or a key that is not listed stops the run at the start. Check at least `assembler`, `strand_specific`, `lineage`, `signalp`, `targetp` and the resources.

| Key | Default | Meaning |
|---|---|---|
| `samples_file` | `samples.txt` | the samples file |
| `trimmomatic_parameters` | `ILLUMINACLIP:{adapters}:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25` | Trimmomatic steps; `{adapters}` is `resources/adapters.fasta` |
| `assembler` | `trinity` | `trinity` or `rnaspades` |
| `trinity_parameters` | `--seqType fq` | extra Trinity options |
| `strand_specific` | `""` | Trinity: `--SS_lib_type RF` (paired RF or FR, single F or R). rnaSPAdes: `--ss rf` or `--ss fr`. Also used for the read quantification |
| `kmers` | `""` | rnaSPAdes k-mer sizes, odd, below 128, ascending, for example `-k 21,33,55`. `""` lets rnaSPAdes choose |
| `rnaspades_filter` | `default` | the rnaSPAdes set that becomes `transcriptome.fasta`: `default`, `hard` (long, well expressed transcripts) or `soft` (also short, low expressed ones) |
| `scratch` | `""` | the environment variable that holds a node-local folder in a cluster job (`$SCRATCH`). The heavy rules work there and copy their files back. In a profile or `--config` write it without the `$` (`SCRATCH`). `""` works in the run folder |
| `annotated_fasta_prefix` | `""` | text put before each name in the annotated files, as `<prefix>\|<name>` |
| `e_value_threshold` | `1e-6` | E-value limit of Pfam, Rfam and BLAST |
| `dispersion` | `0.1` | edgeR dispersion for a pair of conditions where a condition has one replicate |
| `signalp`, `targetp` | `eukarya`, `non-pl` | organism group (`eukarya` or `other`; `pl` or `non-pl`). `""` turns the tool off |
| `signalp_container`, `targetp_container` | `""` | path of a `.sif` image, or `""` for the conda env |
| `lineage` | `viridiplantae_odb12` | BUSCO lineage. `""` lets BUSCO choose |
| `db_dir` | `resources/db` | database folder, see [Annotation databases](#annotation-databases) |
| `pfam_release`, `rfam_release`, `sprot_release` | `38.2`, `15.1`, `2026_03` | database releases; keep the quotes |
| `annotation_chunk_fasta`, `_orfs`, `_pep` | `1000`, `5000`, `5000` | sequences in each annotation job (one job per chunk and tool; a chunk takes up to 2.5 h for blastx, about 1 h for the others) |
| `<rule>_memory`, `_threads`, `_runtime` | see the file | GB, threads and time limit of `rnaspades`, `trinity`, `align_reads`, `kallisto`, `busco`, `transdecoder_predict`, `sprot_blastx`. The defaults suit 18 libraries (about 420 million read pairs). A time limit is minutes or text such as `12h`, `7d` |

The assembly depends on the assembler version, the thread count and the CPU: keep them the same in runs that you compare. For one run, change a limit on the command line: `--set-resources trinity:runtime=10d`. Site flags of Slurm (constraint, account, excluded nodes) go in the profile as `default-resources` or `set-resources`, with two pairs of quotes for text (`slurm_extra: "'--exclude=node1'"`); see `profiles/aurum/config.yaml`.

## Annotation databases

With `db_dir: "resources/db"`, each run downloads Pfam, Rfam and Swiss-Prot (about 6 GB) at the releases of `config.yaml` and checks the md5 of Pfam and Swiss-Prot. BUSCO downloads the latest dataset of its lineage, so BUSCO scores can differ between runs months apart.

For many runs, fill one shared folder once, and never run the target in it a second time:
```bash
mkdir -p /path/to/databases/pfam38.2_rfam15.1_sprot2026_03 && cd /path/to/databases/pfam38.2_rfam15.1_sprot2026_03
cp /path/to/transXpress/config.yaml .
snakemake --profile /path/to/transXpress/profiles/local download_databases
```
`resources/db/versions.txt` lists the releases and the date of the BUSCO lineage. Cite them. Then set `db_dir: "/path/to/databases/.../resources/db"` and a named `lineage` in each run config. A shared folder is only read, BUSCO runs offline, and a missing file stops the run before the first job. In image mode, bind the folder.

## Output

```
runs/<name>/results/
├── summary/                          # the files to use
│   ├── transcriptome.fasta           # assembly (Trinity records sorted by name)
│   ├── transcriptome.gene_trans_map  # gene of each transcript
│   ├── transcriptome_annotated.fasta # transcripts with expression, blastx and rfam in the header
│   ├── transcriptome_{all_orfs,single_best_orf}.pep            # proteins of the two ORF sets
│   ├── transcriptome_{all_orfs,single_best_orf}_annotated.pep  # proteins with TMM, blastp, pfam, rfam, tmhmm, signalp, targetp
│   ├── expr_{tmm,tpm,counts}_{isoform,gene}.tsv  # TMM-scaled TPM, TPM, estimated counts
│   ├── transcriptome_stats.txt       # Trinity statistics
│   ├── busco_report.txt              # BUSCO score
│   └── run_info.txt                  # main-job and assembly node, commit, images, database releases and md5, tool versions
├── annotation/                       # {pfam,sprotblastp}_orfs.out, {rfam,sprotblastx}_fasta.out, {signalp,targetp,tmhmm}_pep.out, chunks/
├── transdecoder/                     # transcriptome_{all_orfs,single_best_orf}.{pep,cds,gff3,bed}, longorfs/
├── assembly/                         # rnaspades/ or trinity.Trinity.fasta (+ .gene_trans_map)
├── expression/                       # kallisto/ (matrices), edger/ (one set per condition pair), bowtie2_rsem/ (only with align_reads)
├── trimmomatic/                      # samples_trimmed.txt, trimmed_reads/
└── qc/                               # fastqc/, multiqc/, assembly/ (ExN50 table and plot, busco/)
```

Notes:
- `all_orfs` keeps every ORF that TransDecoder accepts (a transcript can give several proteins). `single_best_orf` keeps one ORF per transcript, picked by Pfam or Swiss-Prot hit first and by length second.
- Two runs of the same input give the same files, except the ones in `tests/allowed_diffs.tsv`. TransDecoder, `kallisto` and `trinity_stats` set `PERL_HASH_SEED=0` for this.
- rnaSPAdes repeats an assembly at a fixed thread count on the same CPU generation. Trinity does not: two runs differ slightly even with `--inchworm_cpu 1` ([issue 22](https://github.com/titodamiani/transXpress/issues/22)). For a repeatable assembly use rnaSPAdes.
- A run folder that holds its files at other paths runs every step again. Start each run in its own folder.

## Tests

`tests/make_run.sh <folder> <assembler> <db_dir>` makes a T1 run folder (test data, SignalP and TargetP off) directly in the clone. Run each folder twice, then compare:
```bash
tests/make_run.sh t1-a rnaspades /path/to/databases/pfam38.2_rfam15.1_sprot2026_03/resources/db
(cd t1-a && snakemake --profile ../profiles/local --cores 4)     # and the same in t1-b
python tests/compare_runs.py t1-a t1-b       # exit 0 when every file matches or is listed in tests/allowed_diffs.tsv
```
`tests/run_test.sh` runs the test data once with `profiles/local` and `tests/config.yaml` (SignalP and TargetP must be installed). CI (`.github/workflows/test.yml`) runs the paired test for both assemblers and in image mode. Details, and how to add a tool, are in `CONTRIBUTING.md`.
