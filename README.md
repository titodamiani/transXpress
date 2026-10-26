![transXpress](logo/Transxpress_Logo_RGB.png)

transXpress: a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for streamlined de novo transcriptome assembly and annotation

transXpress trims RNA-Seq reads, assembles them with Trinity or rnaSPAdes, finds the open reading frames (ORFs, the stretches of a transcript that code for a protein) with TransDecoder, annotates the transcripts and proteins, and estimates the expression of each transcript.

## Folder layout

| Path | Content |
|---|---|
| `workflow/Snakefile` | the main workflow file. It names the main image (`MAIN_IMAGE`), includes the rule files and holds the default target `all` |
| `workflow/rules/*.smk` | the rules, one file per step: `common.smk` (config checks, `run_info`, `clean`, `delete_intermediates`), `databases.smk`, `trim.smk`, `qc.smk`, `assembly.smk`, `orfs.smk`, `expression.smk`, `annotation.smk` |
| `workflow/scripts/` | the Python scripts of the `script:` rules, and `scratch.sh`, the shell functions of the rules that work on the node's local disk |
| `workflow/envs/` | one conda env file per tool group (`<env>.yaml`), and for each one the pin file `<env>.linux-64.pin.txt`, the exact list of packages for linux-64. `tmhmm.post-deploy.sh` adds `tmhmm.py` from PyPI to the tmhmm env. `default.yaml` is the main env, which holds Snakemake |
| `workflow/schemas/config.schema.yaml` | the type and the default of each config key. Snakemake checks `config.yaml` against it at the start of a run |
| `config.yaml` | the config of a run, with a comment for each key. Copy it into each run folder |
| `profiles/` | the Snakemake options of each place where jobs run: `local`, `slurm`, `aurum` |
| `Dockerfile` | the recipe of the main image |
| `containers/` | `signalp.def` and `targetp.def`, the recipes of the private SignalP and TargetP images |
| `tests/` | the test data, the test config, `make_run.sh`, `compare_runs.py` and `allowed_diffs.tsv` |
| `resources/adapters.fasta` | the adapter sequences for Trimmomatic |
| `.github/workflows/` | `test.yml` (the test runs) and `image.yml` (the pin files and the main image) |

## Dependencies

transXpress needs conda (Miniforge3) and Snakemake. Snakemake installs every tool below, in one of two ways: from the pin files into conda envs, or from the main image (see [Software: conda mode and image mode](#software-conda-mode-and-image-mode)).

| Env file in `workflow/envs/` | Tools |
|---|---|
| `default.yaml` (the main env) | Python 3.11, Snakemake 9.27.0, the executor plugins `slurm` 2.8.0 and `cluster-generic` 1.0.9, Biopython, pandas, matplotlib, seqkit, wget, sra-tools |
| `trimmomatic.yaml` | Trimmomatic 0.41 |
| `qc.yaml` | FastQC 0.12.1, MultiQC 1.35 |
| `rnaspades.yaml` | SPAdes 4.3.0, seqkit 2.14.0 |
| `trinity_utils.yaml` | Trinity 2.13.2, salmon, kallisto 0.48.0, bowtie2, samtools, RSEM, R 4.1 with edgeR and tidyverse |
| `transdecoder.yaml` | TransDecoder 6.0.0 |
| `busco.yaml` | BUSCO 6.1.0 |
| `pfam.yaml` | HMMER 3.4 |
| `rfam.yaml` | Infernal 1.1.5 |
| `blast.yaml` | NCBI BLAST+ 2.17.0 |
| `tmhmm.yaml` | tmhmm.py 1.2.1 |
| `signalp.yaml`, `targetp.yaml` | the base of the named conda envs `signalp` and `targetp`. [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) and [TargetP 2.0](https://services.healthtech.dtu.dk/services/TargetP-2.0/) have a DTU licence, so you install them yourself (see [SignalP and TargetP](#signalp-and-targetp)) |

The pin file of each env gives the exact version and build of every package. `results/run_info.txt` lists them for each run.

## Installation

1. Get the transXpress code:
   ~~~~
   git clone https://github.com/transXpress/transXpress.git
   ~~~~

2. Install [Miniforge3](https://github.com/conda-forge/miniforge):
   ~~~~
   curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
   bash Miniforge3-Linux-x86_64.sh
   rm Miniforge3-Linux-x86_64.sh
   ~~~~
   All commands below use `conda`. conda 23.10 and later use the libmamba solver by default, so the solve speed is the same as with mamba. Snakemake builds the rule environments with `conda`, its default frontend.

3. To make sure that the pipeline uses the R packages of its envs, unset `R_LIBS_SITE`:
   ~~~~
   unset R_LIBS_SITE
   ~~~~

4. Make the main transXpress env:
   ~~~~
   conda config --set channel_priority strict
   conda env create --name transxpress --file workflow/envs/default.yaml
   conda activate transxpress
   ~~~~
   On linux-64, `conda create --name transxpress --file workflow/envs/default.linux-64.pin.txt` gives the exact packages of the pin file.

5. Get the rule envs: pull the main image (image mode) or build the conda envs (conda mode). See [Software: conda mode and image mode](#software-conda-mode-and-image-mode).

6. Install SignalP and TargetP, or turn them off with `signalp: ""` and `targetp: ""`. See [SignalP and TargetP](#signalp-and-targetp).

7. Make a run folder. Put a copy of `config.yaml` there, and a tab-separated file called `samples.txt` that says where to find the raw read FASTQ files:
   ~~~
   cond_A      cond_A_rep1 A_rep1_left.fq    A_rep1_right.fq
   cond_A      cond_A_rep2 A_rep2_left.fq    A_rep2_right.fq
   cond_B      cond_B_rep1 B_rep1_left.fq    B_rep1_right.fq
   cond_B      cond_B_rep2 B_rep2_left.fq    B_rep2_right.fq
   ~~~
   The file names must not hold spaces. [tests/test_samples.txt](./tests/test_samples.txt) is an example.

8. Read `config.yaml` and set the keys for your data. Check at least `assembler`, `strand_specific`, `lineage`, `signalp` and `targetp`. See [Configuration](#configuration).

## Software: conda mode and image mode

Every rule that runs a tool names its env file with `conda:`. Snakemake can get these envs in two ways.

| Mode | Option | Where the rule envs come from |
|---|---|---|
| conda mode | `--software-deployment-method conda` (`profiles/local`, `profiles/slurm`) | Snakemake builds each env in the conda prefix from its pin file `<env>.linux-64.pin.txt`. On another platform it builds the env from the env file |
| image mode | `--software-deployment-method conda apptainer` (`profiles/aurum`) | Apptainer runs each rule in the main image, which holds every rule env, built from the pin files |

On one host, both modes give the same results. The CI job `t1-image` checks this (see [Tests](#tests)).

A rule without a `conda:` env runs in the main env on the host, in both modes. These rules use a Python script of `workflow/scripts/`, Python code in the rule, or a tool of the main env such as seqkit. So the main env (`transxpress`) is necessary in both modes.

### The main image

`MAIN_IMAGE` at the top of `workflow/Snakefile` names the main image:
~~~
MAIN_IMAGE = "docker://ghcr.io/titodamiani/transxpress:1.0.0"
containerized: MAIN_IMAGE
~~~
The image is public on ghcr.io. It holds the 10 rule envs of the `Dockerfile`, each one built from its pin file. Snakemake finds an env in the image at `/conda-envs/<hash>`. The hash covers the env file, its pin file and its post-deploy script, so a change to any of these three files needs another `Dockerfile` and another tag in `MAIN_IMAGE` (see [Adding a tool](#adding-a-tool)).

The workflow `.github/workflows/image.yml` builds the image and pushes it to ghcr.io. `test.yml` calls it on each push to master and on each pull request. It has two jobs:
* `pins` makes the pin file of each env file that has none, with `conda list --explicit --md5` in a clean linux-64 conda. It uploads the pin files as the artifact `pins` and fails, so that you commit them. A manual run of `image.yml` makes the pin files of all env files again and builds no image. An env file with a `name:` line (`signalp.yaml`, `targetp.yaml`) is a named env and gets no pin file.
* `image` checks the env paths of the `Dockerfile` against `snakemake --containerize`, then builds the image and pushes it under the tag of `MAIN_IMAGE`. When ghcr.io holds that tag from the same `Dockerfile`, the job builds nothing. When ghcr.io holds that tag from another `Dockerfile`, the job stops.

### Pull the main image first

Pull the image into an image folder before the first run. Without the file, the main Snakemake process pulls the image at the start of the run. The download is 5.2 GB and the SIF file (the single-file image format of Apptainer) is 3.1 GiB. The pull makes the SIF file with `mksquashfs`, which needed about 27 GB of memory on the Aurum cluster, so on a cluster pull it in a job with enough memory.

Snakemake looks for the image as `<md5 of the MAIN_IMAGE URL>.simg` in the folder of `--apptainer-prefix`. Without that option, the folder is `.snakemake/singularity` in the run folder.
~~~
IMAGES=/path/to/images
URL=docker://ghcr.io/titodamiani/transxpress:1.0.0
apptainer pull "$IMAGES/$(printf %s "$URL" | md5sum | cut -d ' ' -f 1).simg" "$URL"
~~~
For 1.0.0 the file name is `0e79edd4e66970b5267c55f94cb6b128.simg`. Put the image folder in the profile of your site:
~~~
software-deployment-method:
  - conda
  - apptainer
apptainer-prefix: /path/to/images
apptainer-args: "--bind /path/to/transXpress,/path/to/databases,/path/to/reads"
~~~
Apptainer shows a job only some folders of the host, such as `$HOME`, `/tmp`, the run folder and the `--bind` folders. The Apptainer config of a site can add more. A job in the image reads the scripts and the adapters of the transXpress folder, the databases of a shared `db_dir`, the reads and its temporary folder. Bind each of these that is outside the run folder. `profiles/aurum/config.yaml` is an example.

### Conda mode

Choose a folder for the rule envs, so that all your runs share them. Give it to Snakemake with `--conda-prefix /path/to/folder`, or as `conda-prefix: /path/to/folder` in your profile. Without it, Snakemake builds the envs in `.snakemake/conda` in each run folder.

Build the envs before the first run (this takes a while):
~~~~
snakemake --profile /path/to/transXpress/profiles/local --conda-prefix /path/to/folder --conda-create-envs-only
~~~~
On a cluster, build them in a job and not on the shared login node, for example on Slurm:
~~~~
sbatch --time=4:00:00 --mem=24G --cpus-per-task=4 -o envs.%j.out --wrap 'source /path/to/miniforge/etc/profile.d/conda.sh && conda activate transxpress && snakemake --profile /path/to/transXpress/profiles/local --conda-prefix /path/to/folder --conda-create-envs-only'
~~~~
The main Snakemake process builds the envs of the annotation rules that come after the split checkpoints during the first run. When the envs are in place, delete the package caches, which hold a second copy of every package:
~~~~
conda clean --all --yes
pip cache purge
~~~~

## SignalP and TargetP

SignalP 6.0 predicts signal peptides. TargetP 2.0 predicts the targeting peptides of proteins. Both come from DTU Health Tech under a licence that forbids sharing them, so each licensee installs them. Each tool runs in one of two ways:

| Key in `config.yaml` | Value | The rule runs in |
|---|---|---|
| `signalp_container` | `""` | the named conda env `signalp` |
| `signalp_container` | a path to a `.sif` | that Apptainer image, with no conda env |
| `targetp_container` | `""` | the named conda env `targetp` |
| `targetp_container` | a path to a `.sif` | that Apptainer image, with no conda env |

An image path needs `apptainer` in the `software-deployment-method` of the profile. Snakemake stops at the start when the file is missing or apptainer is not a deployment method. `signalp: ""` or `targetp: ""` turns the tool off, and the rule then needs neither an env nor an image.

### The SignalP and TargetP images

`containers/signalp.def` and `containers/targetp.def` are Apptainer recipes. Each one builds a conda env of fixed packages on the base image `condaforge/miniforge3:26.7.2-0` and copies the DTU package into it.

1. Download SignalP 6.0 fast and TargetP 2.0 from DTU, and unpack them.
2. Build the images:
   ~~~
   apptainer build --fakeroot --build-arg SIGNALP_PACKAGE=/path/to/signalp6_fast/signalp-6-package signalp.sif containers/signalp.def
   apptainer build --fakeroot --build-arg TARGETP_DIR=/path/to/targetp-2.0 targetp.sif containers/targetp.def
   ~~~
   The SignalP image is about 4.6 GiB and the TargetP image about 0.8 GiB. The SignalP build needed about 12 GB of memory on the Aurum cluster.
3. Give the paths in `config.yaml`:
   ~~~
   signalp_container: "/path/to/images/signalp.sif"
   targetp_container: "/path/to/images/targetp.sif"
   ~~~

Keep the images private. `results/run_info.txt` records the path, the size, the SIF ID and the build time of each image.

### The named conda envs

1. SignalP 6.0 fast. Download it from https://services.healthtech.dtu.dk/services/SignalP-6.0/ (go to Downloads), then:
   ~~~~
   conda env create -f workflow/envs/signalp.yaml
   conda activate signalp
   tar zxvf signalp-6.0h.fast.tar.gz
   cd signalp6_fast
   pip install signalp-6-package/
   SIGNALP_DIR=$(python -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
   cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/
   conda deactivate
   ~~~~
   Make sure that the conda python runs, or use the full path to python from your conda installation.

2. TargetP 2.0. Download it from https://services.healthtech.dtu.dk/services/TargetP-2.0/, unpack it, and add the targetp `bin` folder to the PATH of the targetp env with an activation script:
   ~~~~
   conda env create -f workflow/envs/targetp.yaml
   conda activate targetp
   tar zxvf targetp-2.0.Linux.tar.gz
   mkdir -p "$CONDA_PREFIX/etc/conda/activate.d" "$CONDA_PREFIX/etc/conda/deactivate.d"
   printf 'export _TARGETP_OLD_PATH="$PATH"\nexport PATH="$PATH:%s/targetp-2.0/bin"\n' "$(pwd)" > "$CONDA_PREFIX/etc/conda/activate.d/targetp.sh"
   printf 'export PATH="$_TARGETP_OLD_PATH"\nunset _TARGETP_OLD_PATH\n' > "$CONDA_PREFIX/etc/conda/deactivate.d/targetp.sh"
   conda deactivate
   ~~~~
   The script holds the full path of the `targetp-2.0` folder. If you move the folder, edit the path in `$CONDA_PREFIX/etc/conda/activate.d/targetp.sh`.

## Configuration

Each run folder holds a `config.yaml`. The `config.yaml` of the repository explains each key in a comment. `workflow/schemas/config.schema.yaml` gives the type of each key and the default of a key that the file leaves out. A value of the wrong type stops the run at the start.

| Key | Value in `config.yaml` | Meaning |
|---|---|---|
| `samples_file` | `"samples.txt"` | the samples file: condition, replicate, left reads, right reads |
| `trimmomatic_parameters` | `"ILLUMINACLIP:{adapters}:2:30:10 ..."` | the Trimmomatic steps. `{adapters}` stands for `resources/adapters.fasta` |
| `assembler` | `"trinity"` | `trinity` or `rnaspades` |
| `trinity_parameters` | `"--seqType fq"` | more options for Trinity |
| `strand_specific` | `""` | the strand option of the assembler, for example `"--SS_lib_type RF"` (Trinity) or `"--ss rf"` (rnaSPAdes). Also used for the read alignment |
| `kmers` | `""` | rnaSPAdes only: the k-mer sizes, for example `"-k 21,33,55"`. `""` lets rnaSPAdes choose |
| `rnaspades_filter` | `"default"` | rnaSPAdes only: the transcript set that becomes `transcriptome.fasta`: `default`, `hard` or `soft` |
| `scratch` | `""` | the environment variable that holds a folder on the node's local disk. See [Working files on the node's local disk](#working-files-on-the-nodes-local-disk) |
| `transdecoder_orfs` | `"all"` | the ORF set of `transcriptome.pep`: `all` or `single_best` |
| `annotated_fasta_prefix` | `""` | a text put before each sequence name in `transcriptome_annotated.*`, as `<prefix>\|<name>`, for example the species |
| `e_value_threshold` | `"1e-6"` | the E-value limit of the Pfam, Rfam and BLAST searches |
| `dispersion` | `"0.1"` | the edgeR dispersion for a pair of conditions that has a condition with one replicate |
| `signalp` | `"eukarya"` | the SignalP organism group: `eukarya` or `other`. `""` turns SignalP off |
| `targetp` | `"non-pl"` | the TargetP organism group: `pl` (plant) or `non-pl`. `""` turns TargetP off |
| `signalp_container`, `targetp_container` | `""` | see [SignalP and TargetP](#signalp-and-targetp) |
| `lineage` | `"viridiplantae_odb12"` | the BUSCO lineage dataset |
| `db_dir` | `"resources/db"` | the database folder. See [Annotation databases](#annotation-databases) |
| `pfam_release`, `rfam_release`, `sprot_release` | `"38.2"`, `"15.1"`, `"2026_03"` | the database releases. Keep the quotes |
| `annotation_chunk_fasta`, `annotation_chunk_orfs`, `annotation_chunk_pep` | `1000`, `5000`, `5000` | the number of sequences in each annotation job. See [Number of jobs](#number-of-jobs) |
| `<rule>_memory`, `<rule>_threads`, `<rule>_runtime` | see the file | the memory in GB, the thread count and the time limit of `rnaspades`, `trinity`, `align_reads`, `kallisto`, `busco`, `transdecoder_predict` and `sprot_blastx` (the rule `sprot_blastx_parallel`) |

The figures of the `<rule>_*` keys suit about 420 million read pairs (18 Illumina libraries). Lower them for a smaller dataset.

## Running transXpress

Run the pipeline in a run folder that holds `config.yaml` and the samples file. A profile is a folder in `profiles/` whose `config.yaml` holds the Snakemake options of the place where the jobs run:

| Profile | Where the jobs run | Rule envs |
|---|---|---|
| `profiles/local` | on this computer, with all its cores | conda mode |
| `profiles/slurm` | as Slurm jobs, through the Snakemake executor plugin `slurm`, at most 50 jobs in the queue | conda mode |
| `profiles/aurum` | as Slurm jobs on the Aurum cluster of IOCB Prague: `profiles/slurm` with the Aurum settings | image mode |

~~~~
conda activate transxpress
cd /path/to/run-folder
snakemake --profile /path/to/transXpress/profiles/local
~~~~
Each profile names the `workflow/Snakefile` of its own copy of transXpress (`snakefile: ../../workflow/Snakefile`), so the command works in any run folder. Options on the command line come after the profile, for example `--conda-prefix /path/to/folder`. A command line option replaces the profile value of the same option.

Each rule gives its thread count and two resources: `mem_mb`, the memory in MB, and `runtime`, the time limit. The rules whose size follows the dataset take these from `config.yaml`. The other rules keep fixed figures in their rule files.

### A profile for your site

1. Copy `profiles/slurm` to `profiles/<site>`. Keep it in `profiles/`, so that `snakefile: ../../workflow/Snakefile` finds the workflow.
2. Add the sbatch flags of your cluster as resources. See [Flags that belong to your cluster](#flags-that-belong-to-your-cluster).
3. For image mode, add `software-deployment-method`, `apptainer-prefix` and `apptainer-args`. See [Pull the main image first](#pull-the-main-image-first).
4. For conda mode, add `conda-prefix` with the folder of the shared rule envs.
5. To work on the node's local disk, add the `scratch` key under `config:`. See [Working files on the node's local disk](#working-files-on-the-nodes-local-disk).

`profiles/aurum/config.yaml` has all of these, each with a comment.

### Other schedulers

These two setups are not tested.

LSF: install the plugin `snakemake-executor-plugin-lsf` in the `transxpress` environment, and make a profile from `profiles/slurm` with these lines in place of `executor: slurm`:
~~~~
executor: lsf
default-resources:
  lsf_project: "'<your LSF project>'"
  lsf_queue: "'<your LSF queue>'"
~~~~

PBS: the plugin `cluster-generic` is in the `transxpress` environment. It submits each job with the command you give. Make a profile from `profiles/slurm` with these lines in place of `executor: slurm`:
~~~~
executor: cluster-generic
cluster-generic-submit-cmd: "qsub -j oe -o logs/ -l select=1:ncpus={threads}:mem={resources.mem_mb}mb -l walltime={resources.runtime}:00"
~~~~
`runtime` is in minutes, so `{resources.runtime}:00` is a PBS time in minutes and seconds. Without a status command (`cluster-generic-status-cmd`), Snakemake learns that a job ended only from a marker file the job writes, so a job that the scheduler stops leaves Snakemake waiting.

### Running the main Snakemake process on a Slurm cluster

The main Snakemake process needs the `transxpress` environment. Run it in its own job and activate the environment first. Replace the paths with yours:
~~~
cd /path/to/run-folder
sbatch --time=7-00:00:00 --mem=24G --cpus-per-task=2 -o main.%j.out --wrap 'source /path/to/miniforge/etc/profile.d/conda.sh && conda activate /path/to/transxpress-env && snakemake --profile /path/to/transXpress/profiles/slurm'
~~~
The time is the limit for the whole run, so it has to be longer than the longest job: `trinity` asks for 7 days. In conda mode, the main process builds the environments of the annotation rules during the run, and a build needs about 10 GB, so the job asks for 24 GB.

Inside a job, the `slurm` plugin prints a warning and removes the `SLURM_` variables of the main job from its own environment, so the jobs it submits get their own settings.

### Running specific steps

You can run specific steps of the pipeline by specifying a rule or a resulting file:
~~~
# run only the rules that make the FastQC and MultiQC report of the input data, to check its quality
snakemake --profile /path/to/transXpress/profiles/slurm results/qc/multiqc/before_trim/multiqc_fastqc.txt
~~~
~~~
# run only the rules to produce the samples_trimmed.txt file
snakemake --profile /path/to/transXpress/profiles/slurm results/trimmomatic/samples_trimmed.txt
~~~

## Output layout

A run writes into `results/`, with one folder per step. Each step folder holds its answer files at the top and its working files in a subfolder. The rule logs go to `logs/`, and the databases of the default `db_dir` go to `resources/`.

~~~
<run folder>/
  config.yaml  samples.txt
  logs/          <rule>.log, and <rule>/<wildcard>.log for a rule with many jobs
  resources/     db/  busco_downloads/
  results/
    run_info.txt
    qc/            fastqc/        before_trim/  after_trim/
                                  WARNING_fastqc_before_trim_overview.txt  FastQC_comparison_after_trim.txt
                   multiqc/       before_trim/  after_trim/  (each holds multiqc_report.html and multiqc_fastqc.txt)
                   busco/  busco_report.txt
                   transcriptome_stats.txt  transcriptome_exN50.tsv  ExN50_plot.pdf
    trimmomatic/   samples_trimmed.txt  samples_trimmed.yaml
                   trimmed_reads/ the trimmed reads
    assembly/      transcriptome.fasta  transcriptome.gene_trans_map
                   rnaspades/     the rnaSPAdes output folder
                   trinity/  trinity.Trinity.fasta  trinity.Trinity.fasta.gene_trans_map
    transdecoder/  transcriptome.orfs  transcriptome.pep
                   transcriptome_all_orfs.{pep,cds,gff3,bed}
                   transcriptome_single_best_orf.{pep,cds,gff3,bed}
                   longorfs/      the TransDecoder.LongOrfs folder
    annotation/    {pfam,sprotblastp}_orfs.out  {rfam,sprotblastx}_fasta.out
                   {signalp,targetp,tmhmm}_pep.out
                   transcriptome_annotated.fasta  transcriptome_annotated.pep  transcriptome_TPM_blast.csv
                   chunks/        the split inputs and the result of each chunk
    expression/    transcriptome_expression_isoform.tsv  transcriptome_expression_gene.tsv
                   kallisto/      kallisto.{isoform,gene}.* matrices, the index, one folder per replicate
                   bowtie2_rsem/  the indexes, one folder per replicate, igv/ (only with align_reads)
                   edger/         the edgeR tests, one set of files for each pair of conditions
~~~

The files most analyses start from:

| File | Content |
|---|---|
| `results/assembly/transcriptome.fasta` | the assembly |
| `results/assembly/transcriptome.gene_trans_map` | the gene of each transcript |
| `results/transdecoder/transcriptome.pep` | the proteins of the ORF set that `transdecoder_orfs` picks |
| `results/annotation/*.out` | the Pfam, Swiss-Prot, Rfam, SignalP, TargetP and TMHMM results |
| `results/annotation/transcriptome_annotated.fasta`, `.pep` | the transcripts and proteins with the annotations in their headers |
| `results/expression/transcriptome_expression_isoform.tsv`, `_gene.tsv` | the TMM-normalized TPM of each transcript and gene, or the TPM when TMM is not possible |
| `results/expression/kallisto/kallisto.{isoform,gene}.counts.matrix` | the estimated read counts |
| `results/qc/busco_report.txt` | the BUSCO score of the assembly |
| `results/run_info.txt` | the host, the processor, the transXpress commit, the deployment methods, the images and the version of each tool |

Two rules delete files, and neither runs unless you name it:

| Rule | Deletes | Keeps |
|---|---|---|
| `clean` | `results/` and `logs/` | `config.yaml`, the samples file, `resources/`, `.snakemake/` |
| `delete_intermediates` | the trimmed reads and the kallisto, bowtie2 and RSEM indexes | every result. A later `snakemake` run in the folder has nothing to do |

~~~
snakemake --profile /path/to/transXpress/profiles/local delete_intermediates
~~~

## Notes on results

* Swiss-Prot, Pfam and Rfam have fixed releases, set in `config.yaml`. See [Annotation databases](#annotation-databases).
* With the default `db_dir: "resources/db"`, BUSCO downloads the latest dataset of its lineage, so BUSCO scores can differ between runs months apart. A shared `db_dir` keeps one lineage dataset for every run.
* The pipeline reads and writes the paths in [Output layout](#output-layout). A run folder that holds its files at other paths runs every step again, so start each run in its own folder.
* The assembly depends on the assembler version, the thread count and the processor. rnaSPAdes repeats an assembly at a fixed thread count, and the SPAdes documentation gives no promise across processor types. The same reads on a different CPU generation, or with a different value of `rnaspades_threads`, can give a different set of transcripts. Two runs to be compared need the same version, the same thread count and the same hardware. `results/run_info.txt` records the processor, the host, the memory, the operating system and the transXpress commit, and the logs of `trimmomatic_parallel`, `rnaspades` and `trinity` name the machine each one ran on.
* TransDecoder and Trinity each have a Perl script that writes records in Perl hash order, which changes in each run. The TransDecoder rules, `kallisto` and `trinity_stats` set `PERL_HASH_SEED=0`, so two runs write the same records in the same order. Upstream reports: [TransDecoder/TransDecoder#212](https://github.com/TransDecoder/TransDecoder/issues/212), [trinityrnaseq/trinityrnaseq#1478](https://github.com/trinityrnaseq/trinityrnaseq/issues/1478).
* Trinity does not repeat an assembly: two runs on the same reads give slightly different transcripts. On one piper09 library, two runs shared 113127 of about 113370 transcripts. `trinity_parameters: "--seqType fq --inchworm_cpu 1"` gives identical runs on the test data, but not on real data. For a repeatable assembly, use rnaSPAdes. Details in [issue 22](https://github.com/titodamiani/transXpress/issues/22).
* Trinity writes the records of `trinity.Trinity.fasta` in another order in each run. `transcriptome_copy` sorts them by name (`seqkit sort --natural-order`) into `transcriptome.fasta` and writes the gene map in the same order.

## Annotation databases

transXpress annotates the transcripts with Swiss-Prot, Pfam and Rfam. BUSCO checks the assembly against a lineage dataset. The three databases take about 6.1 GB on disk: Pfam 4.8 GB, Rfam 0.7 GB, Swiss-Prot 0.6 GB. The download is about 0.6 GB.

`config.yaml` fixes the release of each database. Keep the quotes:
~~~
pfam_release: "38.2"
rfam_release: "15.1"
sprot_release: "2026_03"
~~~
The download checks Pfam and Swiss-Prot against the md5 their servers publish. Rfam publishes no checksum. The BUSCO server keeps only the latest dataset of each lineage, so a key cannot fix the BUSCO lineage. Only a kept copy can.

### One copy in each run folder

This is the default, `db_dir: "resources/db"`. Each run folder downloads its own databases into `resources/db/`. BUSCO downloads its lineage into `resources/busco_downloads/`.

### One shared copy for many run folders

1. Fill a shared folder once. Start from an empty folder, named after the releases:
   ~~~
   mkdir -p /path/to/databases/pfam38.2_rfam15.1_sprot2026_03
   cd /path/to/databases/pfam38.2_rfam15.1_sprot2026_03
   cp /path/to/transXpress/config.yaml .
   snakemake --profile /path/to/transXpress/profiles/local download_databases
   ~~~
   The copied `config.yaml` sets the releases and the BUSCO `lineage`. The databases go to `resources/db/`. `resources/db/versions.txt` lists each release and the date of the BUSCO lineage dataset. Cite these in a paper.

   Do not run `download_databases` in this folder a second time. Snakemake can replace files that other runs are reading.

2. In the `config.yaml` of each run folder, set `db_dir` to that `resources/db/` folder and name the lineage:
   ~~~
   db_dir: "/path/to/databases/pfam38.2_rfam15.1_sprot2026_03/resources/db"
   lineage: "viridiplantae_odb12"
   ~~~
   A run only reads this folder. BUSCO runs offline from `db_dir/busco_downloads`. If a database file is missing, or `lineage` is empty, Snakemake stops before any job starts and names the problem.

In image mode, add the shared folder to `apptainer-args` (`--bind`), so that the jobs in the image can read it.

## Tests

### A quick run

~~~~
cd tests
./run_test.sh
~~~~
`run_test.sh` runs the test data in `tests/` with `profiles/local` on 8 cores. It uses `tests/config.yaml`, which has SignalP and TargetP on, so both tools must be installed. It downloads the databases into `tests/resources/db`. When Graphviz (`dot`) is installed, it also draws the job graph of the run into `tests/dag.svg`.

### The reproducibility test

The reproducibility test runs the test data (T1) twice and checks that the two `results/` folders hold the same bytes. It needs a database folder filled by `download_databases` (see [One shared copy for many run folders](#one-shared-copy-for-many-run-folders)).

1. `tests/make_run.sh <folder> <assembler> <db_dir>` makes a T1 run folder directly inside the repository folder: the test data, `test_samples.txt` and a `config.yaml` from `tests/config.yaml`. In this config SignalP and TargetP are off, and Trinity runs with `--inchworm_cpu 1`.
2. Run the pipeline in each folder.
3. `tests/compare_runs.py <folder A> <folder B>` compares the two `results/` folders file by file.

~~~
tests/make_run.sh t1-a rnaspades /path/to/databases/pfam38.2_rfam15.1_sprot2026_03/resources/db
tests/make_run.sh t1-b rnaspades /path/to/databases/pfam38.2_rfam15.1_sprot2026_03/resources/db
(cd t1-a && snakemake --profile ../profiles/local --cores 4)
(cd t1-b && snakemake --profile ../profiles/local --cores 4)
python tests/compare_runs.py t1-a t1-b
~~~

How `compare_runs.py` compares:
* Every file must exist in both runs. Two files match when their bytes match. A gzip file matches when its decompressed bytes match.
* Before the comparison, it replaces three kinds of paths in the content, because many tools write the folder they ran in: the run folder becomes `@RUN@`, the node-local work folder of a job becomes `@SCRATCH@`, and the folder of a conda env (on the host or in the main image) becomes `@ENV@`.
* A file that differs passes only when a line of `tests/allowed_diffs.tsv` names it. Each line has three tab-separated fields: a path pattern, a mode and the reason. The mode `skip` lets the file differ in any way. The mode `strip:REGEX` compares the file without the lines that REGEX matches.
* The exit code is 0 when every file matches or is allowed to differ, and 1 otherwise.
* `--sequences results/assembly/transcriptome.fasta` compares only the sorted sequences of that file, without the headers.

### The CI workflow

`.github/workflows/test.yml` runs on each push to master, on each pull request and on a manual start. Its jobs:

| Job | What it does |
|---|---|
| `databases` | fills a database folder with `download_databases` at the releases of `test.yml`, and keeps it in the GitHub Actions cache |
| `t1` | for `rnaspades` and for `trinity`: two T1 runs in conda mode, then `compare_runs.py` |
| `image` | calls `image.yml`: the pin files, the env paths of the `Dockerfile`, the main image on ghcr.io |
| `t1-image` | one T1 run with `rnaspades` in conda mode and one in image mode (Apptainer 1.5.2), then `compare_runs.py` |

### Differences that the test does not cover

* TargetP output depends on the number of CPUs of its job: with more CPUs, some probabilities differ in the 6th decimal. `targetp_parallel` asks for 1 thread, and on Slurm each of its jobs gets 1 CPU. The test runs have TargetP off.
* Trinity: `trinity.Trinity.fasta` and its gene map have their records in another order in each run, so `allowed_diffs.tsv` skips them. `transcriptome.fasta` holds the records sorted. On real data the transcripts differ too, also with `--inchworm_cpu 1` (see [Notes on results](#notes-on-results)).
* rnaSPAdes writes the platform of the system it runs on into `results/assembly/rnaspades/params.txt`, with the glibc version. A run in the main image and a run in conda mode on a host with another glibc differ in that line.

## Adding a tool

A tool needs a rule, an env file and the pin file of that env. A tool with a setting also needs a config key.

1. A rule in the rule file of its step, `workflow/rules/<step>.smk`. A rule file for a step of its own also needs an `include:` line in `workflow/Snakefile`. Give the rule `conda: "../envs/<tool>.yaml"`, a `log:`, `threads:` and the resources `mem_mb` and `runtime`. When the output is a final result, add it to the input of `rule all` in `workflow/Snakefile`.
2. An env file `workflow/envs/<tool>.yaml` with the channels `conda-forge` and `bioconda`, a fixed tool version and no `name:` line. A tool from PyPI goes in a post-deploy script, as in `workflow/envs/tmhmm.post-deploy.sh`.
3. The pin file `workflow/envs/<tool>.linux-64.pin.txt`. Push the branch: the `pins` job of `image.yml` makes the pin file, uploads it as the artifact `pins` and fails. Commit the pin file of the artifact.
4. When the tool has a setting, a key in `config.yaml` with a comment, and the same key in `workflow/schemas/config.schema.yaml` with its type and default. The rule reads it as `config["<key>"]`.

For image mode, the main image needs the env too:

1. Run `snakemake --containerize` in a T1 run folder, as the `image` job of `image.yml` does.
2. Copy the block of the env into the `Dockerfile`, and edit it as the other blocks: a `COPY` of the pin file, `conda create --file` of the pin file, and the post-deploy script when there is one.
3. Give `MAIN_IMAGE` in `workflow/Snakefile` another tag, for example `1.1.0`.

The `image` job checks the env paths of the `Dockerfile` against `snakemake --containerize`, then builds the image and pushes it to ghcr.io. A change to an env file, a pin file or a post-deploy script needs the same 3 steps.

## Align reads to the transcriptome assembly and visualize the results in IGV

To align reads to the transcriptome assembly and visualize the results in IGV, use these commands:
~~~~
snakemake --profile /path/to/transXpress/profiles/slurm align_reads
snakemake --profile /path/to/transXpress/profiles/slurm IGV
~~~~

Then load your transcriptome file to IGV: Genomes -> Load Genome from File -> select the file *results/assembly/transcriptome.fasta*

Your sorted .bam files will be in the igv folder:
~~~~
results/expression/bowtie2_rsem/igv/{sample}.sorted.bam
results/expression/bowtie2_rsem/igv/{sample}.sorted.bam.bai
~~~~
Load them to IGV: File -> Load from File -> select the *results/expression/bowtie2_rsem/igv/{sample}.sorted.bam* files

## Possible problems when executing on cluster systems

### Time limit
Slurm jobs need a time limit. Each rule gives it as its `runtime` resource. Most rules ask for 4 hours. `config.yaml` sets the long ones: `rnaspades_runtime` 1 day, `trinity_runtime` 7 days, `align_reads_runtime` 3 days and `sprot_blastx_runtime` 12 hours. A time limit is a number of minutes or a text such as `"30m"`, `"12h"` or `"7d"`.

If a job is cancelled because of the time limit, raise the figure in `config.yaml`, or for one run on the command line, for example `--set-resources trinity:runtime=10d`. The time you need depends on the size of the reads used for the assembly. A job whose time limit reaches into a maintenance of the cluster waits until the maintenance ends, so before a maintenance give the long rules a limit that ends in time.

Many clusters start short jobs sooner. Read the time and the memory a finished job used with `sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS`, and set the limits a margin above them.

See https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Computing-Requirements

### How the pipeline learns that a cluster job stopped
The `slurm` plugin asks Slurm for the state of the jobs of the run with `sacct`. So a job that the node, the kernel or `scancel` stops counts as failed. The first check comes 40 seconds after a submission. While no job ends, the wait between two checks grows from 40 seconds to at most 180 seconds.

`restart-times: 1` in the profile gives a failed job another try. Four rules set `retries: 0` and get no second try: `rnaspades`, `trinity`, `kallisto` and `align_reads`. They fail mostly on memory or time, and a retry asks for the same figures and starts from the beginning.

The plugin writes the Slurm log of each job to `.snakemake/slurm_logs/rule_<rule>/<wildcards>/<job id>.log` in the run folder. It deletes the log of a job that succeeds and keeps the log of a job that fails. `--slurm-keep-successful-logs` keeps every log. The rule logs in `logs/` stay in both cases.

### Number of jobs
The profile keeps at most 50 jobs in the queue (`jobs: 50`) and submits at most one job per second (`max-jobs-per-timespan: "1/1s"`), so a run does not flood the scheduler.

The annotation rules run one job per chunk of sequences and tool. Each job reads the whole database, so larger chunks mean fewer jobs. The chunk sizes are the config keys `annotation_chunk_fasta` (1000 transcripts, for rfam and blastx), `annotation_chunk_orfs` (5000, for pfam and blastp) and `annotation_chunk_pep` (5000, for tmhmm, signalp and targetp). At these sizes a chunk job takes up to about 2.5 hours for blastx and about 1 hour for the others. The merged results follow the order of the sequences at any chunk size.

### Flags that belong to your cluster
Some sbatch flags are site specific: a node to avoid, a processor generation to insist on, an account to charge. They go in the profile of your site, as resources that the `slurm` plugin turns into sbatch flags. `profiles/aurum/config.yaml` is an example.

`default-resources` applies to every job of a run:
~~~~
default-resources:
  constraint: "'gen-d'"
  slurm_account: "'my_account'"
  slurm_extra: "'--exclude=badnode1,badnode2'"
~~~~
A text value needs two pairs of quotes, because Snakemake reads a resource value as a Python expression. `slurm_extra` holds the flags that have no resource of their own. It cannot hold a flag that the plugin sets, such as `--constraint`, `--partition`, `--mem`, `--time` or `--cpus-per-task`.

An assembly that is to be compared with another one needs the same processor on both runs, so a constraint that names the CPU generation goes in `default-resources`. Slurm records no processor model, so a feature such as `gen-d` is as close as a constraint gets.

`set-resources` applies to one rule, for example a rule that has to land on a named node:
~~~~
set-resources:
  rnaspades:
    slurm_partition: "'bigmem'"
    slurm_extra: "'--exclude=badnode1,badnode2 --nodelist=node42'"
~~~~
A rule value replaces the default value, so repeat the flags of `default-resources` that the rule still needs.

A node name and a partition name usually travel together. The partition a cluster picks for the memory of a job does not have to be one that holds the node you named, and sbatch then refuses the job.

### Working files on the node's local disk
On a cluster, the run folder usually sits on network storage that all users share, and many clusters ask that work with many or large files runs on the node's own disk.

The `scratch` key in `config.yaml` names the environment variable that holds a folder on the node's local disk in a cluster job, with its `$`:
~~~~
scratch: "$SCRATCH"
~~~~

Ask your cluster documentation which variable that is. On some clusters `$TMPDIR` is in memory, and a file there counts against the memory of the job. The default `""` runs every rule in the run folder.

In a profile (`config:`) or with `--config` on the command line, write the name without its `$`, for example `scratch: "SCRATCH"`. Snakemake passes these values to each cluster job on the `sbatch` line, and the shell there would put the folder of the main job in place of `$SCRATCH`.

With the key set, every job gets that folder as its temporary folder: Snakemake sets `TMPDIR`, `TMP` and `TEMP` of the job to it. A `tmpdir` value for one rule, for example `--set-resources fastqc:tmpdir=/data/tmp`, replaces it for that rule. `signalp_parallel`, `targetp_parallel` and `tmhmm_parallel` write their work files there. Five rules make their own folder under it, work there, and copy their files to the run folder when they end:

| Rule | Works on the local disk | Copied to the run folder |
|---|---|---|
| `rnaspades` | the whole rnaSPAdes output folder | `results/assembly/rnaspades`, which replaces the folder in the run folder |
| `trinity` | the whole Trinity assembly, with its `read_partitions` folder | `results/assembly/trinity.Trinity.fasta` and its `.gene_trans_map` |
| `kallisto` | the index and the per replicate results | every file, over the files in `results/expression/kallisto` |
| `align_reads` | the indexes and the per replicate results | every file, over the files in `results/expression/bowtie2_rsem` |
| `busco` | the BUSCO output folder, whose per-gene folders it packs into two archives | `results/qc/busco`, which replaces the folder in the run folder |

Each of these logs names its work folder near the top. When the variable is not set in a job, for example on a computer without a cluster, the rule works in the run folder.

`rnaspades`, `trinity`, `kallisto` and `align_reads` stop when the folder is in memory (a `tmpfs` file system), with a message that names the node, because their work files would count against the memory of the job. Exclude such nodes, or set the key to a variable that names a disk.

When `rnaspades`, `kallisto` or `align_reads` fails, for example on an error of the tool or when the kernel stops it for lack of memory, it copies its local files to the place the rule writes them without the local disk, so they can be read to find the cause. `trinity` and `busco` copy nothing on a failure: their work files are the many small files to keep off the shared storage, and Snakemake deletes the output folder of a failed `busco` anyway. A job that Slurm stops at its time limit or with `scancel` copies nothing: Slurm ends the Snakemake process of the job at once, and with it every process the job started. A node that crashes leaves nothing either. In all these cases the rule log `logs/<rule>.log` holds the output of the tool up to the stop, and the Slurm log in `.snakemake/slurm_logs/rule_<rule>/` names a time limit.

The local disk of a node is shared by the jobs on it, and a cluster may not reserve disk space for a job. `rnaspades` needed about 120 GB of working files for 9 libraries of about 23 million read pairs each. `trinity` needed 62 GB and wrote 124,280 files in `read_partitions` for one library of 23.8 million read pairs.

Without the local disk, `trinity` works in `results/assembly/trinity`, and a run after a failure continues from Trinity's own checkpoints there. `--full_cleanup` deletes that folder after a successful assembly.

`rnaspades` runs without `--checkpoints`. On the local disk, the folder is deleted at the end of the job, so a restart point would not survive.

In image mode, the jobs in the image read the folder of the `scratch` key. Bind it with `apptainer-args` when it is outside the run folder, as `profiles/aurum` does with `/scratch`.
