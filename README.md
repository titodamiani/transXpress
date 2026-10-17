![transXpress](logo/Transxpress_Logo_RGB.png)

transXpress: a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for streamlined de novo transcriptome assembly and annotation

## Intro

## Dependencies

transXpress requires:
* snakemake 9.27.0 with the executor plugins `slurm` and `cluster-generic` (install via conda, `workflow/envs/default.yaml`)
* fastqc (install via conda, `envs/qc.yaml`)
* multiqc (install via conda, `envs/qc.yaml`)
* trimmomatic (install via conda, `envs/trimmomatic.yaml`)
* Trinity (install via conda, `trinity_utils.yaml`)
* SPAdes (install via conda, `envs/rnaspades.yaml`)
* TransDecoder (install via conda, `transdecoder.yaml`)
* BioPython (install via conda, `envs/default.yaml`)
* samtools (install via conda, `envs/trinity_utils.yaml`)
* bowtie2 (install via conda, `envs/trinity_utils.yaml`)
* infernal (install via conda, `envs/rfam.yaml`)
* HMMER (install via conda, `envs/pfam.yaml`)
* kallisto (install via conda, `envs/trinity_utils.yaml`)
* NCBI BLAST+ (install via conda, `envs/blast.yaml`)
* R (install via conda, `envs/trinity_utils.yaml`)
* edgeR (install via conda, `envs/trinity_utils.yaml`)
* seqkit (install via conda, `envs/default.yaml`, `envs/rnaspades.yaml`)
* wget (install via conda, `envs/default.yaml`)
* sra-tools (install via conda, `envs/default.yaml`)
* tidyverse (required for Trinity, install via conda, `envs/trinity_utils.yaml`)
* python (install via conda, `workflow/envs/default.yaml`)
* busco 4+ (install via conda, `envs/busco.yaml`)
* rsem (install via conda, `envs/trinity_utils.yaml`)
* [SignalP 6.0](https://services.healthtech.dtu.dk/service.php?SignalP)
* [TargetP 2.0](https://services.healthtech.dtu.dk/service.php?TargetP-2.0)
* tmhmm.py (install via pip, `workflow/envs/tmhmm.yaml`)
* basic Linux utitilies: split, awk, cut, gzip

The conda dependencies are installed in smaller conda environments automatically by transXpress (based on yaml files in the `envs` directory). 

## Installation

1. Checkout the transXpress code into the folder in which you will be performing your assembly:
~~~~
git clone https://github.com/transXpress/transXpress.git
~~~~

2. Install [Miniforge3](https://github.com/conda-forge/miniforge)
~~~~
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh
rm Miniforge3-Linux-x86_64.sh
~~~~
All commands below use `conda`. conda 23.10 and later use the libmamba solver by default, so the solve speed is the same as with mamba. Snakemake builds the rule environments with `conda`, its default frontend.

3. To ensure correct versions of R packages will be used unset R_LIBS_SITE
~~~~
unset R_LIBS_SITE
~~~~

4. Setup main transXpress conda environment:
~~~~
conda config --set channel_priority strict
conda env create --name transxpress --file workflow/envs/default.yaml
conda activate transxpress
~~~~
`workflow/envs/default.yaml` holds Snakemake and its executor plugins for Slurm (`slurm`) and for other schedulers (`cluster-generic`).

5. Choose a folder for the conda environments of the pipeline rules, so that all your runs share them. Give it to Snakemake with `--conda-prefix /path/to/folder`, or as `conda-prefix: /path/to/folder` in your profile (see [Running transXpress](#running-transxpress)). Without it, Snakemake builds the environments in `.snakemake/conda` in each run folder.


6. Create a tab-separated file called *samples.txt* in the assembly directory describing where to find your raw read FASTQ files. Create this file with the following contents:
      ~~~
      cond_A      cond_A_rep1 A_rep1_left.fq    A_rep1_right.fq
      cond_A      cond_A_rep2 A_rep2_left.fq    A_rep2_right.fq
      cond_B      cond_B_rep1 B_rep1_left.fq    B_rep1_right.fq
      cond_B      cond_B_rep2 B_rep2_left.fq    B_rep2_right.fq
      ~~~
    

* You can download reads from SRA with provided script:
    ~~~~
    ./sra_download.sh <SRR0000000> <SRR0000001> <...>
    ~~~~
    where SRR0000000 is an SRA readset ID. E.g., SRR3883773

* See the tests directory for an example of a samples file: [test_samples.txt](./tests/test_samples.txt)

    
7. Take a look at the configuration file *config.yaml* and update as required (you should check: **assembler, targetp, signalp, strand_specific, lineage**).


8. Setup other conda environments (This will take a while):
~~~~
snakemake --profile /path/to/transXpress/profiles/local --conda-prefix /path/to/folder --conda-create-envs-only
~~~~
On a cluster, build them in a job and not on the shared login node, for example on Slurm:
~~~~
sbatch --time=4:00:00 --mem=24G --cpus-per-task=4 -o envs.%j.out --wrap 'source /path/to/miniforge/etc/profile.d/conda.sh && conda activate transxpress && snakemake --profile /path/to/transXpress/profiles/local --conda-prefix /path/to/folder --conda-create-envs-only'
~~~~
The environments of the annotation rules that come after the split checkpoints (rfam and pfam with a shared `db_dir`) are built by the main Snakemake process during the first run.
When the environments are in place, delete the package caches, which hold a second copy of every package:
~~~~
conda clean --all --yes
pip cache purge
~~~~

9. Install SignalP 6.0 (fast):
      * Download SignalP 6.0 fast from https://services.healthtech.dtu.dk/services/SignalP-6.0/ (go to Downloads)
      * Unpack and install signalp:
        ~~~~
         conda env create -f envs/signalp.yaml
         conda activate signalp
         tar zxvf signalp-6.0h.fast.tar.gz
         cd signalp6_fast
         pip install signalp-6-package/
         SIGNALP_DIR=$(python -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
         cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/
         conda deactivate
        ~~~~
        (make sure the conda python is used, or use the full path to python from your conda installation)

10. Install TargetP 2.0:
      * Download TargetP 2.0 from https://services.healthtech.dtu.dk/software.php
      * extract the tarball and add the targetp `bin` folder to the PATH of the targetp environment, with an activation script:
        ~~~~
         conda env create -f envs/targetp.yaml
         conda activate targetp
         tar zxvf targetp-2.0.Linux.tar.gz
         mkdir -p "$CONDA_PREFIX/etc/conda/activate.d" "$CONDA_PREFIX/etc/conda/deactivate.d"
         printf 'export _TARGETP_OLD_PATH="$PATH"\nexport PATH="$PATH:%s/targetp-2.0/bin"\n' "$(pwd)" > "$CONDA_PREFIX/etc/conda/activate.d/targetp.sh"
         printf 'export PATH="$_TARGETP_OLD_PATH"\nunset _TARGETP_OLD_PATH\n' > "$CONDA_PREFIX/etc/conda/deactivate.d/targetp.sh"
         conda deactivate
        ~~~~
        (the script holds the full path of the `targetp-2.0` folder. If you move the folder, edit the path in `$CONDA_PREFIX/etc/conda/activate.d/targetp.sh`.)

## Output layout

A run writes into `results/`, with one folder per step. Each step folder holds its answer files at the top and its working files in a subfolder. The rule logs go to `logs/`, and the databases of the default `db_dir` go to `resources/`.

~~~
<run folder>/
  config.yaml  samples.txt
  logs/          <rule>.log, and <rule>/<wildcard>.log for a rule with many jobs
  resources/     db/  busco_downloads/
  results/
    run_info.txt
    qc/            fastqc_before_trim/  fastqc_after_trim/  multiqc_before_trim/  multiqc_after_trim/
                   multiqc_before_trim.txt  multiqc_after_trim.txt
                   WARNING_fastqc_before_trim_overview.txt  FastQC_comparison_after_trim.txt
                   busco/  busco_report.txt
                   transcriptome_stats.txt  transcriptome_exN50.tsv  ExN50_plot.pdf
    trimmomatic/   samples_trimmed.txt  samples_trimmed.yaml
                   reads/         the trimmed reads
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
* With the default `db_dir: "resources/db"`, BUSCO downloads the newest dataset of its lineage, so BUSCO scores can differ between runs months apart. A shared `db_dir` keeps one lineage dataset for every run.
* The pipeline reads and writes the paths in [Output layout](#output-layout). A run folder that holds its files at other paths runs every step again, so start each run in its own folder.
* The assembly depends on the assembler version, the thread count and the processor. rnaSPAdes repeats an assembly at a fixed thread count, and the SPAdes documentation gives no promise across processor types. The same reads on a different CPU generation, or with a different value of `rnaspades_threads`, can give a different set of transcripts. Two runs to be compared need the same version, the same thread count and the same hardware. `results/run_info.txt` records the processor, the host, the memory, the operating system and the transXpress commit, and the logs of `trimmomatic_parallel`, `rnaspades` and `trinity` name the machine each one ran on.
* TransDecoder and Trinity each have a Perl script that writes records in Perl hash order, which changes in each run. The TransDecoder rules, `kallisto` and `trinity_stats` set `PERL_HASH_SEED=0`, so two runs write the same records in the same order. Upstream reports: [TransDecoder/TransDecoder#212](https://github.com/TransDecoder/TransDecoder/issues/212), [trinityrnaseq/trinityrnaseq#1478](https://github.com/trinityrnaseq/trinityrnaseq/issues/1478).
* Trinity does not repeat an assembly by default: with more than one Inchworm thread, two runs on the same reads give slightly different transcripts. Set `trinity_parameters: "--inchworm_cpu 1"` for a repeatable assembly, or use rnaSPAdes. Tested on the test data only; details in [issue 22](https://github.com/titodamiani/transXpress/issues/22).

## Running transXpress

Run the pipeline in a run folder that holds `config.yaml` and the samples file. A profile is a folder in `profiles/` whose `config.yaml` holds the Snakemake options of the place where the jobs run:

| Profile | Where the jobs run |
|---|---|
| `profiles/local` | on this computer, with all its cores |
| `profiles/slurm` | as Slurm jobs, through the Snakemake executor plugin `slurm` |
| `profiles/aurum` | as Slurm jobs on the Aurum cluster of IOCB Prague: `profiles/slurm` with the Aurum settings |

~~~~
conda activate transxpress
cd /path/to/run-folder
snakemake --profile /path/to/transXpress/profiles/local
~~~~
Each profile names the `workflow/Snakefile` of its own copy of transXpress, so the command works in any run folder. Options on the command line come after the profile, for example `--conda-prefix /path/to/folder`. A command line option replaces the profile value of the same option.

For your own cluster, copy `profiles/slurm` to `profiles/<site>` and add your site flags, see [Flags that belong to your cluster](#flags-that-belong-to-your-cluster).

Each rule gives its thread count and two resources: `mem_mb`, the memory in MB, and `runtime`, the time limit. The rules whose size follows the dataset take these from `config.yaml`, for example `rnaspades_memory` (in GB), `rnaspades_threads` and `rnaspades_runtime`. The other rules keep fixed figures in their rule files.

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
The time is the limit for the whole run, so it has to be longer than the longest job: `trinity` asks for 7 days. The main process builds the environments of the annotation rules during the run, and a build needs about 10 GB, so the job asks for 24 GB.

Inside a job, the `slurm` plugin prints a warning and removes the `SLURM_` variables of the main job from its own environment, so the jobs it submits get their own settings.

### Running specific steps

You can run specific steps of the pipeline by specifying a rule or a resulting file:
~~~
# run only the rules that make the FastQC and MultiQC report of the input data, to check its quality
snakemake --profile /path/to/transXpress/profiles/slurm results/qc/multiqc_before_trim.txt
~~~
~~~
# run only the rules to produce the samples_trimmed.txt file
snakemake --profile /path/to/transXpress/profiles/slurm results/trimmomatic/samples_trimmed.txt
~~~

## Annotation databases

transXpress annotates the transcripts with Swiss-Prot, Pfam and Rfam. BUSCO checks the assembly against a lineage dataset. The three databases take about 6.1 GB on disk: Pfam 4.8 GB, Rfam 0.7 GB, Swiss-Prot 0.6 GB. The download is about 0.6 GB.

`config.yaml` fixes the release of each database. Keep the quotes:
~~~
pfam_release: "38.2"
rfam_release: "15.1"
sprot_release: "2026_03"
~~~
The download checks Pfam and Swiss-Prot against the md5 their servers publish. Rfam publishes no checksum. The BUSCO server keeps only the newest dataset of each lineage, so a key cannot fix the BUSCO lineage. Only a kept copy can.

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

## Running tests
~~~~
cd tests
./run_test.sh
~~~~
`run_test.sh` runs the test data in `tests/` with `profiles/local` on 8 cores.

## Align reads to the transcriptome assembly and visualize the results in IGV
If you want to align reads to the transcriptome assembly and visualize the results in IGV, you can use the following commands:
~~~~
snakemake --profile /path/to/transXpress/profiles/slurm align_reads
snakemake --profile /path/to/transXpress/profiles/slurm IGV
~~~~

Then you can load your transcriptome file to IGV: Genomes -> Load Genome from File -> select the file *results/assembly/transcriptome.fasta*

Your sorted .bam files will be in the igv folder:
~~~~
results/expression/bowtie2_rsem/igv/{sample}.sorted.bam
results/expression/bowtie2_rsem/igv/{sample}.sorted.bam.bai
~~~~
Load them to IGV: File -> Load from File -> select the *results/expression/bowtie2_rsem/igv/{sample}.sorted.bam* files

## The directed acyclic graph (DAG) of the transXpress pipeline execution

![The directed acyclic execution graph](dag.svg )

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

An assembly that is to be compared with another one needs the same processor on both runs, so a constraint that names the CPU generation goes in `default-resources`. Note that Slurm records no processor model, so a feature such as `gen-d` is as close as a constraint gets.

`set-resources` applies to one rule, for example a rule that has to land on a named node:
~~~~
set-resources:
  rnaspades:
    slurm_partition: "'bigmem'"
    slurm_extra: "'--exclude=badnode1,badnode2 --nodelist=node42'"
~~~~
A rule value replaces the default value, so repeat the flags of `default-resources` that the rule still needs.

A node name and a partition name usually travel together. The partition a cluster picks for the memory of a job does not have to be one that holds the node you named, and sbatch then answers `Requested node configuration is not available`.

### Working files on the node's local disk
On a cluster, the run folder usually sits on network storage that all users share, and many clusters ask that work with many or large files runs on the node's own disk.

The `scratch` key in `config.yaml` names the environment variable that holds a folder on the node's local disk in a cluster job, with its `$`:
~~~~
scratch: "$SCRATCH"
~~~~

Ask your cluster documentation which variable that is. On some clusters `$TMPDIR` is in memory, and a file there counts against the memory of the job. The default `""` runs every rule in the run folder.

With the key set, every job gets that folder as its temporary folder: Snakemake sets `TMPDIR`, `TMP` and `TEMP` of the job to it. `signalp_parallel`, `targetp_parallel` and `tmhmm_parallel` write their work files there. Five rules make their own folder under it, work there, and copy their files to the run folder when they end:

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
