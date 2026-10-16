![transXpress](logo/Transxpress_Logo_RGB.png)

transXpress: a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for streamlined de novo transcriptome assembly and annotation

## Intro

## Dependencies

transXpress requires:
* snakemake 5.4.2+ (install via conda, `envs/default.yaml`)
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
* python, numpy, pip (install via conda, `envs/default.yaml`)
* busco 4+ (install via conda, `envs/busco.yaml`)
* rsem (install via conda, `envs/trinity_utils.yaml`)
* [SignalP 6.0](https://services.healthtech.dtu.dk/service.php?SignalP)
* [TargetP 2.0](https://services.healthtech.dtu.dk/service.php?TargetP-2.0)
* tmhmm.py (install via pip, `envs/default.yaml`)
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
All commands below use `conda`. conda 23.10 and later use the libmamba solver by default, so the solve speed is the same as with mamba. Snakemake runs with `--conda-frontend conda`, because mamba is not on the PATH inside cluster jobs.

3. To ensure correct versions of R packages will be used unset R_LIBS_SITE
~~~~
unset R_LIBS_SITE
~~~~

4. Setup main transXpress conda environment:
~~~~
conda config --set channel_priority strict
PIP_NO_BUILD_ISOLATION=0 conda env create --name transxpress --file envs/default.yaml
conda activate transxpress
~~~~
* `PIP_NO_BUILD_ISOLATION=0` lets pip build tmhmm.py with the numpy and Cython of the environment.
* `envs/default.yaml` pins Snakemake 7.32.4, because Snakemake 8 and later have no `--cluster` option.

5. Choose a folder for the conda environments of the pipeline rules and set it in `~/.bashrc`:
~~~~
export TRANSXPRESS_CONDA_PREFIX=/path/to/folder
~~~~
`TRANSXPRESS_CONDA_PREFIX` is a transXpress variable, not a Snakemake one. `transXpress.sh` passes it to the Snakemake option `--conda-prefix`. When you call `snakemake` directly, add `--conda-prefix "$TRANSXPRESS_CONDA_PREFIX"`. If the variable is not set, Snakemake builds the environments in `.snakemake/conda` in the run folder.


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
snakemake --use-conda --conda-frontend conda --conda-prefix "$TRANSXPRESS_CONDA_PREFIX" --conda-create-envs-only --cores 1
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
| `delete_intermediates` | the trimmed reads, the kallisto, bowtie2 and RSEM indexes, and the Trinity working files | every result. A later `snakemake` run in the folder has nothing to do |

~~~
./transXpress.sh delete_intermediates
~~~

## Notes on results

* Swiss-Prot, Pfam and Rfam have fixed releases, set in `config.yaml`. See [Annotation databases](#annotation-databases).
* With the default `db_dir: "resources/db"`, BUSCO downloads the newest dataset of its lineage, so BUSCO scores can differ between runs months apart. A shared `db_dir` keeps one lineage dataset for every run.
* The pipeline reads and writes the paths in [Output layout](#output-layout). A run folder that holds its files at other paths runs every step again, so start each run in its own folder.
* The assembly depends on the assembler version, the thread count and the processor. rnaSPAdes repeats an assembly at a fixed thread count, and the SPAdes documentation gives no promise across processor types. The same reads on a different CPU generation, or with a different value of `rnaspades_threads`, can give a different set of transcripts. Two runs to be compared need the same version, the same thread count and the same hardware. `results/run_info.txt` records the processor, the host, the memory, the operating system and the transXpress commit, and the logs of `trimmomatic_parallel`, `rnaspades`, `trinity_inchworm_chrysalis` and `trinity_final` name the machine each one ran on.
* TransDecoder and Trinity each have a Perl script that writes records in Perl hash order, which changes in each run. The TransDecoder rules, `kallisto` and `trinity_stats` set `PERL_HASH_SEED=0`, so two runs write the same records in the same order. Upstream reports: [TransDecoder/TransDecoder#212](https://github.com/TransDecoder/TransDecoder/issues/212), [trinityrnaseq/trinityrnaseq#1478](https://github.com/trinityrnaseq/trinityrnaseq/issues/1478).
* Trinity does not repeat an assembly, by its own documentation. Its Butterfly stage runs in parallel chunks whose results are collected in the order they finish. Use rnaSPAdes for an assembly that has to be repeatable.

## Running transXpress

There are several options on how to run transXpress

Option 1 - use the provided script:
~~~~
./transXpress.sh
~~~~

Option 2 - run snakemake manually with 10 local threads:
~~~~
snakemake --conda-frontend conda --use-conda --cores 10
~~~~

Option 3 - run snakemake manually on an LSF cluster:
~~~~
snakemake --conda-frontend conda --use-conda --latency-wait 60 --jobs 10000 --cluster 'bsub -n {threads} -R "rusage[mem={params.memory}000] span[hosts=1]" -oo {log}.bsub'
~~~~

Option 4 - define a profile and run snakemake with the profile:
~~~~
snakemake --profile profiles/slurm "$@"
~~~~
The folder `profiles/slurm` holds a profile `config.yaml` for Slurm. You can find more profiles [here](https://github.com/snakemake-profiles/doc)

### Running the main Snakemake process on a Slurm cluster

The main Snakemake process needs the `transxpress` environment. Run it in its own job and activate the environment first. Replace the two paths with yours:
~~~
sbatch --time=5-00:00:00 --mem=4G -o main.%j.out --wrap 'source /path/to/miniforge/etc/profile.d/conda.sh && conda activate /path/to/transxpress-env && bash /path/to/transXpress/transXpress.sh'
~~~
The time is the limit for the whole run.

### Running specific steps

You can run specific steps of the pipeline by specifying a rule or a resulting file:
~~~
# run only the rules until the multiqc_before_trim rule to check quality of the input data
./transXpress.sh multiqc_before_trim
~~~
~~~
# run only the rules to produce the samples_trimmed.txt file
./transXpress.sh results/trimmomatic/samples_trimmed.txt
~~~

or when using profiles:
~~~
# run only the rules until the multiqc_before_trim rule to check quality of the input data
snakemake multiqc_before_trim --profile profiles/slurm "$@"
~~~
~~~
# run only the rules to produce the samples_trimmed.txt file
snakemake results/trimmomatic/samples_trimmed.txt --profile profiles/slurm "$@"
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
   /path/to/transXpress/transXpress.sh --snakefile /path/to/transXpress/Snakefile download_databases
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

## Align reads to the transcriptome assembly and visualize the results in IGV
If you want to align reads to the transcriptome assembly and visualize the results in IGV, you can use the following commands:
~~~~
./transXpress.sh align_reads
./transXpress.sh IGV
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
Slurm jobs need a time limit. The time of each rule is set in `profiles/slurm/config.yaml`: `default-resources` gives 4 hours, and `set-resources` gives longer times to the rules that need them. `transXpress.sh` uses this profile on Slurm clusters.

If a job is cancelled because of the time limit, raise the time of that rule in `set-resources`. The time you need depends on the size of the reads used for the assembly.

See https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Computing-Requirements

### How the pipeline learns that a cluster job stopped
Snakemake asks the cluster for the state of each job it submitted. Without that, Snakemake waits for the job to write a marker file, and a job stopped by the node, by the kernel or by `scancel` writes no marker. Snakemake then waits for that job for ever, with an empty queue and no error message.

On Slurm, `profiles/slurm/config.yaml` sets `cluster-status: slurm-status.sh`. The script maps the Slurm job state to the three words Snakemake expects. It answers most checks from a `squeue` listing it caches in `.snakemake/slurm-status-cache` for 60 seconds, and asks `sacct` only about a job that has left the queue, which keeps the load on the cluster to about one query per minute.

A state the script does not recognise counts as a failure, so an unknown state stops the job instead of holding up the run. `restart-times` in the same profile then gives the job another try. Five rules set `retries: 0` and get no second try: `rnaspades`, `trinity_inchworm_chrysalis`, `trinity_final`, `kallisto` and `align_reads`. They fail mostly on memory or time, and a retry asks for the same figures and starts from the beginning. When Slurm does not answer at all, the script reports the job as still running, so a cluster outage does not end the run.

For the job state names, see https://slurm.schedmd.com/sacct.html

On LSF and PBS, `transXpress.sh` passes `--cluster` without a status command, so a job that the cluster stops can still leave Snakemake waiting. Write a status script for your scheduler and pass it with `--cluster-status`.

### Flags that belong to your cluster
Some sbatch flags are site specific: a node to avoid, a processor generation to insist on, an account to charge. `profiles/slurm/config.yaml` takes them in two places, and both are empty by default.

`TRANSXPRESS_SBATCH_EXTRA` applies to every job of a run. Export it before you start the pipeline:
~~~~
export TRANSXPRESS_SBATCH_EXTRA="--constraint=gen-d --exclude=badnode1"
./transXpress.sh
~~~~

An assembly that is to be compared with another one needs the same processor on both runs, so a constraint that names the CPU generation goes here. Note that Slurm records no processor model, so a feature such as `gen-d` is as close as a constraint gets.

The `extra` resource applies to one rule. Put it in `set-resources` in the same profile, for a rule that has to land on a named node:
~~~~
set-resources:
  - rnaspades:extra=--nodelist=node42 --partition=bigmem
~~~~

Write that value without quotes. A quoted value keeps its quotes, so sbatch reads two flags as one and rejects them with `Invalid specification`. A rule can hold both `extra` and a time limit, as two entries that name the same rule.

A node name and a partition name usually travel together. The partition a cluster picks for the memory of a job does not have to be one that holds the node you named, and sbatch then answers `Requested node configuration is not available`.

Each job needs a time limit, so `--time` stays on the sbatch line and does not belong in either hook.

### Working files on the node's local disk
Four rules write many large or short-lived working files: `rnaspades`, `trinity_butterfly_parallel`, `kallisto` and `align_reads`. On a cluster, the run folder usually sits on network storage that all users share, and many clusters ask that such work runs on the node's own disk.

The `scratch` key in `config.yaml` names the environment variable that holds a folder on the node's local disk in a cluster job, with its `$`:
~~~~
scratch: "$SCRATCH"
~~~~

Ask your cluster documentation which variable that is. On some clusters `$TMPDIR` is in memory, and a file there counts against the memory of the job. The default `""` runs every rule in the run folder.

With the key set, each of the four rules makes its own folder under that variable, works there, and copies its files to the run folder when it ends:

| Rule | Works on the local disk | Copied to the run folder |
|---|---|---|
| `rnaspades` | the whole rnaSPAdes output folder | `results/assembly/rnaspades`, which replaces the folder in the run folder |
| `trinity_butterfly_parallel` | the files of each small assembly, through Trinity's `--workdir` | nothing: Trinity writes the fasta of each assembly to `results/assembly/trinity/read_partitions` itself |
| `kallisto` | the index and the per replicate results | every file, over the files in `results/expression/kallisto` |
| `align_reads` | the indexes and the per replicate results | every file, over the files in `results/expression/bowtie2_rsem` |

When the variable is not set in a job, for example on a computer without a cluster, the rule works in the run folder.

When a rule fails, for example on an error of the tool or when the kernel stops it for lack of memory, it copies its local files to the place the rule writes them without the local disk, so they can be read to find the cause. A job that Slurm stops at its time limit or with `scancel` copies nothing: Slurm ends the Snakemake process of the job at once, and with it every process the job started. A node that crashes leaves nothing either. In all these cases the rule log `logs/<rule>.log` holds the output of the tool up to the stop, and the Slurm log `logs/<rule>.log.slurm` names a time limit.

The local disk of a node is shared by the jobs on it, and a cluster may not reserve disk space for a job. `rnaspades` needed about 120 GB of working files for 9 libraries of about 23 million read pairs each.

`rnaspades` runs without `--checkpoints`. On the local disk, the folder is deleted at the end of the job, so a restart point would not survive.
