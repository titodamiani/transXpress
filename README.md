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

## Notes on results

* `download_sprot` downloads the current Swiss-Prot release, so BLAST hits can differ between runs.
* In a run folder made by an older transXpress version, add `--rerun-triggers mtime` to the Snakemake command to keep the old results. The default triggers also look at code and parameter changes and would run the rules again.

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

### Running specific steps

You can run specific steps of the pipeline by specifying a rule or a resulting file:
~~~
# run only the rules until the multiqc_before_trim rule to check quality of the input data
./transXpress.sh multiqc_before_trim
~~~
~~~
# run only the rules to produce samples_trimmed.txt file
./transXpress.sh samples_trimmed.txt
~~~

or when using profiles:
~~~
# run only the rules until the multiqc_before_trim rule to check quality of the input data
snakemake multiqc_before_trim --profile profiles/slurm "$@"
~~~
~~~
# run only the rules to produce samples_trimmed.txt file
snakemake samples_trimmed.txt --profile profiles/slurm "$@"
~~~

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

Then you can load your transcriptome file to IGV: Genomes -> Load Genome from File -> select the file *transcriptome.fasta*

Your sorted .bam files will be in the bowtie_alignments folder: 
~~~~
bowtie_alignments/{sample}.sorted.bam
bowtie_alignments/{sample}.sorted.bam.bai
~~~~
Load them to IGV: File -> Load from File -> select the *bowtie_alignments/{sample}.sorted.bam* files

## The directed acyclic graph (DAG) of the transXpress pipeline execution

![The directed acyclic execution graph](dag.svg )

## Possible problems when executing on cluster systems

### Time limit
Slurm jobs need a time limit. The time of each rule is set in `profiles/slurm/config.yaml`: `default-resources` gives 4 hours, and `set-resources` gives longer times to the rules that need them. `transXpress.sh` uses this profile on Slurm clusters.

If a job is cancelled because of the time limit, raise the time of that rule in `set-resources`. The time you need depends on the size of the reads used for the assembly.

See https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Computing-Requirements

### Pipeline hangs when cluster cancels the job
It is possible that cluster cancels the job, but pipeline seems to be still running. This can happen because the pipeline does not receive information whether cluster job completed successfully, failed or is still running. You can add `--cluster-status` option and add script which detects the job status. 

See https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#using-cluster-status

Alternatively, you can use snakemake [profiles](https://github.com/Snakemake-Profiles/doc) which also contain status checking script. 

See https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles 