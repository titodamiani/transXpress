checkpoint trimmomatic_split:
  """
  Splits the samples file into one file for each line.
  """
  input:
    samples=checked_samples_file
  output:
    directory("results/trimmomatic/trimmed_reads")
  log:
    "logs/trimmomatic_split.log"
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  shell:
    """
    mkdir -p {output} &> {log}
    split -d -l 1 {input} {output}/sample_ &>> {log}
    """


rule trimmomatic_parallel:
  """
  Trimmomatic on the reads of one line of the samples file.
  """
  input:
    "results/trimmomatic/trimmed_reads/sample_{job_index}"
  output:
    "results/trimmomatic/trimmed_reads/completed_{job_index}"
  log:
    "logs/trimmomatic_parallel/{job_index}.log"
  conda:
    "../envs/trimmomatic.yaml"
  params:
    # a function: in a plain string Snakemake takes {adapters} for a wildcard
    trimmomatic=lambda wildcards: str(config["trimmomatic_parameters"])
  threads:
    8
  resources:
    # Java heap (-Xmx in the shell) is mem_mb - 1024; 16 threads need about 2 GB of heap
    mem_mb=lambda wildcards, threads: 2048 + 256 * threads,
    runtime="4h"
  shell:
    """
    # record the computer this ran on
    source {workflow.basedir}/scripts/tools.sh
    record_host > {log}
    # read returns 1 for a last line without a newline, which split copies as it is. Snakemake made the input.
    read SAMPLE REPLICATE F_READS R_READS < {input} || true
    if [ -z "$REPLICATE" ]; then
      touch {output}
      exit 0
    fi
    # {{adapters}} stands for resources/adapters.fasta; the shell fills in the path, so a run can go on from another clone
    ADAPTERS=$(cd {workflow.basedir}/../resources && pwd)/adapters.fasta
    TRIMMOMATIC_PARAMETERS='{params.trimmomatic}'
    TRIMMOMATIC_PARAMETERS=${{TRIMMOMATIC_PARAMETERS//'{{adapters}}'/$ADAPTERS}}
    OUT=results/trimmomatic/trimmed_reads/{wildcards.job_index}
    # the bioconda wrapper sets a 1 GB heap; -Xmx replaces it
    HEAP=-Xmx$(( {resources.mem_mb} - 1024 ))m
    if [ -n "$R_READS" ]; then
      PARTS="R1-P R1-U R2-P R2-U"
      trimmomatic $HEAP PE -threads {threads} $F_READS $R_READS $OUT.tmp.R1-P.qtrim.fastq.gz $OUT.tmp.R1-U.qtrim.fastq.gz $OUT.tmp.R2-P.qtrim.fastq.gz $OUT.tmp.R2-U.qtrim.fastq.gz $TRIMMOMATIC_PARAMETERS &>> {log}
    else
      PARTS="U"
      trimmomatic $HEAP SE -threads {threads} $F_READS $OUT.tmp.U.qtrim.fastq.gz $TRIMMOMATIC_PARAMETERS &>> {log}
    fi
    # Trimmomatic can report success for a broken gzip file, so each file gets its name only after gzip -t
    for PART in $PARTS; do
      gzip -t $OUT.tmp.$PART.qtrim.fastq.gz &>> {log}
    done
    for PART in $PARTS; do
      mv $OUT.tmp.$PART.qtrim.fastq.gz $OUT.$PART.qtrim.fastq.gz &>> {log}
    done
    if [ -n "$R_READS" ]; then
      echo $SAMPLE $REPLICATE $OUT.R1-P.qtrim.fastq.gz $OUT.R2-P.qtrim.fastq.gz > {output} 2>> {log}
    else
      echo $SAMPLE $REPLICATE $OUT.U.qtrim.fastq.gz > {output} 2>> {log}
    fi
    """


rule trimmomatic_merge:
  """
  Lists the trimmed reads of all lines of the samples file.
  """
  input:
    trimmomatic_completed_parallel_jobs
  output:
    samples_trimmed="results/trimmomatic/samples_trimmed.txt"
  log:
    "logs/trimmomatic_merge.log"
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  shell:
    """
    cat {input} > {output.samples_trimmed} 2> {log}
    """


rule samples_yaml_conversion:
  """
  Converts samples_trimmed.txt to the YAML dataset file of rnaSPAdes.
  """
  input:
    samples_txt="results/trimmomatic/samples_trimmed.txt"
  output:
    samples_yaml="results/trimmomatic/samples_trimmed.yaml"
  log:
    "logs/samples_yaml_conversion.log"
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  script:
    "../scripts/samples_yaml.py"
