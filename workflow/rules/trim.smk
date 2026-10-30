checkpoint trimmomatic_split:
  """
  Splits file with information about all reads files into separate files
  so they can be processed with trimmomatic in parallel.
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
  Processes individual read files with trimmomatic.
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
    # a function, because Snakemake takes {adapters} in a plain string for a wildcard
    trimmomatic=lambda wildcards: str(config["trimmomatic_parameters"])
  threads:
    8
  resources:
    mem_mb=4096,
    runtime="4h"
  shell:
    """
    # record the computer this ran on
    source {workflow.basedir}/scripts/tools.sh
    record_host > {log}
    # read returns 1 for a last line without a newline, which split copies as it is. Snakemake made the input.
    read SAMPLE REPLICATE F_READS R_READS < {input} || true
    # an empty line of the samples file gives an empty marker
    if [ -z "$REPLICATE" ]; then
      touch {output}
      exit 0
    fi
    # {{adapters}} in the trimmomatic parameters stands for the adapter file that ships with transXpress. The
    # shell fills in its path, so the path of the clone is not a param and a run can go on from another clone.
    ADAPTERS=$(cd {workflow.basedir}/../resources && pwd)/adapters.fasta
    TRIMMOMATIC_PARAMETERS='{params.trimmomatic}'
    TRIMMOMATIC_PARAMETERS=${{TRIMMOMATIC_PARAMETERS//'{{adapters}}'/$ADAPTERS}}
    # the output lines hold the fields separated by one space
    OUT=results/trimmomatic/trimmed_reads/{wildcards.job_index}
    if [ -n "$R_READS" ]; then
      trimmomatic PE -threads {threads} $F_READS $R_READS $OUT.R1-P.qtrim.fastq.gz $OUT.R1-U.qtrim.fastq.gz $OUT.R2-P.qtrim.fastq.gz $OUT.R2-U.qtrim.fastq.gz $TRIMMOMATIC_PARAMETERS &>> {log}
      echo $SAMPLE $REPLICATE $OUT.R1-P.qtrim.fastq.gz $OUT.R2-P.qtrim.fastq.gz > {output} 2>> {log}
    else
      trimmomatic SE -threads {threads} $F_READS $OUT.U.qtrim.fastq.gz $TRIMMOMATIC_PARAMETERS &>> {log}
      echo $SAMPLE $REPLICATE $OUT.U.qtrim.fastq.gz > {output} 2>> {log}
    fi
    """


rule trimmomatic_merge:
  """
  Creates a file with information about all reads files processed with trimmomatic.
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
  Converts the file with information about all trimmed reads files into a yaml
  format used by rnaSPAdes.
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
