@functools.lru_cache(maxsize=None)
def check_samples(path):
  """
  Checks the samples file before any job starts. Every line that is not empty has 3 columns (sample,
  replicate and one reads file) or 4 columns (a second reads file), separated by white space. A replicate
  name appears once, and every reads file exists.
  """
  if not os.path.isfile(path):
    raise WorkflowError("The samples file %s does not exist. samples_file in config.yaml names it." % path)
  problems = []
  replicates = {}
  with open(path) as handle:
    for number, line in enumerate(handle, 1):
      fields = line.split()
      if not fields:
        continue
      if len(fields) not in (3, 4):
        problems.append("line %d has %d columns, and it needs 3 or 4" % (number, len(fields)))
        continue
      if fields[1] in replicates:
        problems.append("line %d repeats the replicate name %s of line %d" % (number, fields[1], replicates[fields[1]]))
      replicates.setdefault(fields[1], number)
      for reads in fields[2:]:
        if not os.path.isfile(reads):
          problems.append("line %d: the reads file %s does not exist" % (number, reads))
  if problems:
    raise WorkflowError("The samples file %s has problems:\n  %s" % (path, "\n  ".join(problems)))
  return path


def checked_samples_file(wildcards):
  return check_samples(config["samples_file"])


checkpoint trimmomatic_split:
  """
  Splits file with information about all reads files into separate files
  so they can be processed with trimmomatic in parallel.
  """   
  input:
    samples=checked_samples_file
  output:
    directory("results/trimmomatic/reads")
  log:
    "logs/trimmomatic_split.log"
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
    "results/trimmomatic/reads/sample_{job_index}"
  output:
    "results/trimmomatic/reads/completed_{job_index}"
  log:
    "logs/trimmomatic_parallel/{job_index}.log"
  conda:
    "../envs/trimmomatic.yaml"
  params:
    memory="4",
    trimmomatic=trimmomatic_parameters
  threads:
    8
  shell:
    """
    # record the computer this ran on, because the assembly depends on it
    ( hostname ; grep -m1 "model name" /proc/cpuinfo 2>/dev/null || uname -m ; nproc 2>/dev/null || true ) > {log}
    # note: trimmomatic can use gzipped files directly
    read SAMPLE REPLICATE F_READS R_READS < {input}
    # If the sample line is empty, ignore it
    if [ -z "$REPLICATE" ]; then
      touch {output}
      exit 0
    fi
    if [ ! -z "$R_READS" ]; then
      trimmomatic PE -threads {threads} $F_READS $R_READS results/trimmomatic/reads/{wildcards[job_index]}.R1-P.qtrim.fastq.gz results/trimmomatic/reads/{wildcards[job_index]}.R1-U.qtrim.fastq.gz results/trimmomatic/reads/{wildcards[job_index]}.R2-P.qtrim.fastq.gz results/trimmomatic/reads/{wildcards[job_index]}.R2-U.qtrim.fastq.gz {params.trimmomatic} &>> {log}
      echo $SAMPLE	$REPLICATE	results/trimmomatic/reads/{wildcards[job_index]}.R1-P.qtrim.fastq.gz	results/trimmomatic/reads/{wildcards[job_index]}.R2-P.qtrim.fastq.gz > {output} 2>> {log}
    else
      trimmomatic SE -threads {threads} $F_READS results/trimmomatic/reads/{wildcards[job_index]}.U.qtrim.fastq.gz {params.trimmomatic} &>> {log}
      echo $SAMPLE      $REPLICATE      results/trimmomatic/reads/{wildcards[job_index]}.U.qtrim.fastq.gz > {output} 2>> {log}
    fi
    """


def trimmomatic_completed_parallel_jobs(wildcards):
  """
  Returns names of files with information about files processed with trimmomatic.
  """
  parallel_dir = checkpoints.trimmomatic_split.get(**wildcards).output[0]
  # sorted, because glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "sample_{job_index}")).job_index)
  completed_ids = expand(os.path.join(parallel_dir,"completed_{job_index}"), job_index=job_ids)
  return completed_ids


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
  run:
    import os
    import os.path
    import pprint

    # rnaSPAdes reads a relative path in this file from the folder of the file, so each read path is
    # written relative to that folder
    yaml_dir = os.path.dirname(output["samples_yaml"])
    sample_list = []
    with open(input["samples_txt"], "r") as input_handle:
      for line in input_handle:
        row = re.split("[\t ]+", line)
        if (len(row) > 3): # paired reads
          paired_dict = {}
          paired_dict['orientation'] = 'fr'
          paired_dict['type'] = 'paired-end'
          f_reads = row[2].strip()
          r_reads = row[3].strip()
          assert os.path.isfile(f_reads)
          assert os.path.isfile(r_reads)
          paired_dict['left reads'] = [os.path.relpath(f_reads, yaml_dir)]
          paired_dict['right reads'] = [os.path.relpath(r_reads, yaml_dir)]
          ##Depends on where the unpaired reads are written to
          ##Presumably could get from sample txt as well
          ##Note: Commented out, to make it consistent with Trinity. May change in future
          ##unpaired_reads_f = f[:-16]+"U.qtrim.fastq.gz"
          ##unpaired_reads_r = r[:-16]+"U.qtrim.fastq.gz"
          ##assert os.path.isfile(unpaired_reads_f)
          ##assert os.path.isfile(unpaired_reads_r)
          ##unpaired_reads = [unpaired_reads_f] + [unpaired_reads_r]
          ##if len(unpaired_reads) > 0:
          ##    sample_dict['single reads'] = unpaired_reads
          sample_list.append(paired_dict)

        if (len(row) == 3): # unpaired reads
          unpaired_dict = {}
          unpaired_dict['type'] = 'single'
          u_reads = row[2].strip()
          assert os.path.isfile(u_reads)
          unpaired_dict['single reads'] = [os.path.relpath(u_reads, yaml_dir)]
          sample_list.append(unpaired_dict)

    with open(output["samples_yaml"], "w") as output_handle:
      output_handle.write(pprint.pformat(sample_list))
