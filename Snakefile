import os
import sys
import shutil
import re
import csv
import Bio.SeqIO
import subprocess
import tempfile
from snakemake.utils import min_version
from Bio.Seq import Seq
import pandas as pd

min_version("5.4.1")

configfile: "config.yaml"

# The releases of the annotation databases. The download rules fetch exactly these releases.
PFAM_RELEASE = str(config.get("pfam_release", "38.2"))
RFAM_RELEASE = str(config.get("rfam_release", "15.1"))
SPROT_RELEASE = str(config.get("sprot_release", "2026_03"))
if not re.fullmatch(r"\d{4}_\d{2}", SPROT_RELEASE):
  raise ValueError("sprot_release must look like \"2026_03\", in quotes. The config gives %s" % SPROT_RELEASE)

# The ORF set that transcriptome.pep holds and the downstream rules read. Both sets are always saved.
TRANSDECODER_ORFS = str(config.get("transdecoder_orfs", "all"))
if TRANSDECODER_ORFS not in ["all", "single_best"]:
  raise ValueError("transdecoder_orfs must be \"all\" or \"single_best\". The config gives %s" % TRANSDECODER_ORFS)
TRANSDECODER_ORFS_FILE = {"all": "results/transdecoder/transcriptome_all_orfs.pep", "single_best": "results/transdecoder/transcriptome_single_best_orf.pep"}[TRANSDECODER_ORFS]

# The rnaSPAdes transcript set that becomes transcriptome.fasta. rnaSPAdes always writes all three.
RNASPADES_FILTER = str(config.get("rnaspades_filter", "default"))
if RNASPADES_FILTER not in ["default", "hard", "soft"]:
  raise ValueError("rnaspades_filter must be \"default\", \"hard\" or \"soft\". The config gives %s" % RNASPADES_FILTER)

# The organism group that SignalP predicts signal peptides for.
SIGNALP_ORGANISM = str(config.get("signalp", "eukarya"))
if SIGNALP_ORGANISM not in ["eukarya", "other"]:
  raise ValueError("signalp must be \"eukarya\" or \"other\". The config gives %s" % SIGNALP_ORGANISM)

# The number of sequences in each annotation chunk. Each chunk is one cluster job for each tool that reads it.
ANNOTATION_CHUNK = {}
for chunk_split, chunk_default in [("fasta", 1000), ("orfs", 5000), ("pep", 5000)]:
  chunk_value = config.get("annotation_chunk_" + chunk_split, chunk_default)
  if isinstance(chunk_value, bool) or not isinstance(chunk_value, int) or chunk_value < 1:
    raise ValueError("annotation_chunk_%s must be a whole number above 0. The config gives %s" % (chunk_split, chunk_value))
  ANNOTATION_CHUNK[chunk_split] = chunk_value

# The environment variable that holds a node-local folder in a cluster job, for example "$SCRATCH".
# "" runs every rule in the run folder. The heavy rules read the variable with bash, inside the job,
# because only the job knows its value.
# A key with no value reaches here as None, which means "" too
SCRATCH = str(config.get("scratch") or "")
if SCRATCH and not re.fullmatch(r"\$[A-Za-z_][A-Za-z0-9_]*", SCRATCH):
  raise ValueError("scratch must be \"\" or one environment variable with its $, for example \"$SCRATCH\". The config gives %s" % SCRATCH)
SCRATCH_VAR = SCRATCH[1:]

def job_tmpdir(wildcards):
  """
  The temporary folder of a job: the folder in the variable that the scratch key names, or the
  system default when the key is "" or the variable is not set. Snakemake calls this inside the job.
  """
  return (os.environ.get(SCRATCH_VAR, "") if SCRATCH_VAR else "") or tempfile.gettempdir()

# The folder of the annotation databases. "resources/db" keeps one copy in each run folder, which the
# pipeline downloads. Any other folder is a shared copy that the pipeline only reads.
DB_DIR = str(config.get("db_dir", "resources/db")).rstrip("/")
DB_SHARED = DB_DIR != "resources/db"
LINEAGE = config.get("lineage", "")

# The annotation rules run only after the fasta_split checkpoints, hours into a run, so a
# shared folder is checked here, before any job starts.
if DB_SHARED:
  if not LINEAGE:
    raise ValueError("A shared db_dir needs a named BUSCO lineage, for example lineage: \"viridiplantae_odb12\"")
  db_missing = [path for path in [
      DB_DIR + "/versions.txt",
      DB_DIR + "/uniprot_sprot.fasta",
      DB_DIR + "/Pfam-A.hmm",
      DB_DIR + "/Rfam.cm",
      DB_DIR + "/busco_downloads/lineages/" + LINEAGE]
    if not os.path.exists(path)]
  if db_missing:
    raise ValueError("db_dir %s lacks:\n  %s\nFill the folder with the download_databases target, see README.md."
      % (DB_DIR, "\n  ".join(db_missing)))

# for executable in ["samtools", "bowtie2", "kallisto", "signalp6", "targetp", "Trinity", "blastx", "blastp", "makeblastdb", "cmscan", "hmmscan", "fastqc", "rnaspades.py", "seqkit", "R"]:
#    if not shutil.which(executable):
#        sys.stderr.write("Warning: Cannot find %s in your PATH\n" % (executable))

#TRINITY_EXECUTABLE_PATH=shutil.which("Trinity")
#TRINITY_HOME=os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))) ##Have to resolve the symbolic link that conda makes

# https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinity-Assembly-And-Analysis

# These rules don't need to be sent to a cluster
localrules: all, clean, delete_intermediates, run_info, trimmomatic_split, trimmomatic_merge, samples_yaml_conversion, trinity_butterfly_split, transcriptome_copy, transdecoder_select, compare_qc_after_trim, fastqc_before_trim_warnings, IGV

# A rule with a wildcard keeps its logs in logs/<rule>/, and the wildcard is only in the file name.
# The cluster writes its own log next to the rule log before the job starts, and not every cluster
# makes a missing folder, so every log folder is made here.
onstart:
  for workflow_rule in workflow.rules:
    for log_file in workflow_rule.log:
      log_folder = os.path.dirname(str(log_file))
      if log_folder:
        os.makedirs(log_folder, exist_ok=True)


rule all:
  """
  List of target files of the transxpress pipeline.
  """
  input:
    "results/trimmomatic/samples_trimmed.txt",
    "results/assembly/transcriptome.fasta",
    "results/transdecoder/transcriptome.pep",
    expand("results/transdecoder/transcriptome_all_orfs.{ext}", ext=["pep", "cds", "gff3", "bed"]),
    expand("results/transdecoder/transcriptome_single_best_orf.{ext}", ext=["pep", "cds", "gff3", "bed"]),
    "results/qc/transcriptome_stats.txt",
    "results/qc/ExN50_plot.pdf",
    "results/annotation/transcriptome_annotated.fasta",
    "results/annotation/transcriptome_annotated.pep",
    "results/annotation/transcriptome_TPM_blast.csv",
    "results/qc/busco",
    "results/qc/busco_report.txt",
    "results/qc/fastqc_before_trim",
    "results/qc/multiqc_before_trim",
    "results/qc/fastqc_after_trim",
    "results/qc/multiqc_after_trim",
    "results/qc/multiqc_before_trim.txt",
    "results/qc/multiqc_after_trim.txt",
    "results/qc/WARNING_fastqc_before_trim_overview.txt",
    "results/qc/FastQC_comparison_after_trim.txt",
    "results/expression/edger",
    "results/run_info.txt"

rule clean:
  """
  Not run by default. Deletes the results and the logs of the run, so the next run starts from the
  beginning. Keeps config.yaml, the samples file, resources/ and .snakemake/.
  """
  shell:
    """
    rm -rf results logs
    """


rule delete_intermediates:
  """
  Not run by default. Deletes the large working files that no rule names as its output, so every
  result stays and a later snakemake run in this folder has nothing to do: the trimmed reads, the
  kallisto, bowtie2 and RSEM indexes, and the Trinity working folder apart from the files that the
  Trinity rules name.
  """
  shell:
    """
    rm -f results/trimmomatic/reads/*.fastq.gz
    rm -f results/expression/kallisto/transcriptome.fasta.kallisto_idx
    # the bowtie2.ok marker goes with the index, so a later align_reads run builds the index again
    rm -f results/expression/bowtie2_rsem/transcriptome.fasta.bowtie2.* results/expression/bowtie2_rsem/transcriptome.fasta.RSEM.*
    if [ -d results/assembly/trinity ]; then
      find results/assembly/trinity -mindepth 1 -maxdepth 1 ! -name recursive_trinity.cmds ! -name recursive_trinity.cmds.completed ! -name parallel_jobs -exec rm -rf {{}} +
    fi
    """


rule run_info:
  """
  Records the computer the pipeline ran on and the transXpress commit.
  The assembly depends on the hardware, the assembler version and the thread
  count, so a later run can be compared against this one.
  """
  output:
    "results/run_info.txt"
  params:
    memory="2"
  threads:
    1
  run:
    import platform
    import socket
    import subprocess

    def field(path, key):
      try:
        with open(path) as handle:
          for line in handle:
            if line.startswith(key):
              return line.split(":", 1)[-1].strip()
      except OSError:
        pass
      return ""

    cpu = field("/proc/cpuinfo", "model name") or platform.processor() or platform.machine()
    memory_total = field("/proc/meminfo", "MemTotal") or "unknown"

    os_name = ""
    try:
      with open("/etc/os-release") as handle:
        for line in handle:
          if line.startswith("PRETTY_NAME="):
            os_name = line.split("=", 1)[1].strip().strip('"')
    except OSError:
      pass
    if not os_name:
      os_name = platform.platform()

    try:
      commit = subprocess.check_output(
        ["git", "-C", workflow.basedir, "rev-parse", "HEAD"],
        stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
      commit = "unknown"

    with open(output[0], "w") as out:
      out.write("host\t%s\n" % socket.gethostname())
      out.write("cpu\t%s\n" % cpu)
      out.write("cores\t%s\n" % os.cpu_count())
      out.write("memory_total\t%s\n" % memory_total)
      out.write("os\t%s\n" % os_name)
      out.write("kernel\t%s %s %s\n" % (platform.system(), platform.release(), platform.machine()))
      out.write("python\t%s\n" % platform.python_version())
      out.write("transxpress_commit\t%s\n" % commit)
      out.write("assembler\t%s\n" % config["assembler"])
      out.write("assembler_threads\t%s\n" % (
        config.get("rnaspades_threads", 32) if config["assembler"] == "rnaspades"
        else config.get("trinity_final_threads", 16)))


rule fastqc_before_trim:
  """
  Runs fastQC on individual input reads files.
  """
  input:
    samples=config["samples_file"]
  output:
    directory("results/qc/fastqc_before_trim")
  log:
    "logs/fastqc_before_trim.log"
  conda:
    "envs/qc.yaml"
  params:
    memory="6"
  threads:
    8
  shell:
    """
    # run fastqc on input files
    # making directory for the fastqc summary files
    mkdir {output} &> {log}
    FILES=$(awk '{{ printf "%s\\n%s\\n", $3,$4}}' {input}) &>> {log}
    # running fastqc on the files, {threads} of them at a time
    fastqc -f fastq -t {threads} -o {output} $FILES &>> {log}
    """


rule multiqc_before_trim:
  """
  Creates multiQC report from individual fastQC reports.
  """
  input:
    "results/qc/fastqc_before_trim"
  output:
    out_dir=directory("results/qc/multiqc_before_trim"),
    report="results/qc/multiqc_before_trim.txt"
  log:
    "logs/multiqc_before_trim.log"
  conda:
    "envs/qc.yaml"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # common error resolved by those two export commands
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    mkdir {output.out_dir} &> {log}
    multiqc -o {output.out_dir} {input} &>> {log}
    cp {output.out_dir}/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """

rule fastqc_before_trim_warnings:
  """
  Prints warnings for FastQC modules which produced some warnings or fails of the input data.
  """
  input:
    fastqc_report="results/qc/multiqc_before_trim.txt"
  output:
    warning_file="results/qc/WARNING_fastqc_before_trim_overview.txt"
  log:
    "logs/fastqc_before_trim_warnings.log"
  run:

    mapping = {
    "basic_statistics": "BASIC STATISTICS",
    "per_base_sequence_quality": "PER BASE SEQUENCE QUALITY",
    "per_tile_sequence_quality": "PER TILE SEQUENCE QUALITY",
    "per_sequence_quality_scores": "PER SEQUENCE QUALITY SCORES",
    "per_base_sequence_content": "PER BASE SEQUENCE CONTENT",
    "per_sequence_gc_content": "PER SEQUENCE GC CONTENT",
    "per_base_n_content": "PER BASE N CONTENT",
    "sequence_length_distribution": "SEQUENCE LENGTH DISTRIBUTION",
    "sequence_duplication_levels": "SEQUENCE DUPLICATION LEVELS",
    "overrepresented_sequences": "OVERREPRESENTED SEQUENCES",
    "adapter_content": "ADAPTER CONTENT"
    }

    reasons_mapping = {
    "basic_statistics": "BASIC STATISTICS",
    
    "per_base_sequence_quality":
    """A warning is issued if the lower quartile for any base is less than 10, or if the median for any base is less than 25.
This module will raise a failure if the lower quartile for any base is less than 5 or if the median for any base is less than 20.

Common reasons for warnings:
The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

If the quality of the library falls to a low level then the most common remedy is to perform quality trimming where reads are truncated based on their average quality. For most libraries where this type of degradation has occurred you will often be simultaneously running into the issue of adapter read-through so a combined adapter and quality trimming step is often employed.

Another possibility is that a warn / error is triggered because of a short loss of quality earlier in the run, which then recovers to produce later good quality sequence. This can happen if there is a transient problem with the run (bubbles passing through a flowcell for example). You can normally see this type of error by looking at the per-tile quality plot (if available for your platform). In these cases trimming is not advisable as it will remove later good sequence, but you might want to consider masking bases during subsequent mapping or assembly.

If your library has reads of varying length then you can find a warning or error is triggered from this module because of very low coverage for a given base range. Before committing to any action, check how many sequences were responsible for triggering an error by looking at the sequence length distribution module results.
""",
    
    "per_tile_sequence_quality":
"""This module issues a warning if any tile shows a mean Phred score more than 2 less than the mean for that base across all tiles.
This module will issue a warning if any tile shows a mean Phred score more than 5 less than the mean for that base across all tiles.

Common reasons for warnings:
Whilst warnings in this module can be triggered by individual specific events we have also observed that greater variation in the phred scores attributed to tiles can also appear when a flowcell is generally overloaded. In this case events appear all over the flowcell rather than being confined to a specific area or range of cycles. We would generally ignore errors which mildly affected a small number of tiles for only 1 or 2 cycles, but would pursue larger effects which showed high deviation in scores, or which persisted for several cycles.
""",
    
    "per_sequence_quality_scores": 
"""A warning is raised if the most frequently observed mean quality is below 27 - this equates to a 0.2% error rate.
An error is raised if the most frequently observed mean quality is below 20 - this equates to a 1% error rate.

Common reasons for warnings:
This module is generally fairly robust and errors here usually indicate a general loss of quality within a run. For long runs this may be alleviated through quality trimming. If a bi-modal, or complex distribution is seen then the results should be evaluated in concert with the per-tile qualities (if available) since this might indicate the reason for the loss in quality of a subset of sequences.
""",
    
    "per_base_sequence_content": 
"""Common reasons for warnings:
The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

If the quality of the library falls to a low level then the most common remedy is to perform quality trimming where reads are truncated based on their average quality. For most libraries where this type of degradation has occurred you will often be simultaneously running into the issue of adapter read-through so a combined adapter and quality trimming step is often employed.

Another possibility is that a warn / error is triggered because of a short loss of quality earlier in the run, which then recovers to produce later good quality sequence. This can happen if there is a transient problem with the run (bubbles passing through a flowcell for example). You can normally see this type of error by looking at the per-tile quality plot (if available for your platform). In these cases trimming is not advisable as it will remove later good sequence, but you might want to consider masking bases during subsequent mapping or assembly.

If your library has reads of varying length then you can find a warning or error is triggered from this module because of very low coverage for a given base range. Before committing to any action, check how many sequences were responsible for triggering an error by looking at the sequence length distribution module results.
""",
    
    "per_sequence_gc_content": 
"""A warning is raised if the sum of the deviations from the normal distribution represents more than 15% of the reads.
This module will indicate a failure if the sum of the deviations from the normal distribution represents more than 30% of the reads.

Common reasons for warnings:
Warnings in this module usually indicate a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example), which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species.
""",
    
    "per_base_n_content": 
"""This module raises a warning if any position shows an N content of >5%.
This module will raise an error if any position shows an N content of >20%.

Common reasons for warnings:
The most common reason for the inclusion of significant proportions of Ns is a general loss of quality, so the results of this module should be evaluated in concert with those of the various quality modules. You should check the coverage of a specific bin, since it is possible that the last bin in this analysis could contain very few sequences, and an error could be prematurely triggered in this case.

Another common scenario is the incidence of a high proportions of N at a small number of positions early in the library, against a background of generally good quality. Such deviations can occur when you have very biased sequence composition in the library to the point that base callers can become confused and make poor calls. This type of problem will be apparent when looking at the per-base sequence content results.
""",
    
    "sequence_length_distribution": 
"""This module will raise a warning if all sequences are not the same length.
This module will raise an error if any of the sequences have zero length.

Common reasons for warnings:
For some sequencing platforms it is entirely normal to have different read lengths so warnings here can be ignored.""",
    
    "sequence_duplication_levels": 
"""This module will issue a warning if non-unique sequences make up more than 20% of the total.
This module will issue a error if non-unique sequences make up more than 50% of the total.

Common reasons for warnings:
The underlying assumption of this module is of a diverse unenriched library. Any deviation from this assumption will naturally generate duplicates and can lead to warnings or errors from this module.
        
In general there are two potential types of duplicate in a library, technical duplicates arising from PCR artefacts, or biological duplicates which are natural collisions where different copies of exactly the same sequence are randomly selected. From a sequence level there is no way to distinguish between these two types and both will be reported as duplicates here.
        
A warning or error in this module is simply a statement that you have exhausted the diversity in at least part of your library and are re-sequencing the same sequences. In a supposedly diverse library this would suggest that the diversity has been partially or completely exhausted and that you are therefore wasting sequencing capacity. However in some library types you will naturally tend to over-sequence parts of the library and therefore generate duplication and will therefore expect to see warnings or error from this module.
        
In RNA-Seq libraries sequences from different transcripts will be present at wildly different levels in the starting population. In order to be able to observe lowly expressed transcripts it is therefore common to greatly over-sequence high expressed transcripts, and this will potentially create large set of duplicates. This will result in high overall duplication in this test, and will often produce peaks in the higher duplication bins. This duplication will come from physically connected regions, and an examination of the distribution of duplicates in a specific genomic region will allow the distinction between over-sequencing and general technical duplication, but these distinctions are not possible from raw fastq files. A similar situation can arise in highly enriched ChIP-Seq libraries although the duplication there is less pronounced. Finally, if you have a library where the sequence start points are constrained (a library constructed around restriction sites for example, or an unfragmented small RNA library) then the constrained start sites will generate huge dupliction levels which should not be treated as a problem, nor removed by deduplication. In these types of library you should consider using a system such as random barcoding to allow the distinction of technical and biological duplicates.
""",
    
    "overrepresented_sequences": 
"""This module will issue a warning if any sequence is found to represent more than 0.1% of the total.
This module will issue an error if any sequence is found to represent more than 1% of the total.
        
Common reasons for warnings:
This module will often be triggered when used to analyse small RNA libraries where sequences are not subjected to random fragmentation, and the same sequence may natrually be present in a significant proportion of the library.
""",
    
    "adapter_content": 
"""This module will issue a warning if any sequence is present in more than 5% of all reads.
This module will issue an error if any sequence is present in more than 10% of all reads.
        
Common reasons for warnings:
Any library where a reasonable proportion of the insert sizes are shorter than the read length will trigger this module. This does not indicate a problem as such - just that the sequences will need to be adapter trimmed before proceeding with any downstream analysis.
"""
    }



    with open(output["warning_file"], "w") as out_file:
      df1 = pd.read_csv(input["fastqc_report"], sep='\t')
      df1 = df1.set_index('Sample')

      # filter only columns containing 'pass','warn' or 'fail'
      df1_filtered = df1.loc[:, df1.isin(['pass','warn','fail']).all()]
      df1_filtered = df1_filtered.loc[:, df1_filtered.columns != 'basic_statistics'] # basic statistics module always passes

      n_samples = df1.shape[0]

      for column in df1_filtered.columns:
        if column in mapping.keys():
          out_file.writelines(mapping[column]+'\n')
        else:
          out_file.writelines(column+'\n')
        results = df1_filtered[column].tolist()
        if ('warn' in results) or ('fail' in results):
          n_warn = results.count('warn')
          n_fail = results.count('fail')
          n_warn_perc = round(n_warn/n_samples*100,2)
          n_fail_perc = round(n_fail/n_samples*100,2)
          if n_warn > 0:
            out_file.writelines(f'{n_warn} out of {n_samples} produced "WARN" ({n_warn_perc}%)\n')
          if n_fail > 0:
            out_file.writelines(f'{n_fail} out of {n_samples} produced "FAIL" ({n_fail_perc}%)\n')
          if column in reasons_mapping.keys():
            out_file.writelines('\n' + reasons_mapping[column])
        else:
          out_file.writelines('ok'+'\n')
      out_file.writelines('source: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/\n')
    
    # print to stdout
    with open(output["warning_file"], "r") as out_file:
      for line in out_file.readlines():
        print(line.strip())
      print('See this output in: ' + output["warning_file"])
    
checkpoint trimmomatic_split:
  """
  Splits file with information about all reads files into separate files
  so they can be processed with trimmomatic in parallel.
  """   
  input:
    samples=config["samples_file"]
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
    "envs/trimmomatic.yaml"
  params:
    memory="4"
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
      trimmomatic PE -threads {threads} $F_READS $R_READS results/trimmomatic/reads/{wildcards[job_index]}.R1-P.qtrim.fastq.gz results/trimmomatic/reads/{wildcards[job_index]}.R1-U.qtrim.fastq.gz results/trimmomatic/reads/{wildcards[job_index]}.R2-P.qtrim.fastq.gz results/trimmomatic/reads/{wildcards[job_index]}.R2-U.qtrim.fastq.gz {config[trimmomatic_parameters]} &>> {log}
      echo $SAMPLE	$REPLICATE	results/trimmomatic/reads/{wildcards[job_index]}.R1-P.qtrim.fastq.gz	results/trimmomatic/reads/{wildcards[job_index]}.R2-P.qtrim.fastq.gz > {output} 2>> {log}
    else
      trimmomatic SE -threads {threads} $F_READS results/trimmomatic/reads/{wildcards[job_index]}.U.qtrim.fastq.gz {config[trimmomatic_parameters]} &>> {log}
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

rule fastqc_after_trim:
  """
  Runs fastQC on individual trimmed reads files.
  """
  input:
    samples="results/trimmomatic/samples_trimmed.txt"
  output:
    directory("results/qc/fastqc_after_trim")
  log:
    "logs/fastqc_after_trim.log"
  conda:
    "envs/qc.yaml"
  params:
    memory="6"
  threads:
    8
  shell:
    """
    # run fastqc on input files
    # making directory for the fastqc summary files
    mkdir {output} &> {log}
    FILES=$(awk '{{ printf "%s\\n%s\\n", $3,$4}}' {input}) &>> {log}
    # running fastqc on the files, {threads} of them at a time
    fastqc -f fastq -t {threads} -o {output} $FILES &>> {log}
    """


rule multiqc_after_trim:
  """
  Creates multiQC report from individual fastQC reports after the trimming.
  """
  input:
    "results/qc/fastqc_after_trim"
  output:
    out_directory=directory("results/qc/multiqc_after_trim"),
    report="results/qc/multiqc_after_trim.txt"
  log:
    "logs/multiqc_after_trim.log"
  conda:
    "envs/qc.yaml"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # common error resolved by those two export commands
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    mkdir {output.out_directory} &> {log}
    multiqc -o {output.out_directory} {input} &>> {log}
    cp {output.out_directory}/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """

rule compare_qc_after_trim:
  """
  Compares the quality of reads before and after trimming.
  Compares PER BASE SEQUENCE QUALITY, PER SEQUENCE QUALITY SCORES, 
  OVERREPRESENTED SEQUENCES and ADAPTER CONTENT modules of FastQC.
  """
  input:
    before="results/qc/multiqc_before_trim.txt",
    after="results/qc/multiqc_after_trim.txt"
  output:
    result="results/qc/FastQC_comparison_after_trim.txt"
  log:
    "logs/compare_qc_after_trim.log"
  run:

    mapping = {
      "per_base_sequence_quality": "PER BASE SEQUENCE QUALITY",
      "per_sequence_quality_scores": "PER SEQUENCE QUALITY SCORES",
      "overrepresented_sequences": "OVERREPRESENTED SEQUENCES",
      "adapter_content": "ADAPTER CONTENT"
    }

    def compare_module(module_name, before_trim, after_trim, out_file):
      module_name = mapping[module_name]
      out_file.writelines(f"""**{module_name}**
Before trimming:
pass: {str(before_trim.count('pass'))} out of {str(len(before_trim))} ({str(before_trim.count('pass')/len(before_trim)*100)}%)
warn: {str(before_trim.count('warn'))} out of {str(len(before_trim))} ({str(before_trim.count('warn')/len(before_trim)*100)}%)
fail: {str(before_trim.count('fail'))} out of {str(len(before_trim))} ({str(before_trim.count('fail')/len(before_trim)*100)}%)
      
After trimming:
pass: {str(after_trim.count('pass'))} out of {str(len(after_trim))} ({str(after_trim.count('pass')/len(after_trim)*100)}%)
warn: {str(after_trim.count('warn'))} out of {str(len(after_trim))} ({str(after_trim.count('warn')/len(after_trim)*100)}%)
fail: {str(after_trim.count('fail'))} out of {str(len(after_trim))} ({str(after_trim.count('fail')/len(after_trim)*100)}%)
{'-' * 50}
""")

    before = pd.read_csv(input["before"], sep='\t')
    before = before.set_index('Sample')

    after = pd.read_csv(input["after"], sep='\t')
    after = after.set_index('Sample')

    modules = ['per_base_sequence_quality', 'per_sequence_quality_scores', 'overrepresented_sequences', 'adapter_content']

    with open(output["result"], "w") as out_file:
      out_file.write("FASTQC COMPARISON BEFORE AND AFTER TRIMMING\n")
      out_file.write('-' * 50 + '\n')

      for module in modules:
        before_list = before[module].tolist()
        after_list = after[module].tolist()
        compare_module(module, before_list, after_list, out_file)

    # print the results also in the terminal
    with open(output["result"], "r") as out_file:
      for line in out_file.readlines():
        print(line.strip())
      print('See this output in ' + output["result"])

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


rule trinity_inchworm_chrysalis:
  """
  Runs first two stages of Trinity assembly (Inchworm and Chrysalis).
  """
  input:
    samples="results/trimmomatic/samples_trimmed.txt",
  output:
    "results/assembly/trinity/recursive_trinity.cmds"
  log:
    "logs/trinity_inchworm_chrysalis.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory=config.get("trinity_inchworm_chrysalis_memory", "450")
  threads:
    config.get("trinity_inchworm_chrysalis_threads", 16)
  # The common failures here are memory and time. A retry asks for the same figures and starts from the
  # beginning, so it would fail the same way, after hours on a large node.
  retries: 0
  shell:
    """
    # record the computer this ran on, because the assembly depends on it
    ( hostname ; grep -m1 "model name" /proc/cpuinfo 2>/dev/null || uname -m ; nproc 2>/dev/null || true ) > {log}

    # Trinity ticks off the step that writes {output} with results/assembly/trinity/recursive_trinity.cmds.ok and
    # does not check that the file that step wrote is still in place. trinity_final consumes the file and
    # leaves recursive_trinity.cmds.completed behind, so after a pass the tick is there and the file is
    # not. Trinity then skips the step, exits 0 and writes no output. So the tick goes first.
    rm -f results/assembly/trinity/recursive_trinity.cmds.ok &>> {log}

    # trinity_final deletes read_partitions once the assembly is written, because it holds one folder per
    # read cluster and can reach millions of files. Trinity does not partition the reads again while
    # results/assembly/trinity/partitioned_reads.files.list.ok is in place, and it leaves the top folder of
    # read_partitions behind, so neither the folder nor its emptiness tells the two cases apart. Test the
    # list itself: if the first read file it names is gone, the partitions are gone with it and every
    # Butterfly chunk would stop on a missing input. A crashed run whose partitions are still in place
    # keeps its tick, so its hours of partitioning are not repeated.
    PARTITION_LIST=results/assembly/trinity/partitioned_reads.files.list
    if [ ! -s $PARTITION_LIST ] || [ ! -e "$(head -1 $PARTITION_LIST)" ]; then
      echo "read partitions are missing, so Trinity partitions the reads again" &>> {log}
      rm -f $PARTITION_LIST.ok &>> {log}
    fi

    Trinity --no_distributed_trinity_exec --max_memory {params.memory}G --CPU {threads} --output results/assembly/trinity --samples_file {input} {config[trinity_parameters]} {config[strand_specific]} &>> {log}
    """


checkpoint trinity_butterfly_split:
  """
  Preparation for last stage of Trinity (Butterfly) parallelization by splitting
  the commands into independent parts.
  """
  input:
    "results/assembly/trinity/recursive_trinity.cmds"
  output:
    directory("results/assembly/trinity/parallel_jobs")
  log:
    "logs/trinity_butterfly_split.log"
  params:
    cmds_per_chunk=config.get("trinity_butterfly_cmds_per_chunk", 20),
    max_chunks=config.get("trinity_butterfly_max_chunks", 1000)
  shell:
    """
    mkdir -p {output} &> {log}
    # Each chunk is one Slurm job. A chunk holds at least cmds_per_chunk commands, so small data runs
    # a few jobs, and there are at most max_chunks chunks, so large data does not submit more jobs.
    COMMANDS=$(wc -l < {input})
    CHUNKS=$(( (COMMANDS + {params.cmds_per_chunk} - 1) / {params.cmds_per_chunk} ))
    if [ $CHUNKS -gt {params.max_chunks} ]; then CHUNKS={params.max_chunks}; fi
    if [ $CHUNKS -lt 1 ]; then CHUNKS=1; fi
    echo "$COMMANDS commands in $CHUNKS chunks" &>> {log}
    # Note: it is important to use -d and not --numeric-suffixes, see https://github.com/transXpress/transXpress-snakemake/issues/12
    split -n l/$CHUNKS -e -d {input} {output}/job_ &>> {log}
    """


rule trinity_butterfly_parallel:
  """
  Runs Trinity Butterfly commands (which were split into parts) in parallel.
  """
  input:
    "results/assembly/trinity/parallel_jobs/job_{job_index}"
  output:
    # Keep this inside parallel_jobs. Each split deletes the folder, so every chunk runs again. After
    # trinity_inchworm_chrysalis runs again, the job_ files can hold the same commands, and Snakemake
    # skips a chunk whose job_ file has the same content and whose output sits in another folder.
    "results/assembly/trinity/parallel_jobs/completed_{job_index}"
  log:
    "logs/trinity_butterfly_parallel/{job_index}.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    # A measured chunk job peaks near 550 MB. This rule has one instance per
    # chunk, up to trinity_butterfly_max_chunks of them, so an oversized request here decides how many
    # land on one node and leaves the rest of its cores idle.
    memory="3",
    scratch_var=SCRATCH_VAR
  threads:
    1
  shell:
    """
    # Each line of the chunk is one small Trinity assembly, which writes and deletes many files. When the
    # scratch key names a variable that is set in this job, --workdir puts those files on the node's local
    # disk. Trinity still writes the fasta of each assembly into read_partitions.
    SCRATCH_VAR="{params.scratch_var}"
    LOCAL_ROOT=""
    if [ -n "$SCRATCH_VAR" ]; then LOCAL_ROOT="${{!SCRATCH_VAR:-}}"; fi
    if [ -n "$LOCAL_ROOT" ]; then
      RUN=$(pwd)
      WORK=$(mktemp -d "$LOCAL_ROOT/transxpress.XXXXXX")
      echo "work folder: $WORK" &> {log}
      # On a failure, the work folders go to read_partitions, where they would stand without --workdir.
      # A job that Slurm kills at its time limit or with scancel copies nothing.
      on_exit() {{
        rc=$?
        set +e
        trap - EXIT TERM
        if [ $rc -ne 0 ] && [ -d "$WORK/trinity_work" ]; then
          echo "the chunk stopped with code $rc, so $WORK/trinity_work goes to results/assembly/trinity/read_partitions" &>> {log}
          cp -rp "$WORK/trinity_work/." results/assembly/trinity/read_partitions/ &>> {log}
        fi
        rm -rf "$WORK"
        exit $rc
      }}
      trap on_exit EXIT
      trap "exit 143" TERM
      # Trinity accepts only a work folder with "trinity" in its name
      sed 's|$| --workdir '"$WORK"'/trinity_work|' {input} | bash &>> {log}
    else
      bash {input} &> {log}
    fi
    cp -p {input} {output} &>> {log}
    """


def trinity_completed_parallel_jobs(wildcards):
  """
  Returns filenames of the files processed in parallel.
  """
  parallel_dir = checkpoints.trinity_butterfly_split.get(**wildcards).output[0]
  # sorted, because glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "job_{job_index}")).job_index)
  completed_ids = expand(os.path.join(parallel_dir,"completed_{job_index}"), job_index=job_ids)
  return completed_ids

rule trinity_butterfly_parallel_merge:
  input:
    jobs=trinity_completed_parallel_jobs,
    cmds="results/assembly/trinity/recursive_trinity.cmds"
  output:
    cmds_completed="results/assembly/trinity/recursive_trinity.cmds.completed"
  log:
    "logs/trinity_butterfly_parallel_merge.log"
  params:
    memory="10"
  threads:
    1
  shell:
    """
    # Snakemake runs this script as one string of at most 131071 bytes, so naming every chunk file here
    # stops near 2900 chunks with "Argument list too long", see
    # https://bitbucket.org/snakemake/snakemake/issues/878/errno-7-argument-list-too-long-path-to-bin
    # xargs passes the files to cat in calls that fit the limit, so the chunk count has no limit here.
    # The list comes from the job_ files, as the rule input does, so a stray completed_ file is not read
    # and a missing one stops cat. The chunk suffixes have a fixed width, so sort gives the order of the split.
    {{ find results/assembly/trinity/parallel_jobs -maxdepth 1 -name 'job_*' | LC_ALL=C sort | sed 's|/job_|/completed_|' | xargs -r cat ; }} > {output.cmds_completed} 2> {log}
    """


rule trinity_final:
  """
  Runs the final Trinity assembly.
  """
  input:
    cmds_completed="results/assembly/trinity/recursive_trinity.cmds.completed",
    samples="results/trimmomatic/samples_trimmed.txt"
  output:
    # Trinity writes these two beside its output folder, not inside it
    transcriptome="results/assembly/trinity.Trinity.fasta",
    gene_trans_map="results/assembly/trinity.Trinity.fasta.gene_trans_map"
  log:
    "logs/trinity_final.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory=config.get("trinity_final_memory", "450")
  threads:
    config.get("trinity_final_threads", 16)
  # The common failures here are memory and time. A retry asks for the same figures and starts from the
  # beginning, so it would fail the same way, after hours on a large node.
  retries: 0
  shell:
    """
    # record the computer this ran on, because the assembly depends on it
    ( hostname ; grep -m1 "model name" /proc/cpuinfo 2>/dev/null || uname -m ; nproc 2>/dev/null || true ) > {log}
    # read_partitions is deleted after a successful assembly (see below), and Trinity needs it to collect the transcripts
    # Trinity does not partition the reads again while its own checkpoint files exist, so they must be removed first
    if [ ! -d results/assembly/trinity/read_partitions ]; then
      echo "ERROR: results/assembly/trinity/read_partitions is missing. It was deleted after an earlier assembly." | tee -a {log} >&2
      echo "To assemble again: rm -f results/assembly/trinity/partitioned_reads.files.list.ok results/assembly/trinity/recursive_trinity.cmds.ok" | tee -a {log} >&2
      echo "then run snakemake again with: --forcerun trinity_inchworm_chrysalis" | tee -a {log} >&2
      exit 1
    fi

    # Trinity ticks off its last step with results/assembly/trinity/__salmon_filt.chkpts/salmon.filt.ok and does not
    # check that the file that step wrote is still in place. With the tick of an earlier assembly present,
    # Trinity skips the one command that writes {output.transcriptome}, prints "Trinity.fasta 0 MByte" and
    # exits 0 with no assembly. So the ticks and the half-finished output of an earlier assembly go first.
    rm -rf results/assembly/trinity/__salmon_filt.chkpts results/assembly/trinity/Trinity.tmp.fasta &>> {log}

    Trinity --max_memory {params.memory}G --CPU {threads} --output results/assembly/trinity --samples_file {input.samples} {config[trinity_parameters]} {config[strand_specific]} &>> {log}

    # read_partitions holds one folder per read cluster (can be millions of files) and is not needed after the assembly
    # Delete it only once the assembly is in place. Trinity can exit 0 without writing it, and read_partitions
    # is the only input the assembly can be repeated from, so deleting it on a failed run loses the run.
    if [ -s {output.transcriptome} ] && [ -s {output.gene_trans_map} ]; then
      # rename first, so that a half-deleted folder is never left under its real name
      rm -rf results/assembly/trinity/read_partitions.deleting &>> {log}
      mv results/assembly/trinity/read_partitions results/assembly/trinity/read_partitions.deleting &>> {log}
      rm -rf results/assembly/trinity/read_partitions.deleting &>> {log}
    else
      # fail here, so the log names the cause and Snakemake does not report only a missing output file
      echo "ERROR: Trinity exited 0 but did not write {output.transcriptome} and {output.gene_trans_map}." | tee -a {log} >&2
      echo "results/assembly/trinity/read_partitions is kept, so the assembly can run again. See {log}" | tee -a {log} >&2
      exit 1
    fi
    """


rule rnaspades:
  """
  Runs rnaSPAdes assembly. rnaSPAdes writes three transcript sets, and
  rnaspades_filter in config.yaml picks the one that becomes transcriptome.fasta.
  """
  input:
    samples="results/trimmomatic/samples_trimmed.yaml",
  output:
    default="results/assembly/rnaspades/transcripts.fasta",
    hard="results/assembly/rnaspades/hard_filtered_transcripts.fasta",
    soft="results/assembly/rnaspades/soft_filtered_transcripts.fasta"
    ##config["annotated_fasta_prefix"]+"_spades_out/K"+THEKMER+"/assembly_graph_with_scaffolds.gfa"
  log:
    "logs/rnaspades.log"
  conda:
    "envs/rnaspades.yaml"
  params:
    memory=config.get("rnaspades_memory", "400"),
    scratch_var=SCRATCH_VAR
  threads:
    config.get("rnaspades_threads", 32)
  # The common failures here are memory and time. A retry asks for the same figures and starts from the
  # beginning, so it would fail the same way, after hours on a large node.
  retries: 0
  shell:
    """
    # record the computer this ran on, because the assembly depends on it
    ( hostname ; grep -m1 "model name" /proc/cpuinfo 2>/dev/null || uname -m ; nproc 2>/dev/null || true ) > {log}

    # rnaSPAdes writes about 100 GB of working files on a large dataset. When the scratch key names a
    # variable that is set in this job, it works in a folder on the node's local disk, and its output
    # folder is copied to the run folder at the end. Otherwise it works in the run folder.
    OUT=results/assembly/rnaspades
    SPADES_OUT=$OUT
    SCRATCH_VAR="{params.scratch_var}"
    LOCAL_ROOT=""
    if [ -n "$SCRATCH_VAR" ]; then LOCAL_ROOT="${{!SCRATCH_VAR:-}}"; fi
    # A scratch folder in RAM (tmpfs) counts against the memory of the job, and this rule writes large
    # work files, so it stops there
    if [ -n "$LOCAL_ROOT" ] && [ "$(stat -f -c %T "$LOCAL_ROOT" 2>/dev/null)" = "tmpfs" ]; then
      echo "ERROR: the scratch folder $LOCAL_ROOT is in RAM (tmpfs) on $(hostname). Run this rule on a node with a local disk." | tee -a {log} >&2
      exit 1
    fi
    if [ -n "$LOCAL_ROOT" ]; then
      RUN=$(pwd)
      WORK=$(mktemp -d "$LOCAL_ROOT/transxpress.XXXXXX")
      SPADES_OUT=$WORK/rnaspades
      echo "work folder: $WORK" &>> {log}
      # On a failure, the local folder goes to $OUT, where it would stand without the local disk.
      # A job that Slurm kills at its time limit or with scancel copies nothing. The local folder is deleted
      # in all cases.
      on_exit() {{
        rc=$?
        set +e
        trap - EXIT TERM
        if [ $rc -ne 0 ]; then
          echo "rnaspades stopped with code $rc, so $WORK/rnaspades goes to $OUT" &>> {log}
          mkdir -p "$OUT" &>> {log}
          cp -rp "$WORK/rnaspades/." "$OUT/" &>> {log}
        fi
        rm -rf "$WORK"
        exit $rc
      }}
      trap on_exit EXIT
      trap "exit 143" TERM
    fi

    # No --checkpoints: Slurm deletes the local folder when the job ends, so a restart point would not survive
    rnaspades.py --dataset {input.samples} -t {threads} -m {params.memory} -o "$SPADES_OUT" {config[strand_specific]} --only-assembler {config[kmers]} &>> {log}

    if [ "$SPADES_OUT" != "$OUT" ]; then
      SIZE=$(du -sh "$WORK" | cut -f1)
      START=$(date +%s)
      # replace the whole folder, so no file of an earlier failed run stays beside the new assembly
      rm -rf "$OUT" &>> {log}
      cp -rp "$SPADES_OUT" "$(dirname "$OUT")/" &>> {log}
      echo "copied $SIZE from $WORK to $OUT in $(( $(date +%s) - START )) s" &>> {log}
    fi
    """
 

def transcriptome_copy_input(wildcards):
  if config["assembler"] == "rnaspades":
    return {"transcriptome": getattr(rules.rnaspades.output, RNASPADES_FILTER)}
  return {"transcriptome": rules.trinity_final.output.transcriptome,
          "gene_trans_map": rules.trinity_final.output.gene_trans_map}


rule transcriptome_copy:
  """
  Copies the assembled transcriptome to the transxpress directory, with its
  gene to transcript map. Trinity writes the map. For rnaSPAdes the rule builds
  it from the transcript names: NODE_97_length_6237_cov_11.98_g8_i2 is isoform 2
  of gene g8, so the map line is "g8<TAB>NODE_97_length_6237_cov_11.98_g8_i2".
  A name without _g<N>_i<M> at its end stops the rule.
  """
  input:
    unpack(transcriptome_copy_input)
  output:
    transcriptome="results/assembly/transcriptome.fasta",
    gene_trans_map="results/assembly/transcriptome.gene_trans_map"
  log:
    "logs/transcriptome_copy.log"
  run:
    shell("cp -p {input.transcriptome} {output.transcriptome} &> {log}")
    if config["assembler"] == "rnaspades":
      shell("""seqkit seq -n -i {output.transcriptome} 2>> {log} | awk 'BEGIN {{ OFS = "\\t" }} match($0, /_g[0-9]+_i[0-9]+$/) {{ split(substr($0, RSTART + 1), part, "_"); print part[1], $0; next }} {{ print "ERROR: transcript " $0 " has no _g<N>_i<M> at the end of its name, so its gene is unknown" > "/dev/stderr"; bad = 1 }} END {{ exit bad }}' > {output.gene_trans_map} 2>> {log}""")
    else:
      shell("cp -p {input.gene_trans_map} {output.gene_trans_map} &>> {log}")


rule trinity_stats:
  """
  Runs Trinity script to get statistics about the assembled transcriptome
  (number of transcripts, genes, GC content, EXN50).
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
    expression="results/expression/transcriptome_expression_isoform.tsv"
  output:
    stats="results/qc/transcriptome_stats.txt",
    exN50="results/qc/transcriptome_exN50.tsv",
    exN50plot="results/qc/ExN50_plot.pdf"
  log:
    "logs/trinity_stats.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')
    # contig_ExN50_statistic.pl reads the genes in Perl hash order, which is random in each run, and
    # genes tied on expression keep that order. A fixed Perl hash seed gives the same order in each run.
    export PERL_HASH_SEED=0

    assembler={config[assembler]}
    if [ "$assembler" = 'trinity' ]; then
        $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} > {output.stats} 2> {log}
    else
        $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} | sed -e 's/trinity/rnaspades/g' > {output.stats} 2> {log}
    fi
    # contig_ExN50_statistic.pl writes its .E-inputs table and the R script writes ExN50_plot.pdf into the
    # folder they run in, so both run in the folder of the outputs
    RUN=$(pwd)
    cd $(dirname {output.exN50})
    $TRINITY_HOME/util/misc/contig_ExN50_statistic.pl "$RUN/{input.expression}" "$RUN/{input.transcriptome}" > $(basename {output.exN50}) 2>> "$RUN/{log}"
    $TRINITY_HOME/util/misc/plot_ExN50_statistic.Rscript $(basename {output.exN50}) &>> "$RUN/{log}"
    """

rule busco:
  """
  Runs BUSCO to assess the completeness of the transcriptome.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta"
  output:
    out_directory=directory("results/qc/busco"),
    report="results/qc/busco_report.txt"
  log:
    "logs/busco.log"
  conda:
    "envs/busco.yaml"
  params:
    memory=config.get("busco_memory", "10"),
    lineage=LINEAGE,
    # a shared db_dir holds the lineage, and BUSCO must not download or update anything there
    download_option=("--offline --download_path " + DB_DIR + "/busco_downloads") if DB_SHARED else "--download_path resources/busco_downloads",
    scratch_var=SCRATCH_VAR
  threads:
    config.get("busco_threads", 16)
  shell:
    """
    # BUSCO writes thousands of small files. When the scratch key names a variable that is set in this
    # job, it works in a folder on the node's local disk, packs the per-gene folders there and copies its
    # output folder to the run folder at the end. Otherwise it works in the run folder. A failed run copies
    # nothing, because Snakemake deletes the output folder of a failed job; BUSCO's messages are in {log}.
    OUT={output.out_directory}
    BUSCO_PATH=$(dirname $OUT)
    SCRATCH_VAR="{params.scratch_var}"
    LOCAL_ROOT=""
    if [ -n "$SCRATCH_VAR" ]; then LOCAL_ROOT="${{!SCRATCH_VAR:-}}"; fi
    : > {log}
    if [ -n "$LOCAL_ROOT" ]; then
      WORK=$(mktemp -d "$LOCAL_ROOT/transxpress.XXXXXX")
      BUSCO_PATH=$WORK
      echo "work folder: $WORK" &>> {log}
      trap 'rm -rf "$WORK"' EXIT
    fi

    # an empty or missing lineage lets BUSCO choose the lineage automatically
    if [ -z "{params.lineage}" ]; then
      lineage_option="--auto-lineage"
    else
      lineage_option="-l {params.lineage}"
    fi

    # -o takes a folder name and --out_path the folder that holds it
    busco -m transcriptome -i {input.transcriptome} -o $(basename $OUT) --out_path $BUSCO_PATH $lineage_option -c {threads} {params.download_option} &>> {log}
    echo "BUSCO run completed successfully" &>> {log}

    # pack the per-gene folders (thousands of small files) into one archive each
    for d in $BUSCO_PATH/$(basename $OUT)/run_*/busco_sequences $BUSCO_PATH/$(basename $OUT)/run_*/hmmer_output; do
      if [ -d "$d" ]; then
        tar -czf "$d.tar.gz" -C "$(dirname "$d")" "$(basename "$d")" &>> {log}
        rm -rf "$d" &>> {log}
      fi
    done

    if [ -n "$LOCAL_ROOT" ]; then
      SIZE=$(du -sh "$WORK" | cut -f1)
      START=$(date +%s)
      rm -rf $OUT &>> {log}
      cp -rp $WORK/$(basename $OUT) $(dirname $OUT)/ &>> {log}
      echo "copied $SIZE from $WORK to $OUT in $(( $(date +%s) - START )) s" &>> {log}
    fi

    # with --auto-lineage BUSCO also writes a short_summary.generic.* file, the report is the specific one
    cp $OUT/short_summary.specific.*.txt {output.report} &>> {log}
    """

rule transdecoder_longorfs:
  """
  Runs first stage of Transdecoder extracting the long open reading frames.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
  output:
    orfs="results/transdecoder/transcriptome.orfs"
  log:
    "logs/transdecoder_longorfs.log"
  conda:
    "envs/transdecoder.yaml"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    # LongOrfs names its folder after the file name of the transcriptome
    TD_DIR=results/transdecoder/longorfs/$(basename {input.transcriptome}).transdecoder_dir
    # LongOrfs skips a step when its checkpoint file exists, so a rerun starts from an empty folder
    rm -rf $TD_DIR &> {log}
    # A fixed Perl hash seed gives the same output order in each run
    export PERL_HASH_SEED=0
    # TransDecoder 6 keeps the stage scripts in util/ next to the TransDecoder wrapper, outside bin/
    TD_UTIL=$(dirname $(readlink -f $(command -v TransDecoder)))/util
    $TD_UTIL/TransDecoder.LongOrfs -t {input.transcriptome} --output_dir results/transdecoder/longorfs &>> {log}
    cp -p $TD_DIR/longest_orfs.pep {output.orfs} &>> {log}
    """


rule transdecoder_predict_all:
  """
  Runs second stage of Transdecoder predicting the likely coding regions based
  on Pfam and SwissProt hits. Keeps every ORF that Predict accepts, so a
  transcript can give several proteins. Saves the proteins (.pep), the coding
  sequences (.cds), and the ORF positions on the transcripts (.gff3, .bed).

  Predict skips a step when its checkpoint file exists, and the checkpoint name
  does not include --single_best_only. So each ORF set runs on its own copy of
  the LongOrfs folder, without the checkpoints of an earlier Predict run. The
  copy is deleted when the four files are saved.

  Predict reads the ORFs of one transcript in Perl hash order, which is random
  in each run. A fixed hash seed gives the same output order and the same pick
  among tied ORFs in each run.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
    orfs="results/transdecoder/transcriptome.orfs",
    pfam="results/annotation/pfam_orfs.out",
    blastp="results/annotation/sprotblastp_orfs.out"
  output:
    pep="results/transdecoder/transcriptome_all_orfs.pep",
    cds="results/transdecoder/transcriptome_all_orfs.cds",
    gff3="results/transdecoder/transcriptome_all_orfs.gff3",
    bed="results/transdecoder/transcriptome_all_orfs.bed"
  log:
    "logs/transdecoder_predict_all.log"
  conda:
    "envs/transdecoder.yaml"
  params:
    memory=config.get("transdecoder_predict_memory", "6"),
    workdir="results/transdecoder/predict_all"
  threads:
    config.get("transdecoder_predict_threads", 1)
  shell:
    """
    rm -rf {params.workdir} &> {log}
    mkdir -p {params.workdir} &>> {log}
    # Predict and LongOrfs name their files after the file name of the transcriptome
    NAME=$(basename {input.transcriptome})
    cp -rp results/transdecoder/longorfs/$NAME.transdecoder_dir {params.workdir}/ &>> {log}
    rm -rf {params.workdir}/$NAME.transdecoder_dir/__checkpoints_TDpredict &>> {log}
    export PERL_HASH_SEED=0
    # TransDecoder 6 keeps the stage scripts in util/ next to the TransDecoder wrapper, outside bin/
    TD_UTIL=$(dirname $(readlink -f $(command -v TransDecoder)))/util
    $TD_UTIL/TransDecoder.Predict -t {input.transcriptome} --output_dir {params.workdir} --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.pep {output.pep} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.cds {output.cds} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.gff3 {output.gff3} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.bed {output.bed} &>> {log}
    rm -rf {params.workdir} &>> {log}
    """


rule transdecoder_predict_single_best:
  """
  Runs second stage of Transdecoder like transdecoder_predict_all, with
  --single_best_only. Keeps one ORF per transcript, picked by Pfam or SwissProt
  hit first and by length second.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
    orfs="results/transdecoder/transcriptome.orfs",
    pfam="results/annotation/pfam_orfs.out",
    blastp="results/annotation/sprotblastp_orfs.out"
  output:
    pep="results/transdecoder/transcriptome_single_best_orf.pep",
    cds="results/transdecoder/transcriptome_single_best_orf.cds",
    gff3="results/transdecoder/transcriptome_single_best_orf.gff3",
    bed="results/transdecoder/transcriptome_single_best_orf.bed"
  log:
    "logs/transdecoder_predict_single_best.log"
  conda:
    "envs/transdecoder.yaml"
  params:
    memory=config.get("transdecoder_predict_memory", "6"),
    workdir="results/transdecoder/predict_single_best"
  threads:
    config.get("transdecoder_predict_threads", 1)
  shell:
    """
    rm -rf {params.workdir} &> {log}
    mkdir -p {params.workdir} &>> {log}
    # Predict and LongOrfs name their files after the file name of the transcriptome
    NAME=$(basename {input.transcriptome})
    cp -rp results/transdecoder/longorfs/$NAME.transdecoder_dir {params.workdir}/ &>> {log}
    rm -rf {params.workdir}/$NAME.transdecoder_dir/__checkpoints_TDpredict &>> {log}
    export PERL_HASH_SEED=0
    # TransDecoder 6 keeps the stage scripts in util/ next to the TransDecoder wrapper, outside bin/
    TD_UTIL=$(dirname $(readlink -f $(command -v TransDecoder)))/util
    $TD_UTIL/TransDecoder.Predict -t {input.transcriptome} --output_dir {params.workdir} --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} --single_best_only &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.pep {output.pep} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.cds {output.cds} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.gff3 {output.gff3} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.bed {output.bed} &>> {log}
    rm -rf {params.workdir} &>> {log}
    """


rule transdecoder_select:
  """
  Copies the ORF set that transdecoder_orfs in config.yaml names to
  transcriptome.pep, which the downstream rules read.
  """
  input:
    TRANSDECODER_ORFS_FILE
  output:
    "results/transdecoder/transcriptome.pep"
  log:
    "logs/transdecoder_select.log"
  shell:
    """
    cp -p {input} {output} &> {log}
    """


checkpoint align_reads:
  """
  Not run by default.
  Runs Trinity script aligning the reads to the transcriptome using Bowtie2
  and estimating the abundance with RSEM.
  """
  input:
    samples="results/trimmomatic/samples_trimmed.txt",
    transcriptome="results/assembly/transcriptome.fasta",
    gene_trans_map="results/assembly/transcriptome.gene_trans_map"
  output:
    "results/expression/bowtie2_rsem/transcriptome.fasta.bowtie2.ok"
  log:
    "logs/align_reads.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory=config.get("align_reads_memory", "100"),
    scratch_var=SCRATCH_VAR,
    out_dir="results/expression/bowtie2_rsem"
  threads:
    config.get("align_reads_threads", 16)
  # The common failures here are memory and time. A retry asks for the same figures and starts from the
  # beginning, so it would fail the same way, after hours on a large node.
  retries: 0
  shell:
    """
    # align_and_estimate_abundance.pl writes the indexes beside the transcriptome and one folder per
    # replicate into the folder it runs in, so it runs in a folder that holds links to the transcriptome
    # and the gene map, and a samples file with the full read paths. That folder is the tool folder
    # {params.out_dir}, or a folder on the node's local disk when the scratch key names a variable that is
    # set in this job. At the end, everything but the links and the samples file is copied to the tool folder.
    RUN=$(pwd)
    LOG="$RUN/{log}"
    : > "$LOG"
    OUT="$RUN/{params.out_dir}"
    mkdir -p "$OUT"
    SCRATCH_VAR="{params.scratch_var}"
    LOCAL_ROOT=""
    if [ -n "$SCRATCH_VAR" ]; then LOCAL_ROOT="${{!SCRATCH_VAR:-}}"; fi
    # A scratch folder in RAM (tmpfs) counts against the memory of the job, and this rule writes large
    # work files, so it stops there
    if [ -n "$LOCAL_ROOT" ] && [ "$(stat -f -c %T "$LOCAL_ROOT" 2>/dev/null)" = "tmpfs" ]; then
      echo "ERROR: the scratch folder $LOCAL_ROOT is in RAM (tmpfs) on $(hostname). Run this rule on a node with a local disk." | tee -a "$LOG" >&2
      exit 1
    fi
    if [ -n "$LOCAL_ROOT" ]; then
      WORK=$(mktemp -d "$LOCAL_ROOT/transxpress.XXXXXX")
      echo "work folder: $WORK" &>> "$LOG"
      copy_back() {{
        ( cd "$WORK" && find . -mindepth 1 -maxdepth 1 ! -type l ! -name samples_trimmed.txt -exec cp -rp -t "$OUT" {{}} + )
      }}
      # On a failure, the files go to the tool folder, where they would stand without the local disk.
      # A job that Slurm kills at its time limit or with scancel copies nothing.
      on_exit() {{
        rc=$?
        set +e
        trap - EXIT TERM
        cd "$RUN"
        if [ $rc -ne 0 ]; then
          echo "align_reads stopped with code $rc, so the files in $WORK go to {params.out_dir}" &>> "$LOG"
          copy_back &>> "$LOG"
        fi
        rm -rf "$WORK"
        exit $rc
      }}
      trap on_exit EXIT
      trap "exit 143" TERM
    else
      WORK="$OUT"
    fi
    cd "$WORK"
    ln -sfn "$RUN/{input.transcriptome}" transcriptome.fasta
    ln -sfn "$RUN/{input.gene_trans_map}" transcriptome.gene_trans_map
    awk -v run="$RUN" 'BEGIN {{ OFS = "\t" }} {{ for (i = 3; i <= NF; i++) $i = run "/" $i; print }}' "$RUN/{input.samples}" > samples_trimmed.txt

    # Remove old markers so the aligner runs again on the trimmed reads.
    for rep in $(awk '{{print $2}}' samples_trimmed.txt); do rm -f "$rep/bowtie2.bam.ok"; done

    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')


    assembler="{config[assembler]}"
    strand_specific="{config[strand_specific]}"
    
    if [ $assembler = "rnaspades" ]
    then
      if [[ $strand_specific = "--ss rf" ]]
      then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta --SS_lib_type RF --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"
      elif [[ $strand_specific = "--ss fr" ]]
      then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta --SS_lib_type FR --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"
      else
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"
      fi
    else
      $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta {config[strand_specific]} --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"
    fi

    # for_rsem.bam is a large temporary copy of the alignment. Nothing reads it later.
    for rep in $(awk '{{print $2}}' samples_trimmed.txt); do rm -f "$rep/bowtie2.bam.for_rsem.bam"; done

    # the links and the samples file serve only this rule
    rm -f transcriptome.fasta transcriptome.gene_trans_map samples_trimmed.txt
    if [ -n "$LOCAL_ROOT" ]; then
      SIZE=$(du -sh "$WORK" | cut -f1)
      START=$(date +%s)
      copy_back &>> "$LOG"
      echo "copied $SIZE from $WORK to {params.out_dir} in $(( $(date +%s) - START )) s" &>> "$LOG"
    fi
    """

checkpoint prepare_samples_for_IGV:
  """
  Not run by default.
  Prepares the bam alignment files to be used by IGV by sorting
  them into results/expression/bowtie2_rsem/igv and indexing them.
  """
  input:
    alignment="results/expression/bowtie2_rsem/{sample}/bowtie2.bam"
  output:
    sorted_alignment="results/expression/bowtie2_rsem/igv/{sample}.sorted.bam",
    indexed_alignment="results/expression/bowtie2_rsem/igv/{sample}.sorted.bam.bai"
  log:
    "logs/prepare_samples_for_IGV/{sample}.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="10"
  threads:
    4
  shell:
    """
    samtools sort --threads {threads} -m 2G -o {output.sorted_alignment} {input.alignment} &> {log}
    samtools index {output.sorted_alignment} &>> {log}
    """ 

def get_igv(wildcards):
  # sorted, because glob_wildcards returns the folder order, which differs between runs
  indeces = sorted(glob_wildcards("results/expression/bowtie2_rsem/{sample}/bowtie2.bam").sample)
  completed = expand("results/expression/bowtie2_rsem/igv/{sample}.sorted.bam.bai", sample=indeces)
  return completed

rule IGV:
  """
  Not run by default.
  Rule to be called to get all bowtie alignments ready for IGV.
  """
  input:
    get_igv,
    "results/expression/bowtie2_rsem/transcriptome.fasta.bowtie2.ok"
  output:
    "results/expression/bowtie2_rsem/igv/igv.ok"
  log:
    "logs/IGV.log"
  shell:
    """
    touch {output}
    """

rule trinity_DE:
  """
  Runs Trinity script to perform differential expression analysis using edgeR.
  """
  input:
    samples=config["samples_file"],
    expression="results/expression/kallisto/kallisto.gene.counts.matrix"
  output:
    directory("results/expression/edger")
  log:
    "logs/trinity_DE.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')


    num_replicates=`awk '{{print $2}}' {input.samples} | sort | uniq | wc -l` &> {log}
    num_samples=`awk '{{print $1}}' {input.samples} | sort | uniq | wc -l` &>> {log}
    num_replicates_minus_samples=$((num_replicates - num_samples)) &>> {log}
    if [ $num_replicates_minus_samples -gt 1 ] &>> {log}
    then &>> {log}
	    $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input.expression} --method edgeR --output {output} --samples_file {input.samples} &>> {log}
    else &>> {log}
	    echo "No biological replicates to run proper differential expression analysis, last-resorting to edgeR with --dispersion 0.1" &>> {log}
      $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input.expression} --method edgeR --output {output} --samples_file {input.samples} --dispersion {config[dispersion]} &>> {log}
    fi &>> {log}
    """


# a more elegant way to do this is:
# seqkit seq -i {input} | seqkit replace -s -p "\*" -r "" | seqkit split -f -s 1 -O {output}
# however, there seems to be a problem in seqkit: https://github.com/shenwei356/seqkit/issues/65
checkpoint fasta_split_fasta:
  """
  Splits the transcriptome.fasta into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "results/assembly/transcriptome.fasta"
  output:
    directory("results/annotation/chunks/fasta")
  log:
    "logs/fasta_split_fasta.log"
  params:
    memory="2",
    chunk_size=ANNOTATION_CHUNK["fasta"]
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % params.chunk_size == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / params.chunk_size) + 1);
          filename = os.path.join(output[0], fileindex + ".fasta")
          output_handle = open(filename, "w");
        
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        

checkpoint fasta_split_orfs:
  """
  Splits the transcriptome.orfs into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "results/transdecoder/transcriptome.orfs"
  output:
    directory("results/annotation/chunks/orfs")
  log:
    "logs/fasta_split_orfs.log"
  params:
    memory="2",
    chunk_size=ANNOTATION_CHUNK["orfs"]
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % params.chunk_size == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / params.chunk_size) + 1);
          filename = os.path.join(output[0], fileindex + ".orfs")
          output_handle = open(filename, "w");
        
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        
checkpoint fasta_split_pep:
  """
  Splits the transcriptome.pep into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "results/transdecoder/transcriptome.pep"
  output:
    directory("results/annotation/chunks/pep")
  log:
    "logs/fasta_split_pep.log"
  params:
    memory="2",
    chunk_size=ANNOTATION_CHUNK["pep"]
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % params.chunk_size == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / params.chunk_size) + 1);
          filename = os.path.join(output[0], fileindex + ".pep")
          output_handle = open(filename, "w");
        # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan) 
        #if wildcards["extension"] == "pep":
        record.seq = record.seq.strip("*")
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        

def parallel_annotation_tasks_fasta(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_fasta.get(**wildcards).output[0]
  # sorted by number, so the merged file follows the order of the sequences at any chunk size.
  # glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "{index}.fasta")).index, key=int)
  completed_files = expand("results/annotation/chunks/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files

def parallel_annotation_tasks_orfs(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_orfs.get(**wildcards).output[0]
  # sorted by number, so the merged file follows the order of the sequences at any chunk size.
  # glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "{index}.orfs")).index, key=int)
  completed_files = expand("results/annotation/chunks/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files

def parallel_annotation_tasks_pep(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_pep.get(**wildcards).output[0]
  # sorted by number, so the merged file follows the order of the sequences at any chunk size.
  # glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "{index}.pep")).index, key=int)
  completed_files = expand("results/annotation/chunks/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files


rule annotation_merge_fasta:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_fasta
  output:
    "results/annotation/{task}_fasta.out"
  log:
    "logs/annotation_merge_fasta/{task}.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """

rule annotation_merge_orfs:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_orfs
  output:
    "results/annotation/{task}_orfs.out"
  log:
    "logs/annotation_merge_orfs/{task}.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """


rule annotation_merge_pep:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_pep
  output:
    "results/annotation/{task}_pep.out"
  log:
    "logs/annotation_merge_pep/{task}.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """


rule rfam_parallel:
  """
  Runs cmscan on Rfam database on smaller nucleotide files in parallel.
  """
  input:
    fasta="results/annotation/chunks/fasta/{index}.fasta",
    db=DB_DIR + "/Rfam.cm"
  output:
    "results/annotation/chunks/rfam/{index}.out"
  log:
    "logs/rfam_parallel/{index}.log"
  conda:
    "envs/rfam.yaml"
  params:
    memory="5"
  threads:
    4
  shell:
    """
    cmscan -E {config[e_value_threshold]} --rfam --cpu {threads} --tblout {output} {input[db]} {input[fasta]} &> {log}
    """


rule pfam_parallel:
  """
  Runs hmmscan on Pfam database on smaller protein files in parallel.
  """
  input:
    fasta="results/annotation/chunks/orfs/{index}.orfs",
    db=DB_DIR + "/Pfam-A.hmm"
  output:
    "results/annotation/chunks/pfam/{index}.out"
  log:
    "logs/pfam_parallel/{index}.log"
  conda:
    "envs/pfam.yaml"
  params:
    memory="2"
  threads:
    2
  shell:
    """
    # Transdecoder requires --domtblout output
    hmmscan -E {config[e_value_threshold]} --cpu {threads} --domtblout {output} {input[db]} {input[fasta]} &> {log}
    """


rule sprot_blastp_parallel:
  """
  Runs blast search on SwissProt database on smaller protein files in parallel.
  """
  input:
    fasta="results/annotation/chunks/orfs/{index}.orfs",
    db=DB_DIR + "/uniprot_sprot.fasta"
  output:
    "results/annotation/chunks/sprotblastp/{index}.out"
  log:
    "logs/sprot_blastp_parallel/{index}.log"
  conda:
    "envs/blast.yaml"
  params:
    memory="2"
  threads:
    2
  shell:
    """
    blastp -query {input[fasta]} -db {input[db]} -num_threads {threads} -evalue {config[e_value_threshold]} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output} &> {log}
    """


rule sprot_blastx_parallel:
  """
  Runs blast search on SwissProt database on smaller nucleotide files in parallel.
  """
  input:
    fasta="results/annotation/chunks/fasta/{index}.fasta",
    db=DB_DIR + "/uniprot_sprot.fasta"
  output:
    "results/annotation/chunks/sprotblastx/{index}.out"
  log:
    "logs/sprot_blastx_parallel/{index}.log"
  conda:
    "envs/blast.yaml"
  params:
    memory="2"
  threads:
    4
  shell:
    """
    blastx -query {input[fasta]} -db {input[db]} -num_threads {threads} -evalue {config[e_value_threshold]} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output} &> {log}
    """


rule tmhmm_parallel:
  """
  Runs tmhmm.py on smaller protein files in parallel to predict transmembrane
  topology.
  """
  input:
    "results/annotation/chunks/pep/{index}.pep"
  output:
    "results/annotation/chunks/tmhmm/{index}.out"
  log:
    "logs/tmhmm_parallel/{index}.log"
  params:
    memory="2"
  threads:
    1
  run:
    # open the input file, output file and log file
    # tmhmm.py writes its result files into the current directory, so each chunk works in its own temporary directory
    # in the job's temporary folder (see job_tmpdir). The directory and the per-protein files in it are deleted when
    # the chunk is done
    with open(input[0], "r") as input_handle, open(output[0], "a+") as output_handle, open(log[0], "a+") as log_handle, tempfile.TemporaryDirectory(dir=job_tmpdir(wildcards)) as tmp_dir:
      log_handle.write("work folder: " + tmp_dir + "\n")
      log_handle.flush()

      # iterate through individual sequences in input file, tmhmm.py can be executed with just one sequence at a time
      for record in Bio.SeqIO.parse(input_handle, "fasta"):

        # individual fasta headers must have an ID and description, separated by a space, otherwise error occurs
        # see: https://github.com/dansondergaard/tmhmm.py/issues/16
        # adding artificial "description" in header
        record.description = " x"

        # tmhmm.py does not handle 'X' in the sequence (happens when using RNASpades as an assembler and 'N' occurs in the transcript sequence)
        # see: https://github.com/dansondergaard/tmhmm.py/issues/9
        # removing 'X' from the sequence
        record.seq = Seq(re.sub('X',"",str(record.seq)))

        # saving this fasta file
        fasta_file =  os.path.join(tmp_dir, record.id + ".fasta")
        Bio.SeqIO.write(record, fasta_file, "fasta")

        # executing the tmhmm.py on created fasta
        cmd_tmhmm = ['tmhmm', '-f', fasta_file]
        subprocess.run(cmd_tmhmm, stdout = log_handle, stderr = log_handle, cwd = tmp_dir)

        # 3 files are created afterwards in tmp_dir
        # record.id.summary     record.id.annotation    record.id.plot

        summary_file = os.path.join(tmp_dir, record.id + ".summary")
        annotation_file = os.path.join(tmp_dir, record.id + ".annotation")
        plot_file = os.path.join(tmp_dir, record.id + ".plot")

        # creating shorter output from summary file
        with open(summary_file, "r") as inp:
          final_topology = ''
          num_of_helices = 0

          for line in inp:
            # split the line of format "start end topology" into list [start, end, topology]
            split_line = re.split(r'\s+', line, 2)
            start = split_line[0].strip()
            end = split_line[1].strip()
            topology = split_line[2].strip()

            # update final string based on the topology
            if (topology == 'inside'):
              final_topology += 'i'
            if (topology == 'outside'):
              final_topology += 'o'
            if ('transmembrane helix' in topology):
              final_topology += start + '-' + end
              num_of_helices += 1

          output_handle.write(record.id + "\tPredHel=" + str(num_of_helices) + "\tTopology=" + final_topology + "\n")
        os.remove(summary_file)
        os.remove(annotation_file)
        os.remove(plot_file)

rule signalp_parallel:
  """
  Runs signalp on smaller protein files in parallel to predict signal peptides.
  """
  input:
    "results/annotation/chunks/pep/{index}.pep"
  output:
    "results/annotation/chunks/signalp/{index}.out"
  log:
    "logs/signalp_parallel/{index}.log"
  conda:
    "signalp"
  params:
    memory="4",
    organism=SIGNALP_ORGANISM
  threads:
    1
  shell:
    """
    # Snakemake sets TMPDIR to the job's temporary folder, see job_tmpdir
    TMP=$(mktemp -d "${{TMPDIR:-/tmp}}/signalp.XXXXXX")
    echo "work folder: $TMP" &> {log}
    signalp6 --fastafile {input} --organism {params.organism} --output_dir $TMP/out --format none --mode fast &>> {log}
    mv $TMP/out/prediction_results.txt {output}
    rm -r $TMP
    """

rule targetp_parallel:
  """
  Runs TargetP on smaller protein files to predict presence and type of
  targeting peptide.
  """
  input:
    "results/annotation/chunks/pep/{index}.pep"
  output:
    "results/annotation/chunks/targetp/{index}.out"
  log:
    "logs/targetp_parallel/{index}.log"
  conda:
    "targetp"
  params:
    memory="3"
  threads:
    1
  shell:
    """
    # targetp writes its result into the folder it runs in, which is the job's temporary folder
    # (Snakemake sets TMPDIR, see job_tmpdir)
    RUN=$(pwd)
    TMP=$(mktemp -d "${{TMPDIR:-/tmp}}/targetp.XXXXXX")
    echo "work folder: $TMP" &> "$RUN/{log}"
    cd $TMP
    targetp -fasta "$RUN/{input}" -format short -org {config[targetp]} -prefix targetp_{wildcards.index} &>> "$RUN/{log}"
    mv targetp_{wildcards.index}_summary.targetp2 "$RUN/{output}" &>> "$RUN/{log}"
    cd "$RUN"
    rm -r $TMP
    """


rule kallisto:
  """
  Runs Trinity script to perform transcript expression quantification 
  using Kallisto.
  """
  input:
    samples="results/trimmomatic/samples_trimmed.txt",
    transcriptome="results/assembly/transcriptome.fasta",
    gene_trans_map="results/assembly/transcriptome.gene_trans_map"
  output:
    "results/expression/transcriptome_expression_isoform.tsv",
    "results/expression/transcriptome_expression_gene.tsv",
    "results/expression/kallisto/kallisto.gene.counts.matrix"
  log:
    "logs/kallisto.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory=config.get("kallisto_memory", "20"),
    scratch_var=SCRATCH_VAR,
    out_dir="results/expression/kallisto"
  threads:
    config.get("kallisto_threads", 16)
  # The common failures here are memory and time. A retry asks for the same figures and starts from the
  # beginning, so it would fail the same way, after hours on a large node.
  retries: 0
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')
    # abundance_estimates_to_matrix.pl writes the matrix rows in Perl hash order, which is random
    # in each run. A fixed Perl hash seed gives the same row order in each run.
    export PERL_HASH_SEED=0

    # align_and_estimate_abundance.pl writes the indexes beside the transcriptome and one folder per
    # replicate into the folder it runs in, so it runs in a folder that holds links to the transcriptome
    # and the gene map, and a samples file with the full read paths. That folder is the tool folder
    # {params.out_dir}, or a folder on the node's local disk when the scratch key names a variable that is
    # set in this job. At the end, everything but the links and the samples file is copied to the tool folder.
    RUN=$(pwd)
    LOG="$RUN/{log}"
    : > "$LOG"
    OUT="$RUN/{params.out_dir}"
    mkdir -p "$OUT"
    SCRATCH_VAR="{params.scratch_var}"
    LOCAL_ROOT=""
    if [ -n "$SCRATCH_VAR" ]; then LOCAL_ROOT="${{!SCRATCH_VAR:-}}"; fi
    # A scratch folder in RAM (tmpfs) counts against the memory of the job, and this rule writes large
    # work files, so it stops there
    if [ -n "$LOCAL_ROOT" ] && [ "$(stat -f -c %T "$LOCAL_ROOT" 2>/dev/null)" = "tmpfs" ]; then
      echo "ERROR: the scratch folder $LOCAL_ROOT is in RAM (tmpfs) on $(hostname). Run this rule on a node with a local disk." | tee -a "$LOG" >&2
      exit 1
    fi
    if [ -n "$LOCAL_ROOT" ]; then
      WORK=$(mktemp -d "$LOCAL_ROOT/transxpress.XXXXXX")
      echo "work folder: $WORK" &>> "$LOG"
      copy_back() {{
        ( cd "$WORK" && find . -mindepth 1 -maxdepth 1 ! -type l ! -name samples_trimmed.txt -exec cp -rp -t "$OUT" {{}} + )
      }}
      # On a failure, the files go to the tool folder, where they would stand without the local disk.
      # A job that Slurm kills at its time limit or with scancel copies nothing.
      on_exit() {{
        rc=$?
        set +e
        trap - EXIT TERM
        cd "$RUN"
        if [ $rc -ne 0 ]; then
          echo "kallisto stopped with code $rc, so the files in $WORK go to {params.out_dir}" &>> "$LOG"
          copy_back &>> "$LOG"
        fi
        rm -rf "$WORK"
        exit $rc
      }}
      trap on_exit EXIT
      trap "exit 143" TERM
    else
      WORK="$OUT"
    fi
    cd "$WORK"
    ln -sfn "$RUN/{input.transcriptome}" transcriptome.fasta
    ln -sfn "$RUN/{input.gene_trans_map}" transcriptome.gene_trans_map
    awk -v run="$RUN" 'BEGIN {{ OFS = "\t" }} {{ for (i = 3; i <= NF; i++) $i = run "/" $i; print }}' "$RUN/{input.samples}" > samples_trimmed.txt

    assembler="{config[assembler]}"
    strand_specific="{config[strand_specific]}"

    if [ $assembler = "rnaspades" ]
    then
      if [[ $strand_specific = "--ss rf" ]]
      then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta --SS_lib_type RF --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"
      elif [[ $strand_specific = "--ss fr" ]]
      then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta --SS_lib_type FR --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"
      else
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"
      fi
    else
      $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta {config[strand_specific]} --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"
    fi
    
    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map transcriptome.gene_trans_map */abundance.tsv &>> "$LOG"
    if [ -f kallisto.isoform.TMM.EXPR.matrix ]; then
      cp -p kallisto.isoform.TMM.EXPR.matrix "$RUN/{output[0]}" &>> "$LOG"
    elif [ -f kallisto.isoform.TPM.not_cross_norm ]; then
      cp -p kallisto.isoform.TPM.not_cross_norm "$RUN/{output[0]}" &>> "$LOG"
    else
      echo Neither kallisto.isoform.TMM.EXPR.matrix or kallisto.isoform.TPM.not_cross_norm were produced
      exit 1
    fi
    if [ -f kallisto.gene.TMM.EXPR.matrix ]; then
      cp -p kallisto.gene.TMM.EXPR.matrix "$RUN/{output[1]}" &>> "$LOG"
    elif [ -f kallisto.gene.TPM.not_cross_norm ]; then
      cp -p kallisto.gene.TPM.not_cross_norm "$RUN/{output[1]}" &>> "$LOG"
    else
      echo Neither kallisto.gene.TMM.EXPR.matrix or kallisto.gene.TPM.not_cross_norm were produced
      exit 1
    fi

    # the links and the samples file serve only this rule
    rm -f transcriptome.fasta transcriptome.gene_trans_map samples_trimmed.txt
    if [ -n "$LOCAL_ROOT" ]; then
      SIZE=$(du -sh "$WORK" | cut -f1)
      START=$(date +%s)
      copy_back &>> "$LOG"
      echo "copied $SIZE from $WORK to {params.out_dir} in $(( $(date +%s) - START )) s" &>> "$LOG"
    fi
    """

 
# The download rules exist only when db_dir is "resources/db", so a run never writes to a shared
# database folder. A shared folder is filled once with the download_databases target.
if not DB_SHARED:

  localrules: download_databases

  rule download_sprot:
    """
    Downloads the Swiss-Prot release named by sprot_release, checks its md5 and prepares it for BLAST.
    UniProt serves only the newest release as a plain fasta file. An older release comes from the
    release archive, which holds the same fasta file.
    """
    output:
      fasta="resources/db/uniprot_sprot.fasta",
      version="resources/db/uniprot_sprot.version"
    log:
      "logs/download_sprot.log"
    conda:
      "envs/blast.yaml"
    params:
      memory="2",
      release=SPROT_RELEASE,
      current_url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete",
      archive_url="https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-" + SPROT_RELEASE + "/knowledgebase"
    threads:
      1
    shell:
      """
      db=$(dirname {output.fasta})

      # the release that UniProt serves as current
      current_release() {{
        wget -q -O - "{params.current_url}/reldate.txt" | awk '$1 == "UniProtKB/Swiss-Prot" && $2 == "Release" {{ print $3 }}'
      }}

      # the md5 that the RELEASE.metalink file in folder $1 gives for file $2
      metalink_md5() {{
        wget -q -O - "$1/RELEASE.metalink" | awk -v name="$2" '/<file name=/ {{ f = index($0, "\\"" name "\\"") > 0 }} f && !found && /type="md5"/ {{ split($0, a, /[<>]/); print a[3]; found = 1 }}'
      }}

      before=$(current_release)
      echo "Swiss-Prot release asked: {params.release}, current release: $before" > {log}

      if [ "$before" = "{params.release}" ]; then
        download="$db/uniprot_sprot.fasta.gz"
        wget -O "$download" "{params.current_url}/uniprot_sprot.fasta.gz" &>> {log}
        expected=$(metalink_md5 "{params.current_url}" uniprot_sprot.fasta.gz)
      else
        archive="uniprot_sprot-only{params.release}.tar.gz"
        download="$db/$archive"
        if ! wget -O "$download" "{params.archive_url}/$archive" &>> {log}; then
          echo "Swiss-Prot {params.release} is not the current release ($before) and has no archive at {params.archive_url}" >> {log}
          exit 1
        fi
        expected=$(metalink_md5 "{params.archive_url}" "$archive")
      fi

      actual=$(md5sum "$download" | cut -d " " -f 1)
      echo "md5 of $download: expected $expected, actual $actual" >> {log}
      if [ -z "$expected" ] || [ "$expected" != "$actual" ]; then
        echo "The md5 does not match RELEASE.metalink" >> {log}
        exit 1
      fi

      if [ "$before" = "{params.release}" ]; then
        # UniProt could switch to the next release while the file downloads
        after=$(current_release)
        if [ "$after" != "$before" ]; then
          echo "The current release changed from $before to $after during the download. Run the rule again." >> {log}
          exit 1
        fi
      else
        tar -xzf "$download" -C "$db" uniprot_sprot.fasta.gz &>> {log}
        rm "$download"
      fi

      gunzip -f "$db/uniprot_sprot.fasta.gz" &>> {log}
      makeblastdb -in {output.fasta} -dbtype prot &>> {log}
      printf "Swiss-Prot\\t%s\\n" "{params.release}" > {output.version}
      """


  rule download_pfam:
    """
    Downloads the Pfam release named by pfam_release, checks its md5 and prepares it for hmmscan.
    """
    output:
      hmm="resources/db/Pfam-A.hmm",
      version="resources/db/Pfam-A.version"
    log:
      "logs/download_pfam.log"
    conda:
      "envs/pfam.yaml"
    params:
      memory="2",
      release=PFAM_RELEASE,
      url="https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam" + PFAM_RELEASE
    threads:
      1
    shell:
      """
      wget -O {output.hmm}.gz "{params.url}/Pfam-A.hmm.gz" &> {log}
      expected=$(wget -q -O - "{params.url}/md5_checksums" | awk '$2 == "Pfam-A.hmm.gz" {{ print $1 }}')
      actual=$(md5sum {output.hmm}.gz | cut -d " " -f 1)
      echo "md5 of Pfam-A.hmm.gz: expected $expected, actual $actual" >> {log}
      if [ -z "$expected" ] || [ "$expected" != "$actual" ]; then
        echo "The md5 does not match md5_checksums of Pfam {params.release}" >> {log}
        exit 1
      fi
      gunzip -f {output.hmm}.gz &>> {log}
      # -f replaces index files left by an earlier attempt
      hmmpress -f {output.hmm} &>> {log}
      printf "Pfam\\t%s\\n" "{params.release}" > {output.version}
      """


  rule download_rfam:
    """
    Downloads the Rfam release named by rfam_release and prepares it for cmscan.
    Rfam publishes no checksum, so the file is not checked.
    """
    output:
      cm="resources/db/Rfam.cm",
      version="resources/db/Rfam.version"
    log:
      "logs/download_rfam.log"
    conda:
      "envs/rfam.yaml"
    params:
      memory="2",
      release=RFAM_RELEASE,
      url="https://ftp.ebi.ac.uk/pub/databases/Rfam/" + RFAM_RELEASE
    threads:
      1
    shell:
      """
      wget -O {output.cm}.gz "{params.url}/Rfam.cm.gz" &> {log}
      gunzip -f {output.cm}.gz &>> {log}
      # -F replaces index files left by an earlier attempt
      cmpress -F {output.cm} &>> {log}
      printf "Rfam\\t%s\\n" "{params.release}" > {output.version}
      """


  rule download_busco_lineage:
    """
    Downloads the BUSCO lineage named by lineage into resources/db/busco_downloads, for a shared database
    folder. A run with db_dir "resources/db" does not read this copy: BUSCO downloads its own lineage
    into resources/busco_downloads.
    """
    output:
      directory("resources/db/busco_downloads/lineages/" + LINEAGE)
    log:
      "logs/download_busco_lineage.log"
    conda:
      "envs/busco.yaml"
    params:
      memory="2",
      lineage=LINEAGE
    threads:
      1
    shell:
      """
      # --download has to come last on the busco command line
      busco --download_path resources/db/busco_downloads --download {params.lineage} &> {log}
      """


  rule download_databases:
    """
    Downloads the annotation databases and the BUSCO lineage into resources/db/ and lists their releases
    in resources/db/versions.txt. Run it once in an empty folder to fill a shared database folder.
    """
    input:
      versions=["resources/db/uniprot_sprot.version", "resources/db/Pfam-A.version", "resources/db/Rfam.version"],
      lineage=["resources/db/busco_downloads/lineages/" + LINEAGE] if LINEAGE else []
    output:
      "resources/db/versions.txt"
    params:
      memory="2"
    threads:
      1
    run:
      with open(output[0], "w") as out:
        for path in input.versions:
          with open(path) as handle:
            out.write(handle.read())
        for path in input.lineage:
          fields = {}
          with open(os.path.join(path, "dataset.cfg")) as handle:
            for line in handle:
              if "=" in line:
                key, value = line.strip().split("=", 1)
                fields[key] = value
          out.write("BUSCO\t%s\t%s\n" % (fields.get("name", os.path.basename(path)), fields.get("creation_date", "unknown")))


rule annotated_fasta:
  """
  Puts annotations in the headers of transcripts/proteins in the 
  .fasta/.pep transcriptome files.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
    proteome="results/transdecoder/transcriptome.pep",
    expression="results/expression/transcriptome_expression_isoform.tsv",
    blastx_results="results/annotation/sprotblastx_fasta.out",
    rfam_results="results/annotation/rfam_fasta.out",
    blastp_results="results/annotation/sprotblastp_orfs.out",
    pfam_results="results/annotation/pfam_orfs.out",
    tmhmm_results="results/annotation/tmhmm_pep.out",
    signalp_results="results/annotation/signalp_pep.out",
    targetp_results="results/annotation/targetp_pep.out"
  output:
    transcriptome_annotated="results/annotation/transcriptome_annotated.fasta",
    proteome_annotated="results/annotation/transcriptome_annotated.pep",
    tpm_blast_table="results/annotation/transcriptome_TPM_blast.csv"
  log:
    "logs/annotated_fasta.log"
  params:
    memory="4"
  threads:
    1
  run:
    # Open log file
    with open(log[0], "w") as log_handle:
      ## Annotation map: transcript id -> description
      expression_annotations = {}
      blastx_annotations = {}
      blastp_annotations = {}
      pfam_annotations = {}
      rfam_annotations = {}
      tmhmm_annotations = {}
      signalp_annotations = {}
      targetp_annotations = {}
  
      ## Load kallisto results
      print ("Loading expression values from", input["expression"], file=log_handle)
      with open(input["expression"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter='\t')
        columns = next(csv_reader)
        for row in csv_reader:
          expression_annotations[row[0]] = columns[1] + "=" + str(row[1])
          for i in range(2, len(columns)):
            expression_annotations[row[0]] += " " + columns[i] + "=" + str(row[i])

      ## Load blastx results
      print ("Loading blastx results from", input["blastx_results"], file=log_handle)
      with open(input["blastx_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 13): continue
          blastx_annotations[row[0]] = row[12] + " E=" + str(row[10])
  
      ## Load blastp results
      print ("Loading blastp results from", input["blastp_results"], file=log_handle)
      with open(input["blastp_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 13): continue
          blastp_annotations[row[0]] = row[12] + " E=" + str(row[10])

      ## Load pfam results
      print ("Loading pfam predictions from", input["pfam_results"], file=log_handle)
      with open(input["pfam_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = re.split(" +", line, 22)
          if (len(row) < 23): continue
          if row[3] not in pfam_annotations:
            pfam_annotations[row[3]] = row[1] + " " + row[22] + "E=" + str(row[6])

      ## Load rfam results
      print ("Loading rfam predictions from", input["rfam_results"], file=log_handle)
      with open(input["rfam_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = re.split(" +", line, 17)
          if (len(row) < 18): continue
          rfam_annotations[row[2]] = row[1] + " " + row[17] + "E=" + str(row[15])
  
      ## Load tmhmm results
      print ("Loading tmhmm predictions from", input["tmhmm_results"], file=log_handle)
      with open(input["tmhmm_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) > 2):
            tmhmm_annotations[row[0]] = row[1] + " " + row[2]

      ## Load signalp results
      print ("Loading signalp predictions from", input["signalp_results"], file=log_handle)
      with open(input["signalp_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 2): continue
          prediction = str(row[1])
          if prediction == "OTHER":
            signalp_annotations[row[0]] = "OTHER (No SP)"
          elif prediction == "SP":
            signalp_annotations[row[0]] = "standard secretory signal peptides (Sec/SPI)"
          elif prediction == "LIPO":
            signalp_annotations[row[0]] = "lipoprotein signal peptides (Sec/SPII)"
          elif prediction == "TAT":
            signalp_annotations[row[0]] = "Tat signal peptides (Tat/SPI)"
          elif prediction == "TATLIPO":
            signalp_annotations[row[0]] = "Tat lipoprotein signal peptides (Tat/SPII)"
          elif prediction == "PILIN":
            signalp_annotations[row[0]] = "Pilin and pilin-like signal peptides (Sec/SPIII)"
          else:
            signalp_annotations[row[0]] = str(row[1])

      ## Load targetp results
      translation =    {"SP": "Signal peptide",
                        "mTP": "Mitochondrial transit peptide",
                        "cTP": "chloroplast transit peptide",
                        "luTP": "thylakoidal lumen composite transit peptide",
                        "noTP": "no targeting peptide" }
      print("Loading targetp predictions from", input["targetp_results"], file=log_handle)
      with open(input["targetp_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = line.split()
          if (len(row) >= 8):
            targetp_annotations[row[0]] = translation[row[1]] + ", " + " ".join(row[7:])
          elif (len(row) >= 2):
            targetp_annotations[row[0]] = translation[row[1]]
      
      ## Do the work
      print ("Annotating FASTA file", input["transcriptome"], "to", output["transcriptome_annotated"], file=log_handle)
      with open(input["transcriptome"], "r") as input_fasta_handle, open(output["transcriptome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = record.id
          record.description = "TPM: " + expression_annotations.get(transcript_id)
          if transcript_id in blastx_annotations:
            record.description += "; blastx: " + blastx_annotations.get(transcript_id)
          if transcript_id in rfam_annotations:
            record.description += "; rfam: " + rfam_annotations.get(transcript_id)
          Bio.SeqIO.write(record, output_fasta_handle, "fasta")
      
      print ("Annotating FASTA file", input["proteome"], "to", output["proteome_annotated"], file=log_handle)
      with open(input["proteome"], "r") as input_fasta_handle, open(output["proteome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = re.sub("\\.p[0-9]+$", "", record.id)
          record.description = "transdecoder: " + re.search("ORF type:([^,]+,score=[^,]+)", record.description).group(1)
          if transcript_id in expression_annotations:
            record.description += "; TPM: " + expression_annotations.get(transcript_id)
          if record.id in blastp_annotations:
            record.description += "; blastp: " + blastp_annotations.get(record.id)
          if record.id in pfam_annotations:
            record.description += "; pfam: " + pfam_annotations.get(record.id)
          if transcript_id in rfam_annotations:
            record.description += "; rfam: " + rfam_annotations.get(transcript_id)
          if record.id in tmhmm_annotations:
            record.description += "; tmhmm: " + tmhmm_annotations.get(record.id)
          if record.id in signalp_annotations:
            record.description += "; signalp: " + signalp_annotations.get(record.id)
          if record.id in targetp_annotations:
            record.description += "; targetp: " + targetp_annotations.get(record.id)
          # Add sequence ID prefix from configuration
          if config["annotated_fasta_prefix"]:
            record.id = config["annotated_fasta_prefix"] + "|" + record.id
          Bio.SeqIO.write(record, output_fasta_handle, "fasta")

      print ("Generating transcriptome_TPM_blast.csv table", file=log_handle)
      with open(input["expression"], "r") as input_csv_handle, open(output["tpm_blast_table"], "w") as output_csv_handle:
        csv_reader = csv.reader(input_csv_handle, delimiter="\t")
        csv_writer = csv.writer(output_csv_handle, delimiter=",")
        csv_columns = next(csv_reader)
        csv_columns[0] = "transcript"
        for i in range(1, len(csv_columns)):
          csv_columns[i] = "TPM(" + csv_columns[i] + ")"
        csv_columns.append("blastx")
        csv_writer.writerow(csv_columns)
        for row in csv_reader:
          row.append(blastx_annotations.get(row[0], ""))
          csv_writer.writerow(row)


# Every job gets its temporary folder (TMPDIR, TMP and TEMP) from job_tmpdir, so the temporary files of
# the tools go to the node's local disk when the scratch key is set. Snakemake evaluates the tmpdir
# resource inside the job, where the variable is set.
for workflow_rule in workflow.rules:
  workflow_rule.resources["tmpdir"] = job_tmpdir
