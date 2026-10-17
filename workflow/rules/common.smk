import functools
import os
import sys
import shutil
import re
import csv
import Bio.SeqIO
import subprocess
import tempfile
from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version, validate
from Bio.Seq import Seq
import pandas as pd

min_version("9.27.0")

configfile: "config.yaml"

# tmhmm.py 1.2.1 imports numpy in its setup.py, so pip builds it in envs/tmhmm.yaml against the numpy of
# that env, without an isolated build env. Snakemake passes this variable to conda when it builds the env.
os.environ.setdefault("PIP_NO_BUILD_ISOLATION", "0")

validate(config, os.path.join(workflow.basedir, "schemas", "config.schema.yaml"), set_default=True)

def gb_to_mb(gb):
  """
  The mem_mb resource for a memory figure of config.yaml, which is in GB. Slurm reads the --mem value
  in MB of 1024 KB, so 1 GB is 1024 MB here.
  """
  return int(float(gb) * 1024)

# The releases of the annotation databases. The download rules fetch exactly these releases.
PFAM_RELEASE = config["pfam_release"]
RFAM_RELEASE = config["rfam_release"]
SPROT_RELEASE = config["sprot_release"]

# The adapter file that ships with transXpress. {adapters} in trimmomatic_parameters stands for its path.
ADAPTERS = os.path.normpath(os.path.join(workflow.basedir, "..", "resources", "adapters.fasta"))

def trimmomatic_parameters(wildcards):
  return str(config["trimmomatic_parameters"]).replace("{adapters}", ADAPTERS)

# The ORF set that transcriptome.pep holds and the downstream rules read. Both sets are always saved.
TRANSDECODER_ORFS = config["transdecoder_orfs"]
TRANSDECODER_ORFS_FILE = {"all": "results/transdecoder/transcriptome_all_orfs.pep", "single_best": "results/transdecoder/transcriptome_single_best_orf.pep"}[TRANSDECODER_ORFS]

# The rnaSPAdes transcript set that becomes transcriptome.fasta. rnaSPAdes always writes all three.
RNASPADES_FILTER = config["rnaspades_filter"]

# The organism groups of SignalP and TargetP. "" turns the tool off.
# A key with no value reaches here as None, which means "" too
SIGNALP_ORGANISM = config["signalp"] or ""
TARGETP_ORGANISM = config["targetp"] or ""

# The number of sequences in each annotation chunk. Each chunk is one cluster job for each tool that reads it.
ANNOTATION_CHUNK = {split: config["annotation_chunk_" + split] for split in ["fasta", "orfs", "pep"]}

# The environment variable that holds a node-local folder in a cluster job, for example "$SCRATCH".
# "" runs every rule in the run folder. The heavy rules read the variable with bash, inside the job,
# because only the job knows its value.
SCRATCH = config["scratch"] or ""
SCRATCH_VAR = SCRATCH[1:]
# The shell functions that the heavy rules use to work in that folder
SCRATCH_SH = os.path.join(workflow.basedir, "scripts", "scratch.sh")

def job_tmpdir(wildcards):
  """
  The temporary folder of a job: the folder in the variable that the scratch key names, or the
  system default when the key is "" or the variable is not set. Snakemake calls this inside the job.
  """
  return (os.environ.get(SCRATCH_VAR, "") if SCRATCH_VAR else "") or tempfile.gettempdir()

# The folder of the annotation databases. "resources/db" keeps one copy in each run folder, which the
# pipeline downloads. Any other folder is a shared copy that the pipeline only reads.
DB_DIR = config["db_dir"].rstrip("/")
DB_SHARED = DB_DIR != "resources/db"
LINEAGE = config["lineage"] or ""

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
localrules: all, clean, delete_intermediates, run_info, trimmomatic_split, trimmomatic_merge, samples_yaml_conversion, transcriptome_copy, transdecoder_select, compare_qc_after_trim, fastqc_before_trim_warnings, IGV

# A rule with a wildcard keeps its logs in logs/<rule>/, and the wildcard is only in the file name.
# The cluster writes its own log next to the rule log before the job starts, and not every cluster
# makes a missing folder, so every log folder is made here.
onstart:
  for workflow_rule in workflow.rules:
    for log_file in workflow_rule.log:
      log_folder = os.path.dirname(str(log_file))
      if log_folder:
        os.makedirs(log_folder, exist_ok=True)

rule clean:
  """
  Not run by default. Deletes the results and the logs of the run, so the next run starts from the
  beginning. Keeps config.yaml, the samples file, resources/ and .snakemake/.
  """
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  shell:
    """
    rm -rf results logs
    """


rule delete_intermediates:
  """
  Not run by default. Deletes the large working files that no rule names as its output, so every
  result stays and a later snakemake run in this folder has nothing to do: the trimmed reads and the
  kallisto, bowtie2 and RSEM indexes.
  """
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  shell:
    """
    rm -f results/trimmomatic/reads/*.fastq.gz
    rm -f results/expression/kallisto/transcriptome.fasta.kallisto_idx
    # the bowtie2.ok marker goes with the index, so a later align_reads run builds the index again
    rm -f results/expression/bowtie2_rsem/transcriptome.fasta.bowtie2.* results/expression/bowtie2_rsem/transcriptome.fasta.RSEM.*
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
    basedir=workflow.basedir
  threads:
    1
  resources:
    mem_mb=2048,
    runtime="4h"
  script:
    "../scripts/run_info.py"
