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

# The adapter file that ships with transXpress. {adapters} in trimmomatic_parameters stands for its path.
ADAPTERS = os.path.normpath(os.path.join(workflow.basedir, "..", "resources", "adapters.fasta"))
TRIMMOMATIC_PARAMETERS = str(config["trimmomatic_parameters"]).replace("{adapters}", ADAPTERS)

# The ORF set that transcriptome.pep holds and the downstream rules read. Both sets are always saved.
TRANSDECODER_ORFS = str(config.get("transdecoder_orfs", "all"))
if TRANSDECODER_ORFS not in ["all", "single_best"]:
  raise ValueError("transdecoder_orfs must be \"all\" or \"single_best\". The config gives %s" % TRANSDECODER_ORFS)
TRANSDECODER_ORFS_FILE = {"all": "results/transdecoder/transcriptome_all_orfs.pep", "single_best": "results/transdecoder/transcriptome_single_best_orf.pep"}[TRANSDECODER_ORFS]

# The rnaSPAdes transcript set that becomes transcriptome.fasta. rnaSPAdes always writes all three.
RNASPADES_FILTER = str(config.get("rnaspades_filter", "default"))
if RNASPADES_FILTER not in ["default", "hard", "soft"]:
  raise ValueError("rnaspades_filter must be \"default\", \"hard\" or \"soft\". The config gives %s" % RNASPADES_FILTER)

# The organism group that SignalP predicts signal peptides for. "" turns SignalP off.
# A key with no value reaches here as None, which means "" too
SIGNALP_ORGANISM = str(config.get("signalp", "eukarya") or "")
if SIGNALP_ORGANISM not in ["", "eukarya", "other"]:
  raise ValueError("signalp must be \"eukarya\", \"other\" or \"\". The config gives %s" % SIGNALP_ORGANISM)

# The organism group that TargetP predicts targeting peptides for. "" turns TargetP off.
TARGETP_ORGANISM = str(config.get("targetp", "non-pl") or "")
if TARGETP_ORGANISM not in ["", "pl", "non-pl"]:
  raise ValueError("targetp must be \"pl\", \"non-pl\" or \"\". The config gives %s" % TARGETP_ORGANISM)

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
  result stays and a later snakemake run in this folder has nothing to do: the trimmed reads and the
  kallisto, bowtie2 and RSEM indexes.
  """
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
    memory="2",
    basedir=workflow.basedir
  threads:
    1
  script:
    "../scripts/run_info.py"
