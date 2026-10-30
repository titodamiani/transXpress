import functools
import hashlib
import os
import tempfile
from snakemake.exceptions import WorkflowError
from snakemake.settings.types import DeploymentMethod
from snakemake.utils import min_version, validate

min_version("9.27.0")

configfile: "config.yaml"

# tmhmm.py 1.2.1 imports numpy in its setup.py, so pip must build it without an isolated build env.
os.environ.setdefault("PIP_NO_BUILD_ISOLATION", "0")

validate(config, os.path.join(workflow.basedir, "schemas", "config.schema.yaml"), set_default=True)

def gb_to_mb(gb):
  """
  The mem_mb resource for a memory figure of config.yaml in GB (1 GB is 1024 MB, as Slurm reads --mem).
  """
  return int(float(gb) * 1024)

PFAM_RELEASE = config["pfam_release"]
RFAM_RELEASE = config["rfam_release"]
SPROT_RELEASE = config["sprot_release"]

# The TMM-scaled TPM of each transcript, which ExN50 and the annotated FASTA headers read.
TMM_ISOFORM = "results/summary/expr_tmm_isoform.tsv"

RNASPADES_FILTER = config["rnaspades_filter"]

# Organism groups of SignalP and TargetP. "" or an empty key turns the tool off.
SIGNALP_ORGANISM = config["signalp"] or ""
TARGETP_ORGANISM = config["targetp"] or ""

# The image of SignalP or TargetP, "" for the conda env. A rule with an image has no conda env.
def tool_container(key, organism):
  path = config[key] or ""
  if not path or not organism:
    return ""
  path = os.path.abspath(path)
  if not os.path.isfile(path):
    raise ValueError("%s: no such file: %s" % (key, path))
  if DeploymentMethod.APPTAINER not in workflow.deployment_settings.deployment_method:
    raise ValueError("%s is set, so the run needs --software-deployment-method conda apptainer" % key)
  return path

SIGNALP_CONTAINER = tool_container("signalp_container", SIGNALP_ORGANISM)
TARGETP_CONTAINER = tool_container("targetp_container", TARGETP_ORGANISM)

ANNOTATION_CHUNK = {split: config[f"annotation_chunk_{split}"] for split in ["fasta", "orfs", "pep"]}

# The environment variable that holds the node-local folder. The heavy rules read it with bash inside the job.
SCRATCH = config["scratch"] or ""
SCRATCH_VAR = SCRATCH.lstrip("$")

def job_tmpdir(wildcards):
  """
  The temporary folder of a job: the folder in the scratch variable, or the system default. Snakemake calls this inside the job.
  """
  return os.environ.get(SCRATCH_VAR) or tempfile.gettempdir()

# "resources/db" is downloaded into the run folder. Any other db_dir is a shared copy that is only read.
DB_DIR = config["db_dir"].rstrip("/")
DB_SHARED = DB_DIR != "resources/db"
LINEAGE = config["lineage"] or ""

# Annotation starts hours into a run, so a shared folder (with the BLAST, hmmpress and cmpress index files) is checked before the first job.
if DB_SHARED:
  if not LINEAGE:
    raise ValueError("A shared db_dir needs a named BUSCO lineage, for example lineage: \"viridiplantae_odb12\"")
  db_missing = [path for path in [
      f"{DB_DIR}/versions.txt",
      f"{DB_DIR}/uniprot_sprot.fasta",
      f"{DB_DIR}/uniprot_sprot.fasta.pin",
      f"{DB_DIR}/Pfam-A.hmm",
      f"{DB_DIR}/Pfam-A.hmm.h3i",
      f"{DB_DIR}/Rfam.cm",
      f"{DB_DIR}/Rfam.cm.i1i",
      f"{DB_DIR}/busco_downloads/lineages/{LINEAGE}"]
    if not os.path.exists(path)]
  if db_missing:
    raise ValueError("db_dir %s lacks:\n  %s\nFill the folder with the download_databases target, see README.md."
      % (DB_DIR, "\n  ".join(db_missing)))

@functools.lru_cache(maxsize=None)
def check_samples(path):
  """
  Checks the samples file before any job starts: 3 or 4 columns on each non-empty line, unique replicate
  names, existing reads files, two or more libraries (TMM needs them).
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
  if len(replicates) < 2:
    problems.append("it has %d libraries, and it needs 2 or more" % len(replicates))
  if problems:
    raise WorkflowError("The samples file %s has problems:\n  %s" % (path, "\n  ".join(problems)))
  return path


def checked_samples_file(wildcards):
  return check_samples(config["samples_file"])


def trimmomatic_completed_parallel_jobs(wildcards):
  """
  The marker files of the trimmomatic_parallel jobs, one for each line of the samples file.
  """
  parallel_dir = checkpoints.trimmomatic_split.get(**wildcards).output[0]
  # sorted: glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "sample_{job_index}")).job_index)
  return expand(os.path.join(parallel_dir, "completed_{job_index}"), job_index=job_ids)


def transcriptome_copy_input(wildcards):
  if config["assembler"] == "rnaspades":
    return {"transcriptome": getattr(rules.rnaspades.output, RNASPADES_FILTER)}
  return {"transcriptome": rules.trinity.output.transcriptome,
          "gene_trans_map": rules.trinity.output.gene_trans_map}


def abundance_strand_option():
  """
  The strand option of align_and_estimate_abundance.pl. Trinity: strand_specific as it is. rnaSPAdes:
  "--ss rf" and "--ss fr" become "--SS_lib_type RF" and "FR", "" gives none, anything else stops the run.
  """
  if config["assembler"] == "rnaspades":
    options = {"": "", "--ss rf": "--SS_lib_type RF", "--ss fr": "--SS_lib_type FR"}
    if config["strand_specific"] not in options:
      raise ValueError('With rnaSPAdes, strand_specific is "", "--ss rf" or "--ss fr", not "%s"' % config["strand_specific"])
    return options[config["strand_specific"]]
  return config["strand_specific"]


def get_igv(wildcards):
  # one indexed alignment per replicate, the second column of the samples file
  with open(check_samples(config["samples_file"])) as handle:
    replicates = [line.split()[1] for line in handle if line.split()]
  return expand("results/expression/bowtie2_rsem/igv/{sample}.sorted.bam.bai", sample=replicates)


def parallel_annotation_tasks(wildcards):
  """
  The chunk outputs of one annotation tool for one sequence set.
  """
  parallel_dir = checkpoints.fasta_split.get(split=wildcards["split"]).output[0]
  # sorted by number, so the merged file follows the order of the sequences at any chunk size.
  # glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "{index}.%s" % wildcards["split"])).index, key=int)
  return expand("results/annotation/chunks/{task}/{index}.out", index=job_ids, task=wildcards["task"])

# These rules run in the main job, not as cluster jobs
localrules: all, clean, delete_intermediates, run_info, trimmomatic_split, trimmomatic_merge, samples_yaml_conversion, transcriptome_copy, summary_peps, orf_sets_pep, compare_qc_after_trim, fastqc_before_trim_warnings, IGV

# The cluster writes its log into the log folder before the job starts, and not every cluster makes the folder.
onstart:
  for workflow_rule in workflow.rules:
    for log_file in workflow_rule.log:
      log_folder = os.path.dirname(str(log_file))
      if log_folder:
        os.makedirs(log_folder, exist_ok=True)

rule clean:
  """
  Not run by default. Deletes results/ and logs/. Keeps config.yaml, the samples file, resources/ and .snakemake/.
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
  Not run by default. Deletes the trimmed reads and the kallisto, bowtie2 and RSEM indexes, which no rule names
  as output, so a later run has nothing to do.
  """
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  shell:
    """
    rm -f results/trimmomatic/trimmed_reads/*.fastq.gz
    rm -f results/expression/kallisto/transcriptome.fasta.kallisto_idx
    # the bowtie2.ok marker goes with the index
    rm -f results/expression/bowtie2_rsem/transcriptome.fasta.bowtie2.* results/expression/bowtie2_rsem/transcriptome.fasta.RSEM.*
    """


DEPLOYMENT = " ".join(sorted(method.name.lower().replace("_", "-")
  for method in workflow.deployment_settings.deployment_method))

def main_image_file():
  """
  The image file of MAIN_IMAGE: <md5 of the URL>.simg in the apptainer prefix.
  """
  prefix = workflow.deployment_settings.apptainer_prefix or os.path.join(".snakemake", "singularity")
  return os.path.join(str(prefix), hashlib.md5(MAIN_IMAGE.encode()).hexdigest() + ".simg")

rule run_info:
  """
  Records the computer of the main job and of the assembly job, the transXpress commit, the images, the
  database releases and md5 sums and the version of each tool. It runs after the assembly and BUSCO, which it reads.
  """
  input:
    transcriptome="results/summary/transcriptome.fasta",
    busco="results/summary/busco_report.txt",
    pfam=f"{DB_DIR}/Pfam-A.hmm",
    rfam=f"{DB_DIR}/Rfam.cm",
    sprot=f"{DB_DIR}/uniprot_sprot.fasta"
  output:
    "results/summary/run_info.txt"
  params:
    deployment=DEPLOYMENT,
    # the rule envs come from the image only with apptainer
    image=MAIN_IMAGE if DeploymentMethod.APPTAINER in workflow.deployment_settings.deployment_method else "",
    image_file=main_image_file(),
    signalp_container=SIGNALP_CONTAINER,
    targetp_container=TARGETP_CONTAINER
  threads:
    1
  resources:
    mem_mb=2048,
    runtime="4h"
  script:
    "../scripts/run_info.py"
