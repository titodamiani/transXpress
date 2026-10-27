import functools
import hashlib
import os
import tempfile
from snakemake.exceptions import WorkflowError
from snakemake.settings.types import DeploymentMethod
from snakemake.utils import min_version, validate

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

# The ORF set that transcriptome.pep holds and the downstream rules read. Both sets are always saved.
TRANSDECODER_ORFS = config["transdecoder_orfs"]
TRANSDECODER_ORFS_FILE = {"all": "results/transdecoder/transcriptome_all_orfs.pep", "single_best": "results/transdecoder/transcriptome_single_best_orf.pep"}[TRANSDECODER_ORFS]

# The isoform expression matrix that ExN50 and the annotated FASTA headers read: the TPM of each transcript,
# scaled across the libraries with TMM. The ExN50 page of the Trinity wiki uses the TMM matrix.
TMM_ISOFORM = "results/summary/expr_tmm_isoform.tsv"

# The rnaSPAdes transcript set that becomes transcriptome.fasta. rnaSPAdes always writes all three.
RNASPADES_FILTER = config["rnaspades_filter"]

# The organism groups of SignalP and TargetP. "" turns the tool off.
# A key with no value reaches here as None, which means "" too
SIGNALP_ORGANISM = config["signalp"] or ""
TARGETP_ORGANISM = config["targetp"] or ""

# The Apptainer images of SignalP and TargetP. "" runs the tool in its conda env. A rule with an image
# has no conda env, because Snakemake would activate the conda env inside the image.
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

# The number of sequences in each annotation chunk. Each chunk is one cluster job for each tool that reads it.
ANNOTATION_CHUNK = {split: config[f"annotation_chunk_{split}"] for split in ["fasta", "orfs", "pep"]}

# The environment variable that holds a node-local folder in a cluster job, for example "$SCRATCH" or
# "SCRATCH". "" runs every rule in the run folder. The heavy rules read the variable with bash, inside
# the job, because only the job knows its value.
SCRATCH = config["scratch"] or ""
SCRATCH_VAR = SCRATCH.lstrip("$")

def job_tmpdir(wildcards):
  """
  The temporary folder of a job: the folder in the variable that the scratch key names, or the
  system default when the key is "" or the variable is not set. Snakemake calls this inside the job.
  """
  return os.environ.get(SCRATCH_VAR) or tempfile.gettempdir()

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
      f"{DB_DIR}/versions.txt",
      f"{DB_DIR}/uniprot_sprot.fasta",
      f"{DB_DIR}/Pfam-A.hmm",
      f"{DB_DIR}/Rfam.cm",
      f"{DB_DIR}/busco_downloads/lineages/{LINEAGE}"]
    if not os.path.exists(path)]
  if db_missing:
    raise ValueError("db_dir %s lacks:\n  %s\nFill the folder with the download_databases target, see README.md."
      % (DB_DIR, "\n  ".join(db_missing)))

@functools.lru_cache(maxsize=None)
def check_samples(path):
  """
  Checks the samples file before any job starts. Every line that is not empty has 3 columns (sample,
  replicate and one reads file) or 4 columns (a second reads file), separated by white space. A replicate
  name appears once, and every reads file exists. The file has two or more libraries, because the TMM
  normalization of the expression needs them.
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
  # The TMM normalization of the expression needs two or more libraries
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
  # sorted, because glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "sample_{job_index}")).job_index)
  return expand(os.path.join(parallel_dir, "completed_{job_index}"), job_index=job_ids)


def transcriptome_copy_input(wildcards):
  if config["assembler"] == "rnaspades":
    return {"transcriptome": getattr(rules.rnaspades.output, RNASPADES_FILTER)}
  return {"transcriptome": rules.trinity.output.transcriptome,
          "gene_trans_map": rules.trinity.output.gene_trans_map}


def abundance_strand_option():
  """
  The strand option of align_and_estimate_abundance.pl. With Trinity, strand_specific is already a
  Trinity option, for example "--SS_lib_type RF". With rnaSPAdes, "--ss rf" and "--ss fr" become
  "--SS_lib_type RF" and "--SS_lib_type FR", and "" (unstranded reads) gives no option. Any other value
  stops the run before the first job: rnaSPAdes takes only rf and fr, and a value that is not one of
  them would make the expression estimates run as if the reads had no strand.
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
    rm -f results/trimmomatic/trimmed_reads/*.fastq.gz
    rm -f results/expression/kallisto/transcriptome.fasta.kallisto_idx
    # the bowtie2.ok marker goes with the index, so a later align_reads run builds the index again
    rm -f results/expression/bowtie2_rsem/transcriptome.fasta.bowtie2.* results/expression/bowtie2_rsem/transcriptome.fasta.RSEM.*
    """


# The software deployment methods of the run, for example "apptainer conda"
DEPLOYMENT = " ".join(sorted(method.name.lower().replace("_", "-")
  for method in workflow.deployment_settings.deployment_method))

def main_image_file():
  """
  The image file of MAIN_IMAGE that apptainer runs: Snakemake keeps it as <md5 of the URL>.simg in the
  apptainer prefix, .snakemake/singularity by default.
  """
  prefix = workflow.deployment_settings.apptainer_prefix or os.path.join(".snakemake", "singularity")
  return os.path.join(str(prefix), hashlib.md5(MAIN_IMAGE.encode()).hexdigest() + ".simg")

rule run_info:
  """
  Records the computer the pipeline ran on, the transXpress commit, the images and the version of each
  tool. The assembly depends on the hardware, the assembler version and the thread count, so a later run
  can be compared against this one.
  """
  output:
    "results/summary/run_info.txt"
  params:
    deployment=DEPLOYMENT,
    # the rule envs come from the main image only when apptainer is a deployment method
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
