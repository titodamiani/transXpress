# The three sequence sets that the annotation tools read: the transcripts (fasta), the long ORFs (orfs)
# and the proteins (pep). The extension of a chunk file is the name of its set.
ANNOTATION_SPLIT_INPUT = {
  "fasta": "results/assembly/transcriptome.fasta",
  "orfs": "results/transdecoder/transcriptome.orfs",
  "pep": "results/transdecoder/transcriptome.pep"}

wildcard_constraints:
  split="fasta|orfs|pep"

checkpoint fasta_split:
  """
  Splits the transcriptome.fasta, the transcriptome.orfs or the transcriptome.pep into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    lambda wildcards: ANNOTATION_SPLIT_INPUT[wildcards.split]
  output:
    directory("results/annotation/chunks/{split}")
  log:
    "logs/fasta_split_{split}.log"
  params:
    chunk_size=lambda wildcards: ANNOTATION_CHUNK[wildcards.split],
    extension=lambda wildcards: wildcards.split,
    strip_stop=lambda wildcards: wildcards.split == "pep"
  threads:
    1
  resources:
    mem_mb=2048,
    runtime="4h"
  script:
    "../scripts/split_fasta.py"


rule annotation_merge:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks
  output:
    "results/annotation/{task}_{split}.out"
  log:
    "logs/annotation_merge_{split}/{task}.log"
  threads:
    1
  resources:
    mem_mb=2048,
    runtime="4h"
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
    db=f"{DB_DIR}/Rfam.cm"
  output:
    "results/annotation/chunks/rfam/{index}.out"
  log:
    "logs/rfam_parallel/{index}.log"
  conda:
    "../envs/rfam.yaml"
  params:
    evalue=config["e_value_threshold"]
  threads:
    4
  resources:
    mem_mb=5120,
    runtime="4h"
  shell:
    """
    # The tool writes {output}.part. The chunk gets its name only when the file ends with the "# [ok]" line
    # of a finished cmscan, so a job that is killed or hangs leaves no cut file under the name Snakemake trusts.
    rm -f {output}.part
    cmscan -E {params.evalue} --rfam --cpu {threads} --tblout {output}.part {input[db]} {input[fasta]} &> {log}
    tail -n 1 {output}.part | grep -qx '# \\[ok\\]' || {{ echo "{output}.part does not end with # [ok]" >> {log}; exit 1; }}
    mv {output}.part {output}
    """


rule pfam_parallel:
  """
  Runs hmmscan on Pfam database on smaller protein files in parallel.
  """
  input:
    fasta="results/annotation/chunks/orfs/{index}.orfs",
    db=f"{DB_DIR}/Pfam-A.hmm"
  output:
    "results/annotation/chunks/pfam/{index}.out"
  log:
    "logs/pfam_parallel/{index}.log"
  conda:
    "../envs/pfam.yaml"
  params:
    evalue=config["e_value_threshold"]
  threads:
    2
  resources:
    mem_mb=2048,
    runtime="4h"
  shell:
    """
    # The tool writes {output}.part. The chunk gets its name only when the file ends with the "# [ok]" line
    # of a finished hmmscan, so a job that is killed or hangs leaves no cut file under the name Snakemake trusts.
    rm -f {output}.part
    # Transdecoder requires --domtblout output
    hmmscan -E {params.evalue} --cpu {threads} --domtblout {output}.part {input[db]} {input[fasta]} &> {log}
    tail -n 1 {output}.part | grep -qx '# \\[ok\\]' || {{ echo "{output}.part does not end with # [ok]" >> {log}; exit 1; }}
    mv {output}.part {output}
    """


rule sprot_blastp_parallel:
  """
  Runs blast search on SwissProt database on smaller protein files in parallel.
  """
  input:
    fasta="results/annotation/chunks/orfs/{index}.orfs",
    db=f"{DB_DIR}/uniprot_sprot.fasta"
  output:
    "results/annotation/chunks/sprotblastp/{index}.out"
  log:
    "logs/sprot_blastp_parallel/{index}.log"
  conda:
    "../envs/blast.yaml"
  params:
    evalue=config["e_value_threshold"]
  threads:
    2
  resources:
    mem_mb=2048,
    runtime="4h"
  shell:
    """
    # The tool writes {output}.part. The chunk gets its name only when blastp has ended without an error
    # (-outfmt 6 has no end line to check), so a job that is killed or hangs leaves no cut file under the name
    # Snakemake trusts.
    rm -f {output}.part
    blastp -query {input[fasta]} -db {input[db]} -num_threads {threads} -evalue {params.evalue} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output}.part &> {log}
    mv {output}.part {output}
    """


rule sprot_blastx_parallel:
  """
  Runs blast search on SwissProt database on smaller nucleotide files in parallel.
  """
  input:
    fasta="results/annotation/chunks/fasta/{index}.fasta",
    db=f"{DB_DIR}/uniprot_sprot.fasta"
  output:
    "results/annotation/chunks/sprotblastx/{index}.out"
  log:
    "logs/sprot_blastx_parallel/{index}.log"
  conda:
    "../envs/blast.yaml"
  params:
    evalue=config["e_value_threshold"]
  threads:
    config["sprot_blastx_threads"]
  resources:
    mem_mb=gb_to_mb(config["sprot_blastx_memory"]),
    runtime=config["sprot_blastx_runtime"]
  shell:
    """
    # The tool writes {output}.part. The chunk gets its name only when blastx has ended without an error
    # (-outfmt 6 has no end line to check), so a job that is killed or hangs leaves no cut file under the name
    # Snakemake trusts.
    rm -f {output}.part
    blastx -query {input[fasta]} -db {input[db]} -num_threads {threads} -evalue {params.evalue} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output}.part &> {log}
    mv {output}.part {output}
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
  conda:
    "../envs/tmhmm.yaml"
  threads:
    1
  resources:
    mem_mb=2048,
    runtime="4h"
  script:
    "../scripts/tmhmm_chunk.py"


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
  # the image of the signalp_container key, or the conda env "signalp" when the key is "". No env when the tool
  # is off, so snakemake --containerize, which cannot put a named env in an image, takes the rule.
  conda:
    None if SIGNALP_CONTAINER or not SIGNALP_ORGANISM else "signalp"
  container:
    SIGNALP_CONTAINER or None
  params:
    organism=SIGNALP_ORGANISM
  threads:
    1
  resources:
    mem_mb=4096,
    runtime="4h"
  shell:
    """
    # Snakemake sets TMPDIR to the job's temporary folder, see job_tmpdir
    TMP=$(mktemp -d "${{TMPDIR:-/tmp}}/signalp.XXXXXX")
    echo "work folder: $TMP" &> {log}
    signalp6 --fastafile {input} --organism {params.organism} --output_dir $TMP/out --format none --mode fast &>> {log}
    # The temporary folder can be on another disk than the run folder, so the move there is a copy. The copy goes
    # to {output}.part, and the rename gives the chunk its name only when the copy is whole.
    rm -f {output}.part
    mv $TMP/out/prediction_results.txt {output}.part
    mv {output}.part {output}
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
  # the image of the targetp_container key, or the conda env "targetp" when the key is "". No env when the tool
  # is off, so snakemake --containerize, which cannot put a named env in an image, takes the rule.
  conda:
    None if TARGETP_CONTAINER or not TARGETP_ORGANISM else "targetp"
  container:
    TARGETP_CONTAINER or None
  params:
    organism=TARGETP_ORGANISM
  threads:
    1
  resources:
    mem_mb=3072,
    runtime="4h"
  shell:
    """
    # targetp writes its result into the folder it runs in, which is the job's temporary folder
    # (Snakemake sets TMPDIR, see job_tmpdir)
    RUN=$(pwd)
    TMP=$(mktemp -d "${{TMPDIR:-/tmp}}/targetp.XXXXXX")
    echo "work folder: $TMP" &> "$RUN/{log}"
    cd $TMP
    targetp -fasta "$RUN/{input}" -format short -org {params.organism} -prefix targetp_{wildcards.index} &>> "$RUN/{log}"
    # The temporary folder can be on another disk than the run folder, so the move there is a copy. The copy goes
    # to {output}.part, and the rename gives the chunk its name only when the copy is whole.
    rm -f "$RUN/{output}.part"
    mv targetp_{wildcards.index}_summary.targetp2 "$RUN/{output}.part" &>> "$RUN/{log}"
    mv "$RUN/{output}.part" "$RUN/{output}"
    cd "$RUN"
    rm -r $TMP
    """


rule write_summaries:
  """
  Puts annotations in the headers of transcripts/proteins in the
  .fasta/.pep transcriptome files.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
    proteome="results/transdecoder/transcriptome.pep",
    expression=TMM_ISOFORM,
    blastx_results="results/annotation/sprotblastx_fasta.out",
    rfam_results="results/annotation/rfam_fasta.out",
    blastp_results="results/annotation/sprotblastp_orfs.out",
    pfam_results="results/annotation/pfam_orfs.out",
    tmhmm_results="results/annotation/tmhmm_pep.out",
    signalp_results=["results/annotation/signalp_pep.out"] if SIGNALP_ORGANISM else [],
    targetp_results=["results/annotation/targetp_pep.out"] if TARGETP_ORGANISM else []
  output:
    transcriptome_annotated="results/summary/transcriptome_annotated.fasta",
    proteome_annotated="results/summary/transcriptome_annotated.pep"
  log:
    "logs/write_summaries.log"
  threads:
    1
  resources:
    mem_mb=4096,
    runtime="4h"
  script:
    "../scripts/write_summaries.py"
