# The sequence sets the annotation tools read: transcripts (fasta), long ORFs (orfs), proteins of both ORF sets (pep).
# The extension of a chunk file is the name of its set.
ANNOTATION_SPLIT_INPUT = {
  "fasta": "results/summary/transcriptome.fasta",
  "orfs": "results/transdecoder/transcriptome.orfs",
  "pep": "results/transdecoder/transcriptome_both_sets.pep"}

wildcard_constraints:
  split="fasta|orfs|pep"

checkpoint fasta_split:
  """
  Splits a sequence set into chunks of annotation_chunk_<set> sequences.
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
    # {output}.part gets the chunk name only when cmscan finished (AGENTS.md "Traps")
    rm -f {output}.part
    cmscan -E {params.evalue} --rfam --cpu {threads} --tblout {output}.part {input[db]} {input[fasta]} &> {log}
    source {workflow.basedir}/scripts/tools.sh
    publish_part {output} {log}
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
    1
  resources:
    mem_mb=2048,
    runtime="4h"
  shell:
    """
    # {output}.part gets the chunk name only when hmmscan finished (AGENTS.md "Traps")
    rm -f {output}.part
    # TransDecoder reads --domtblout. One thread: with more, the width of the name column changes from run to run,
    # and the chunk takes as long (cpu 1 / cpu 2 = 0.99 to 1.02).
    hmmscan -E {params.evalue} --cpu {threads} --domtblout {output}.part {input[db]} {input[fasta]} &> {log}
    source {workflow.basedir}/scripts/tools.sh
    publish_part {output} {log}
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
    # {output}.part gets the chunk name only after a clean exit (-outfmt 6 has no end line to check)
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
    # {output}.part gets the chunk name only after a clean exit (-outfmt 6 has no end line to check)
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
  # the image of signalp_container, or the conda env "signalp". No env when the tool is off (--containerize cannot take a named env).
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
    TMP=$(mktemp -d "${{TMPDIR:-/tmp}}/signalp.XXXXXX")
    trap 'rm -rf "$TMP"' EXIT
    echo "work folder: $TMP" &> {log}
    signalp6 --fastafile {input} --organism {params.organism} --output_dir $TMP/out --format none --mode fast &>> {log}
    # the move from the temporary folder can be a copy, so the chunk gets its name only after it
    rm -f {output}.part
    mv $TMP/out/prediction_results.txt {output}.part
    mv {output}.part {output}
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
  # the image of targetp_container, or the conda env "targetp". No env when the tool is off (--containerize cannot take a named env).
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
    # targetp writes its result into the folder it runs in
    RUN=$(pwd)
    TMP=$(mktemp -d "${{TMPDIR:-/tmp}}/targetp.XXXXXX")
    trap 'rm -rf "$TMP"' EXIT
    echo "work folder: $TMP" &> "$RUN/{log}"
    cd $TMP
    targetp -fasta "$RUN/{input}" -format short -org {params.organism} -prefix targetp_{wildcards.index} &>> "$RUN/{log}"
    # the move from the temporary folder can be a copy, so the chunk gets its name only after it
    rm -f "$RUN/{output}.part"
    mv targetp_{wildcards.index}_summary.targetp2 "$RUN/{output}.part" &>> "$RUN/{log}"
    mv "$RUN/{output}.part" "$RUN/{output}"
    """


rule write_summaries:
  """
  Puts the annotations in the headers of the transcripts and of the proteins of each ORF set.
  """
  input:
    transcriptome="results/summary/transcriptome.fasta",
    proteome_all_orfs="results/transdecoder/transcriptome_all_orfs.pep",
    proteome_single_best_orf="results/transdecoder/transcriptome_single_best_orf.pep",
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
    proteome_all_orfs_annotated="results/summary/transcriptome_all_orfs_annotated.pep",
    proteome_single_best_orf_annotated="results/summary/transcriptome_single_best_orf_annotated.pep"
  log:
    "logs/write_summaries.log"
  threads:
    1
  resources:
    mem_mb=4096,
    runtime="4h"
  script:
    "../scripts/write_summaries.py"
