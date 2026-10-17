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
    chunk_size=ANNOTATION_CHUNK["fasta"],
    extension="fasta",
    strip_stop=False
  threads:
    1
  script:
    "../scripts/split_fasta.py"
        


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
    chunk_size=ANNOTATION_CHUNK["orfs"],
    extension="orfs",
    strip_stop=False
  threads:
    1
  script:
    "../scripts/split_fasta.py"
        


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
    chunk_size=ANNOTATION_CHUNK["pep"],
    extension="pep",
    strip_stop=True
  threads:
    1
  script:
    "../scripts/split_fasta.py"
        


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
    "../envs/rfam.yaml"
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
    "../envs/pfam.yaml"
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
    "../envs/blast.yaml"
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
    "../envs/blast.yaml"
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
    signalp_results=["results/annotation/signalp_pep.out"] if SIGNALP_ORGANISM else [],
    targetp_results=["results/annotation/targetp_pep.out"] if TARGETP_ORGANISM else []
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
  script:
    "../scripts/annotated_fasta.py"
