rule align_reads:
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
    "../envs/trinity_utils.yaml"
  params:
    scratch_var=SCRATCH_VAR,
    scratch_sh=SCRATCH_SH,
    out_dir=lambda wildcards, output: os.path.dirname(output[0]),
    strand=abundance_strand_option()
  threads:
    config["align_reads_threads"]
  resources:
    mem_mb=gb_to_mb(config["align_reads_memory"]),
    runtime=config["align_reads_runtime"]
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
    source {params.scratch_sh}
    copy_back() {{
      ( cd "$WORK" && find . -mindepth 1 -maxdepth 1 ! -type l ! -name samples_trimmed.txt -exec cp -rp -t "$OUT" {{}} + )
    }}
    # On a failure, the files go to the tool folder, where they would stand without the local disk.
    on_failure() {{
      cd "$RUN"
      echo "align_reads stopped with code $rc, so the files in $WORK go to {params.out_dir}" &>> "$LOG"
      copy_back &>> "$LOG"
    }}
    scratch_start "{params.scratch_var}" "$LOG" stop on_failure
    if [ -z "$LOCAL_ROOT" ]; then WORK="$OUT"; fi
    cd "$WORK"
    ln -sfn "$RUN/{input.transcriptome}" transcriptome.fasta
    ln -sfn "$RUN/{input.gene_trans_map}" transcriptome.gene_trans_map
    awk -v run="$RUN" 'BEGIN {{ OFS = "\t" }} {{ for (i = 3; i <= NF; i++) $i = run "/" $i; print }}' "$RUN/{input.samples}" > samples_trimmed.txt

    # Remove old markers so the aligner runs again on the trimmed reads.
    for rep in $(awk '{{print $2}}' samples_trimmed.txt); do rm -f "$rep/bowtie2.bam.ok"; done

    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')
    $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta {params.strand} --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"

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


rule prepare_samples_for_IGV:
  """
  Not run by default.
  Prepares the bam alignment files to be used by IGV by sorting
  them into results/expression/bowtie2_rsem/igv and indexing them.
  """
  input:
    # align_reads writes one bowtie2.bam per replicate, and names only its marker as the output
    "results/expression/bowtie2_rsem/transcriptome.fasta.bowtie2.ok"
  params:
    alignment="results/expression/bowtie2_rsem/{sample}/bowtie2.bam"
  output:
    sorted_alignment="results/expression/bowtie2_rsem/igv/{sample}.sorted.bam",
    indexed_alignment="results/expression/bowtie2_rsem/igv/{sample}.sorted.bam.bai"
  log:
    "logs/prepare_samples_for_IGV/{sample}.log"
  conda:
    "../envs/trinity_utils.yaml"
  threads:
    4
  resources:
    mem_mb=10240,
    runtime="4h"
  shell:
    """
    samtools sort --threads {threads} -m 2G -o {output.sorted_alignment} {params.alignment} &> {log}
    samtools index {output.sorted_alignment} &>> {log}
    """


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
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
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
    "../envs/trinity_utils.yaml"
  params:
    dispersion=config["dispersion"]
  threads:
    1
  resources:
    mem_mb=2048,
    runtime="4h"
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')

    # run_DE_analysis.pl estimates the dispersion for each pair of samples in which both samples have
    # more than one replicate, and uses --dispersion for a pair in which a sample has one replicate
    $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input.expression} --method edgeR --output {output} --samples_file {input.samples} --dispersion {params.dispersion} &> {log}
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
    isoform="results/expression/transcriptome_expression_isoform.tsv",
    gene="results/expression/transcriptome_expression_gene.tsv",
    counts="results/expression/kallisto/kallisto.gene.counts.matrix"
  log:
    "logs/kallisto.log"
  conda:
    "../envs/trinity_utils.yaml"
  params:
    scratch_var=SCRATCH_VAR,
    scratch_sh=SCRATCH_SH,
    out_dir=lambda wildcards, output: os.path.dirname(output.counts),
    strand=abundance_strand_option()
  threads:
    config["kallisto_threads"]
  resources:
    mem_mb=gb_to_mb(config["kallisto_memory"]),
    runtime=config["kallisto_runtime"]
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
    source {params.scratch_sh}
    copy_back() {{
      ( cd "$WORK" && find . -mindepth 1 -maxdepth 1 ! -type l ! -name samples_trimmed.txt -exec cp -rp -t "$OUT" {{}} + )
    }}
    # On a failure, the files go to the tool folder, where they would stand without the local disk.
    on_failure() {{
      cd "$RUN"
      echo "kallisto stopped with code $rc, so the files in $WORK go to {params.out_dir}" &>> "$LOG"
      copy_back &>> "$LOG"
    }}
    scratch_start "{params.scratch_var}" "$LOG" stop on_failure
    if [ -z "$LOCAL_ROOT" ]; then WORK="$OUT"; fi
    cd "$WORK"
    ln -sfn "$RUN/{input.transcriptome}" transcriptome.fasta
    ln -sfn "$RUN/{input.gene_trans_map}" transcriptome.gene_trans_map
    awk -v run="$RUN" 'BEGIN {{ OFS = "\t" }} {{ for (i = 3; i <= NF; i++) $i = run "/" $i; print }}' "$RUN/{input.samples}" > samples_trimmed.txt

    $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta {params.strand} --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"

    # the matrix columns follow the order of the files on the command line, here the order of the samples file
    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map transcriptome.gene_trans_map $(awk '{{print $2 "/abundance.tsv"}}' samples_trimmed.txt) &>> "$LOG"
    if [ -f kallisto.isoform.TMM.EXPR.matrix ]; then
      cp -p kallisto.isoform.TMM.EXPR.matrix "$RUN/{output.isoform}" &>> "$LOG"
    elif [ -f kallisto.isoform.TPM.not_cross_norm ]; then
      cp -p kallisto.isoform.TPM.not_cross_norm "$RUN/{output.isoform}" &>> "$LOG"
    else
      echo "Neither kallisto.isoform.TMM.EXPR.matrix nor kallisto.isoform.TPM.not_cross_norm was produced" &>> "$LOG"
      exit 1
    fi
    if [ -f kallisto.gene.TMM.EXPR.matrix ]; then
      cp -p kallisto.gene.TMM.EXPR.matrix "$RUN/{output.gene}" &>> "$LOG"
    elif [ -f kallisto.gene.TPM.not_cross_norm ]; then
      cp -p kallisto.gene.TPM.not_cross_norm "$RUN/{output.gene}" &>> "$LOG"
    else
      echo "Neither kallisto.gene.TMM.EXPR.matrix nor kallisto.gene.TPM.not_cross_norm was produced" &>> "$LOG"
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
