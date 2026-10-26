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
    # abundance_start (scripts/scratch.sh) works in the tool folder {params.out_dir}, or on the node's
    # local disk when the scratch key names a variable that is set in this job
    source {workflow.basedir}/scripts/scratch.sh
    abundance_start align_reads "{params.scratch_var}" {log} {params.out_dir} {input.transcriptome} {input.gene_trans_map} {input.samples}

    # Remove old markers so the aligner runs again on the trimmed reads.
    for rep in $(awk '{{print $2}}' samples_trimmed.txt); do rm -f "$rep/bowtie2.bam.ok"; done

    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')
    $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta {params.strand} --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"

    # for_rsem.bam is a large temporary copy of the alignment. Nothing reads it later.
    for rep in $(awk '{{print $2}}' samples_trimmed.txt); do rm -f "$rep/bowtie2.bam.for_rsem.bam"; done

    abundance_finish
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
    tmm_isoform="results/expression/transcriptome_tmm_isoform.tsv",
    tpm_isoform="results/expression/transcriptome_tpm_isoform.tsv",
    counts_isoform="results/expression/transcriptome_counts_isoform.tsv",
    tmm_gene="results/expression/transcriptome_tmm_gene.tsv",
    tpm_gene="results/expression/transcriptome_tpm_gene.tsv",
    counts_gene="results/expression/transcriptome_counts_gene.tsv",
    counts="results/expression/kallisto/kallisto.gene.counts.matrix"
  log:
    "logs/kallisto.log"
  conda:
    "../envs/trinity_utils.yaml"
  params:
    scratch_var=SCRATCH_VAR,
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

    # abundance_start (scripts/scratch.sh) works in the tool folder {params.out_dir}, or on the node's
    # local disk when the scratch key names a variable that is set in this job
    source {workflow.basedir}/scripts/scratch.sh
    abundance_start kallisto "{params.scratch_var}" {log} {params.out_dir} {input.transcriptome} {input.gene_trans_map} {input.samples}

    $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts transcriptome.fasta {params.strand} --seqType fq --samples_file samples_trimmed.txt --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map transcriptome.gene_trans_map &>> "$LOG"

    # the matrix columns follow the order of the files on the command line, here the order of the samples file
    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map transcriptome.gene_trans_map $(awk '{{print $2 "/abundance.tsv"}}' samples_trimmed.txt) &>> "$LOG"
    # TMM needs two or more libraries. check_samples stops a run with fewer before the first job, so a
    # missing matrix here is an error.
    for f in kallisto.{{isoform,gene}}.TMM.EXPR.matrix kallisto.{{isoform,gene}}.TPM.not_cross_norm kallisto.{{isoform,gene}}.counts.matrix; do
      if [ ! -f $f ]; then
        echo "$f was not produced" &>> "$LOG"
        exit 1
      fi
    done
    cp -p kallisto.isoform.TMM.EXPR.matrix "$RUN/{output.tmm_isoform}" &>> "$LOG"
    cp -p kallisto.isoform.TPM.not_cross_norm "$RUN/{output.tpm_isoform}" &>> "$LOG"
    cp -p kallisto.isoform.counts.matrix "$RUN/{output.counts_isoform}" &>> "$LOG"
    cp -p kallisto.gene.TMM.EXPR.matrix "$RUN/{output.tmm_gene}" &>> "$LOG"
    cp -p kallisto.gene.TPM.not_cross_norm "$RUN/{output.tpm_gene}" &>> "$LOG"
    cp -p kallisto.gene.counts.matrix "$RUN/{output.counts_gene}" &>> "$LOG"

    abundance_finish
    """
