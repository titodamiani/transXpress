rule transdecoder_longorfs:
  """
  TransDecoder.LongOrfs: the long ORFs.
  """
  input:
    transcriptome="results/summary/transcriptome.fasta",
  output:
    orfs="results/transdecoder/transcriptome.orfs"
  log:
    "logs/transdecoder_longorfs.log"
  conda:
    "../envs/transdecoder.yaml"
  threads:
    1
  resources:
    mem_mb=2048,
    runtime="4h"
  shell:
    """
    TD_DIR=results/transdecoder/longorfs/$(basename {input.transcriptome}).transdecoder_dir
    # LongOrfs skips a step whose checkpoint file exists, so a rerun starts from an empty folder
    rm -rf $TD_DIR &> {log}
    # a fixed Perl hash seed makes the output order repeatable
    export PERL_HASH_SEED=0
    # the stage scripts are in util/ beside the TransDecoder wrapper
    TD_UTIL=$(dirname $(readlink -f $(command -v TransDecoder)))/util
    $TD_UTIL/TransDecoder.LongOrfs -t {input.transcriptome} --output_dir results/transdecoder/longorfs &>> {log}
    cp -p $TD_DIR/longest_orfs.pep {output.orfs} &>> {log}
    """


rule transdecoder_predict:
  """
  TransDecoder.Predict with the Pfam and Swiss-Prot hits: .pep, .cds, .gff3, .bed for two ORF sets.
  all_orfs keeps every accepted ORF. single_best_orf (--single_best_only) keeps one ORF per transcript.

  The checkpoint name does not include --single_best_only, so each set runs on its own copy of the
  LongOrfs folder without Predict checkpoints. A fixed Perl hash seed makes the output order and the pick
  among tied ORFs repeatable.
  """
  input:
    transcriptome="results/summary/transcriptome.fasta",
    orfs="results/transdecoder/transcriptome.orfs",
    pfam="results/annotation/pfam_orfs.out",
    blastp="results/annotation/sprotblastp_orfs.out"
  output:
    pep="results/transdecoder/transcriptome_{orf_set}.pep",
    cds="results/transdecoder/transcriptome_{orf_set}.cds",
    gff3="results/transdecoder/transcriptome_{orf_set}.gff3",
    bed="results/transdecoder/transcriptome_{orf_set}.bed"
  wildcard_constraints:
    orf_set="all_orfs|single_best_orf"
  log:
    "logs/transdecoder_predict_{orf_set}.log"
  conda:
    "../envs/transdecoder.yaml"
  params:
    workdir="results/transdecoder/predict_{orf_set}",
    single_best=lambda wildcards: "--single_best_only" if wildcards.orf_set == "single_best_orf" else ""
  threads:
    config["transdecoder_predict_threads"]
  resources:
    mem_mb=gb_to_mb(config["transdecoder_predict_memory"]),
    runtime=config["transdecoder_predict_runtime"]
  shell:
    """
    rm -rf {params.workdir} &> {log}
    mkdir -p {params.workdir} &>> {log}
    # both stages name their files after the transcriptome file
    NAME=$(basename {input.transcriptome})
    cp -rp results/transdecoder/longorfs/$NAME.transdecoder_dir {params.workdir}/ &>> {log}
    rm -rf {params.workdir}/$NAME.transdecoder_dir/__checkpoints_TDpredict &>> {log}
    export PERL_HASH_SEED=0
    # the stage scripts are in util/ beside the TransDecoder wrapper
    TD_UTIL=$(dirname $(readlink -f $(command -v TransDecoder)))/util
    $TD_UTIL/TransDecoder.Predict -t {input.transcriptome} --output_dir {params.workdir} --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} {params.single_best} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.pep {output.pep} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.cds {output.cds} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.gff3 {output.gff3} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.bed {output.bed} &>> {log}
    rm -rf {params.workdir} &>> {log}
    """


rule summary_peps:
  """
  Copies the proteins of both ORF sets to results/summary.
  """
  input:
    expand("results/transdecoder/transcriptome_{orf_set}.pep", orf_set=["all_orfs", "single_best_orf"])
  output:
    expand("results/summary/transcriptome_{orf_set}.pep", orf_set=["all_orfs", "single_best_orf"])
  log:
    "logs/summary_peps.log"
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  shell:
    """
    cp -p {input[0]} {output[0]} &> {log}
    cp -p {input[1]} {output[1]} &>> {log}
    """


rule orf_sets_pep:
  """
  Joins the proteins of both ORF sets, so TMHMM, SignalP and TargetP read each protein once.
  A protein ID with two different sequences in the sets stops the rule.
  """
  input:
    all_orfs="results/transdecoder/transcriptome_all_orfs.pep",
    single_best_orf="results/transdecoder/transcriptome_single_best_orf.pep"
  output:
    "results/transdecoder/transcriptome_both_sets.pep"
  log:
    "logs/orf_sets_pep.log"
  threads:
    1
  resources:
    mem_mb=2048,
    runtime="4h"
  script:
    "../scripts/orf_sets_pep.py"
