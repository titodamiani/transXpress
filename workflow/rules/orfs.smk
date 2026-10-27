rule transdecoder_longorfs:
  """
  Runs first stage of Transdecoder extracting the long open reading frames.
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
    # LongOrfs names its folder after the file name of the transcriptome
    TD_DIR=results/transdecoder/longorfs/$(basename {input.transcriptome}).transdecoder_dir
    # LongOrfs skips a step when its checkpoint file exists, so a rerun starts from an empty folder
    rm -rf $TD_DIR &> {log}
    # A fixed Perl hash seed gives the same output order in each run
    export PERL_HASH_SEED=0
    # TransDecoder 6 keeps the stage scripts in util/ next to the TransDecoder wrapper, outside bin/
    TD_UTIL=$(dirname $(readlink -f $(command -v TransDecoder)))/util
    $TD_UTIL/TransDecoder.LongOrfs -t {input.transcriptome} --output_dir results/transdecoder/longorfs &>> {log}
    cp -p $TD_DIR/longest_orfs.pep {output.orfs} &>> {log}
    """


rule transdecoder_predict:
  """
  Runs second stage of Transdecoder predicting the likely coding regions based
  on Pfam and SwissProt hits. Saves the proteins (.pep), the coding sequences
  (.cds), and the ORF positions on the transcripts (.gff3, .bed). It runs for
  two ORF sets:
  all_orfs keeps every ORF that Predict accepts, so a transcript can give
  several proteins.
  single_best_orf runs with --single_best_only. It keeps one ORF per transcript,
  picked by Pfam or SwissProt hit first and by length second.

  Predict skips a step when its checkpoint file exists, and the checkpoint name
  does not include --single_best_only. So each ORF set runs on its own copy of
  the LongOrfs folder, without the checkpoints of an earlier Predict run. The
  copy is deleted when the four files are saved.

  Predict reads the ORFs of one transcript in Perl hash order, which is random
  in each run. A fixed hash seed gives the same output order and the same pick
  among tied ORFs in each run.
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
    # Predict and LongOrfs name their files after the file name of the transcriptome
    NAME=$(basename {input.transcriptome})
    cp -rp results/transdecoder/longorfs/$NAME.transdecoder_dir {params.workdir}/ &>> {log}
    rm -rf {params.workdir}/$NAME.transdecoder_dir/__checkpoints_TDpredict &>> {log}
    export PERL_HASH_SEED=0
    # TransDecoder 6 keeps the stage scripts in util/ next to the TransDecoder wrapper, outside bin/
    TD_UTIL=$(dirname $(readlink -f $(command -v TransDecoder)))/util
    $TD_UTIL/TransDecoder.Predict -t {input.transcriptome} --output_dir {params.workdir} --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} {params.single_best} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.pep {output.pep} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.cds {output.cds} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.gff3 {output.gff3} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.bed {output.bed} &>> {log}
    rm -rf {params.workdir} &>> {log}
    """


rule transdecoder_select:
  """
  Copies the ORF set that transdecoder_orfs in config.yaml names to
  transcriptome.pep, which the downstream rules read.
  """
  input:
    TRANSDECODER_ORFS_FILE
  output:
    "results/summary/transcriptome.pep"
  log:
    "logs/transdecoder_select.log"
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  shell:
    """
    cp -p {input} {output} &> {log}
    """
