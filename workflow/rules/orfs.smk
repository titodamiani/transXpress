rule transdecoder_longorfs:
  """
  Runs first stage of Transdecoder extracting the long open reading frames.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
  output:
    orfs="results/transdecoder/transcriptome.orfs"
  log:
    "logs/transdecoder_longorfs.log"
  conda:
    "../envs/transdecoder.yaml"
  params:
    memory="2"
  threads:
    1
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


rule transdecoder_predict_all:
  """
  Runs second stage of Transdecoder predicting the likely coding regions based
  on Pfam and SwissProt hits. Keeps every ORF that Predict accepts, so a
  transcript can give several proteins. Saves the proteins (.pep), the coding
  sequences (.cds), and the ORF positions on the transcripts (.gff3, .bed).

  Predict skips a step when its checkpoint file exists, and the checkpoint name
  does not include --single_best_only. So each ORF set runs on its own copy of
  the LongOrfs folder, without the checkpoints of an earlier Predict run. The
  copy is deleted when the four files are saved.

  Predict reads the ORFs of one transcript in Perl hash order, which is random
  in each run. A fixed hash seed gives the same output order and the same pick
  among tied ORFs in each run.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
    orfs="results/transdecoder/transcriptome.orfs",
    pfam="results/annotation/pfam_orfs.out",
    blastp="results/annotation/sprotblastp_orfs.out"
  output:
    pep="results/transdecoder/transcriptome_all_orfs.pep",
    cds="results/transdecoder/transcriptome_all_orfs.cds",
    gff3="results/transdecoder/transcriptome_all_orfs.gff3",
    bed="results/transdecoder/transcriptome_all_orfs.bed"
  log:
    "logs/transdecoder_predict_all.log"
  conda:
    "../envs/transdecoder.yaml"
  params:
    memory=config.get("transdecoder_predict_memory", "6"),
    workdir="results/transdecoder/predict_all"
  threads:
    config.get("transdecoder_predict_threads", 1)
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
    $TD_UTIL/TransDecoder.Predict -t {input.transcriptome} --output_dir {params.workdir} --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.pep {output.pep} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.cds {output.cds} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.gff3 {output.gff3} &>> {log}
    cp -p {params.workdir}/$NAME.transdecoder.bed {output.bed} &>> {log}
    rm -rf {params.workdir} &>> {log}
    """


rule transdecoder_predict_single_best:
  """
  Runs second stage of Transdecoder like transdecoder_predict_all, with
  --single_best_only. Keeps one ORF per transcript, picked by Pfam or SwissProt
  hit first and by length second.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
    orfs="results/transdecoder/transcriptome.orfs",
    pfam="results/annotation/pfam_orfs.out",
    blastp="results/annotation/sprotblastp_orfs.out"
  output:
    pep="results/transdecoder/transcriptome_single_best_orf.pep",
    cds="results/transdecoder/transcriptome_single_best_orf.cds",
    gff3="results/transdecoder/transcriptome_single_best_orf.gff3",
    bed="results/transdecoder/transcriptome_single_best_orf.bed"
  log:
    "logs/transdecoder_predict_single_best.log"
  conda:
    "../envs/transdecoder.yaml"
  params:
    memory=config.get("transdecoder_predict_memory", "6"),
    workdir="results/transdecoder/predict_single_best"
  threads:
    config.get("transdecoder_predict_threads", 1)
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
    $TD_UTIL/TransDecoder.Predict -t {input.transcriptome} --output_dir {params.workdir} --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} --single_best_only &>> {log}
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
    "results/transdecoder/transcriptome.pep"
  log:
    "logs/transdecoder_select.log"
  shell:
    """
    cp -p {input} {output} &> {log}
    """
