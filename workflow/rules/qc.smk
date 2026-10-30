rule fastqc:
  """
  Runs fastQC on individual reads files: the input reads (before_trim) or the trimmed reads (after_trim).
  """
  input:
    samples=lambda wildcards: checked_samples_file(wildcards) if wildcards.stage == "before" else "results/trimmomatic/samples_trimmed.txt"
  output:
    directory("results/qc/fastqc/{stage}_trim")
  wildcard_constraints:
    stage="before|after"
  log:
    "logs/fastqc_{stage}_trim.log"
  conda:
    "../envs/qc.yaml"
  threads:
    8
  resources:
    mem_mb=6144,
    runtime="4h"
  shell:
    """
    mkdir {output} &> {log}
    # the reads files, columns 3 and 4 of the samples file
    FILES=$(awk '{{ printf "%s\\n%s\\n", $3,$4}}' {input})
    # FastQC draws the text of its graphs with the fonts that fontconfig finds. A fonts.conf with only the
    # fonts of the env makes the graphs the same on every host and in the main image, which has no system
    # fonts. Java keeps the fonts it got in <user.home>/.java/fonts, so user.home is the job's folder.
    # (Snakemake sets TMPDIR, see job_tmpdir)
    TMP=$(mktemp -d "${{TMPDIR:-/tmp}}/fastqc.XXXXXX")
    trap 'rm -rf "$TMP"' EXIT
    ENV=$(dirname $(dirname $(command -v fastqc)))
    printf '<?xml version="1.0"?>\\n<!DOCTYPE fontconfig SYSTEM "urn:fontconfig:fonts.dtd">\\n<fontconfig>\\n<dir>%s/fonts</dir>\\n<dir>%s/share/fonts</dir>\\n<include ignore_missing="yes">%s/etc/fonts/conf.d</include>\\n<cachedir>%s/cache</cachedir>\\n</fontconfig>\\n' \
      "$ENV" "$ENV" "$ENV" "$TMP" > "$TMP/fonts.conf"
    export FONTCONFIG_FILE="$TMP/fonts.conf" JAVA_TOOL_OPTIONS="-Duser.home=$TMP"
    # {threads} files at a time
    fastqc -f fastq -t {threads} -o {output} $FILES &>> {log}
    """


rule multiqc:
  """
  Creates multiQC report from individual fastQC reports, before or after the trimming.
  """
  input:
    "results/qc/fastqc/{stage}_trim"
  output:
    out_dir=directory("results/qc/multiqc/{stage}_trim"),
    report="results/qc/multiqc/{stage}_trim/multiqc_fastqc.txt"
  wildcard_constraints:
    stage="before|after"
  log:
    "logs/multiqc_{stage}_trim.log"
  conda:
    "../envs/qc.yaml"
  threads:
    1
  resources:
    mem_mb=4096,
    runtime="4h"
  shell:
    """
    # MultiQC stops when the locale is not UTF-8
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    # Snakemake makes this folder when it makes the folder of the report in it
    mkdir -p {output.out_dir} &> {log}
    multiqc -o {output.out_dir} {input} &>> {log}
    cp {output.out_dir}/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """


rule fastqc_before_trim_warnings:
  """
  Prints warnings for FastQC modules which produced some warnings or fails of the input data.
  """
  input:
    fastqc_report="results/qc/multiqc/before_trim/multiqc_fastqc.txt"
  output:
    warning_file="results/qc/fastqc/WARNING_fastqc_before_trim_overview.txt"
  log:
    "logs/fastqc_before_trim_warnings.log"
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  script:
    "../scripts/fastqc_warnings.py"


rule compare_qc_after_trim:
  """
  Compares the quality of reads before and after trimming.
  Compares PER BASE SEQUENCE QUALITY, PER SEQUENCE QUALITY SCORES,
  OVERREPRESENTED SEQUENCES and ADAPTER CONTENT modules of FastQC.
  """
  input:
    before="results/qc/multiqc/before_trim/multiqc_fastqc.txt",
    after="results/qc/multiqc/after_trim/multiqc_fastqc.txt"
  output:
    result="results/qc/fastqc/FastQC_comparison_after_trim.txt"
  log:
    "logs/compare_qc_after_trim.log"
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  script:
    "../scripts/fastqc_compare.py"


rule trinity_stats:
  """
  Runs Trinity script to get statistics about the assembled transcriptome
  (number of transcripts, genes, GC content, EXN50).
  """
  input:
    transcriptome="results/summary/transcriptome.fasta",
    expression=TMM_ISOFORM
  output:
    stats="results/summary/transcriptome_stats.txt",
    exN50="results/qc/assembly/transcriptome_exN50_tmm.tsv",
    exN50plot="results/qc/assembly/transcriptome_exN50_tmm.pdf"
  log:
    "logs/trinity_stats.log"
  conda:
    "../envs/trinity_utils.yaml"
  params:
    assembler=config["assembler"]
  threads:
    1
  resources:
    mem_mb=2048,
    runtime="4h"
  shell:
    """
    source {workflow.basedir}/scripts/tools.sh
    TRINITY_HOME=$(trinity_home)
    # contig_ExN50_statistic.pl reads the genes in Perl hash order, which is random in each run, and
    # genes tied on expression keep that order. A fixed Perl hash seed gives the same order in each run.
    export PERL_HASH_SEED=0

    # TrinityStats.pl names the assembler "trinity" in its output
    $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} 2> {log} | sed -e 's/trinity/{params.assembler}/g' > {output.stats}
    # contig_ExN50_statistic.pl writes its .E-inputs table and the R script writes ExN50_plot.pdf into the
    # folder they run in, so both run in the folder of the outputs. The R script gives the plot a fixed name.
    RUN=$(pwd)
    cd $(dirname {output.exN50})
    $TRINITY_HOME/util/misc/contig_ExN50_statistic.pl "$RUN/{input.expression}" "$RUN/{input.transcriptome}" > $(basename {output.exN50}) 2>> "$RUN/{log}"
    $TRINITY_HOME/util/misc/plot_ExN50_statistic.Rscript $(basename {output.exN50}) &>> "$RUN/{log}"
    mv ExN50_plot.pdf $(basename {output.exN50plot})
    """


rule busco:
  """
  Runs BUSCO to assess the completeness of the transcriptome.
  """
  input:
    transcriptome="results/summary/transcriptome.fasta"
  output:
    out_directory=directory("results/qc/assembly/busco"),
    report="results/summary/busco_report.txt"
  log:
    "logs/busco.log"
  conda:
    "../envs/busco.yaml"
  params:
    lineage=LINEAGE,
    # a shared db_dir holds the lineage, and BUSCO must not download or update anything there
    download_option=f"--offline --download_path {DB_DIR}/busco_downloads" if DB_SHARED else "--download_path resources/busco_downloads",
    scratch_var=SCRATCH_VAR
  threads:
    config["busco_threads"]
  resources:
    mem_mb=gb_to_mb(config["busco_memory"]),
    runtime=config["busco_runtime"]
  shell:
    """
    # BUSCO writes thousands of small files. When the scratch key names a variable that is set in this
    # job, it works in a folder on the node's local disk, packs the per-gene folders there and copies its
    # output folder to the run folder at the end. Otherwise it works in the run folder. A failed run copies
    # nothing, because Snakemake deletes the output folder of a failed job; BUSCO's messages are in {log}.
    OUT={output.out_directory}
    BUSCO_PATH=$(dirname $OUT)
    : > {log}
    source {workflow.basedir}/scripts/scratch.sh
    scratch_start "{params.scratch_var}" {log} allow
    if [ -n "$LOCAL_ROOT" ]; then BUSCO_PATH=$WORK; fi
    # Snakemake makes the output folder when it makes the folder of the report in it, and BUSCO makes the folder itself
    rm -rf $OUT

    # an empty or missing lineage lets BUSCO choose the lineage automatically
    if [ -z "{params.lineage}" ]; then
      lineage_option="--auto-lineage"
    else
      lineage_option="-l {params.lineage}"
    fi

    # -o takes a folder name and --out_path the folder that holds it
    busco -m transcriptome -i {input.transcriptome} -o $(basename $OUT) --out_path $BUSCO_PATH $lineage_option -c {threads} {params.download_option} &>> {log}
    echo "BUSCO run completed successfully" &>> {log}

    # pack the per-gene folders (thousands of small files) into one archive each. The archive holds
    # the files in name order with a fixed time and owner, and gzip -n leaves out the name and time
    # of the archive, so the same files give the same archive bytes.
    for d in $BUSCO_PATH/$(basename $OUT)/run_*/busco_sequences $BUSCO_PATH/$(basename $OUT)/run_*/hmmer_output; do
      if [ -d "$d" ]; then
        ( tar --format=gnu --sort=name --mtime=@0 --owner=0 --group=0 --numeric-owner -cf - -C "$(dirname "$d")" "$(basename "$d")" | gzip -n > "$d.tar.gz" ) &>> {log}
        rm -rf "$d" &>> {log}
      fi
    done

    if [ -n "$LOCAL_ROOT" ]; then
      SIZE=$(du -sh "$WORK" | cut -f1)
      START=$(date +%s)
      rm -rf $OUT &>> {log}
      cp -rp $WORK/$(basename $OUT) $(dirname $OUT)/ &>> {log}
      echo "copied $SIZE from $WORK to $OUT in $(( $(date +%s) - START )) s" &>> {log}
    fi

    # with --auto-lineage BUSCO also writes a short_summary.generic.* file, the report is the specific one
    cp $OUT/short_summary.specific.*.txt {output.report} &>> {log}
    """
