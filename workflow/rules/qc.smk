rule fastqc:
  """
  FastQC on the input reads (before_trim) or the trimmed reads (after_trim).
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
    FILES=$(awk '{{ printf "%s\\n%s\\n", $3,$4}}' {input})
    # A fonts.conf with only the fonts of the env makes the graphs the same on every host and in the image (no system fonts).
    # Java caches fonts in <user.home>/.java/fonts, so user.home is the job's folder.
    TMP=$(mktemp -d "${{TMPDIR:-/tmp}}/fastqc.XXXXXX")
    trap 'rm -rf "$TMP"' EXIT
    ENV=$(dirname $(dirname $(command -v fastqc)))
    printf '<?xml version="1.0"?>\\n<!DOCTYPE fontconfig SYSTEM "urn:fontconfig:fonts.dtd">\\n<fontconfig>\\n<dir>%s/fonts</dir>\\n<dir>%s/share/fonts</dir>\\n<include ignore_missing="yes">%s/etc/fonts/conf.d</include>\\n<cachedir>%s/cache</cachedir>\\n</fontconfig>\\n' \
      "$ENV" "$ENV" "$ENV" "$TMP" > "$TMP/fonts.conf"
    export FONTCONFIG_FILE="$TMP/fonts.conf" JAVA_TOOL_OPTIONS="-Duser.home=$TMP"
    fastqc -f fastq -t {threads} -o {output} $FILES &>> {log}
    """


rule multiqc:
  """
  MultiQC report of the FastQC reports before or after trimming.
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
    # Snakemake makes the folder of the report
    mkdir -p {output.out_dir} &> {log}
    multiqc -o {output.out_dir} {input} &>> {log}
    cp {output.out_dir}/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """


rule fastqc_before_trim_warnings:
  """
  Lists the FastQC modules that warn or fail on the input reads.
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
  Compares four FastQC modules (base quality, sequence quality scores, overrepresented sequences, adapter content) before and after trimming.
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
  Statistics of the transcriptome: counts, GC content, ExN50.
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
    # contig_ExN50_statistic.pl keeps Perl hash order among genes tied on expression; a fixed seed makes it repeatable
    export PERL_HASH_SEED=0

    # TrinityStats.pl names the assembler "trinity" in its output
    $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} 2> {log} | sed -e 's/trinity/{params.assembler}/g' > {output.stats}
    # both tools write into the folder they run in, so they run in the folder of the outputs; the plot has a fixed name
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
    # BUSCO writes thousands of small files: it works on the local disk when the scratch variable is set
    # and copies its output folder back at the end. A failed run copies nothing.
    OUT={output.out_directory}
    BUSCO_PATH=$(dirname $OUT)
    : > {log}
    source {workflow.basedir}/scripts/scratch.sh
    scratch_start "{params.scratch_var}" {log} allow
    if [ -n "$LOCAL_ROOT" ]; then BUSCO_PATH=$WORK; fi
    # BUSCO makes its output folder itself
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

    # one archive for each per-gene folder; name order, fixed time and owner, gzip -n: the same files give the same bytes
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

    # --auto-lineage also writes short_summary.generic.*; the report is the specific one
    cp $OUT/short_summary.specific.*.txt {output.report} &>> {log}
    """
