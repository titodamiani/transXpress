rule fastqc:
  """
  Runs fastQC on individual reads files: the input reads (before_trim) or the trimmed reads (after_trim).
  """
  input:
    samples=lambda wildcards: config["samples_file"] if wildcards.stage == "before" else "results/trimmomatic/samples_trimmed.txt"
  output:
    directory("results/qc/fastqc_{stage}_trim")
  wildcard_constraints:
    stage="before|after"
  log:
    "logs/fastqc_{stage}_trim.log"
  conda:
    "../envs/qc.yaml"
  params:
    memory="6"
  threads:
    8
  shell:
    """
    # run fastqc on input files
    # making directory for the fastqc summary files
    mkdir {output} &> {log}
    FILES=$(awk '{{ printf "%s\\n%s\\n", $3,$4}}' {input}) &>> {log}
    # running fastqc on the files, {threads} of them at a time
    fastqc -f fastq -t {threads} -o {output} $FILES &>> {log}
    """


rule multiqc:
  """
  Creates multiQC report from individual fastQC reports, before or after the trimming.
  """
  input:
    "results/qc/fastqc_{stage}_trim"
  output:
    out_dir=directory("results/qc/multiqc_{stage}_trim"),
    report="results/qc/multiqc_{stage}_trim.txt"
  wildcard_constraints:
    stage="before|after"
  log:
    "logs/multiqc_{stage}_trim.log"
  conda:
    "../envs/qc.yaml"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # common error resolved by those two export commands
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    mkdir {output.out_dir} &> {log}
    multiqc -o {output.out_dir} {input} &>> {log}
    cp {output.out_dir}/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """


rule fastqc_before_trim_warnings:
  """
  Prints warnings for FastQC modules which produced some warnings or fails of the input data.
  """
  input:
    fastqc_report="results/qc/multiqc_before_trim.txt"
  output:
    warning_file="results/qc/WARNING_fastqc_before_trim_overview.txt"
  log:
    "logs/fastqc_before_trim_warnings.log"
  script:
    "../scripts/fastqc_warnings.py"
    


rule compare_qc_after_trim:
  """
  Compares the quality of reads before and after trimming.
  Compares PER BASE SEQUENCE QUALITY, PER SEQUENCE QUALITY SCORES, 
  OVERREPRESENTED SEQUENCES and ADAPTER CONTENT modules of FastQC.
  """
  input:
    before="results/qc/multiqc_before_trim.txt",
    after="results/qc/multiqc_after_trim.txt"
  output:
    result="results/qc/FastQC_comparison_after_trim.txt"
  log:
    "logs/compare_qc_after_trim.log"
  run:

    mapping = {
      "per_base_sequence_quality": "PER BASE SEQUENCE QUALITY",
      "per_sequence_quality_scores": "PER SEQUENCE QUALITY SCORES",
      "overrepresented_sequences": "OVERREPRESENTED SEQUENCES",
      "adapter_content": "ADAPTER CONTENT"
    }

    def compare_module(module_name, before_trim, after_trim, out_file):
      module_name = mapping[module_name]
      out_file.writelines(f"""**{module_name}**
Before trimming:
pass: {str(before_trim.count('pass'))} out of {str(len(before_trim))} ({str(before_trim.count('pass')/len(before_trim)*100)}%)
warn: {str(before_trim.count('warn'))} out of {str(len(before_trim))} ({str(before_trim.count('warn')/len(before_trim)*100)}%)
fail: {str(before_trim.count('fail'))} out of {str(len(before_trim))} ({str(before_trim.count('fail')/len(before_trim)*100)}%)
      
After trimming:
pass: {str(after_trim.count('pass'))} out of {str(len(after_trim))} ({str(after_trim.count('pass')/len(after_trim)*100)}%)
warn: {str(after_trim.count('warn'))} out of {str(len(after_trim))} ({str(after_trim.count('warn')/len(after_trim)*100)}%)
fail: {str(after_trim.count('fail'))} out of {str(len(after_trim))} ({str(after_trim.count('fail')/len(after_trim)*100)}%)
{'-' * 50}
""")

    before = pd.read_csv(input["before"], sep='\t')
    before = before.set_index('Sample')

    after = pd.read_csv(input["after"], sep='\t')
    after = after.set_index('Sample')

    modules = ['per_base_sequence_quality', 'per_sequence_quality_scores', 'overrepresented_sequences', 'adapter_content']

    with open(output["result"], "w") as out_file:
      out_file.write("FASTQC COMPARISON BEFORE AND AFTER TRIMMING\n")
      out_file.write('-' * 50 + '\n')

      for module in modules:
        before_list = before[module].tolist()
        after_list = after[module].tolist()
        compare_module(module, before_list, after_list, out_file)

    # print the results also in the terminal
    with open(output["result"], "r") as out_file:
      for line in out_file.readlines():
        print(line.strip())
      print('See this output in ' + output["result"])


rule trinity_stats:
  """
  Runs Trinity script to get statistics about the assembled transcriptome
  (number of transcripts, genes, GC content, EXN50).
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
    expression="results/expression/transcriptome_expression_isoform.tsv"
  output:
    stats="results/qc/transcriptome_stats.txt",
    exN50="results/qc/transcriptome_exN50.tsv",
    exN50plot="results/qc/ExN50_plot.pdf"
  log:
    "logs/trinity_stats.log"
  conda:
    "../envs/trinity_utils.yaml"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')
    # contig_ExN50_statistic.pl reads the genes in Perl hash order, which is random in each run, and
    # genes tied on expression keep that order. A fixed Perl hash seed gives the same order in each run.
    export PERL_HASH_SEED=0

    assembler={config[assembler]}
    if [ "$assembler" = 'trinity' ]; then
        $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} > {output.stats} 2> {log}
    else
        $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} | sed -e 's/trinity/rnaspades/g' > {output.stats} 2> {log}
    fi
    # contig_ExN50_statistic.pl writes its .E-inputs table and the R script writes ExN50_plot.pdf into the
    # folder they run in, so both run in the folder of the outputs
    RUN=$(pwd)
    cd $(dirname {output.exN50})
    $TRINITY_HOME/util/misc/contig_ExN50_statistic.pl "$RUN/{input.expression}" "$RUN/{input.transcriptome}" > $(basename {output.exN50}) 2>> "$RUN/{log}"
    $TRINITY_HOME/util/misc/plot_ExN50_statistic.Rscript $(basename {output.exN50}) &>> "$RUN/{log}"
    """


rule busco:
  """
  Runs BUSCO to assess the completeness of the transcriptome.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta"
  output:
    out_directory=directory("results/qc/busco"),
    report="results/qc/busco_report.txt"
  log:
    "logs/busco.log"
  conda:
    "../envs/busco.yaml"
  params:
    memory=config.get("busco_memory", "10"),
    lineage=LINEAGE,
    # a shared db_dir holds the lineage, and BUSCO must not download or update anything there
    download_option=("--offline --download_path " + DB_DIR + "/busco_downloads") if DB_SHARED else "--download_path resources/busco_downloads",
    scratch_var=SCRATCH_VAR,
    scratch_sh=SCRATCH_SH
  threads:
    config.get("busco_threads", 16)
  shell:
    """
    # BUSCO writes thousands of small files. When the scratch key names a variable that is set in this
    # job, it works in a folder on the node's local disk, packs the per-gene folders there and copies its
    # output folder to the run folder at the end. Otherwise it works in the run folder. A failed run copies
    # nothing, because Snakemake deletes the output folder of a failed job; BUSCO's messages are in {log}.
    OUT={output.out_directory}
    BUSCO_PATH=$(dirname $OUT)
    : > {log}
    source {params.scratch_sh}
    scratch_start "{params.scratch_var}" {log} allow
    if [ -n "$LOCAL_ROOT" ]; then BUSCO_PATH=$WORK; fi

    # an empty or missing lineage lets BUSCO choose the lineage automatically
    if [ -z "{params.lineage}" ]; then
      lineage_option="--auto-lineage"
    else
      lineage_option="-l {params.lineage}"
    fi

    # -o takes a folder name and --out_path the folder that holds it
    busco -m transcriptome -i {input.transcriptome} -o $(basename $OUT) --out_path $BUSCO_PATH $lineage_option -c {threads} {params.download_option} &>> {log}
    echo "BUSCO run completed successfully" &>> {log}

    # pack the per-gene folders (thousands of small files) into one archive each
    for d in $BUSCO_PATH/$(basename $OUT)/run_*/busco_sequences $BUSCO_PATH/$(basename $OUT)/run_*/hmmer_output; do
      if [ -d "$d" ]; then
        tar -czf "$d.tar.gz" -C "$(dirname "$d")" "$(basename "$d")" &>> {log}
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
