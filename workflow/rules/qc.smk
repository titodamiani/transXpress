rule fastqc_before_trim:
  """
  Runs fastQC on individual input reads files.
  """
  input:
    samples=config["samples_file"]
  output:
    directory("results/qc/fastqc_before_trim")
  log:
    "logs/fastqc_before_trim.log"
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


rule multiqc_before_trim:
  """
  Creates multiQC report from individual fastQC reports.
  """
  input:
    "results/qc/fastqc_before_trim"
  output:
    out_dir=directory("results/qc/multiqc_before_trim"),
    report="results/qc/multiqc_before_trim.txt"
  log:
    "logs/multiqc_before_trim.log"
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
  run:

    mapping = {
    "basic_statistics": "BASIC STATISTICS",
    "per_base_sequence_quality": "PER BASE SEQUENCE QUALITY",
    "per_tile_sequence_quality": "PER TILE SEQUENCE QUALITY",
    "per_sequence_quality_scores": "PER SEQUENCE QUALITY SCORES",
    "per_base_sequence_content": "PER BASE SEQUENCE CONTENT",
    "per_sequence_gc_content": "PER SEQUENCE GC CONTENT",
    "per_base_n_content": "PER BASE N CONTENT",
    "sequence_length_distribution": "SEQUENCE LENGTH DISTRIBUTION",
    "sequence_duplication_levels": "SEQUENCE DUPLICATION LEVELS",
    "overrepresented_sequences": "OVERREPRESENTED SEQUENCES",
    "adapter_content": "ADAPTER CONTENT"
    }

    reasons_mapping = {
    "basic_statistics": "BASIC STATISTICS",
    
    "per_base_sequence_quality":
    """A warning is issued if the lower quartile for any base is less than 10, or if the median for any base is less than 25.
This module will raise a failure if the lower quartile for any base is less than 5 or if the median for any base is less than 20.

Common reasons for warnings:
The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

If the quality of the library falls to a low level then the most common remedy is to perform quality trimming where reads are truncated based on their average quality. For most libraries where this type of degradation has occurred you will often be simultaneously running into the issue of adapter read-through so a combined adapter and quality trimming step is often employed.

Another possibility is that a warn / error is triggered because of a short loss of quality earlier in the run, which then recovers to produce later good quality sequence. This can happen if there is a transient problem with the run (bubbles passing through a flowcell for example). You can normally see this type of error by looking at the per-tile quality plot (if available for your platform). In these cases trimming is not advisable as it will remove later good sequence, but you might want to consider masking bases during subsequent mapping or assembly.

If your library has reads of varying length then you can find a warning or error is triggered from this module because of very low coverage for a given base range. Before committing to any action, check how many sequences were responsible for triggering an error by looking at the sequence length distribution module results.
""",
    
    "per_tile_sequence_quality":
"""This module issues a warning if any tile shows a mean Phred score more than 2 less than the mean for that base across all tiles.
This module will issue a warning if any tile shows a mean Phred score more than 5 less than the mean for that base across all tiles.

Common reasons for warnings:
Whilst warnings in this module can be triggered by individual specific events we have also observed that greater variation in the phred scores attributed to tiles can also appear when a flowcell is generally overloaded. In this case events appear all over the flowcell rather than being confined to a specific area or range of cycles. We would generally ignore errors which mildly affected a small number of tiles for only 1 or 2 cycles, but would pursue larger effects which showed high deviation in scores, or which persisted for several cycles.
""",
    
    "per_sequence_quality_scores": 
"""A warning is raised if the most frequently observed mean quality is below 27 - this equates to a 0.2% error rate.
An error is raised if the most frequently observed mean quality is below 20 - this equates to a 1% error rate.

Common reasons for warnings:
This module is generally fairly robust and errors here usually indicate a general loss of quality within a run. For long runs this may be alleviated through quality trimming. If a bi-modal, or complex distribution is seen then the results should be evaluated in concert with the per-tile qualities (if available) since this might indicate the reason for the loss in quality of a subset of sequences.
""",
    
    "per_base_sequence_content": 
"""Common reasons for warnings:
The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

If the quality of the library falls to a low level then the most common remedy is to perform quality trimming where reads are truncated based on their average quality. For most libraries where this type of degradation has occurred you will often be simultaneously running into the issue of adapter read-through so a combined adapter and quality trimming step is often employed.

Another possibility is that a warn / error is triggered because of a short loss of quality earlier in the run, which then recovers to produce later good quality sequence. This can happen if there is a transient problem with the run (bubbles passing through a flowcell for example). You can normally see this type of error by looking at the per-tile quality plot (if available for your platform). In these cases trimming is not advisable as it will remove later good sequence, but you might want to consider masking bases during subsequent mapping or assembly.

If your library has reads of varying length then you can find a warning or error is triggered from this module because of very low coverage for a given base range. Before committing to any action, check how many sequences were responsible for triggering an error by looking at the sequence length distribution module results.
""",
    
    "per_sequence_gc_content": 
"""A warning is raised if the sum of the deviations from the normal distribution represents more than 15% of the reads.
This module will indicate a failure if the sum of the deviations from the normal distribution represents more than 30% of the reads.

Common reasons for warnings:
Warnings in this module usually indicate a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example), which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species.
""",
    
    "per_base_n_content": 
"""This module raises a warning if any position shows an N content of >5%.
This module will raise an error if any position shows an N content of >20%.

Common reasons for warnings:
The most common reason for the inclusion of significant proportions of Ns is a general loss of quality, so the results of this module should be evaluated in concert with those of the various quality modules. You should check the coverage of a specific bin, since it is possible that the last bin in this analysis could contain very few sequences, and an error could be prematurely triggered in this case.

Another common scenario is the incidence of a high proportions of N at a small number of positions early in the library, against a background of generally good quality. Such deviations can occur when you have very biased sequence composition in the library to the point that base callers can become confused and make poor calls. This type of problem will be apparent when looking at the per-base sequence content results.
""",
    
    "sequence_length_distribution": 
"""This module will raise a warning if all sequences are not the same length.
This module will raise an error if any of the sequences have zero length.

Common reasons for warnings:
For some sequencing platforms it is entirely normal to have different read lengths so warnings here can be ignored.""",
    
    "sequence_duplication_levels": 
"""This module will issue a warning if non-unique sequences make up more than 20% of the total.
This module will issue a error if non-unique sequences make up more than 50% of the total.

Common reasons for warnings:
The underlying assumption of this module is of a diverse unenriched library. Any deviation from this assumption will naturally generate duplicates and can lead to warnings or errors from this module.
        
In general there are two potential types of duplicate in a library, technical duplicates arising from PCR artefacts, or biological duplicates which are natural collisions where different copies of exactly the same sequence are randomly selected. From a sequence level there is no way to distinguish between these two types and both will be reported as duplicates here.
        
A warning or error in this module is simply a statement that you have exhausted the diversity in at least part of your library and are re-sequencing the same sequences. In a supposedly diverse library this would suggest that the diversity has been partially or completely exhausted and that you are therefore wasting sequencing capacity. However in some library types you will naturally tend to over-sequence parts of the library and therefore generate duplication and will therefore expect to see warnings or error from this module.
        
In RNA-Seq libraries sequences from different transcripts will be present at wildly different levels in the starting population. In order to be able to observe lowly expressed transcripts it is therefore common to greatly over-sequence high expressed transcripts, and this will potentially create large set of duplicates. This will result in high overall duplication in this test, and will often produce peaks in the higher duplication bins. This duplication will come from physically connected regions, and an examination of the distribution of duplicates in a specific genomic region will allow the distinction between over-sequencing and general technical duplication, but these distinctions are not possible from raw fastq files. A similar situation can arise in highly enriched ChIP-Seq libraries although the duplication there is less pronounced. Finally, if you have a library where the sequence start points are constrained (a library constructed around restriction sites for example, or an unfragmented small RNA library) then the constrained start sites will generate huge dupliction levels which should not be treated as a problem, nor removed by deduplication. In these types of library you should consider using a system such as random barcoding to allow the distinction of technical and biological duplicates.
""",
    
    "overrepresented_sequences": 
"""This module will issue a warning if any sequence is found to represent more than 0.1% of the total.
This module will issue an error if any sequence is found to represent more than 1% of the total.
        
Common reasons for warnings:
This module will often be triggered when used to analyse small RNA libraries where sequences are not subjected to random fragmentation, and the same sequence may natrually be present in a significant proportion of the library.
""",
    
    "adapter_content": 
"""This module will issue a warning if any sequence is present in more than 5% of all reads.
This module will issue an error if any sequence is present in more than 10% of all reads.
        
Common reasons for warnings:
Any library where a reasonable proportion of the insert sizes are shorter than the read length will trigger this module. This does not indicate a problem as such - just that the sequences will need to be adapter trimmed before proceeding with any downstream analysis.
"""
    }



    with open(output["warning_file"], "w") as out_file:
      df1 = pd.read_csv(input["fastqc_report"], sep='\t')
      df1 = df1.set_index('Sample')

      # filter only columns containing 'pass','warn' or 'fail'
      df1_filtered = df1.loc[:, df1.isin(['pass','warn','fail']).all()]
      df1_filtered = df1_filtered.loc[:, df1_filtered.columns != 'basic_statistics'] # basic statistics module always passes

      n_samples = df1.shape[0]

      for column in df1_filtered.columns:
        if column in mapping.keys():
          out_file.writelines(mapping[column]+'\n')
        else:
          out_file.writelines(column+'\n')
        results = df1_filtered[column].tolist()
        if ('warn' in results) or ('fail' in results):
          n_warn = results.count('warn')
          n_fail = results.count('fail')
          n_warn_perc = round(n_warn/n_samples*100,2)
          n_fail_perc = round(n_fail/n_samples*100,2)
          if n_warn > 0:
            out_file.writelines(f'{n_warn} out of {n_samples} produced "WARN" ({n_warn_perc}%)\n')
          if n_fail > 0:
            out_file.writelines(f'{n_fail} out of {n_samples} produced "FAIL" ({n_fail_perc}%)\n')
          if column in reasons_mapping.keys():
            out_file.writelines('\n' + reasons_mapping[column])
        else:
          out_file.writelines('ok'+'\n')
      out_file.writelines('source: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/\n')
    
    # print to stdout
    with open(output["warning_file"], "r") as out_file:
      for line in out_file.readlines():
        print(line.strip())
      print('See this output in: ' + output["warning_file"])
    


rule fastqc_after_trim:
  """
  Runs fastQC on individual trimmed reads files.
  """
  input:
    samples="results/trimmomatic/samples_trimmed.txt"
  output:
    directory("results/qc/fastqc_after_trim")
  log:
    "logs/fastqc_after_trim.log"
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


rule multiqc_after_trim:
  """
  Creates multiQC report from individual fastQC reports after the trimming.
  """
  input:
    "results/qc/fastqc_after_trim"
  output:
    out_directory=directory("results/qc/multiqc_after_trim"),
    report="results/qc/multiqc_after_trim.txt"
  log:
    "logs/multiqc_after_trim.log"
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
    mkdir {output.out_directory} &> {log}
    multiqc -o {output.out_directory} {input} &>> {log}
    cp {output.out_directory}/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """


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
    scratch_var=SCRATCH_VAR
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
    SCRATCH_VAR="{params.scratch_var}"
    LOCAL_ROOT=""
    if [ -n "$SCRATCH_VAR" ]; then LOCAL_ROOT="${{!SCRATCH_VAR:-}}"; fi
    : > {log}
    if [ -n "$LOCAL_ROOT" ]; then
      WORK=$(mktemp -d "$LOCAL_ROOT/transxpress.XXXXXX")
      BUSCO_PATH=$WORK
      echo "work folder: $WORK" &>> {log}
      trap 'rm -rf "$WORK"' EXIT
    fi

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
