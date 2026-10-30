rule trinity:
  """
  Runs the Trinity assembly in one job, on the node-local disk when the scratch variable is set.
  """
  input:
    samples="results/trimmomatic/samples_trimmed.txt"
  output:
    # Trinity writes these two beside its output folder
    transcriptome="results/assembly/trinity.Trinity.fasta",
    gene_trans_map="results/assembly/trinity.Trinity.fasta.gene_trans_map"
  log:
    "logs/trinity.log"
  conda:
    "../envs/trinity_utils.yaml"
  params:
    scratch_var=SCRATCH_VAR,
    trinity_parameters=config["trinity_parameters"],
    strand_specific=config["strand_specific"]
  threads:
    config["trinity_threads"]
  resources:
    mem_mb=gb_to_mb(config["trinity_memory"]),
    runtime=config["trinity_runtime"]
  # a retry asks for the same memory and time and fails the same way
  retries: 0
  shell:
    """
    # record the computer this ran on, because the assembly depends on it
    source {workflow.basedir}/scripts/tools.sh
    record_host > {log}

    # read_partitions can hold millions of files, so Trinity works on the local disk and a failed run copies nothing back.
    # Without it, Trinity works in results/assembly/trinity and a rerun continues from its checkpoints.
    # --full_cleanup deletes the work folder after a successful assembly.
    TRINITY_OUT=results/assembly/trinity
    source {workflow.basedir}/scripts/scratch.sh
    scratch_start "{params.scratch_var}" {log} stop
    if [ -n "$LOCAL_ROOT" ]; then
      # Trinity wants "trinity" in the name of the output folder
      TRINITY_OUT=$WORK/trinity
    fi

    Trinity --max_memory $(( {resources.mem_mb} / 1024 ))G --CPU {threads} --output "$TRINITY_OUT" --full_cleanup --samples_file {input.samples} {params.trinity_parameters} {params.strand_specific} &>> {log}

    if [ -n "$LOCAL_ROOT" ]; then
      START=$(date +%s)
      cp -p "$TRINITY_OUT.Trinity.fasta" {output.transcriptome} &>> {log}
      cp -p "$TRINITY_OUT.Trinity.fasta.gene_trans_map" {output.gene_trans_map} &>> {log}
      echo "copied the assembly and the gene map from $WORK to results/assembly in $(( $(date +%s) - START )) s" &>> {log}
    fi
    """


rule rnaspades:
  """
  Runs rnaSPAdes. It writes three transcript sets, and rnaspades_filter picks the one that becomes transcriptome.fasta.
  """
  input:
    samples="results/trimmomatic/samples_trimmed.yaml",
  output:
    default="results/assembly/rnaspades/transcripts.fasta",
    hard="results/assembly/rnaspades/hard_filtered_transcripts.fasta",
    soft="results/assembly/rnaspades/soft_filtered_transcripts.fasta"
  log:
    "logs/rnaspades.log"
  conda:
    "../envs/rnaspades.yaml"
  params:
    scratch_var=SCRATCH_VAR,
    strand_specific=config["strand_specific"],
    kmers=config["kmers"]
  threads:
    config["rnaspades_threads"]
  resources:
    mem_mb=gb_to_mb(config["rnaspades_memory"]),
    runtime=config["rnaspades_runtime"]
  # a retry asks for the same memory and time and fails the same way
  retries: 0
  shell:
    """
    # record the computer this ran on, because the assembly depends on it
    source {workflow.basedir}/scripts/tools.sh
    record_host > {log}

    # about 100 GB of working files on a large dataset: work on the local disk and copy the output folder back
    OUT=results/assembly/rnaspades
    SPADES_OUT=$OUT
    source {workflow.basedir}/scripts/scratch.sh
    on_failure() {{
      echo "rnaspades stopped with code $rc, so $WORK/rnaspades goes to $OUT" &>> {log}
      mkdir -p "$OUT" &>> {log}
      cp -rp "$WORK/rnaspades/." "$OUT/" &>> {log}
    }}
    scratch_start "{params.scratch_var}" {log} stop on_failure
    if [ -n "$LOCAL_ROOT" ]; then
      SPADES_OUT=$WORK/rnaspades
    fi

    # No --checkpoints: Slurm deletes the local folder at the end of the job
    rnaspades.py --dataset {input.samples} -t {threads} -m $(( {resources.mem_mb} / 1024 )) -o "$SPADES_OUT" {params.strand_specific} --only-assembler {params.kmers} &>> {log}

    if [ "$SPADES_OUT" != "$OUT" ]; then
      SIZE=$(du -sh "$WORK" | cut -f1)
      START=$(date +%s)
      # replace the whole folder, so no file of an earlier failed run stays
      rm -rf "$OUT" &>> {log}
      cp -rp "$SPADES_OUT" "$(dirname "$OUT")/" &>> {log}
      echo "copied $SIZE from $WORK to $OUT in $(( $(date +%s) - START )) s" &>> {log}
    fi
    """


rule transcriptome_copy:
  """
  Copies the assembly to results/summary with its gene to transcript map. Trinity writes the map. For
  rnaSPAdes the map comes from the names: NODE_97_length_6237_cov_11.98_g8_i2 is gene g8; a name without
  _g<N>_i<M> at its end stops the rule. The Trinity copy is sorted by name in natural order (TRINITY_DN2
  before TRINITY_DN10), because Trinity writes the records in another order in each run.

  A rule of its own, so a failed name check does not delete an assembly (Snakemake deletes the outputs of a
  failed rule) and a change of rnaspades_filter reruns only this rule.
  """
  input:
    unpack(transcriptome_copy_input)
  output:
    transcriptome="results/summary/transcriptome.fasta",
    gene_trans_map="results/summary/transcriptome.gene_trans_map"
  log:
    "logs/transcriptome_copy.log"
  threads:
    1
  resources:
    mem_mb=1024,
    runtime="4h"
  run:
    if config["assembler"] == "rnaspades":
      shell("cp -p {input.transcriptome} {output.transcriptome} &> {log}")
      shell("""seqkit seq -n -i {output.transcriptome} 2>> {log} | awk 'BEGIN {{ OFS = "\\t" }} match($0, /_g[0-9]+_i[0-9]+$/) {{ split(substr($0, RSTART + 1), part, "_"); print part[1], $0; next }} {{ print "ERROR: transcript " $0 " has no _g<N>_i<M> at the end of its name, so its gene is unknown" > "/dev/stderr"; bad = 1 }} END {{ exit bad }}' > {output.gene_trans_map} 2>> {log}""")
    else:
      shell("seqkit sort --natural-order {input.transcriptome} > {output.transcriptome} 2> {log}")
      shell("""seqkit seq -n -i {output.transcriptome} 2>> {log} | awk 'BEGIN {{ OFS = "\\t" }} NR == FNR {{ gene[$2] = $1; next }} $1 in gene {{ print gene[$1], $1; next }} {{ print "ERROR: transcript " $1 " is not in the Trinity gene map" > "/dev/stderr"; bad = 1 }} END {{ exit bad }}' {input.gene_trans_map} - > {output.gene_trans_map} 2>> {log}""")
