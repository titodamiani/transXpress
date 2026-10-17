rule trinity:
  """
  Runs the whole Trinity assembly in one job. When the scratch key names a variable that is set in
  this job, Trinity works in a folder on the node's local disk, and only the assembly and its gene
  map come to the run folder.
  """
  input:
    samples="results/trimmomatic/samples_trimmed.txt"
  output:
    # Trinity writes these two beside its output folder, not inside it
    transcriptome="results/assembly/trinity.Trinity.fasta",
    gene_trans_map="results/assembly/trinity.Trinity.fasta.gene_trans_map"
  log:
    "logs/trinity.log"
  conda:
    "../envs/trinity_utils.yaml"
  params:
    memory=config.get("trinity_memory", "450"),
    scratch_var=SCRATCH_VAR
  threads:
    config.get("trinity_threads", 32)
  # The common failures here are memory and time. A retry asks for the same figures and starts from the
  # beginning, so it would fail the same way, after hours on a large node.
  retries: 0
  shell:
    """
    # record the computer this ran on, because the assembly depends on it
    ( hostname ; grep -m1 "model name" /proc/cpuinfo 2>/dev/null || uname -m ; nproc 2>/dev/null || true ) > {log}

    # Trinity writes one folder per read cluster in read_partitions, which can reach millions of files.
    # On the local disk they never reach the run folder. A failed run copies nothing back: the files are
    # the ones to keep off the shared disk, and Trinity's messages are in {log}. Without the local disk,
    # Trinity works in results/assembly/trinity, and a rerun after a failure continues from Trinity's own
    # checkpoints there. --full_cleanup deletes the work folder after a successful assembly.
    TRINITY_OUT=results/assembly/trinity
    SCRATCH_VAR="{params.scratch_var}"
    LOCAL_ROOT=""
    if [ -n "$SCRATCH_VAR" ]; then LOCAL_ROOT="${{!SCRATCH_VAR:-}}"; fi
    # A scratch folder in RAM (tmpfs) counts against the memory of the job, and this rule writes large
    # work files, so it stops there
    if [ -n "$LOCAL_ROOT" ] && [ "$(stat -f -c %T "$LOCAL_ROOT" 2>/dev/null)" = "tmpfs" ]; then
      echo "ERROR: the scratch folder $LOCAL_ROOT is in RAM (tmpfs) on $(hostname). Run this rule on a node with a local disk." | tee -a {log} >&2
      exit 1
    fi
    if [ -n "$LOCAL_ROOT" ]; then
      WORK=$(mktemp -d "$LOCAL_ROOT/transxpress.XXXXXX")
      # Trinity accepts only an output folder with "trinity" in its name
      TRINITY_OUT=$WORK/trinity
      echo "work folder: $WORK" &>> {log}
      trap 'rm -rf "$WORK"' EXIT
      trap "exit 143" TERM
    fi

    Trinity --max_memory {params.memory}G --CPU {threads} --output "$TRINITY_OUT" --full_cleanup --samples_file {input.samples} {config[trinity_parameters]} {config[strand_specific]} &>> {log}

    if [ -n "$LOCAL_ROOT" ]; then
      START=$(date +%s)
      cp -p "$TRINITY_OUT.Trinity.fasta" {output.transcriptome} &>> {log}
      cp -p "$TRINITY_OUT.Trinity.fasta.gene_trans_map" {output.gene_trans_map} &>> {log}
      echo "copied the assembly and the gene map from $WORK to results/assembly in $(( $(date +%s) - START )) s" &>> {log}
    fi
    """


rule rnaspades:
  """
  Runs rnaSPAdes assembly. rnaSPAdes writes three transcript sets, and
  rnaspades_filter in config.yaml picks the one that becomes transcriptome.fasta.
  """
  input:
    samples="results/trimmomatic/samples_trimmed.yaml",
  output:
    default="results/assembly/rnaspades/transcripts.fasta",
    hard="results/assembly/rnaspades/hard_filtered_transcripts.fasta",
    soft="results/assembly/rnaspades/soft_filtered_transcripts.fasta"
    ##config["annotated_fasta_prefix"]+"_spades_out/K"+THEKMER+"/assembly_graph_with_scaffolds.gfa"
  log:
    "logs/rnaspades.log"
  conda:
    "../envs/rnaspades.yaml"
  params:
    memory=config.get("rnaspades_memory", "400"),
    scratch_var=SCRATCH_VAR
  threads:
    config.get("rnaspades_threads", 32)
  # The common failures here are memory and time. A retry asks for the same figures and starts from the
  # beginning, so it would fail the same way, after hours on a large node.
  retries: 0
  shell:
    """
    # record the computer this ran on, because the assembly depends on it
    ( hostname ; grep -m1 "model name" /proc/cpuinfo 2>/dev/null || uname -m ; nproc 2>/dev/null || true ) > {log}

    # rnaSPAdes writes about 100 GB of working files on a large dataset. When the scratch key names a
    # variable that is set in this job, it works in a folder on the node's local disk, and its output
    # folder is copied to the run folder at the end. Otherwise it works in the run folder.
    OUT=results/assembly/rnaspades
    SPADES_OUT=$OUT
    SCRATCH_VAR="{params.scratch_var}"
    LOCAL_ROOT=""
    if [ -n "$SCRATCH_VAR" ]; then LOCAL_ROOT="${{!SCRATCH_VAR:-}}"; fi
    # A scratch folder in RAM (tmpfs) counts against the memory of the job, and this rule writes large
    # work files, so it stops there
    if [ -n "$LOCAL_ROOT" ] && [ "$(stat -f -c %T "$LOCAL_ROOT" 2>/dev/null)" = "tmpfs" ]; then
      echo "ERROR: the scratch folder $LOCAL_ROOT is in RAM (tmpfs) on $(hostname). Run this rule on a node with a local disk." | tee -a {log} >&2
      exit 1
    fi
    if [ -n "$LOCAL_ROOT" ]; then
      RUN=$(pwd)
      WORK=$(mktemp -d "$LOCAL_ROOT/transxpress.XXXXXX")
      SPADES_OUT=$WORK/rnaspades
      echo "work folder: $WORK" &>> {log}
      # On a failure, the local folder goes to $OUT, where it would stand without the local disk.
      # A job that Slurm kills at its time limit or with scancel copies nothing. The local folder is deleted
      # in all cases.
      on_exit() {{
        rc=$?
        set +e
        trap - EXIT TERM
        if [ $rc -ne 0 ]; then
          echo "rnaspades stopped with code $rc, so $WORK/rnaspades goes to $OUT" &>> {log}
          mkdir -p "$OUT" &>> {log}
          cp -rp "$WORK/rnaspades/." "$OUT/" &>> {log}
        fi
        rm -rf "$WORK"
        exit $rc
      }}
      trap on_exit EXIT
      trap "exit 143" TERM
    fi

    # No --checkpoints: Slurm deletes the local folder when the job ends, so a restart point would not survive
    rnaspades.py --dataset {input.samples} -t {threads} -m {params.memory} -o "$SPADES_OUT" {config[strand_specific]} --only-assembler {config[kmers]} &>> {log}

    if [ "$SPADES_OUT" != "$OUT" ]; then
      SIZE=$(du -sh "$WORK" | cut -f1)
      START=$(date +%s)
      # replace the whole folder, so no file of an earlier failed run stays beside the new assembly
      rm -rf "$OUT" &>> {log}
      cp -rp "$SPADES_OUT" "$(dirname "$OUT")/" &>> {log}
      echo "copied $SIZE from $WORK to $OUT in $(( $(date +%s) - START )) s" &>> {log}
    fi
    """
 


def transcriptome_copy_input(wildcards):
  if config["assembler"] == "rnaspades":
    return {"transcriptome": getattr(rules.rnaspades.output, RNASPADES_FILTER)}
  return {"transcriptome": rules.trinity.output.transcriptome,
          "gene_trans_map": rules.trinity.output.gene_trans_map}


rule transcriptome_copy:
  """
  Copies the assembled transcriptome to the transxpress directory, with its
  gene to transcript map. Trinity writes the map. For rnaSPAdes the rule builds
  it from the transcript names: NODE_97_length_6237_cov_11.98_g8_i2 is isoform 2
  of gene g8, so the map line is "g8<TAB>NODE_97_length_6237_cov_11.98_g8_i2".
  A name without _g<N>_i<M> at its end stops the rule.
  """
  input:
    unpack(transcriptome_copy_input)
  output:
    transcriptome="results/assembly/transcriptome.fasta",
    gene_trans_map="results/assembly/transcriptome.gene_trans_map"
  log:
    "logs/transcriptome_copy.log"
  run:
    shell("cp -p {input.transcriptome} {output.transcriptome} &> {log}")
    if config["assembler"] == "rnaspades":
      shell("""seqkit seq -n -i {output.transcriptome} 2>> {log} | awk 'BEGIN {{ OFS = "\\t" }} match($0, /_g[0-9]+_i[0-9]+$/) {{ split(substr($0, RSTART + 1), part, "_"); print part[1], $0; next }} {{ print "ERROR: transcript " $0 " has no _g<N>_i<M> at the end of its name, so its gene is unknown" > "/dev/stderr"; bad = 1 }} END {{ exit bad }}' > {output.gene_trans_map} 2>> {log}""")
    else:
      shell("cp -p {input.gene_trans_map} {output.gene_trans_map} &>> {log}")
