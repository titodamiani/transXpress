# a more elegant way to do this is:
# seqkit seq -i {input} | seqkit replace -s -p "\*" -r "" | seqkit split -f -s 1 -O {output}
# however, there seems to be a problem in seqkit: https://github.com/shenwei356/seqkit/issues/65
checkpoint fasta_split_fasta:
  """
  Splits the transcriptome.fasta into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "results/assembly/transcriptome.fasta"
  output:
    directory("results/annotation/chunks/fasta")
  log:
    "logs/fasta_split_fasta.log"
  params:
    memory="2",
    chunk_size=ANNOTATION_CHUNK["fasta"]
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % params.chunk_size == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / params.chunk_size) + 1);
          filename = os.path.join(output[0], fileindex + ".fasta")
          output_handle = open(filename, "w");
        
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        


checkpoint fasta_split_orfs:
  """
  Splits the transcriptome.orfs into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "results/transdecoder/transcriptome.orfs"
  output:
    directory("results/annotation/chunks/orfs")
  log:
    "logs/fasta_split_orfs.log"
  params:
    memory="2",
    chunk_size=ANNOTATION_CHUNK["orfs"]
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % params.chunk_size == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / params.chunk_size) + 1);
          filename = os.path.join(output[0], fileindex + ".orfs")
          output_handle = open(filename, "w");
        
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        


checkpoint fasta_split_pep:
  """
  Splits the transcriptome.pep into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "results/transdecoder/transcriptome.pep"
  output:
    directory("results/annotation/chunks/pep")
  log:
    "logs/fasta_split_pep.log"
  params:
    memory="2",
    chunk_size=ANNOTATION_CHUNK["pep"]
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % params.chunk_size == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / params.chunk_size) + 1);
          filename = os.path.join(output[0], fileindex + ".pep")
          output_handle = open(filename, "w");
        # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan) 
        #if wildcards["extension"] == "pep":
        record.seq = record.seq.strip("*")
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        


def parallel_annotation_tasks_fasta(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_fasta.get(**wildcards).output[0]
  # sorted by number, so the merged file follows the order of the sequences at any chunk size.
  # glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "{index}.fasta")).index, key=int)
  completed_files = expand("results/annotation/chunks/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files


def parallel_annotation_tasks_orfs(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_orfs.get(**wildcards).output[0]
  # sorted by number, so the merged file follows the order of the sequences at any chunk size.
  # glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "{index}.orfs")).index, key=int)
  completed_files = expand("results/annotation/chunks/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files


def parallel_annotation_tasks_pep(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_pep.get(**wildcards).output[0]
  # sorted by number, so the merged file follows the order of the sequences at any chunk size.
  # glob_wildcards returns the folder order, which differs between runs
  job_ids = sorted(glob_wildcards(os.path.join(parallel_dir, "{index}.pep")).index, key=int)
  completed_files = expand("results/annotation/chunks/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files


rule annotation_merge_fasta:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_fasta
  output:
    "results/annotation/{task}_fasta.out"
  log:
    "logs/annotation_merge_fasta/{task}.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """


rule annotation_merge_orfs:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_orfs
  output:
    "results/annotation/{task}_orfs.out"
  log:
    "logs/annotation_merge_orfs/{task}.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """


rule annotation_merge_pep:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_pep
  output:
    "results/annotation/{task}_pep.out"
  log:
    "logs/annotation_merge_pep/{task}.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """


rule rfam_parallel:
  """
  Runs cmscan on Rfam database on smaller nucleotide files in parallel.
  """
  input:
    fasta="results/annotation/chunks/fasta/{index}.fasta",
    db=DB_DIR + "/Rfam.cm"
  output:
    "results/annotation/chunks/rfam/{index}.out"
  log:
    "logs/rfam_parallel/{index}.log"
  conda:
    "../envs/rfam.yaml"
  params:
    memory="5"
  threads:
    4
  shell:
    """
    cmscan -E {config[e_value_threshold]} --rfam --cpu {threads} --tblout {output} {input[db]} {input[fasta]} &> {log}
    """


rule pfam_parallel:
  """
  Runs hmmscan on Pfam database on smaller protein files in parallel.
  """
  input:
    fasta="results/annotation/chunks/orfs/{index}.orfs",
    db=DB_DIR + "/Pfam-A.hmm"
  output:
    "results/annotation/chunks/pfam/{index}.out"
  log:
    "logs/pfam_parallel/{index}.log"
  conda:
    "../envs/pfam.yaml"
  params:
    memory="2"
  threads:
    2
  shell:
    """
    # Transdecoder requires --domtblout output
    hmmscan -E {config[e_value_threshold]} --cpu {threads} --domtblout {output} {input[db]} {input[fasta]} &> {log}
    """


rule sprot_blastp_parallel:
  """
  Runs blast search on SwissProt database on smaller protein files in parallel.
  """
  input:
    fasta="results/annotation/chunks/orfs/{index}.orfs",
    db=DB_DIR + "/uniprot_sprot.fasta"
  output:
    "results/annotation/chunks/sprotblastp/{index}.out"
  log:
    "logs/sprot_blastp_parallel/{index}.log"
  conda:
    "../envs/blast.yaml"
  params:
    memory="2"
  threads:
    2
  shell:
    """
    blastp -query {input[fasta]} -db {input[db]} -num_threads {threads} -evalue {config[e_value_threshold]} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output} &> {log}
    """


rule sprot_blastx_parallel:
  """
  Runs blast search on SwissProt database on smaller nucleotide files in parallel.
  """
  input:
    fasta="results/annotation/chunks/fasta/{index}.fasta",
    db=DB_DIR + "/uniprot_sprot.fasta"
  output:
    "results/annotation/chunks/sprotblastx/{index}.out"
  log:
    "logs/sprot_blastx_parallel/{index}.log"
  conda:
    "../envs/blast.yaml"
  params:
    memory="2"
  threads:
    4
  shell:
    """
    blastx -query {input[fasta]} -db {input[db]} -num_threads {threads} -evalue {config[e_value_threshold]} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output} &> {log}
    """


rule tmhmm_parallel:
  """
  Runs tmhmm.py on smaller protein files in parallel to predict transmembrane
  topology.
  """
  input:
    "results/annotation/chunks/pep/{index}.pep"
  output:
    "results/annotation/chunks/tmhmm/{index}.out"
  log:
    "logs/tmhmm_parallel/{index}.log"
  params:
    memory="2"
  threads:
    1
  run:
    # open the input file, output file and log file
    # tmhmm.py writes its result files into the current directory, so each chunk works in its own temporary directory
    # in the job's temporary folder (see job_tmpdir). The directory and the per-protein files in it are deleted when
    # the chunk is done
    with open(input[0], "r") as input_handle, open(output[0], "a+") as output_handle, open(log[0], "a+") as log_handle, tempfile.TemporaryDirectory(dir=job_tmpdir(wildcards)) as tmp_dir:
      log_handle.write("work folder: " + tmp_dir + "\n")
      log_handle.flush()

      # iterate through individual sequences in input file, tmhmm.py can be executed with just one sequence at a time
      for record in Bio.SeqIO.parse(input_handle, "fasta"):

        # individual fasta headers must have an ID and description, separated by a space, otherwise error occurs
        # see: https://github.com/dansondergaard/tmhmm.py/issues/16
        # adding artificial "description" in header
        record.description = " x"

        # tmhmm.py does not handle 'X' in the sequence (happens when using RNASpades as an assembler and 'N' occurs in the transcript sequence)
        # see: https://github.com/dansondergaard/tmhmm.py/issues/9
        # removing 'X' from the sequence
        record.seq = Seq(re.sub('X',"",str(record.seq)))

        # saving this fasta file
        fasta_file =  os.path.join(tmp_dir, record.id + ".fasta")
        Bio.SeqIO.write(record, fasta_file, "fasta")

        # executing the tmhmm.py on created fasta
        cmd_tmhmm = ['tmhmm', '-f', fasta_file]
        subprocess.run(cmd_tmhmm, stdout = log_handle, stderr = log_handle, cwd = tmp_dir)

        # 3 files are created afterwards in tmp_dir
        # record.id.summary     record.id.annotation    record.id.plot

        summary_file = os.path.join(tmp_dir, record.id + ".summary")
        annotation_file = os.path.join(tmp_dir, record.id + ".annotation")
        plot_file = os.path.join(tmp_dir, record.id + ".plot")

        # creating shorter output from summary file
        with open(summary_file, "r") as inp:
          final_topology = ''
          num_of_helices = 0

          for line in inp:
            # split the line of format "start end topology" into list [start, end, topology]
            split_line = re.split(r'\s+', line, 2)
            start = split_line[0].strip()
            end = split_line[1].strip()
            topology = split_line[2].strip()

            # update final string based on the topology
            if (topology == 'inside'):
              final_topology += 'i'
            if (topology == 'outside'):
              final_topology += 'o'
            if ('transmembrane helix' in topology):
              final_topology += start + '-' + end
              num_of_helices += 1

          output_handle.write(record.id + "\tPredHel=" + str(num_of_helices) + "\tTopology=" + final_topology + "\n")
        os.remove(summary_file)
        os.remove(annotation_file)
        os.remove(plot_file)


rule signalp_parallel:
  """
  Runs signalp on smaller protein files in parallel to predict signal peptides.
  """
  input:
    "results/annotation/chunks/pep/{index}.pep"
  output:
    "results/annotation/chunks/signalp/{index}.out"
  log:
    "logs/signalp_parallel/{index}.log"
  conda:
    "signalp"
  params:
    memory="4",
    organism=SIGNALP_ORGANISM
  threads:
    1
  shell:
    """
    # Snakemake sets TMPDIR to the job's temporary folder, see job_tmpdir
    TMP=$(mktemp -d "${{TMPDIR:-/tmp}}/signalp.XXXXXX")
    echo "work folder: $TMP" &> {log}
    signalp6 --fastafile {input} --organism {params.organism} --output_dir $TMP/out --format none --mode fast &>> {log}
    mv $TMP/out/prediction_results.txt {output}
    rm -r $TMP
    """


rule targetp_parallel:
  """
  Runs TargetP on smaller protein files to predict presence and type of
  targeting peptide.
  """
  input:
    "results/annotation/chunks/pep/{index}.pep"
  output:
    "results/annotation/chunks/targetp/{index}.out"
  log:
    "logs/targetp_parallel/{index}.log"
  conda:
    "targetp"
  params:
    memory="3"
  threads:
    1
  shell:
    """
    # targetp writes its result into the folder it runs in, which is the job's temporary folder
    # (Snakemake sets TMPDIR, see job_tmpdir)
    RUN=$(pwd)
    TMP=$(mktemp -d "${{TMPDIR:-/tmp}}/targetp.XXXXXX")
    echo "work folder: $TMP" &> "$RUN/{log}"
    cd $TMP
    targetp -fasta "$RUN/{input}" -format short -org {config[targetp]} -prefix targetp_{wildcards.index} &>> "$RUN/{log}"
    mv targetp_{wildcards.index}_summary.targetp2 "$RUN/{output}" &>> "$RUN/{log}"
    cd "$RUN"
    rm -r $TMP
    """


rule annotated_fasta:
  """
  Puts annotations in the headers of transcripts/proteins in the 
  .fasta/.pep transcriptome files.
  """
  input:
    transcriptome="results/assembly/transcriptome.fasta",
    proteome="results/transdecoder/transcriptome.pep",
    expression="results/expression/transcriptome_expression_isoform.tsv",
    blastx_results="results/annotation/sprotblastx_fasta.out",
    rfam_results="results/annotation/rfam_fasta.out",
    blastp_results="results/annotation/sprotblastp_orfs.out",
    pfam_results="results/annotation/pfam_orfs.out",
    tmhmm_results="results/annotation/tmhmm_pep.out",
    signalp_results=["results/annotation/signalp_pep.out"] if SIGNALP_ORGANISM else [],
    targetp_results=["results/annotation/targetp_pep.out"] if TARGETP_ORGANISM else []
  output:
    transcriptome_annotated="results/annotation/transcriptome_annotated.fasta",
    proteome_annotated="results/annotation/transcriptome_annotated.pep",
    tpm_blast_table="results/annotation/transcriptome_TPM_blast.csv"
  log:
    "logs/annotated_fasta.log"
  params:
    memory="4"
  threads:
    1
  run:
    # Open log file
    with open(log[0], "w") as log_handle:
      ## Annotation map: transcript id -> description
      expression_annotations = {}
      blastx_annotations = {}
      blastp_annotations = {}
      pfam_annotations = {}
      rfam_annotations = {}
      tmhmm_annotations = {}
      signalp_annotations = {}
      targetp_annotations = {}
  
      ## Load kallisto results
      print ("Loading expression values from", input["expression"], file=log_handle)
      with open(input["expression"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter='\t')
        columns = next(csv_reader)
        for row in csv_reader:
          expression_annotations[row[0]] = columns[1] + "=" + str(row[1])
          for i in range(2, len(columns)):
            expression_annotations[row[0]] += " " + columns[i] + "=" + str(row[i])

      ## Load blastx results
      print ("Loading blastx results from", input["blastx_results"], file=log_handle)
      with open(input["blastx_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 13): continue
          blastx_annotations[row[0]] = row[12] + " E=" + str(row[10])
  
      ## Load blastp results
      print ("Loading blastp results from", input["blastp_results"], file=log_handle)
      with open(input["blastp_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 13): continue
          blastp_annotations[row[0]] = row[12] + " E=" + str(row[10])

      ## Load pfam results
      print ("Loading pfam predictions from", input["pfam_results"], file=log_handle)
      with open(input["pfam_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = re.split(" +", line, 22)
          if (len(row) < 23): continue
          if row[3] not in pfam_annotations:
            pfam_annotations[row[3]] = row[1] + " " + row[22] + "E=" + str(row[6])

      ## Load rfam results
      print ("Loading rfam predictions from", input["rfam_results"], file=log_handle)
      with open(input["rfam_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = re.split(" +", line, 17)
          if (len(row) < 18): continue
          rfam_annotations[row[2]] = row[1] + " " + row[17] + "E=" + str(row[15])
  
      ## Load tmhmm results
      print ("Loading tmhmm predictions from", input["tmhmm_results"], file=log_handle)
      with open(input["tmhmm_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) > 2):
            tmhmm_annotations[row[0]] = row[1] + " " + row[2]

      ## Load signalp results
      if input["signalp_results"]:
        print ("Loading signalp predictions from", input["signalp_results"], file=log_handle)
        with open(input["signalp_results"]) as input_handle:
          csv_reader = csv.reader(input_handle, delimiter="\t")
          for row in csv_reader:
            if (len(row) < 2): continue
            prediction = str(row[1])
            if prediction == "OTHER":
              signalp_annotations[row[0]] = "OTHER (No SP)"
            elif prediction == "SP":
              signalp_annotations[row[0]] = "standard secretory signal peptides (Sec/SPI)"
            elif prediction == "LIPO":
              signalp_annotations[row[0]] = "lipoprotein signal peptides (Sec/SPII)"
            elif prediction == "TAT":
              signalp_annotations[row[0]] = "Tat signal peptides (Tat/SPI)"
            elif prediction == "TATLIPO":
              signalp_annotations[row[0]] = "Tat lipoprotein signal peptides (Tat/SPII)"
            elif prediction == "PILIN":
              signalp_annotations[row[0]] = "Pilin and pilin-like signal peptides (Sec/SPIII)"
            else:
              signalp_annotations[row[0]] = str(row[1])

      ## Load targetp results
      if input["targetp_results"]:
        translation =    {"SP": "Signal peptide",
                          "mTP": "Mitochondrial transit peptide",
                          "cTP": "chloroplast transit peptide",
                          "luTP": "thylakoidal lumen composite transit peptide",
                          "noTP": "no targeting peptide" }
        print("Loading targetp predictions from", input["targetp_results"], file=log_handle)
        with open(input["targetp_results"]) as input_handle:
          for line in input_handle:
            if (line.startswith("#")): continue
            row = line.split()
            if (len(row) >= 8):
              targetp_annotations[row[0]] = translation[row[1]] + ", " + " ".join(row[7:])
            elif (len(row) >= 2):
              targetp_annotations[row[0]] = translation[row[1]]

      ## Do the work
      print ("Annotating FASTA file", input["transcriptome"], "to", output["transcriptome_annotated"], file=log_handle)
      with open(input["transcriptome"], "r") as input_fasta_handle, open(output["transcriptome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = record.id
          record.description = "TPM: " + expression_annotations.get(transcript_id)
          if transcript_id in blastx_annotations:
            record.description += "; blastx: " + blastx_annotations.get(transcript_id)
          if transcript_id in rfam_annotations:
            record.description += "; rfam: " + rfam_annotations.get(transcript_id)
          Bio.SeqIO.write(record, output_fasta_handle, "fasta")
      
      print ("Annotating FASTA file", input["proteome"], "to", output["proteome_annotated"], file=log_handle)
      with open(input["proteome"], "r") as input_fasta_handle, open(output["proteome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = re.sub("\\.p[0-9]+$", "", record.id)
          record.description = "transdecoder: " + re.search("ORF type:([^,]+,score=[^,]+)", record.description).group(1)
          if transcript_id in expression_annotations:
            record.description += "; TPM: " + expression_annotations.get(transcript_id)
          if record.id in blastp_annotations:
            record.description += "; blastp: " + blastp_annotations.get(record.id)
          if record.id in pfam_annotations:
            record.description += "; pfam: " + pfam_annotations.get(record.id)
          if transcript_id in rfam_annotations:
            record.description += "; rfam: " + rfam_annotations.get(transcript_id)
          if record.id in tmhmm_annotations:
            record.description += "; tmhmm: " + tmhmm_annotations.get(record.id)
          if record.id in signalp_annotations:
            record.description += "; signalp: " + signalp_annotations.get(record.id)
          if record.id in targetp_annotations:
            record.description += "; targetp: " + targetp_annotations.get(record.id)
          # Add sequence ID prefix from configuration
          if config["annotated_fasta_prefix"]:
            record.id = config["annotated_fasta_prefix"] + "|" + record.id
          Bio.SeqIO.write(record, output_fasta_handle, "fasta")

      print ("Generating transcriptome_TPM_blast.csv table", file=log_handle)
      with open(input["expression"], "r") as input_csv_handle, open(output["tpm_blast_table"], "w") as output_csv_handle:
        csv_reader = csv.reader(input_csv_handle, delimiter="\t")
        csv_writer = csv.writer(output_csv_handle, delimiter=",")
        csv_columns = next(csv_reader)
        csv_columns[0] = "transcript"
        for i in range(1, len(csv_columns)):
          csv_columns[i] = "TPM(" + csv_columns[i] + ")"
        csv_columns.append("blastx")
        csv_writer.writerow(csv_columns)
        for row in csv_reader:
          row.append(blastx_annotations.get(row[0], ""))
          csv_writer.writerow(row)
