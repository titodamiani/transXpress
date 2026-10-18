import Bio.SeqIO
import re
import csv

input, output, log, config = snakemake.input, snakemake.output, snakemake.log, snakemake.config

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
      # cmscan lists the hits of each transcript from the best one down, so the first hit is kept
      if row[2] not in rfam_annotations:
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
    print ("Loading signalp predictions from", input["signalp_results"][0], file=log_handle)
    with open(input["signalp_results"][0]) as input_handle:
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
    print("Loading targetp predictions from", input["targetp_results"][0], file=log_handle)
    with open(input["targetp_results"][0]) as input_handle:
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
