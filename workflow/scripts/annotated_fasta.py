import Bio.SeqIO
import re
import csv

input, output, log, config = snakemake.input, snakemake.output, snakemake.log, snakemake.config

SIGNALP_NAMES = {
  "OTHER": "OTHER (No SP)",
  "SP": "standard secretory signal peptides (Sec/SPI)",
  "LIPO": "lipoprotein signal peptides (Sec/SPII)",
  "TAT": "Tat signal peptides (Tat/SPI)",
  "TATLIPO": "Tat lipoprotein signal peptides (Tat/SPII)",
  "PILIN": "Pilin and pilin-like signal peptides (Sec/SPIII)"}

TARGETP_NAMES = {
  "SP": "Signal peptide",
  "mTP": "Mitochondrial transit peptide",
  "cTP": "chloroplast transit peptide",
  "luTP": "thylakoidal lumen composite transit peptide",
  "noTP": "no targeting peptide"}

def blast_hits(path):
  """
  Sequence ID to "<subject title> E=<e-value>", from a blast table with the columns "6 std stitle".
  """
  hits = {}
  with open(path) as handle:
    for row in csv.reader(handle, delimiter="\t"):
      if len(row) >= 13:
        hits[row[0]] = row[12] + " E=" + row[10]
  return hits

def hmmer_hits(path, columns, id_column, name_column, evalue_column):
  """
  Sequence ID to "<target name> <description>E=<e-value>" of the first hit of each sequence, from a
  table of hmmscan or cmscan, whose description is the last of its columns.
  """
  hits = {}
  with open(path) as handle:
    for line in handle:
      if line.startswith("#"):
        continue
      row = re.split(" +", line, columns - 1)
      # the tools list the hits of each sequence from the best one down
      if len(row) == columns and row[id_column] not in hits:
        hits[row[id_column]] = row[name_column] + " " + row[columns - 1] + "E=" + row[evalue_column]
  return hits

with open(log[0], "w") as log_handle:
  # transcript ID to "<sample>=<TMM> <sample>=<TMM> ...", the TPM scaled across the libraries with TMM
  print("Loading expression values from", input["expression"], file=log_handle)
  expression_annotations = {}
  with open(input["expression"]) as input_handle:
    csv_reader = csv.reader(input_handle, delimiter="\t")
    columns = next(csv_reader)
    for row in csv_reader:
      expression_annotations[row[0]] = " ".join(columns[i] + "=" + row[i] for i in range(1, len(columns)))

  print("Loading blastx results from", input["blastx_results"], file=log_handle)
  blastx_annotations = blast_hits(input["blastx_results"])

  print("Loading blastp results from", input["blastp_results"], file=log_handle)
  blastp_annotations = blast_hits(input["blastp_results"])

  print("Loading pfam predictions from", input["pfam_results"], file=log_handle)
  pfam_annotations = hmmer_hits(input["pfam_results"], 23, 3, 1, 6)

  print("Loading rfam predictions from", input["rfam_results"], file=log_handle)
  rfam_annotations = hmmer_hits(input["rfam_results"], 18, 2, 1, 15)

  print("Loading tmhmm predictions from", input["tmhmm_results"], file=log_handle)
  tmhmm_annotations = {}
  with open(input["tmhmm_results"]) as input_handle:
    for row in csv.reader(input_handle, delimiter="\t"):
      if len(row) > 2:
        tmhmm_annotations[row[0]] = row[1] + " " + row[2]

  signalp_annotations = {}
  if input["signalp_results"]:
    print("Loading signalp predictions from", input["signalp_results"][0], file=log_handle)
    with open(input["signalp_results"][0]) as input_handle:
      for row in csv.reader(input_handle, delimiter="\t"):
        if len(row) >= 2:
          signalp_annotations[row[0]] = SIGNALP_NAMES.get(row[1], row[1])

  targetp_annotations = {}
  if input["targetp_results"]:
    print("Loading targetp predictions from", input["targetp_results"][0], file=log_handle)
    with open(input["targetp_results"][0]) as input_handle:
      for line in input_handle:
        if line.startswith("#"):
          continue
        row = line.split()
        if len(row) >= 8:
          targetp_annotations[row[0]] = TARGETP_NAMES[row[1]] + ", " + " ".join(row[7:])
        elif len(row) >= 2:
          targetp_annotations[row[0]] = TARGETP_NAMES[row[1]]

  print("Annotating FASTA file", input["transcriptome"], "to", output["transcriptome_annotated"], file=log_handle)
  with open(input["transcriptome"], "r") as input_fasta_handle, open(output["transcriptome_annotated"], "w") as output_fasta_handle:
    for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
      transcript_id = record.id
      record.description = "TMM: " + expression_annotations[transcript_id]
      if transcript_id in blastx_annotations:
        record.description += "; blastx: " + blastx_annotations[transcript_id]
      if transcript_id in rfam_annotations:
        record.description += "; rfam: " + rfam_annotations[transcript_id]
      Bio.SeqIO.write(record, output_fasta_handle, "fasta")

  print("Annotating FASTA file", input["proteome"], "to", output["proteome_annotated"], file=log_handle)
  with open(input["proteome"], "r") as input_fasta_handle, open(output["proteome_annotated"], "w") as output_fasta_handle:
    for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
      # a TransDecoder protein ID is the transcript ID with .p<N> added
      transcript_id = re.sub("\\.p[0-9]+$", "", record.id)
      record.description = "transdecoder: " + re.search("ORF type:([^,]+,score=[^,]+)", record.description).group(1)
      for name, annotations, key in [
          ("TMM", expression_annotations, transcript_id),
          ("blastp", blastp_annotations, record.id),
          ("pfam", pfam_annotations, record.id),
          ("rfam", rfam_annotations, transcript_id),
          ("tmhmm", tmhmm_annotations, record.id),
          ("signalp", signalp_annotations, record.id),
          ("targetp", targetp_annotations, record.id)]:
        if key in annotations:
          record.description += "; %s: %s" % (name, annotations[key])
      if config["annotated_fasta_prefix"]:
        record.id = config["annotated_fasta_prefix"] + "|" + record.id
      Bio.SeqIO.write(record, output_fasta_handle, "fasta")
