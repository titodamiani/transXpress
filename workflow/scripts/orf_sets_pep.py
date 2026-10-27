import Bio.SeqIO

input, output, log = snakemake.input, snakemake.output, snakemake.log

# the protein ID to its sequence without the stop sign
sequences = {}
extra = 0
with open(output[0], "w") as output_handle:
  for orf_set in ["all_orfs", "single_best_orf"]:
    for record in Bio.SeqIO.parse(input[orf_set], "fasta"):
      sequence = str(record.seq).strip("*")
      if record.id in sequences:
        # the same protein in both sets
        assert sequences[record.id] == sequence, "%s has two different sequences in the two ORF sets" % record.id
        continue
      sequences[record.id] = sequence
      if orf_set == "single_best_orf":
        extra += 1
      Bio.SeqIO.write(record, output_handle, "fasta")

with open(log[0], "w") as log_handle:
  print(len(sequences), "proteins,", extra, "of them only in single_best_orf", file=log_handle)
