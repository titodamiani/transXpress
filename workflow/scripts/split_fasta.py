import os

import Bio.SeqIO

input, output, params = snakemake.input, snakemake.output, snakemake.params

os.makedirs(output[0], exist_ok=True)
with open(input[0], "r") as input_handle:
  output_handle = None
  count = 0
  for record in Bio.SeqIO.parse(input_handle, "fasta"):
    if count % params.chunk_size == 0:
      if output_handle is not None:
        output_handle.close()
      filename = os.path.join(output[0], "%d.%s" % (count // params.chunk_size + 1, params.extension))
      output_handle = open(filename, "w")

    # the stop codons (*) at the ends of the proteins go, because some annotation tools reject them
    if params.strip_stop:
      record.seq = record.seq.strip("*")
    # the headers keep only the sequence ID, because some tools copy the description into their output
    record.description = ""
    Bio.SeqIO.write(record, output_handle, "fasta")
    count += 1
  if output_handle is not None:
    output_handle.close()
