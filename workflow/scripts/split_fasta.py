import os

import Bio.SeqIO

input, output, params = snakemake.input, snakemake.output, snakemake.params

os.makedirs(output[0], exist_ok=True)
with open(input[0], "r") as input_handle:
  output_handle = None
  count = 0
  for record in Bio.SeqIO.parse(input_handle, "fasta"):
    if (count % params.chunk_size == 0):
      if output_handle is not None:
        output_handle.close()
      fileindex = str(int(count / params.chunk_size) + 1)
      filename = os.path.join(output[0], fileindex + "." + params.extension)
      output_handle = open(filename, "w")

    # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan)
    if params.strip_stop:
      record.seq = record.seq.strip("*")
    # Strip the description, because some tools (e.g. deeploc) include it in their output
    record.description = ""
    Bio.SeqIO.write(record, output_handle, "fasta")
    count += 1
  if output_handle is not None:
    output_handle.close()
