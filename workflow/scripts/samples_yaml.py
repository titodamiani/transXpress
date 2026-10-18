import os
import pprint
import re

input, output = snakemake.input, snakemake.output

# rnaSPAdes reads a relative path in this file from the folder of the file, so each read path is
# written relative to that folder
yaml_dir = os.path.dirname(output["samples_yaml"])
sample_list = []
with open(input["samples_txt"], "r") as input_handle:
  for line in input_handle:
    row = re.split("[\t ]+", line)
    if len(row) > 3: # paired reads
      f_reads = row[2].strip()
      r_reads = row[3].strip()
      assert os.path.isfile(f_reads)
      assert os.path.isfile(r_reads)
      sample_list.append({
        "orientation": "fr",
        "type": "paired-end",
        "left reads": [os.path.relpath(f_reads, yaml_dir)],
        "right reads": [os.path.relpath(r_reads, yaml_dir)]})
    if len(row) == 3: # unpaired reads
      u_reads = row[2].strip()
      assert os.path.isfile(u_reads)
      sample_list.append({
        "type": "single",
        "single reads": [os.path.relpath(u_reads, yaml_dir)]})

with open(output["samples_yaml"], "w") as output_handle:
  output_handle.write(pprint.pformat(sample_list))
