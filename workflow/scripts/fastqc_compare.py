import pandas as pd

input, output = snakemake.input, snakemake.output

# The FastQC modules that the trimming can change, and their titles in the report
MODULES = {
  "per_base_sequence_quality": "PER BASE SEQUENCE QUALITY",
  "per_sequence_quality_scores": "PER SEQUENCE QUALITY SCORES",
  "overrepresented_sequences": "OVERREPRESENTED SEQUENCES",
  "adapter_content": "ADAPTER CONTENT"}

def counts(results):
  """
  One line for each FastQC state: the number of reads files in that state, and their percentage.
  """
  return "".join(f"{state}: {results.count(state)} out of {len(results)} ({results.count(state) / len(results) * 100}%)\n"
    for state in ["pass", "warn", "fail"])

before = pd.read_csv(input["before"], sep="\t").set_index("Sample")
after = pd.read_csv(input["after"], sep="\t").set_index("Sample")

with open(output["result"], "w") as out_file:
  out_file.write("FASTQC COMPARISON BEFORE AND AFTER TRIMMING\n")
  out_file.write("-" * 50 + "\n")
  for module, title in MODULES.items():
    # a line of 6 spaces separates the two stages in the report
    out_file.write(f"**{title}**\nBefore trimming:\n{counts(before[module].tolist())}      \n"
      f"After trimming:\n{counts(after[module].tolist())}{'-' * 50}\n")

# print the results also in the terminal
with open(output["result"], "r") as out_file:
  for line in out_file.readlines():
    print(line.strip())
  print("See this output in " + output["result"])
