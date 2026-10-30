import Bio.SeqIO
import os
import re
import subprocess
import tempfile
from Bio.Seq import Seq

input, output, log = snakemake.input, snakemake.output, snakemake.log

# tmhmm.py takes one sequence at a time and writes its result files into the current directory, so each
# chunk works in its own temporary directory in the job's temporary folder (see job_tmpdir). The directory
# and the files in it are deleted when the chunk is done.
# The result goes to <chunk>.out.part and gets the name of the chunk only when every sequence is done, so a job
# that is killed or hangs leaves no cut file under the name Snakemake trusts.
part_file = output[0] + ".part"
with open(input[0], "r") as input_handle, open(part_file, "w") as output_handle, open(log[0], "w") as log_handle, tempfile.TemporaryDirectory(dir=snakemake.resources["tmpdir"]) as tmp_dir:
  log_handle.write("work folder: " + tmp_dir + "\n")
  log_handle.flush()

  for record in Bio.SeqIO.parse(input_handle, "fasta"):
    # tmhmm.py needs a description after the ID in the header
    # (https://github.com/dansondergaard/tmhmm.py/issues/16)
    record.description = " x"
    # tmhmm.py does not accept X in a sequence, which an N in an rnaSPAdes transcript gives
    # (https://github.com/dansondergaard/tmhmm.py/issues/9)
    record.seq = Seq(str(record.seq).replace("X", ""))

    fasta_file = os.path.join(tmp_dir, record.id + ".fasta")
    Bio.SeqIO.write(record, fasta_file, "fasta")
    # a failure of tmhmm.py stops the job
    subprocess.run(["tmhmm", "-f", fasta_file], stdout=log_handle, stderr=log_handle, cwd=tmp_dir, check=True)

    # tmhmm.py writes <ID>.summary, <ID>.annotation and <ID>.plot. Each line of the summary is
    # "<start> <end> <region>": the region is inside, outside or a transmembrane helix.
    summary_file = os.path.join(tmp_dir, record.id + ".summary")
    topology = ""
    helices = 0
    with open(summary_file, "r") as summary:
      for line in summary:
        start, end, region = [field.strip() for field in re.split(r"\s+", line, 2)]
        if region == "inside":
          topology += "i"
        if region == "outside":
          topology += "o"
        if "transmembrane helix" in region:
          topology += start + "-" + end
          helices += 1
    output_handle.write(record.id + "\tPredHel=" + str(helices) + "\tTopology=" + topology + "\n")
    for extension in ["summary", "annotation", "plot"]:
      os.remove(os.path.join(tmp_dir, record.id + "." + extension))

os.replace(part_file, output[0])
