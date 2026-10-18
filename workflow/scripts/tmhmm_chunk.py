import Bio.SeqIO
import os
import re
import subprocess
import tempfile
from Bio.Seq import Seq

input, output, log = snakemake.input, snakemake.output, snakemake.log

# open the input file, output file and log file
# tmhmm.py writes its result files into the current directory, so each chunk works in its own temporary directory
# in the job's temporary folder (see job_tmpdir). The directory and the per-protein files in it are deleted when
# the chunk is done
with open(input[0], "r") as input_handle, open(output[0], "w") as output_handle, open(log[0], "a+") as log_handle, tempfile.TemporaryDirectory(dir=snakemake.resources["tmpdir"]) as tmp_dir:
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

    # executing the tmhmm.py on created fasta, a failure of tmhmm.py stops the job
    cmd_tmhmm = ['tmhmm', '-f', fasta_file]
    subprocess.run(cmd_tmhmm, stdout = log_handle, stderr = log_handle, cwd = tmp_dir, check = True)

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
