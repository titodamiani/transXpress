import os
import subprocess

output, config = snakemake.output, snakemake.config

import platform
import socket
import subprocess

def field(path, key):
  try:
    with open(path) as handle:
      for line in handle:
        if line.startswith(key):
          return line.split(":", 1)[-1].strip()
  except OSError:
    pass
  return ""

cpu = field("/proc/cpuinfo", "model name") or platform.processor() or platform.machine()
memory_total = field("/proc/meminfo", "MemTotal") or "unknown"

os_name = ""
try:
  with open("/etc/os-release") as handle:
    for line in handle:
      if line.startswith("PRETTY_NAME="):
        os_name = line.split("=", 1)[1].strip().strip('"')
except OSError:
  pass
if not os_name:
  os_name = platform.platform()

try:
  commit = subprocess.check_output(
    ["git", "-C", snakemake.params.basedir, "rev-parse", "HEAD"],
    stderr=subprocess.DEVNULL).decode().strip()
except Exception:
  commit = "unknown"

with open(output[0], "w") as out:
  out.write("host\t%s\n" % socket.gethostname())
  out.write("cpu\t%s\n" % cpu)
  out.write("cores\t%s\n" % os.cpu_count())
  out.write("memory_total\t%s\n" % memory_total)
  out.write("os\t%s\n" % os_name)
  out.write("kernel\t%s %s %s\n" % (platform.system(), platform.release(), platform.machine()))
  out.write("python\t%s\n" % platform.python_version())
  out.write("transxpress_commit\t%s\n" % commit)
  out.write("assembler\t%s\n" % config["assembler"])
  out.write("assembler_threads\t%s\n" % (
    config.get("rnaspades_threads", 32) if config["assembler"] == "rnaspades"
    else config.get("trinity_threads", 32)))
