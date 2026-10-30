import datetime
import glob
import json
import os
import platform
import re
import socket
import struct
import subprocess
import sys
import uuid

output, config, params = snakemake.output, snakemake.config, snakemake.params
# the workflow folder of this clone: the folder above scripts/
basedir = os.path.dirname(snakemake.scriptdir)

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
    ["git", "-C", basedir, "rev-parse", "HEAD"],
    stderr=subprocess.DEVNULL).decode().strip()
except Exception:
  commit = "unknown"

def image_info(path):
  """
  An image file as "<path>, <size> bytes, sif id <ID>, built <UTC time>". The ID and the build time come
  from the global header of a SIF file (https://github.com/apptainer/sif): 32 bytes of launch script,
  the magic "SIF_MAGIC", 3 bytes of version, 3 bytes of architecture, the 16-byte ID and the build time
  in seconds as a little-endian int64.
  """
  if not path:
    return ""
  if not os.path.isfile(path):
    return "%s, no such file" % path
  info = "%s, %d bytes" % (path, os.path.getsize(path))
  with open(path, "rb") as handle:
    header = handle.read(72)
  if len(header) == 72 and header[32:41] == b"SIF_MAGIC":
    built = datetime.datetime.fromtimestamp(struct.unpack("<q", header[64:72])[0], datetime.timezone.utc)
    info += ", sif id %s, built %s" % (uuid.UUID(bytes=header[48:64]), built.strftime("%Y-%m-%d %H:%M:%S UTC"))
  return info

def env_packages(env_file):
  """
  The package names of an env file, and the exact pip packages (name==version) of its pip section.
  """
  names, pip = [], []
  section = ""
  in_pip = False
  with open(env_file) as handle:
    for line in handle:
      line = line.split("#", 1)[0].rstrip()
      key = re.match(r"^([\w-]+):", line)
      if key:
        section = key.group(1)
        continue
      item = re.match(r"^(\s*)- (.+)$", line)
      if section != "dependencies" or not item:
        continue
      spec = item.group(2).strip().strip("\"'")
      if spec == "pip:":
        in_pip = True
      elif in_pip and len(item.group(1)) > 2:
        pip.append(spec)
      else:
        in_pip = False
        names.append(re.split(r"[=<>!\[ ]", spec, 1)[0])
  return names, pip

def pinned_versions(pin_file):
  """
  Package name to "<version>=<build>", from the URLs of a pin file (conda list --explicit).
  """
  versions = {}
  with open(pin_file) as handle:
    for line in handle:
      line = line.strip()
      if not line or line.startswith(("#", "@")):
        continue
      archive = line.split("#", 1)[0].rsplit("/", 1)[-1]
      archive = re.sub(r"\.(conda|tar\.bz2)$", "", archive)
      name, version, build = archive.rsplit("-", 2)
      versions[name] = "%s=%s" % (version, build)
  return versions

def main_env_versions():
  """
  Package name to "<version>=<build>" of the env this script runs in, the main env, from its conda-meta.
  """
  versions = {}
  for record in glob.glob(os.path.join(sys.prefix, "conda-meta", "*.json")):
    with open(record) as handle:
      meta = json.load(handle)
    versions[meta["name"]] = "%s=%s" % (meta["version"], meta["build"])
  return versions

# The tools of each env file at the versions of its pin file. default.yaml is the main env. An env file with
# a name: line (signalp, targetp) is a named env that the run does not build.
platform_name = {("Linux", "x86_64"): "linux-64", ("Linux", "aarch64"): "linux-aarch64",
  ("Darwin", "x86_64"): "osx-64", ("Darwin", "arm64"): "osx-arm64"}.get((platform.system(), platform.machine()), "")
envs_dir = os.path.join(basedir, "envs")
tools = []
for env_file in sorted(glob.glob(os.path.join(envs_dir, "*.yaml"))):
  with open(env_file) as handle:
    if re.search(r"^name:", handle.read(), re.MULTILINE):
      continue
  env = os.path.basename(env_file)
  names, pip = env_packages(env_file)
  if env == "default.yaml":
    versions, missing = main_env_versions(), "not in the main env"
  else:
    pin_file = os.path.join(envs_dir, env.replace(".yaml", ".%s.pin.txt" % platform_name))
    versions = pinned_versions(pin_file) if os.path.isfile(pin_file) else {}
    missing = "not pinned for " + (platform_name or platform.platform())
  for name in names:
    tools.append((env, name, versions.get(name, missing)))
  for spec in pip:
    name, _, version = spec.partition("==")
    tools.append((env, name, version + "=pip" if version else "not pinned"))

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
  out.write("assembler_threads\t%s\n" % config["%s_threads" % config["assembler"]])
  out.write("deployment\t%s\n" % params.deployment)
  out.write("image\t%s\n" % params.image)
  out.write("image_file\t%s\n" % (image_info(params.image_file) if params.image else ""))
  out.write("signalp_container\t%s\n" % image_info(params.signalp_container))
  out.write("targetp_container\t%s\n" % image_info(params.targetp_container))
  for env, name, version in tools:
    out.write("tool\t%s\t%s\t%s\n" % (env, name, version))
