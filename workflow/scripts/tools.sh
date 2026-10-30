# Shell functions that several rules share. A rule sources this file in its shell block.

# The computer a job ran on: host name, processor and CPU count.
record_host() {
  hostname
  grep -m1 "model name" /proc/cpuinfo 2>/dev/null || uname -m
  nproc 2>/dev/null || true
}

# The folder of the Trinity install. The Trinity command in the env is a link into it.
trinity_home() {
  python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))'
}

# publish_part FILE LOG: FILE.part gets the name FILE when its last line is "# [ok]", the end of a finished
# hmmscan or cmscan. Otherwise the job fails.
publish_part() {
  tail -n 1 "$1.part" | grep -qx '# \[ok\]' || { echo "$1.part does not end with # [ok]" >> "$2"; exit 1; }
  mv "$1.part" "$1"
}
