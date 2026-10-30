#!/usr/bin/env python3
"""Compare the results/ folders of two transXpress runs, file by file.

Usage: compare_runs.py RUN_A RUN_B [--allowed tests/allowed_diffs.tsv] [--quiet]
       compare_runs.py RUN_A RUN_B --sequences results/summary/transcriptome.fasta

Two files match when their bytes match (a gzip file: its decompressed bytes), after the run folder, the
node-local work folder and the conda env folder in the content become @RUN@, @SCRATCH@ and @ENV@.
A file that differs passes only when a line of the allowed-diffs file names it: pattern, mode (skip,
strip:REGEX, sorted:REGEX), reason, separated by tabs. See CONTRIBUTING.md, section "Tests".
With --sequences only the sorted sequences of that fasta file are compared.
Exit code 0: every file matches or is allowed to differ. 1: a file differs, or is in one run only.
"""
import argparse
import collections
import concurrent.futures
import fnmatch
import gzip
import hashlib
import re
import sys
import time
from pathlib import Path


SCRATCH_FOLDER = re.compile(rb"/scratch/slurm_jobs/[^/\s]+/\d+/transxpress\.\w+")
# a path that ends in the folder of a Snakemake conda env: 32 hex digits, and "_" in some Snakemake versions.
# A match starts only at the start of a word, so a long line (a sequence) takes linear time.
ENV_FOLDER = re.compile(rb"(?<![^\s\"'()])[^\s\"'()]*/[0-9a-f]{32}_?(?=/)")


def read_allowed(path):
    rules = []
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 3 or not (fields[1] == "skip" or fields[1].startswith(("strip:", "sorted:"))):
            sys.exit("%s line %d: expected pattern, tab, mode (skip, strip:REGEX or sorted:REGEX), tab, reason"
                     % (path, number))
        rules.append((fields[0], fields[1], fields[2]))
    return rules


def content(path, run):
    """The bytes of a file. A gzip file gives its decompressed bytes. The run folder path becomes @RUN@,
    a scratch folder @SCRATCH@ and an env folder @ENV@."""
    data = path.read_bytes()
    if data[:2] == b"\x1f\x8b":
        try:
            data = gzip.decompress(data)
        except (OSError, EOFError):
            pass
    data = data.replace(str(Path(run).resolve()).encode(), b"@RUN@")
    data = SCRATCH_FOLDER.sub(b"@SCRATCH@", data)
    return ENV_FOLDER.sub(b"@ENV@", data)


def raw_equal(path_a, path_b):
    """True when both files have the same bytes. It reads 1 MiB at a time, so a big file stays out of memory."""
    if path_a.stat().st_size != path_b.stat().st_size:
        return False
    with open(path_a, "rb") as fa, open(path_b, "rb") as fb:
        while True:
            chunk = fa.read(1 << 20)
            if chunk != fb.read(1 << 20):
                return False
            if not chunk:
                return True


def md5(data):
    return hashlib.md5(data).hexdigest()


def stripped(data, mode):
    """The content a strip: or sorted: rule compares."""
    regex = mode.split(":", 1)[1]
    pattern = re.compile(regex.encode())
    lines = [line for line in data.split(b"\n") if not pattern.search(line)]
    if mode.startswith("sorted:"):
        lines = sorted(re.sub(b" +", b" ", line) for line in lines)
    return b"\n".join(lines)


def first_difference(a, b):
    a_lines = a.split(b"\n")
    b_lines = b.split(b"\n")
    for number, (x, y) in enumerate(zip(a_lines, b_lines), 1):
        if x != y:
            # from 20 bytes before the first byte that differs
            start = max(0, next(i for i in range(len(x) + 1) if x[i:i + 1] != y[i:i + 1]) - 20)
            return "line %d, byte %d: %r != %r" % (number, start + 1, x[start:start + 80], y[start:start + 80])
    return "line counts %d and %d" % (len(a_lines), len(b_lines))


def results_files(run):
    root = Path(run) / "results"
    if not root.is_dir():
        sys.exit("%s has no results folder" % run)
    return {str(p.relative_to(run)): p for p in sorted(root.rglob("*")) if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_a")
    parser.add_argument("run_b")
    parser.add_argument("--allowed", default=str(Path(__file__).with_name("allowed_diffs.tsv")))
    parser.add_argument("--quiet", action="store_true", help="print only the files that fail")
    parser.add_argument("--sequences", metavar="FASTA", help="compare only the sorted sequences of this fasta file")
    args = parser.parse_args()

    if args.sequences:
        sets = []
        for run in (args.run_a, args.run_b):
            data = (Path(run) / args.sequences).read_bytes()
            assert data.startswith(b">") or not data, "%s does not start with a header" % (Path(run) / args.sequences)
            # a record starts with ">" at the start of a line; a ">" inside a header does not cut it
            lines = data[1:].split(b"\n>") if data else []
            sets.append(sorted(b"".join(record.split(b"\n")[1:]) for record in lines))
        print("sequences: %d and %d" % (len(sets[0]), len(sets[1])))
        if sets[0] == sets[1]:
            print("OK")
            return 0
        print("FAILED: the sets of sequences differ (%d sequences are in both)"
              % sum((collections.Counter(sets[0]) & collections.Counter(sets[1])).values()))
        return 1

    rules = read_allowed(args.allowed) if Path(args.allowed).is_file() else []
    files_a = results_files(args.run_a)
    files_b = results_files(args.run_b)
    same = 0
    allowed = {}
    failed = []

    def compare(name):
        """('same',), ('allowed', pattern) or ('failed', message) for one file."""
        rule = next((r for r in rules if fnmatch.fnmatchcase(name, r[0])), None)
        if name not in files_a or name not in files_b:
            # a skip rule also lets a file exist in one run only
            if rule is not None and rule[1] == "skip":
                return ("allowed", rule[0])
            return ("failed", "only in %s: %s" % (args.run_a if name in files_a else args.run_b, name))
        if rule is not None and rule[1] == "skip":
            return ("allowed", rule[0])
        # Equal bytes are equal after the replacements too, except when a file holds the run folder path
        # in both runs: the paths differ, so the raw bytes differ and the slow path below decides.
        if raw_equal(files_a[name], files_b[name]):
            return ("same",)
        a = content(files_a[name], args.run_a)
        b = content(files_b[name], args.run_b)
        if a == b:
            return ("same",)
        if rule is None:
            return ("failed", "DIFF %s (%s)" % (name, first_difference(a, b)))
        if stripped(a, rule[1]) == stripped(b, rule[1]):
            return ("allowed", rule[0])
        return ("failed", "DIFF %s (%s; allowed only for mode %s)"
                % (name, first_difference(stripped(a, rule[1]), stripped(b, rule[1])), rule[1]))

    names = sorted(set(files_a) | set(files_b))
    start = time.time()
    # 4 threads: the files are on a shared disk. map keeps the order of names, so the output does not change.
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for number, (name, result) in enumerate(zip(names, pool.map(compare, names)), 1):
            if result[0] == "same":
                same += 1
            elif result[0] == "allowed":
                allowed.setdefault(result[1], []).append(name)
            else:
                failed.append(result[1])
            if number % 100 == 0:
                print("%d of %d files, %d s" % (number, len(names), time.time() - start), file=sys.stderr, flush=True)

    print("same content: %d files" % same)
    if not args.quiet:
        for pattern, reason in [(r[0], r[2]) for r in rules]:
            if pattern in allowed:
                print("allowed to differ: %d files, %s (%s)" % (len(allowed[pattern]), pattern, reason))
    for line in failed:
        print(line)
    if failed:
        print("FAILED: %d files" % len(failed))
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
