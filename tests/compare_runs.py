#!/usr/bin/env python3
"""Compare the results/ folders of two transXpress runs, file by file.

Usage: compare_runs.py RUN_A RUN_B [--allowed tests/allowed_diffs.tsv] [--quiet]
       compare_runs.py RUN_A RUN_B --sequences results/assembly/transcriptome.fasta

Every file under RUN_A/results must exist in RUN_B/results, and the other way round.
Two files match when their bytes match. A gzip file matches when its decompressed bytes match.
Before the comparison, two paths are replaced in the content, because many tools write the folder
they ran in: the absolute path of each run folder becomes "@RUN@", and the node-local work folder of
a job (/scratch/slurm_jobs/<user>/<job id>/transxpress.<random>) becomes "@SCRATCH@".

A file that differs passes only when a line of the allowed-diffs file names it. Each line has
three tab-separated fields: a path pattern (fnmatch, relative to the run folder, "*" also
matches "/"), a mode, and the reason the file may differ. The mode is one of:
  skip         the file may differ in any way, or exist in one run only
  strip:REGEX  the file matches after the lines that REGEX matches are dropped from both
Lines that start with "#" and empty lines are comments.

Exit code 0: every file matches or is allowed to differ. 1: a file differs, or is in one run only.

With --sequences FASTA, only that fasta file of each run is compared: the sorted sequences, without
the headers. Exit code 0 when the two sets of sequences are the same.
"""
import argparse
import collections
import fnmatch
import gzip
import hashlib
import re
import sys
from pathlib import Path


SCRATCH_FOLDER = re.compile(rb"/scratch/slurm_jobs/[^/\s]+/\d+/transxpress\.\w+")


def read_allowed(path):
    rules = []
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 3 or not (fields[1] == "skip" or fields[1].startswith("strip:")):
            sys.exit("%s line %d: expected pattern, tab, mode (skip or strip:REGEX), tab, reason" % (path, number))
        rules.append((fields[0], fields[1], fields[2]))
    return rules


def content(path, run):
    """The bytes of a file. A gzip file gives its decompressed bytes. The run folder path becomes @RUN@."""
    data = path.read_bytes()
    if data[:2] == b"\x1f\x8b":
        try:
            data = gzip.decompress(data)
        except (OSError, EOFError):
            pass
    data = data.replace(str(Path(run).resolve()).encode(), b"@RUN@")
    return SCRATCH_FOLDER.sub(b"@SCRATCH@", data)


def md5(data):
    return hashlib.md5(data).hexdigest()


def stripped(data, regex):
    pattern = re.compile(regex.encode())
    return b"\n".join(line for line in data.split(b"\n") if not pattern.search(line))


def first_difference(a, b):
    a_lines = a.split(b"\n")
    b_lines = b.split(b"\n")
    for number, (x, y) in enumerate(zip(a_lines, b_lines), 1):
        if x != y:
            return "line %d: %r != %r" % (number, x[:80], y[:80])
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
            lines = (Path(run) / args.sequences).read_bytes().split(b">")[1:]
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

    for name in sorted(set(files_a) | set(files_b)):
        rule = next((r for r in rules if fnmatch.fnmatchcase(name, r[0])), None)
        if name not in files_a or name not in files_b:
            # a skip rule also lets a file exist in one run only
            if rule is not None and rule[1] == "skip":
                allowed.setdefault(rule[0], []).append(name)
            else:
                failed.append("only in %s: %s" % (args.run_a if name in files_a else args.run_b, name))
            continue
        a = content(files_a[name], args.run_a)
        b = content(files_b[name], args.run_b)
        if a == b:
            same += 1
            continue
        if rule is None:
            failed.append("DIFF %s (%s)" % (name, first_difference(a, b)))
        elif rule[1] == "skip" or stripped(a, rule[1][6:]) == stripped(b, rule[1][6:]):
            allowed.setdefault(rule[0], []).append(name)
        else:
            failed.append("DIFF %s (%s; allowed only after dropping lines of %s)"
                          % (name, first_difference(stripped(a, rule[1][6:]), stripped(b, rule[1][6:])), rule[1][6:]))

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
