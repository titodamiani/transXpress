#!/usr/bin/env python3
"""Check the strand flag that the kallisto rule passes to align_and_estimate_abundance.pl.

Usage: check_kallisto_strand.py   (run from the clone)
Takes abundance_strand_option() and kallisto_strand_option() from workflow/rules/common.smk, runs them for each
assembler and strand_specific value, and checks that the kallisto shell block passes the result. Exit code 0: all
as expected. 1: a value differs.
"""
import re
import sys
from pathlib import Path

NAMES = {"abundance_strand_option", "kallisto_strand_option"}
RF = '--kallisto_add_opts "--rf-stranded"'
FR = '--kallisto_add_opts "--fr-stranded"'
CASES = [
    ("trinity", "", ""),
    ("trinity", "--SS_lib_type RF", RF),
    ("trinity", "--SS_lib_type FR", FR),
    ("trinity", "--SS_lib_type R", RF),
    ("trinity", "--SS_lib_type F", FR),
    ("rnaspades", "", ""),
    ("rnaspades", "--ss rf", RF),
    ("rnaspades", "--ss fr", FR),
]


def main():
    # common.smk is not valid Python (it has Snakemake rules): cut out the two functions by text
    text = Path("workflow/rules/common.smk").read_text()
    code = ""
    for name in NAMES:
        match = re.search(r"^def %s\(.*?(?=^\S)" % name, text, re.S | re.M)
        assert match, name
        code += match.group(0)
    config = {}
    scope = {"config": config}
    exec(code, scope)
    status = 0
    for assembler, strand, expected in CASES:
        config.update(assembler=assembler, strand_specific=strand)
        found = scope["kallisto_strand_option"]()
        if found != expected:
            print("%s, '%s': expected [%s], got [%s]" % (assembler, strand, expected, found))
            status = 1
    config.update(assembler="rnaspades", strand_specific="--ss xx")
    try:
        scope["kallisto_strand_option"]()
        print("rnaSPAdes with '--ss xx' did not stop")
        status = 1
    except ValueError:
        pass
    rule = Path("workflow/rules/expression.smk").read_text().split("rule kallisto:")[1].split("\nrule ")[0]
    if "{params.kallisto_strand}" not in rule or "kallisto_strand=kallisto_strand_option()" not in rule:
        print("rule kallisto does not pass kallisto_strand")
        status = 1
    print("OK: %d cases" % (len(CASES) + 1) if status == 0 else "FAILED")
    return status


if __name__ == "__main__":
    sys.exit(main())
