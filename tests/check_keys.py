#!/usr/bin/env python3
"""Check that config.yaml, workflow/schemas/config.schema.yaml and the Configuration table of README.md name the same keys.

Usage: check_keys.py   (run from the clone)
README.md rows: `a`, `b` give two keys; `_x` after `annotation_chunk_fasta` is annotation_chunk_x; a row that
starts with `<rule>_memory` gives <rule>_memory, <rule>_threads and <rule>_runtime for each rule named in its
last cell. Exit code 0: the three lists are equal. 1: a key is missing somewhere.
"""
import re
import sys
from pathlib import Path

import yaml


def readme_keys(path):
    keys = set()
    in_table = False
    for line in Path(path).read_text().splitlines():
        if line.startswith("## "):
            in_table = line == "## Configuration"
        if not in_table or not line.startswith("| `"):
            continue
        cells = re.split(r"(?<!\\) \| ", line.strip("| "))
        tokens = re.findall(r"`([^`]+)`", cells[0])
        if tokens[0].startswith("<rule>"):
            rules = [token for token in re.findall(r"`([^`]+)`", cells[-1]) if re.fullmatch(r"[a-z_]+", token)]
            keys.update(rule + token[len("<rule>"):] if token.startswith("<rule>") else rule + token
                        for rule in rules for token in tokens)
            continue
        previous = None
        for token in tokens:
            if token.startswith("_"):
                token = previous.rsplit("_", 1)[0] + token
            keys.add(token)
            previous = token
    return keys


def main():
    config = set(yaml.safe_load(Path("config.yaml").read_text()))
    schema = set(yaml.safe_load(Path("workflow/schemas/config.schema.yaml").read_text())["properties"])
    readme = readme_keys("README.md")
    status = 0
    for name, keys in [("workflow/schemas/config.schema.yaml", schema), ("README.md", readme)]:
        for key in sorted(config - keys):
            print("%s lacks %s (in config.yaml)" % (name, key))
            status = 1
        for key in sorted(keys - config):
            print("config.yaml lacks %s (in %s)" % (key, name))
            status = 1
    print("OK: %d keys" % len(config) if status == 0 else "FAILED")
    return status


if __name__ == "__main__":
    sys.exit(main())
