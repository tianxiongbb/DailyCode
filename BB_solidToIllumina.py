#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import re
import bb_basic as bb


def main():
    if len(sys.argv) < 2:
        bb.fun_print_help("in.fastq", "out.fastq")
    global RECODE
    RECODE = {"0": {"A": "A", "C": "C", "G": "G", "T": "T"},
            "1": {"A": "C", "C": "A", "G": "T", "T": "G"},
            "2": {"A": "G", "C": "T", "G": "A", "T": "C"},
            "3": {"A": "T", "C": "G", "G": "C", "T": "A"}}
    file_solid = bb.fun_open_file(sys.argv[1])
    file_illumina = bb.fun_open_file(sys.argv[2], "w")
    count = 0
    for l in file_solid:
        count += 1
        if count%100000 == 0:
            sys.stdout.write(".")
            sys.stdout.flush()
        if l.startswith("@") or l.startswith("+"):
            file_illumina.write(l)
        elif l.startswith("!"):
            file_illumina.write(l[1:])
        else:
            string_solid = l.strip()
            string_illumina = decode(string_solid)
            file_illumina.write(string_illumina + "\n")
    file_illumina.close()

# --------functions--------
def decode(string_solid):
    string_illumina = RECODE[string_solid[1]][string_solid[0]]
    for i in range(2, len(string_solid)):
        string_illumina += RECODE[string_solid[i]][string_illumina[-1]]
    return string_illumina

# --------process--------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        bb.fun_print_error("user interrupted, abort!")
        sys.exit(0)
