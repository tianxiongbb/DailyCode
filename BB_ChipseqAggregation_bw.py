#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess
import sys
import re
import os
import argparse
import bb_basic as bb


def main():
    print_help()
    matrix = {}
    strand = {}
    for l in bb.fun_open_file(sys.argv[2]):
	line = l.strip().split()
        strand[line[3]] = line[5]
        cmd = "bigWigSummary {0} {1} {2} {3} {4}".format(sys.argv[1], line[0], line[1], line[2], sys.argv[4])
        res = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        sig = res.stdout.readline().strip().split()
	if sig == []:
		matrix[line[3]] = []
		for j in range(int(sys.argv[4])):
			matrix[line[3]].append("0")
	else:
        	matrix[line[3]] = sig
    file_out = bb.fun_open_file(sys.argv[3], "w")
    for gn in matrix:
        file_out.write(gn)
        list_out = []
        for sig in matrix[gn]:
            if sig == "n/a":
                list_out.append("0")
            else:
                list_out.append(sig)
        if strand[gn] == "-":
            list_out = list_out[::-1]
        file_out.write("\t" + "\t".join(list_out) + "\n")
    file_out.close()
    #args = get_args()

# --------functions--------
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    parser.add_argument("args", help="", nargs="*")
    args = parser.parse_args()
    return args


def print_help():
    if len(sys.argv) < 2:
        bb.fun_print_help("in.bw", "gene.bed", "out.mat", "dataPoints")

# --------process--------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        bb.fun_print_error("user interrupted, abort!")
        sys.exit(0)
