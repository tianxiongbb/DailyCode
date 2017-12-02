#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess
import sys
import re
import os
import argparse
import bb_basic as bb


def main():
    args = get_args()
    args.align = prepare_align_file(args)
    list_gene = read_gene_list(args)
    dict_matrix = {}
    for gn in list_gene:
        dict_matrix[gn] = []
    cmd = ["bedtools", "coverage", "-d", "-a", args.gene_position, "-b", args.align]
    bb.fun_print("run bedtools coverage", "green", "black", 1)
    ret = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    factor = 0
    for l in bb.fun_open_file(args.align):
        factor += 1
    print "unique mapped reads number: %s"%factor
    c = 0
    r = 0
    name = "0"
    bb.fun_print("make matrix", "green", "black", 1)
    tt = 0
    for l in ret.stdout:
        tt += 1
        if tt % 1000000 == 0:
            bb.fun_print("%s million lines processed......"%(tt/1000000))
        line = l.strip().split()
        if line[3] != name and name != "0" and r != 0:
            dict_matrix[name].append(str(float(c)/args.resolution*1000000/factor))
            r = 0
            c = 0
        name = line[3]
        r += 1
        c += int(line[7])
        if r == args.resolution:
            dict_matrix[name].append(str(float(c)/args.resolution*1000000/factor))
            r = 0
            c = 0
    file_matrix = bb.fun_open_file(args.output, "w")
    for gn in list_gene:
        file_matrix.write(gn)
        if DICT_STRAND[gn] == "+":
            for dp in dict_matrix[gn]:
                file_matrix.write("\t" + dp)
            file_matrix.write("\n")
        else:
            for dp in dict_matrix[gn][::-1]:
                file_matrix.write("\t" + dp)
            file_matrix.write("\n")
    file_matrix.close()
    if args.format in ["sam", "bam"]:
        subprocess.check_call("rm %s_temp*"%args.align)

# --------functions--------
def read_gene_list(args):
    list_gene = []
    global DICT_STRAND
    DICT_STRAND = {}
    for l in bb.fun_open_file(args.gene_position):
        if l.split()[3] in list_gene:
            bb.fun_print_error("there is repetitive names in gene position file, please remove the repeat name")
        list_gene.append(l.split()[3])
        DICT_STRAND[l.split()[3]] = l.split()[5]
    return list_gene


def prepare_align_file(args):
    if args.format == "bed":
        align_file = args.align
    elif args.format == "bam":
        if args.sorted:
            cmd = ["bedtools", "bamtobed", "-i", args.align, ">", "%s_temp.bed"%args.align]
            subprocess.check_call(cmd)
            align_file = "%s_temp.bed"%args.align
        else:
            cmd = ["samtools", "sort", args.align, "%s_temp.bam"%args.align]
            subprocess.check_call(cmd)
            cmd = ["bedtools", "bamtobed", "-i","%s_temp.bam"%args.align, ">", "%s_temp.bed"%args.align]
            subprocess.check_call(cmd)
            align_file = "%s_temp.bed"%args.align
    elif args.format == "sam":
        cmd = ["samtools", "view", "-bhS", args.align, ">", "%s_temp.unsort.bam"%args.align]
        subprocess.check_call(cmd)
        cmd = ["samtools", "sort", "%s_temp.unsort.bam"%args.align, "%s_temp.bam"%args.align]
        subprocess.check_call(cmd)
        cmd = ["bedtools", "bamtobed", "-i","%s_temp.bam"%args.align, ">", "%s_temp.bed"%args.align]
        subprocess.check_call(cmd)
        align_file = "%s_temp.bed"%args.align
    return align_file
        
            
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    parser.add_argument("-f", "--format", help="chip-seq align file format, sam, bam or bed", default="bed")
    parser.add_argument("-s", "--sorted", help="is align file sorted", action='store_false', default=True)
    parser.add_argument("-a", "--align", help="chip-seq align file")
    parser.add_argument("-g", "--gene_position", help="bed file for each gene region to calculate")
    parser.add_argument("-o", "--output", help="output depth matrix")
    parser.add_argument("-r", "--resolution", help="bin (bp) used for depth calculation", default=10, type=int)
    args = parser.parse_args()
    check_args(args)
    return args

def check_args(args):
    if not os.path.exists(args.align):
        bb.fun_print_error("wrong input pathway: %s, please re-define the right pathway for align file"%args.align)
    if not os.path.exists(args.gene_position):
        bb.fun_print_error("wrong input pathway: %s, please re-define the right pathway for gene position file"%args.gene_position)
    if args.format not in ["sam", "bam", "bed"]:
        bb.fun_print_error("wrong format, the program only support sam, bam or bed")


def print_help():
    if len(sys.argv) < 2:
        bb.fun_print_help(args)

# --------process--------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        bb.fun_print_error("user interrupted, abort!")
        sys.exit(0)
