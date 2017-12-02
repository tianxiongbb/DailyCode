#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import re
import bb_basic as bb


def main():
    if len(sys.argv) < 2:
        bb.fun_print_help("in.bed", "in.bed2", "out.rpm")
    judge_format(sys.argv[1])
    intersectBed(sys.argv[1], sys.argv[2], sys.argv[3])
    dict_bed = calculate_sig(sys.argv[3])
    calculate_rpm(dict_bed, sys.argv[2], sys.argv[3])
    os.system("rm {0}.temp*".format(sys.argv[3]))

# --------functions--------
def judge_format(path):
    file_bed = bb.fun_open_file(path)
    while True:
        line = file_bed.readline().split()
        if line[0][0] != "#":
            global N
            N = len(line)
            bb.fun_print("input bed format is bed" + str(N), "green", "black", 1)
            break


def calculate_rpm(dict_bed, path_bed2, out):
    factor = 0
    file_bed2 = bb.fun_open_file(path_bed2)
    for l in file_bed2:
        line = l.strip().split()
        if line[4] == "1":
            factor += float(line[3])/float(line[4])
    factor = factor/1000000
    file_rpm = bb.fun_open_file(out, "w")
    for pi in dict_bed:
        dict_bed[pi][N-1] = float(dict_bed[pi][N-1]) / factor
        for i in range(3):
            file_rpm.write("%s\t"%dict_bed[pi][i])
        file_rpm.write(pi)
        for i in range(3, N):
            file_rpm.write("\t%s"%dict_bed[pi][i])
        file_rpm.write("\n")
    file_rpm.close()

        
def calculate_sig(out):
    file_bed = bb.fun_open_file("{0}.temp.bed".format(out))
    dict_bed = {}
    for l in file_bed:
        line = l.strip().split()
        if line[3] not in dict_bed:
            dict_bed[line[3]] = line[0:3] + line[4:N] + [0]
            dict_bed[line[3]][1:3] = map(int, dict_bed[line[3]][1:3])
        else:
            dict_bed[line[3]][1] = min(dict_bed[line[3]][1], int(line[1]))
            dict_bed[line[3]][2] = max(dict_bed[line[3]][2], int(line[2]))
    file_inter = bb.fun_open_file("{0}.temp.inter".format(out))
    for l in file_inter:
        line = l.strip().split()
        signal = float(line[N+3]) / float(line[N+4])
        dict_bed[line[3]][N-1] += signal
    return dict_bed


def intersectBed(path1, path2, out):
    command = "sort -k1,1 -k2,2n {0} > {1}.temp.bed".format(path1, out)
    os.system(command)
    #already sorted
    #command = "sort -k1,1 -k2,2n {0} > {1}.temp.bed2".format(path2, out)
    #os.system(command)
    command = "bedtools intersect -nonamecheck -sorted -s -wo -F 0.5 -a {0}.temp.bed -b {1} > {0}.temp.inter".format(out, path2)
    if os.system(command) != 0:
        bb.fun_print_error("intersectBed Error")


# --------process--------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        bb.fun_print_error("user interrupted, abort!")
        sys.exit(0)
