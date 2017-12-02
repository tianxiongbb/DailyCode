#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import re
import bb_basic as bb


def main():
    if len(sys.argv) < 2:
        bb.fun_print_help("species_query", "species_target", "path_out", "chain_files(eg: rn5ToRn6,rn6ToMm10)","[query.ns.bed]",
                "[query.bed]")
    global_para()
    if not os.path.exists(PATH_OUT):
        bb.fun_print("create out path......", "red", "black", 1)
        os.mkdir(PATH_OUT)
    dict_pi = get_pi()
    psl_map()
    dict_exon_intron = read_exon_intron()
    write_exon_intron_length(dict_exon_intron)
    dict_cover = calculate_coverage(dict_exon_intron)
    list_ortho = find_ortho_pi(dict_cover)
    dict_synteny = get_synteny(dict_pi, list_ortho)
    run_lastz(dict_synteny, dict_pi, dict_cover)
    for i in dict_pi:
        if i not in list_ortho and i not in dict_synteny:
            LIST_OUT.append([i, 0, "none", "0", "0", "0", "chrNA", 1, 1, "+"])
    file_out = bb.fun_open_file("{0}/{1}_{2}.conserve.tab".format(PATH_OUT, SPECIES1, SPECIES2), "w")
    for i in LIST_OUT:
        file_out.write(i[0])
        for j in range(1, len(i)):
            file_out.write("\t" + str(i[j]))
        file_out.write("\n")
    file_out.close()

# --------functions--------
def write_exon_intron_length(dict_exon_intron):
    file_out = bb.fun_open_file("{0}/{1}.exon_intron.len".format(PATH_OUT, SPECIES1), "w")
    for pi in dict_exon_intron:
        file_out.write(pi)
        exon_length = 0
        intron_length = 0
        for i in dict_exon_intron[pi]["exon"]:
            exon_length += (i[1] - i[0] + 1)
        for i in dict_exon_intron[pi]["intron"]:
            intron_length += (i[1] - i[0] + 1)
        file_out.write("\t{0}\t{1}\n".format(exon_length, max(0, intron_length)))
    file_out.close()


def calculate_coverage(dict_exon_intron):
    file_psl = bb.fun_open_file("{0}/{1}_{2}.psl".format(PATH_OUT, SPECIES1, SPECIES2))
    dict_cover = {}
    for l in file_psl:
        line = l.strip().split()
        pi = line[9]
        align = []
        map_length = line[18].strip(",").split(",")
        map_start = line[19].strip(",").split(",")
        for i in range(len(map_length)):
            align.append([int(map_start[i]), int(map_start[i]) + int(map_length[i]) - 1])
        exon_length = 0
        intron_length = 0
        for i in dict_exon_intron[pi]["exon"]:
            exon_length += (i[1] - i[0] + 1)
        for i in dict_exon_intron[pi]["intron"]:
            intron_length += (i[1] - i[0] + 1)
        exon_over = 0
        intron_over = 0
        for i in align:
            for j in dict_exon_intron[pi]["exon"]:
                if not (i[1] < j[0] or i[0] > j[1]):
                    exon_over += (min(i[1], j[1]) - max(i[0], j[0]) + 1)
            for j in dict_exon_intron[pi]["intron"]:
                if not (i[1] < j[0] or i[0] > j[1]):
                    intron_over += (min(i[1], j[1]) - max(i[0], j[0]) + 1)
        exon_cover = float(exon_over) / float(exon_length)
        if intron_length == 0:
            intron_cover = -1
        else:
            intron_cover = float(intron_over) / float(intron_length)
        dict_cover[pi] = [str(exon_cover), str(intron_cover)]
    return dict_cover


def read_exon_intron():
    file_pi = bb.fun_open_file(PATH_QUERY_BED)
    file_psl = bb.fun_open_file("{0}/{1}_{2}.psl".format(PATH_OUT, SPECIES1, SPECIES2))
    dict_exon_intron = {}
    dict_min = {}
    for l in file_pi:
        line = l.strip().split()
        line[1] = int(line[1])
        if line[3] not in dict_min:
            dict_min[line[3]] = line[1]
        else:
            if line[1] < dict_min[line[3]]:
                dict_min[line[3]] = line[1]
    file_pi.seek(0)
    for l in file_pi:
        line = l.strip().split()
        line[1:3] = map(int, line[1:3])
        if line[3] not in dict_exon_intron:
            dict_exon_intron[line[3]] = {"exon": [], "intron": []}
        start = dict_min[line[3]]
        dict_exon_intron[line[3]]["exon"].append([line[1] - start, line[2] - start])
    for pi in dict_exon_intron:
        exons = dict_exon_intron[pi]["exon"]
        exons.sort()
        start = exons[0][1] + 1
        for i in range(1, len(exons)):
            end = exons[i][0] - 1
            dict_exon_intron[pi]["intron"].append([start, end])
            start = exons[i][1] + 1
    return dict_exon_intron


def run_lastz(dict_synteny, dict_pi, dict_cover):
    for pi in dict_synteny:
        bb.fun_print("start check conservation for " + pi, "green", "black", 1)
        bb.fun_quick_write("%s\t%s\t%s\t%s\t0\t%s\n"%(dict_synteny[pi][0][0], dict_synteny[pi][0][1], dict_synteny[pi][0][2], pi, dict_synteny[pi][0][3]), "%s/temp1.bed"%PATH_OUT)
        bb.fun_quick_write("%s\t%s\t%s\t%s\t0\t%s\n"%(dict_pi[pi][0], max(1, dict_pi[pi][1] - 10000), min(DICT_CHROM1[dict_pi[pi][0]], dict_pi[pi][2] + 10000), pi, dict_pi[pi][4]), "%s/temp.shuffle1.bed"%PATH_OUT)
        bb.fun_quick_write("%s\t%s\t%s\t%s\t0\t%s\n"%(dict_synteny[pi][1][0], dict_synteny[pi][1][1], dict_synteny[pi][1][2], pi, dict_synteny[pi][1][3]), "%s/temp2.bed"%PATH_OUT)
        command = "bedtools getfasta -fi /data/tongji2/Annotation/Fasta/{0}.fa \
                -fo {1}/temp1.fa -bed {1}/temp1.bed -name -s".format(DICT_SP[SPECIES1], PATH_OUT)
        if os.system(command) != 0:
            bb.fun_print_error("getfasta Error")
        command = "bedtools getfasta -fi /data/tongji2/Annotation/Fasta/{0}.fa \
                -fo {1}/temp.shuffle1.fa -bed {1}/temp.shuffle1.bed -name -s".format(DICT_SP[SPECIES1], PATH_OUT)
        if os.system(command) != 0:
            bb.fun_print_error("getfasta Error")
        command = "bedtools getfasta -fi /data/tongji2/Annotation/Fasta/{0}.fa \
                -fo {1}/temp2.fa -bed {1}/temp2.bed -name -s".format(DICT_SP[SPECIES2], PATH_OUT)
        if os.system(command) != 0:
            bb.fun_print_error("getfasta Error")
        # shuffle
        bb.fun_print("start shuffling......", "blue", "black", font=1)
        shuffle_scores = [0]
        for i in range(0):
            command = "bedtools shuffle -i {0}/temp.shuffle1.bed -g /data/tongji2/Annotation/ChromSize/{1}.chrom.size > {0}/temp.shuffle2.bed".format(PATH_OUT, DICT_SP[SPECIES2])
            if os.system(command) != 0:
                bb.fun_print_error("bedtools shuffle Error")
            command = "bedtools getfasta -fi /data/tongji2/Annotation/Fasta/{0}.fa \
                    -fo {1}/temp.shuffle2.fa -bed {1}/temp.shuffle2.bed -name -s".format(DICT_SP[SPECIES2], PATH_OUT)
            if os.system(command) != 0:
                bb.fun_print_error("getfasta Error")
            command = "lastz {0}/temp.shuffle1.fa {0}/temp.shuffle2.fa --strand=plus --chain \
                    --output={0}/temp.shuffle.lastz --format=general:score".format(PATH_OUT)
            os.system(command)
            file_shuffle = bb.fun_open_file("{0}/temp.shuffle.lastz".format(PATH_OUT))
            lines = file_shuffle.readlines()
            if len(lines) == 1:
                shuffle_scores.append(0)
            else:
                for record in lines:
                    if record[0] != "#":
                        shuffle_scores.append(float(record.strip()))
        shuffle_scores.sort()
        cutoff = shuffle_scores[int(len(shuffle_scores)*0.95)]  #p-value=0.05
        bb.fun_print("p 0.05 cutoff mapping score: " + str(cutoff), "blue", "black", font=1)
        # end shuffle
        command = "lastz {0}/temp1.fa {0}/temp2.fa --strand=plus --chain \
                --rdotplot={0}/temp.rplot --output={0}/temp.lastz \
                --format=general:name1,start1,end1,name2,start2,end2,nmatch,identity,score".format(PATH_OUT)
        if os.path.getsize("{0}/temp1.fa".format(PATH_OUT)) and os.path.getsize("{0}/temp1.fa".format(PATH_OUT)):
            if os.system(command) != 0:
                bb.fun_print_error("lastz Error")
        else:
            LIST_OUT.append([pi, 1, "none", "0", "0", "0", "chrNA", 1, 1, "+"])
            continue
        os.system("sort -k2,2n {0}/temp.lastz > {0}/temp.sort.lastz".format(PATH_OUT))
        file_in = bb.fun_open_file("%s/temp.lastz"%PATH_OUT)
        list_lastz = []
        start = dict_synteny[pi][0][1]
        end = dict_synteny[pi][0][2]
        k = 0
        for l in file_in:
            if l[0] == "#":
                continue
            line = l.strip().split()
            if float(line[-1]) > cutoff:
                list_lastz.append([int(line[1]), int(line[2])])
        if dict_synteny[pi][0][3] == "+":
            pi_start = dict_pi[pi][1] - start
            pi_end = dict_pi[pi][2] - start
        else:
            pi_start = end - dict_pi[pi][2]
            pi_end = end - dict_pi[pi][1]
        for i in list_lastz:
            if pi_end < i[0] or pi_start > i[1]:
                continue
            else:
                LIST_OUT.append([pi, 3, "none",
                    dict_synteny[pi][1][4], dict_cover[pi][0], dict_cover[pi][1],
                    dict_synteny[pi][1][0], dict_synteny[pi][1][1], dict_synteny[pi][1][2],
                    dict_synteny[pi][1][3]])
                k = 1
                break
        if k==0:
            LIST_OUT.append([pi, 2, "none", "0", "0", "0", "chrNA", 1, 1, "+"])


def get_synteny(dict_pi, list_ortho):
    out = {}
    file_in = bb.fun_open_file("{0}/{1}_{2}.bed".format(PATH_OUT, SPECIES1, SPECIES2))
    for l in file_in:
        line = l.strip().split()
        if line[3] not in list_ortho:
            out[line[3]] = []
            out[line[3]].append([dict_pi[line[3]][0],
                max(1, dict_pi[line[3]][1] - 150000),
                min(DICT_CHROM1[dict_pi[line[3]][0]], dict_pi[line[3]][2] + 150000),
                dict_pi[line[3]][4]]) 
            out[line[3]].append([line[0],
                max(1, int(line[1]) - 150000),
                min(DICT_CHROM2[line[0]], int(line[2]) + 150000),
                line[5], line[4]])
    return out


def get_pi():
    out = {}
    file_in = bb.fun_open_file(PATH_QUERY_NS_BED)
    for l in file_in:
        line = l.strip().split()
        out[line[3]] = line[:3] + line[4:]
        out[line[3]][1:3] = map(int, out[line[3]][1:3])
    return out


def find_ortho_pi(dict_cover):
    os.system("sort -k1,1 -k2,2n {0}/{1}_{2}.bed > {0}/t && mv {0}/t {0}/{1}_{2}.bed && sort -k1,1 -k2,2n {3}/{2}/{2}.piRNA.ns.bed > {0}/t && mv {0}/t {3}/{2}/{2}.piRNA.ns.bed".format(PATH_OUT, SPECIES1, SPECIES2, PATH_PI))
    command = "bedtools intersect -sorted -wo -s -a {0}/{1}_{2}.bed -b {3}/{2}/{2}.piRNA.ns.bed > {0}/{1}_{2}.ortho".format(PATH_OUT, SPECIES1, SPECIES2, PATH_PI)
    if os.system(command) != 0:
        bb.fun_print_error("intersectBed Error")
    out = []
    file_in = bb.fun_open_file("{0}/{1}_{2}.ortho".format(PATH_OUT, SPECIES1, SPECIES2))
    for l in file_in:
        line = l.strip().split()
        LIST_OUT.append([line[3], 4, line[9],
            line[4], dict_cover[line[3]][0], dict_cover[line[3]][1],
            line[6], line[7], line[8], line[11]])
        out.append(line[3])
    return out

    
def psl_map():
    TEMP_PSL_IN = PATH_QUERY_NS_BED
    for chains in PATH_CHAIN:
        TEMP_CHROM_SIZE = chains.split("To")[0] + ".chrom.size"
        command = "BB_PslMap {0} {1} {2} {3}/temp_pslmap_out 0.1".format(TEMP_PSL_IN, TEMP_CHROM_SIZE, chains, PATH_OUT)
        if os.system(command) != 0:
            bb.fun_print(PATH_CHAIN)
            bb.fun_print_error("BB_PslMap Error, please check chain files")
        os.system("cp {0}/temp_pslmap_out.psl {0}/temp_pslmap_in.psl".format(PATH_OUT))
        TEMP_PSL_IN = "{0}/temp_pslmap_in.psl".format(PATH_OUT)
    os.system("cp {0}/temp_pslmap_out.psl {0}/{1}_{2}.psl && cp {0}/temp_pslmap_out.bed {0}/{1}_{2}.bed".format(PATH_OUT, SPECIES1, SPECIES2))


def global_para():
    global SPECIES1, SPECIES2, PATH_PI, PATH_OUT, DICT_SP, DICT_CHROM1, DICT_CHROM2, LIST_OUT
    SPECIES1, SPECIES2 = sys.argv[1], sys.argv[2]
    PATH_PI = "/data/tongji2/piRNA/OutputConservation/piFinder"
    PATH_OUT = sys.argv[3].rstrip("/")
    DICT_SP = {"mouse":"mm10", "rat":"rn6", "rabbit":"oryCun2", "human":"hg38",
            "rhesus":"rheMac8", "marmoset":"calJac3", "cow":"bosTau8",
            "pig":"susScr3", "platypus":"ornAna1", "opossum":"monDom5",
            "chicken":"galGal4"}
    DICT_CHROM1 = read_chrom("/data/tongji2/Annotation/ChromSize/%s.chrom.size"%DICT_SP[SPECIES1])
    DICT_CHROM2 = read_chrom("/data/tongji2/Annotation/ChromSize/%s.chrom.size"%DICT_SP[SPECIES2])
    LIST_OUT = []
    global PATH_QUERY_BED, PATH_QUERY_NS_BED
    PATH_QUERY_BED = "0"
    PATH_QUERY_NS_BED = "0"
    if len(sys.argv) > 5:
        PATH_QUERY_NS_BED = sys.argv[5]
        PATH_QUERY_BED = sys.argv[6]
    else:
        PATH_QUERY_NS_BED = "{0}/{1}/{1}.piRNA.ns.bed".format(PATH_PI, SPECIES1)
        PATH_QUERY_BED = "{0}/{1}/{1}.piRNA.bed".format(PATH_PI, SPECIES1)
    global PATH_CHAIN
    PATH_CHAIN = sys.argv[4].split(",")

    
def read_chrom(path):
    out = {}
    file_in = bb.fun_open_file(path)
    for l in file_in:
        line = l.strip().split()
        out[line[0]] = int(line[1])
    return out

# --------process--------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        bb.fun_print_error("user interrupted, abort!")
        sys.exit(0)
