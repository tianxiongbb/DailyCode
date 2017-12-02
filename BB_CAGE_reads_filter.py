#!/usr/bin/env python
import sys
import os

if len(sys.argv)<2:
	print sys.argv[0] + " in.fq out.fq"
	os._exit(0)


f_in = open(sys.argv[1], "r")
list_reads = []
k = 0
for l in f_in:
	if k%4 == 0:
		list_reads.append([])
	list_reads[-1].append(l)
	k += 1

print list_reads[0]
f_out = open(sys.argv[2], "w")
for i in range(len(list_reads)):
	if list_reads[i][1][0:2] == "NG" or list_reads[i][1][0:2] == "GG":
		f_out.write(list_reads[i][0])
		f_out.write(list_reads[i][1][2:])
		f_out.write(list_reads[i][2])
		f_out.write(list_reads[i][3][2:])
f_out.close()

		


