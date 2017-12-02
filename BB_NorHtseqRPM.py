#!/usr/bin/env python
# -*- coding: utf-8 -*-

##----INTRO-----------##
# Name = 'Mar08.NormaReadCount.py'
# Date = 'Mar 08, 2016'
# Update = ', 2016'
#-Kaili Fan-#

# Purpose = 'This program is for Normalize Read Counts out from Htseq.'

##----START-----------##
import time
########################
#print "\033[1;31;38m"
#print " @o@ "*20
#print ' <BEGIN@'+ time.strftime("%Y-%m-%d %X", time.localtime())+ '>'
#print " @o@ "*20
#print "\033[0m"
########################

##----PREPARE---------##
import os, sys, re
#os.chdir('')

###--help--###
if len(sys.argv)<2:
	print "\033[1;31;38m"
	print "usage:"
	print "python Mar08.NormaReadCount.py <option>* [inputFile] [outputFile]"
	print ""
	print " inputFile --The input_file to be processed."
	print " outputFile --The output_file."
	print ""
	print "This file is for ."
	print "\033[0m"
	os._exit(0)

###--Function--###
#--common--
def FuN_MakeDir(Path_):
	Path_ = Path_.strip(); Path_ = Path_.rstrip("\\")
	if not os.path.exists(Path_):
		os.makedirs(Path_)

def FuN_Write(File_, File_D_):
	Out_ = open(File_D_, 'w')
	Out_.writelines(File_)
	Out_.close()

# search Str_ in List_, return location
def ListSearch(List_,Str_):
	Result_= -1
	for i in range(len(List_)):
		if List_[i].find(Str_)!= -1:
			Result_=i
			break
	return Result_

# Print Dic(key+value)
def PrintDic(Dic_):
	Out_=[]
	for key,value in Dic_.items():
		Out_.append(key+'\t'+value+'\n')
	return Out_

# Print Dic(value)
def PrintDic2(Dic_):
	Out_=[]
	for key,value in Dic_.items():
		Out_.append(value)
	return Out_

#--special-


##----MAIN------------##
#-Parameter--#
inputFile=sys.argv[1]
outputFile=sys.argv[2]
#-Data Input--#
RawData=open(inputFile).readlines()
#-Process--#
# get total number of reads
N=0
for line in RawData:
	if line.rstrip().split("\t")[0][0:7]=="Warning":
		continue
	num=int(line.rstrip().split("\t")[1])
	if line.startswith("_"):
		if line.startswith("__no_feature"):
			N+=num
	else:
		N+=num
# normalize
Out=[]
for line in RawData:
	if line.rstrip().split("\t")[0][0:7]=="Warning":
		continue
	name=line.rstrip().split("\t")[0]
	num=int(line.rstrip().split("\t")[1])
	if not line.startswith("_"):
		norNum=round(num*1000000.0/N,2)
		outLine=name+'\t'+str(norNum)+'\n'
		Out.append(outLine)
#Output
FuN_Write(Out,outputFile)



		
##----CODA------------##
########################
#print "\033[1;31;38m"
#print " ^m^ "*20
#print ' <END  @'+ time.strftime("%Y-%m-%d %X", time.localtime())+ '>'
#print " ^m^ "*20
#print "\033[0m"
########################
##----TEST------------##
#FileName = os.path.splitext(os.path.basename(Dir))[1]
