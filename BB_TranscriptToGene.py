#! /usr/bin/env python
# -*- coding: utf-8 -*-
#########readme##############


#########module##############

import os,sys,time,re

#########information#########

pro_name=sys.argv[0]
pro_author="Tianxiong Yu"
pro_date="Oct 22, 2015"
pro_purpose=""
pro_begin_time=time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
begin_time=time.clock()

print "\033[1;31;38m"
print "*"*50
print "ProgramName:\t"+pro_name
print "Author:\t"+pro_author
print "ProgramDate:\t"+pro_date
print "Purpose:\t"+pro_purpose
print "Begin time:\t"+pro_begin_time
print "="*50
print "\033[0m"

#########prepare#############

#os.chdir("/Users/duan/Desktop")
	##set your work dir

#########function############
def input_file(input,n):
	file_list=[]
	if len(input)!=n+1:
		print "The input file's number is wrong!!!"
		return
	else:
		for i in range(1,len(input)):
			file_list.append(open("%s"%input[i],"r"))
		return file_list
		
def delete(li,index):
	li=li[:index]+li[index+1:]
	return li

def unique_list(li):
	new_li=list(set(li))
	new_li.sort(li.index)
	return new_li

def read_file_into_list(file_in):
	l=[]
	while True:
		line=file_in.readline().split()
		if not line:
			break
		l.append(line)
	return l

def read_file_into_list_split_tab(file_in):
	l=[]
	while True:
		line=file_in.readline().split("\t")
		if not line:
			break
		l.append(line)
	return l

def read_file_into_dict(file_in,key_num):
	###key_num is start with 0
	d={}
	while True:
		line=file_in.readline().split()
		if not line:
			break
		key=line[key_num]
		value=delete(line,key_num)
		d[key]=value
	return d

def write_list_into_file(li,file_out):
	###the value in li must be string not number
	if len(li)==1:
		content=li[0]
		file_out.write(content)
	if len(li)>1:
		if isinstance(li[0],list):
			for i in range(len(li)):
				li[i]="\t".join(li[i])
			content="\n".join(li)
			file_out.write(content)
		else:
			content="\n".join(li)
			file_out.write(content)
	file_out.close()

#########code################

##help document##
if len(sys.argv)<2:
	print "usage:"
	print "python BB_TranscriptToGene.py in.bed out.bed"
	print ""
	print "\033[1;31;38m"
	print "="*50
	print "\033[0m"
	os._exit(0)

##run##
print "The program is running!!!"

file_pirn_tran=open(sys.argv[1],"r")
file_pirn_gene=open(sys.argv[2],"w")

list_pirn_tran=read_file_into_list(file_pirn_tran)
dict_pirn_gene={}

for i in list_pirn_tran:
	if not dict_pirn_gene.has_key(i[7]):
		dict_pirn_gene[i[7]]=[i[0],i[1],i[2],i[11],i[5]]
	else:
		if int(i[1])<int(dict_pirn_gene[i[7]][1]):
			dict_pirn_gene[i[7]][1]=i[1]
		if int(i[2])>int(dict_pirn_gene[i[7]][2]):
			dict_pirn_gene[i[7]][2]=i[2]

for key in dict_pirn_gene.keys():
	if re.search("^Mir",key):
		continue
	if re.search("rRNA$",key):
		continue
	if re.search("^U\d+",key):
		continue
	if re.search("^SNO",key):
		continue
	if re.search("Y_RNA$",key):
		continue
	if re.search("^sno",key):
		continue
	if int(dict_pirn_gene[key][2])-int(dict_pirn_gene[key][1])<1000:
		continue
	file_pirn_gene.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(dict_pirn_gene[key][0],dict_pirn_gene[key][1],dict_pirn_gene[key][2],key,dict_pirn_gene[key][3],dict_pirn_gene[key][4]))

file_pirn_gene.close()

#########end#################

pro_end_time=time.strftime('%Y-%m-%d %H:%M:%S',\
	time.localtime(time.time()))
end_time=time.clock()
pro_time=end_time-begin_time

print "\033[1;31;38m"
print "="*50
print "End time:\t"+pro_end_time
print "Program run time:\t%.03f seconds"%pro_time
print "*"*50
print "\033[0m"
