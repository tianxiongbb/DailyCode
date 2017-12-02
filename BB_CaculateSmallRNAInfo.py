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
		if line==[""]:
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
	print "BB_CaculateSmallRNASInfo.py in.fq out.prefix"
	print "\033[1;31;38m"
	print "="*50
	print "\033[0m"
	os._exit(0)

##run##
f_in=open(sys.argv[1],"r")
f_lendist=open("%s_lendist.txt"%sys.argv[2],"w")
f_firstnucldist=open("%s_firstnucldist.txt"%sys.argv[2],"w")
f_tenthnucldist=open("%s_tenthnucldist.txt"%sys.argv[2],"w")

d_len={}
for i in range(10,50):
	d_len[str(i)]=0

d_first_nucl={"A":0,"T":0,"C":0,"G":0,"N":0}
d_tenth_nucl={"A":0,"T":0,"C":0,"G":0,"N":0}

a=0
while True:
	line=f_in.readline().split()
	if not line:
		break
	if line[0][0]=="@":
		a=1
		continue
	if a==1:
		l=len(line[0])
		d_len[str(l)]+=1
		d_first_nucl[line[0][0]]+=1
		d_tenth_nucl[line[0][9]]+=1
		a=0

for i in range(10,50):
	f_lendist.write("%s\t%s\n"%(i,d_len[str(i)]))

for i in ["A","T","C","G","N"]:
	f_firstnucldist.write("%s\t%s\n"%(i,d_first_nucl[i]))
	f_tenthnucldist.write("%s\t%s\n"%(i,d_tenth_nucl[i]))

f_lendist.close()
f_firstnucldist.close()
f_tenthnucldist.close()

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
