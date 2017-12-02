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
# if len(sys.argv)<2:
# 	print "usage:"
# 	print ""
# 	print "\033[1;31;38m"
# 	print "="*50
# 	print "\033[0m"
# 	os._exit(0)
cls=["103,0,31","178,24,43","214,96,77","244,165,130","253,219,199","209,229,240","146,197,222","67,147,195","33,102,172","5,48,97"]
def colorselect(diver):
	if int(diver)<100:
		c=cls[0]
	elif int(diver)<125:
		c=cls[1]
	elif int(diver)<150:
		c=cls[2]
	elif int(diver)<175:
		c=cls[3]
	elif int(diver)<200:
		c=cls[4]
	elif int(diver)<225:
		c=cls[5]
	elif int(diver)<250:
		c=cls[6]
	elif int(diver)<275:
		c=cls[7]
	elif int(diver)<300:
		c=cls[8]
	else:
		c=cls[9]
	return c

##run##
f_in=open("/data/tongji2/piRNA/Output/RMSK_Annotation/mm10/rmsk.txt","r")
f_out={}
f_out["LINE_-"]=open("/data/tongji2/piRNA/Output/RMSK_Annotation/mm10/rmsk_LINE_anti.bed","w")
f_out["LINE_+"]=open("/data/tongji2/piRNA/Output/RMSK_Annotation/mm10/rmsk_LINE_sense.bed","w")
f_out["SINE_-"]=open("/data/tongji2/piRNA/Output/RMSK_Annotation/mm10/rmsk_SINE_anti.bed","w")
f_out["SINE_+"]=open("/data/tongji2/piRNA/Output/RMSK_Annotation/mm10/rmsk_SINE_sense.bed","w")
f_out["LTR_-"]=open("/data/tongji2/piRNA/Output/RMSK_Annotation/mm10/rmsk_LTR_anti.bed","w")
f_out["LTR_+"]=open("/data/tongji2/piRNA/Output/RMSK_Annotation/mm10/rmsk_LTR_sense.bed","w")
f_out["DNA_-"]=open("/data/tongji2/piRNA/Output/RMSK_Annotation/mm10/rmsk_DNA_anti.bed","w")
f_out["DNA_+"]=open("/data/tongji2/piRNA/Output/RMSK_Annotation/mm10/rmsk_DNA_anti.bed","w")
while True:
	l=f_in.readline().split()
	if not l:
		break
	if l[11] not in ["LINE","SINE","LTR","DNA"]:
		continue
	f_out["_".join([l[11],l[9]])].write("%s\t%s\t%s\t%s\t0\t%s\t%s\t%s\t%s\n"%(l[5],l[6],l[7],l[10],l[9],l[6],l[7],colorselect(l[2])))

f_out["LINE_-"].close()
f_out["LINE_+"].close()
f_out["SINE_-"].close()
f_out["SINE_+"].close()
f_out["LTR_-"].close()
f_out["LTR_+"].close()
f_out["DNA_-"].close()
f_out["DNA_+"].close()


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
