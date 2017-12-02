#1/bin/bash

if [ $# -lt 2 ];then
	echo0 1 "BB_trackUrl.sh in.bed out.html"
	exit 0
fi


awk 'BEGIN{FS=OFS="\t";print "<html>\n<table>\n"} {print "<tr>\n<td>"$4"</td>\n<td><a href=\"http://genome.cse.ucsc.edu/cgi-bin/hgTracks? db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtM ode=0&nonVirtPosition=&position="$1"%3A"($2-int(0.2*($3-$2)))"-"($3+int(0.2*($3- $2)))"&hgsid=568800783_Rsol4nbYzx6CYpMyjYyr1A3OxbrS\">Genome Browser</a></td>\n</tr>"} END{print "</table>\n</html>"}' $1 > $2


