#!/bin/bash

echo -e "\n"
if [ $# -eq 0 ]; then
  echo "No argument given, using default size threshold 1000 bps to find deletions"
  #exit 1
elif [ $# -eq 1 ]; then
	echo -e "Finding deletions of size larger than $1 bps\n"
elif [ $# -eq 2 ]; then
	if [ $2 == "dup" ]; then
		echo -e "Finding duplications of size larger than $1 bps\n"
	elif [ $2 == "del" ]; then
		echo -e "Finding deletions of size larger than $1 bps\n"
	fi
else
   echo -e "Something is wrong with your input. It should be like:\nFor deletions: ./svcalls.sh 1000\nFor duplications: ./svcalls.sh 1000 dup\n"
   exit 1
fi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

gunzip ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

if [ $# -eq 1 ] || ([ $# -eq 2 ] && [ $2 == "del" ]); then
	cat ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf |grep "SVTYPE=DEL"|grep PASS|grep -v "DEL_ALU\|DEL_HERV\|DEL_LINE1\|DEL_SVA"|cut -f 1,2,8|sed "s/;END=/\t/"|cut -f 1,2,4|sed "s/;/\t/"|cut -f 1,2,3| awk '{if (($3-$2)>'$1') print}'|grep -v "X\|Y" >1K_dels.bed
fi

if [ $# -eq 2 ] && [ $2 == "dup" ]; then
	cat ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf |grep "SVTYPE=DUP"|grep PASS|cut -f 1,2,8|sed "s/;END=/\t/"|cut -f 1,2,4|sed "s/;/\t/"|cut -f 1,2,3| awk '{if (($3-$2)>'$1') print}'|grep -v "X\|Y" >1K_dups.bed
fi

rm ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf

echo -e "\n"
