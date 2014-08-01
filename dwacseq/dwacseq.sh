#!/bin/bash

#	dwacseq.sh 

if [ $# -ne 11 ]
then
	echo "unexpected number of arguments ($#), expected 11. Exiting"
	exit
fi


testbam=$1
refbam=$2
chromosome=$3
chrprefix=$4
format=$5
amb=$6
scale=$7
finetuning=$8
gwn=$9
removeclonal=${10}
scriptdir=${11}

echo "command: $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}"
echo "gwn: $gwn"


# index bam files
samtools index $testbam
samtools index $refbam


# deal with "chr" prefix Y/N
if [ $chrprefix == "N" ]
then
	chromosome=${chromosome:3}
fi 


# deal with "None" arguments
if [ $gwn == "None" ]
then
	gwn=""
fi

if [ $finetuning == "None" ]
then
	finetuning=""
fi

if [ $removeclonal == "None" ]
then
	removeclonal=""
fi

# call dwacseq
${scriptdir}/dwac-seq.pl --test $testbam --ref $refbam --chr $chromosome --output outfile_ $gwn $finetuning --amb $amb --scale $scale --output_format $format $removeclonal     2>&1


# move output to expected names
mv outfile_*.win windows_output
mv outfile_*.avg avg4hilbert_output
mv outfile_*.seg segmented_output
mv outfile_*.txt tuned_output
