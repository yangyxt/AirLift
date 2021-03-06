#!/bin/sh

GAP_ALN=$1
GAP_BED=$(echo ${GAP_ALN} | sed s/.bam/.bed/g)
OUT_FOLDER=$2
OUT_BEDS="${OUT_FOLDER}/*.bed"
convert2bed --input=bam < ${GAP_ALN} > ${GAP_BED}
awk -v out_fol="$OUT_FOLDER" '{split($4, array, "_"); fname=array[1] "_" array[2] "_" array[3] "_" array[4] ".bed"; print $1 "\t" $2 "\t" $3 >> out_fol "/" fname;}' ${GAP_BED}
for i in `echo ${OUT_BEDS}`; do fname=`basename $i`; mergeBed -i $i > "${OUT_FOLDER}/merged_$fname"; mv "${OUT_FOLDER}/merged_$fname" $i; done

