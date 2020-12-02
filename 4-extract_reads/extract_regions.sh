#!/bin/sh
SCRIPTPATH=$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )
SCRIPTNAME=$(basename "${SCRIPTPATH}")
GAP_ALN=$1
GAP_BED=$(echo ${GAP_ALN} | sed s/.bam/.bed/g)
OUT_FOLDER=$2
OUT_BEDS=${OUT_FOLDER}/*.bed
convert2bed --input=bam < ${GAP_ALN} > ${GAP_BED}
awk -v out_fol="$OUT_FOLDER" '{split($4, array, "_"); fname=array[1]"_"array[2]"_"array[3]"_"array[4]".bed"; printf $1"\t"$2"\t"$3"\n" >> out_fol"/"fname;}' < ${GAP_BED}
echo "Before ls bed files, take a look at the out_folder "${OUT_FOLDER}
for i in `ls ${OUT_FOLDER}/*.bed`; do
    fname=`basename $i`
    sortBed -i $i | mergeBed -i - > ${OUT_FOLDER}/merged_$fname
    mv ${OUT_FOLDER}/merged_$fname $i
done

