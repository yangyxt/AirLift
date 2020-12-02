#!/bin/bash
SCRIPTPATH=$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )
parent_folder_path=$(echo ${SCRIPTPATH} | awk -F '/' '{for(i=2;i<NF;i++) printf "/"$i;}')
. ${parent_folder_path}/common_bash_utils.sh

SRCFOLDER=$1
#ref folders containing ref file:
OLDREF=$2
NEWREF=$3
SEQ_FILE_EXT=$4
READSIZE=$5
#make sure the bam files are indexed (i.e., samtools index)
READ_BAM=$6
FIRST_PAIR=$7
SECOND_PAIR=$8
OUTPUT=$9
THREAD=${10}


function generate_chain(){
    bash "${SRCFOLDER}/1-generate_chain/chain_install.sh" "${OUTPUT}"
    bash "${SRCFOLDER}/1-generate_chain/sample_run.sh" \
    "$(get_current_dir ${OLDREF})" \
    "$(get_current_dir ${NEWREF})" \
    "${OUTPUT}/" \
    "${SRCFOLDER}/1-generate_chain/chain_generate.sh" \
    ${SEQ_FILE_EXT} \
    "${OUTPUT}/"
}

function generate_gaps(){
    for i in `ls ${OUTPUT}/*.chain`; do 
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/${chain_basename}_extracted_gaps.fa.time" \
        python3 "${SRCFOLDER}/2-generate_gaps/extract_gaps.py" \
        "$i" \
        "${NEWREF}" \
        ${READSIZE} \
        "${OUTPUT}/${chain_basename}_extracted_gaps.fa"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/${chain_basename}_kmer_gaps.fasta.time" \
        python3 "${SRCFOLDER}/2-generate_gaps/gaps_to_fasta.py" \
        "${OUTPUT}/${chain_basename}_extracted_gaps.fa" \
        ${READSIZE} \
        "${OUTPUT}/${chain_basename}_kmer_gaps.fasta" \
        10
    done
}

function align_gaps(){
    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        bash "${SRCFOLDER}/3-align_gaps/align_gaps.sh" \
        "${OLDREF}" \
        "${OUTPUT}/${chain_basename}_kmer_gaps.fasta" \
        ${THREAD} \
        "${OUTPUT}/${chain_basename}_aligned_gaps"
    done
}

function extract_reads(){
    mkdir "${OUTPUT}/bedfiles"
    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        mkdir -p "${OUTPUT}/bedfiles/${chain_basename}/"; \
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/extract_regions.time" \
        bash "${SRCFOLDER}/4-extract_reads/extract_regions.sh" \
        "${OUTPUT}/${chain_basename}_aligned_gaps.bam" \
        "${OUTPUT}/bedfiles/${chain_basename}"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        local num_of_files=`ls -l ${OUTPUT}/bedfiles/${chain_basename}/*.bed | wc -l`
        local sed_bed_files=$(echo ${OUTPUT}\/bedfiles\/${chain_basename}/ | sed 's/\//\\\//g')
        local sed_read_scripts=$(echo ${SRCFOLDER}\/4-extract_reads\/extract_reads.sh | sed 's/\//\\\//g')
        local sed_read_bam=$(echo $READ_BAM | sed 's/\//\\\//g')
        cat "${SRCFOLDER}/4-extract_reads/slurm_extract_reads.sh" | \
        sed "s/{NUM_FILES}/$num_of_files/g" | \
        sed "s/{THREAD}/$THREAD/g" | \
        sed "s/{BED_FILES_FOLDER}/$sed_bed_files/g" | \
        sed "s/{EXTRACT_READS_SCRIPT}/$sed_read_scripts/g" | \
        sed "s/{READ_BAM_FILE}/$sed_read_bam/g" | \
        sed "s/{READSIZE}/$READSIZE/g" > "${OUTPUT}/${chain_basename}_samplejob.sh"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        sbatch --wait "${OUTPUT}/${chain_basename}_samplejob.sh"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/merged_${chain_basename}.bed.time" \
        cat "${OUTPUT}/bedfiles/${chain_basename}/"*.bed | sortBed -i - | mergeBed -i - > "${OUTPUT}/bedfiles/${chain_basename}/merged_${chain_basename}.bed"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        mkdir -p "${OUTPUT}/bedfiles/${chain_basename}/retired_bed/"
        mkdir -p "${OUTPUT}/bedfiles/${chain_basename}/constant_bed/"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_regions.bed.time" \
        python3 "${SRCFOLDER}/4-extract_reads/get_retired_regions.py" \
        ${READSIZE} \
        "${OUTPUT}/bedfiles/${chain_basename}/merged_${chain_basename}.bed" \
        $i \
        "${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_regions.bed"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_reads.bed.time" \
        bash "${SRCFOLDER}/4-extract_reads/extract_reads_noprune.sh" \
        ${READ_BAM} \
        "${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_regions.bed" > "${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_reads.bed"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        mkdir -p "${OUTPUT}/bedfiles/${chain_basename}/reads/"
        cat "${OUTPUT}/bedfiles/${chain_basename}/"*.reads "${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_reads.bed" > "${OUTPUT}/bedfiles/${chain_basename}/updated_and_retired_reads.bed"
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/extract_sequences_full.time" \
        bash "${SRCFOLDER}/4-extract_reads/extract_sequence.sh" \
        ${FIRST_PAIR} \
        ${SECOND_PAIR} \
        "${OUTPUT}/bedfiles/${chain_basename}/updated_and_retired_reads.bed" \
        "${OUTPUT}/bedfiles/${chain_basename}/reads/"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/reads/paired_reads.time" \
        bwa mem -M -t $THREAD "${NEWREF}" "${OUTPUT}/bedfiles/${chain_basename}/reads/reads_1.fastq" "${OUTPUT}/bedfiles/${chain_basename}/reads/reads_2.fastq" | \
        samtools view -h -F4 | \
        samtools sort -m 16g -l0 > "${OUTPUT}/bedfiles/${chain_basename}/reads/paired_reads.bam"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/reads/singletons_reads.time" \
        bwa mem -M -t $THREAD "${NEWREF}" "${OUTPUT}/bedfiles/${chain_basename}/reads/singletons.fastq" | \
        samtools view -h -F4 | \
        samtools sort -m 16g -l0 > "${OUTPUT}/bedfiles/${chain_basename}/reads/singletons_reads.bam"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        bash "${SRCFOLDER}/0-align_reads.sh" \
        "${NEWREF}" \
        "${OUTPUT}/bedfiles/${chain_basename}/reads/reads" \
        "${OUTPUT}/bedfiles/${chain_basename}/reads/paired" \
        $THREAD
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        bash "${SRCFOLDER}/0-align_singletons.sh" \
        "${NEWREF}" \
        "${OUTPUT}/bedfiles/${chain_basename}/reads/singletons.fastq" \
        "${OUTPUT}/bedfiles/${chain_basename}/reads/singletons" \
        $THREAD
    done
}

function extract_reads_with_parallel(){
    mkdir "${OUTPUT}/bedfiles" 2> /dev/null
    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        mkdir -p "${OUTPUT}/bedfiles/${chain_basename}" 2> /dev/null
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/extract_regions.time" \
        bash ${SRCFOLDER}/4-extract_reads/extract_regions.sh \
        ${OUTPUT}/${chain_basename}_aligned_gaps.bam \
        ${OUTPUT}/bedfiles/${chain_basename}
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        local num_of_files=`ls -l ${OUTPUT}/bedfiles/${chain_basename}/*.bed | wc -l`
        local sed_bed_files=$(echo ${OUTPUT}\/bedfiles\/${chain_basename} | sed 's/\//\\\//g')
        local sed_read_scripts=$(echo ${SRCFOLDER}\/4-extract_reads\/extract_reads.sh | sed 's/\//\\\//g')
        local sed_read_bam=$(echo $READ_BAM | sed 's/\//\\\//g')
        cat "${SRCFOLDER}/4-extract_reads/parallel_extract_reads.sh" | \
        sed "s/{NUM_FILES}/$num_of_files/g" | \
        sed "s/{THREAD}/$THREAD/g" | \
        sed "s/{BED_FILES_FOLDER}/$sed_bed_files/g" | \
        sed "s/{EXTRACT_READS_SCRIPT}/$sed_read_scripts/g" | \
        sed "s/{READ_BAM_FILE}/$sed_read_bam/g" | \
        sed "s/{READSIZE}/$READSIZE/g" > ${OUTPUT}/${chain_basename}_samplejob.sh
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        bash ${OUTPUT}/${chain_basename}_samplejob.sh
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        local -a beds=($(ls ${OUTPUT}/bedfiles/${chain_basename}/*.bed | awk '{printf $1" ";}'))
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/merged_${chain_basename}.bed.time" \
        cat ${beds[@]} | sortBed -i - | mergeBed -i - > ${OUTPUT}/bedfiles/${chain_basename}/merged_${chain_basename}.bed
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        mkdir -p ${OUTPUT}/bedfiles/${chain_basename}/retired_bed
        mkdir -p ${OUTPUT}/bedfiles/${chain_basename}/constant_bed
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_regions.bed.time" \
        python3 ${SRCFOLDER}/4-extract_reads/get_retired_regions.py \
        ${READSIZE} \
        ${OUTPUT}/bedfiles/${chain_basename}/merged_${chain_basename}.bed \
        $i \
        ${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_regions.bed
    done

    # constant regions:
    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/constant_bed/constant_regions.bed.time" \
        python3 \
        "${SRCFOLDER}/4-extract_reads/get_constant_regions.py" \
        ${READSIZE} \
        $i \
        "${OUTPUT}/bedfiles/${chain_basename}/constant_bed/constant_regions.bed"
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/constant_bed/constant_reads.bed.time" \
        bash ${SRCFOLDER}/4-extract_reads/extract_reads_noprune.sh \
        ${READ_BAM} \
        ${OUTPUT}/bedfiles/${chain_basename}/constant_bed/constant_regions.bed > "${OUTPUT}/bedfiles/${chain_basename}/constant_bed/constant_reads.bed"
    done

    local SAMPLE=$(echo ${READ_BAM} | awk -F '/' '{printf $NF}' | awk -F '.' '{printf $1}')
    echo SAMPLE_ID is ${SAMPLE}
    for i in `ls ${OUTPUT}/*.chain`; do 
        local chain_basename=`basename $i | sed s/.chain//`
        cat <(bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:${SAMPLE}" ${NEWREF} ${FIRST_PAIR} ${SECOND_PAIR} | \
        samtools view -H -) <(/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/constant_bed/constant_liftOver.time" liftOver -minMatch=1 <(samtools view -h -L "${OUTPUT}/bedfiles/${chain_basename}/constant_bed/constant_regions.bed" ${READ_BAM} | bamToBed -i -) ${i} >("${SRCFOLDER}/5-merge/liftBedToSam" <(samtools view -L "${OUTPUT}/bedfiles/${chain_basename}/constant_bed/constant_regions.bed" ${READ_BAM}) - 4 3,4 1,2) "${OUTPUT}/bedfiles/${chain_basename}/constant_bed/constant_unmapped.bed") | \
        samtools sort -l5 > ${OUTPUT}/bedfiles/${chain_basename}/constant_bed/constant_lifted.bam
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_reads.bed.time" \
        bash ${SRCFOLDER}/4-extract_reads/extract_reads_noprune.sh \
        ${READ_BAM} \
        ${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_regions.bed > ${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_reads.bed
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        mkdir -p ${OUTPUT}/bedfiles/${chain_basename}/reads
        local -a reads_files=($(ls ${OUTPUT}/bedfiles/${chain_basename}/*.reads | awk '{printf $1" ";}'))
        cat ${reads_files[@]} ${OUTPUT}/bedfiles/${chain_basename}/retired_bed/retired_reads.bed > ${OUTPUT}/bedfiles/${chain_basename}/updated_and_retired_reads.bed
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/extract_sequences_full.time" \
        bash ${SRCFOLDER}/4-extract_reads/extract_sequence.sh \
        ${FIRST_PAIR} \
        ${SECOND_PAIR} \
        ${OUTPUT}/bedfiles/${chain_basename}/updated_and_retired_reads.bed \
        ${OUTPUT}/bedfiles/${chain_basename}/reads
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/reads/paired_reads.time" \
        bwa mem -M -t $THREAD "${NEWREF}" "${OUTPUT}/bedfiles/${chain_basename}/reads/reads_1.fastq" "${OUTPUT}/bedfiles/${chain_basename}/reads/reads_2.fastq" | \
        samtools view -h -F4 | \
        samtools sort -m 16g -l0 > ${OUTPUT}/bedfiles/${chain_basename}/reads/paired_reads.bam
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chain_basename}/reads/singletons_reads.time" \
        bwa mem -M -t $THREAD "${NEWREF}" "${OUTPUT}/bedfiles/${chain_basename}/reads/singletons.fastq" | \
        samtools view -h -F4 | \
        samtools sort -m 16g -l0 > ${OUTPUT}/bedfiles/${chain_basename}/reads/singletons_reads.bam
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        bash ${SRCFOLDER}/0-align_reads.sh \
        "${NEWREF}" \
        "${OUTPUT}/bedfiles/${chain_basename}/reads/reads" \
        "${OUTPUT}/bedfiles/${chain_basename}/reads/paired" \
        $THREAD
    done

    for i in `ls ${OUTPUT}/*.chain`; do
        local chain_basename=`basename $i | sed s/.chain//`
        bash ${SRCFOLDER}/0-align_singletons.sh \
        "${NEWREF}" \
        "${OUTPUT}/bedfiles/${chain_basename}/reads/singletons.fastq" \
        "${OUTPUT}/bedfiles/${chain_basename}/reads/singletons" \
        $THREAD
    done
}


function main(){
    mkdir -p "${OUTPUT}/"

    local chain_existed=$(ls ${OUTPUT}/*.chain 2> /dev/null)
    
    if [[ ${#chain_existed} -gt 0 ]]; then
        echo "Line "${LINENO}": In function "${FUNCNAME}": Use user offered chain "${chain_existed}
    else
        echo "Line "${LINENO}": In function "${FUNCNAME}": User does not offer chain file. Generating one."
        generate_chain
    fi

    generate_gaps

    align_gaps

    extract_reads_with_parallel
    
}


if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]
then
    main "$@"
fi