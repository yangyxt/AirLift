#!/bin/bash
# Need GNU parallel

module load parallel
source $(which env_parallel.bash)

#variables to change are:
#{BED_FILES_FOLDER}
#{READ_BAM_FILE}
#{EXTRACT_READS_SCRIPT}
#{NUM_FILES}
#{THREAD}

FILES=($(ls {BED_FILES_FOLDER}/*.bed | awk '{printf $1" ";}'))
BED_FILE=${FILES[${SLURM_ARRAY_TASK_ID}]}
OUT_READS=$(echo ${BED_FILE}.reads)
tmp_dir=/paedyl01/disk1/yangyxt/test_tmp

# Check the parallel command first 
env_parallel -j3 --tmpdir ${tmp_dir} -k --dryrun /usr/bin/time -v -p -o "${OUT_READS}.time" bash {EXTRACT_READS_SCRIPT} {READ_BAM_FILE} {1} {READSIZE} '>' {2} ::: $(echo ${FILES[@]}) ::: $(echo ${FILES[@]} | awk '{for(i=1;i<=NF;i++) printf $i".reads ";}')
env_parallel -j3 --tmpdir ${tmp_dir} -k --progress /usr/bin/time -v -p -o "${OUT_READS}.time" bash {EXTRACT_READS_SCRIPT} {READ_BAM_FILE} {1} {READSIZE} '>' {2} ::: $(echo ${FILES[@]}) ::: $(echo ${FILES[@]} | awk '{for(i=1;i<=NF;i++) printf $i".reads ";}')

module unload parallel
