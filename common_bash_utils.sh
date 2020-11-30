#!/bin/bash

function timestamp(){
	date
}

echo_line_no(){
    grep -n "$1" $0 |  sed "s/echo_line_no//" 
    # grep the line(s) containing input $1 with line numbers
    # replace the function name with nothing 
}

function display_table(){
	if [[ -z $2 ]]; then local rows=4; else local rows=${2}; fi
	if [[ -z $3 ]]; then local delimiter="\t"; else local delimiter=${3}; fi

	python3 /paedyl01/disk1/yangyxt/ngs_scripts/display_table.py \
	-df $1 \
	-nr ${rows}
}

# Define a fucntion to judge if some element in an array
function containsElement () {
	local e match="$1"
	shift
	for e; do [[ "$e" == "$match" ]] && return 0; done
	return 1
}

function randomID(){
    dd bs=24 count=1 status=none if=/dev/urandom | base64 | tr +/ _
}

function build_genomicsDB_from_scratch(){
    local genomicDB=${1}
    local chr=${2}
    local sample_map=${3}
    local gatk=/home/yangyxt/software/gatk-4.1.8.1/gatk

    if [[ -d ${genomicDB} ]]; then rm -R -f ${genomicDB}; fi

    local batch_size=$(wc -l ${sample_map} | awk '{printf $1}')
    echo "Line "${LINENO}": In function "${FUNCNAME}: The batch size for probe ${probe} is ${batch_size}

    # After running the script above, all gvcfs corresponding to each sample should be already there. Now we generate a script, which is used to generate the script for consolidated genotyping
    time ${gatk} --java-options "-Xmx8g -Xms2g" GenomicsDBImport \
        --tmp-dir /paedyl01/disk1/yangyxt/test_tmp \
        --genomicsdb-workspace-path ${genomicDB} \
        -R ${ref_gen}/ucsc.hg19.fasta \
        --batch-size ${batch_size} \
        --sample-name-map ${sample_map} \
        --reader-threads 5 \
        --intervals chr${chr}
    check_return_code
}

function check_return_code(){
	local return_code=$(echo $?)
	echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": The last step's return code is "${return_code}
	if [[ ${return_code} -gt 0 ]]; then echo "The last step does not finish in normal way. Exiting the whole script now." && exit 1; else echo "The last step finished properly. Continue"; fi
}

function get_current_dir(){
    local input=${1}
    if [[ -f ${input} ]]; then
        if [[ ${input} =~ /\// ]]; then
            cd=$(echo ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        else
            # The input does not have path info
            cd=$(readlink -f ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        fi
    elif [[ -d ${input} ]]; then
        if [[ ${input} =~ /\// ]]; then
            cd=$(echo ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        else
            # The input does not have path info
            cd=$(readlink -f ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        fi
    elif [[ ${input} =~ \.[a-z]+$ ]]; then
        cd=$(echo ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
        echo ${cd}
    else
        if [[ ${input} =~ /\// ]]; then
            cd=$(echo ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        else
            # The input does not have path info
            cd=$(readlink -f ${input} | awk -F '/' '{for(i=1;i<NF;i++) printf $i"/";}' | awk '{gsub(/\/$/, ""); print}')
            echo ${cd}
        fi
    fi
}

function gnws(){
    # gnws stands for get name of the file without suffix. (still with the path)
    local input=${1}
    echo ${input} | awk -F '.' '{for(i=1;i<NF;i++) printf $i".";}' | awk '{gsub(/\.$/, ""); print}'
}


function get_name(){
    # Get name of the file without suffix and without path
    local input=${1}

    if [[ -d ${input} ]]; then
        local last=$(echo ${input} | awk -F '/' '{printf $NF;}')
        if [[ ${#last} -eq 0 ]]; then
            echo ${input} | awk -F '/' '{printf $(NF-1);}'
        else
            echo ${input} | awk -F '/' '{printf $NF;}'
        fi
    elif [[ -f ${input} ]]; then
        echo ${input} | awk -F '/' '{printf $NF;}' | awk -F '.' '{for(i=1;i<NF;i++) printf $i".";}' | awk '{gsub(/\.$/, ""); print}'
    fi
}

function join_by(){ local IFS="$1"; shift; echo "$*"; }

# Define a fucntion to judge if some element in an array
function containsElement () {
	local e match="$1"
	shift
	for e; do [[ "$e" == "$match" ]] && return 0; done
	return 1
}

function drop_duplicated_rows () {
	export TMPDIR=/paedyl01/disk1/yangyxt/test_tmp/
	head -1 $1 > $1.tmp
	tail -n +2 $1 | sort - | uniq - >> $1.tmp
	mv $1.tmp $1
}

# Define a function that check the size of a table once it below 2k the whole script exit 1.
function check_result_size() {
	local size=$(ls -l ${1} | awk '{gsub(/ +/, " "); print}' | awk '{print $5}')
	local lines=$(wc -l ${1} | awk '{print $1}')
	if [[ ${size} -le 2800 ]] || [[ ${lines} -le 1 ]]
	then 
		echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "The table "$1" after processing is just too oddly small, call it an end for now."
		ls -lh $1
		exit 1
	else
		echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "The table "$1" after processing is not too small to only fit a header line. Continue"
		ls -lh $1
		if [[ $1 =~ (\.bed$|\.txt$|\.ped$|\.tsv$|\.csv$|\.table$|\.tmp$|\.sam$) ]]; then display_table $1; fi
	fi
}

function check_odd_tab() {
    local -a odd_tabs=($(bash /paedyl01/disk1/yangyxt/ngs_scripts/check_odd_tab_rows.sh $1))
    echo "Line "${LINENO}": In function "${FUNCNAME}: odd_tabs = ${odd_tabs[@]}
    if [[ ${#odd_tabs[@]} -ge 2 ]]
    then
        sed -i 's/\t*$//' $1
        local -a second_odd=($(bash /paedyl01/disk1/yangyxt/ngs_scripts/check_odd_tab_rows.sh $1))
        echo "Line "${LINENO}": In function "${FUNCNAME}: second_odd = ${second_odd[@]}
        if [[ ${#second_odd[@]} -ge 2 ]]
        then
            echo "Line "${LINENO}": In function "${FUNCNAME}: "The table "$1" has some odd rows with errorenous tabs."
            echo "Line "${LINENO}": In function "${FUNCNAME}: "The rows having odd tabs are ${second_odd[@]}"
            exit 1
        else
            echo "Line "${LINENO}": In function "${FUNCNAME}: "The table"$1" is normal and continue."
            ls -lh $1
        fi
    else
        echo "Line "${LINENO}": In function "${FUNCNAME}: "The table"$1" is normal and continue."
        ls -lh $1
    fi 
}

function drop_duplicated_cols(){
	local table=$1
	local -a x_cols=($(awk -F '\t' 'NR == 1{ for(i=1;i<=NF;i++) {if ($i ~ /_x$/) printf i" ";} }' < ${table}))
	local -a y_cols=($(awk -F '\t' 'NR == 1{ for(i=1;i<=NF;i++) {if ($i ~ /_y$/) printf i" ";} }' < ${table}))
    local -a extra_cols=($(awk -F '\t' 'NR == 1{ for(i=1;i<=NF;i++) {if ($i ~ /\.1$/) printf i" ";} }' < ${table}))
    local -a tmp_cols=($(awk -F '\t' 'NR == 1{ for(i=1;i<=NF;i++) {if ($i ~ /^uniq_ID/) printf i" ";} }' < ${table}))

	if [[ ${#x_cols[@]} -eq ${#y_cols[@]} ]]; then
		if [[ ${#x_cols[@]} -gt 0 ]]; then
			local tobe_dropped_cols=$(echo ${y_cols[@]} | awk '{for(i=1;i<NF;i++) {printf $i",";} printf $NF}')
			cut --complement -f ${tobe_dropped_cols} < ${table} > ${table}.tmp && \
			awk -F '\t' '{ if(NR==1) {for (i=1;i<NF;i++) {gsub(/_x$/, "", $i); printf $i"\t";} printf $NF"\n";} else print; }' < ${table}.tmp > ${table}
		else
			echo "Line "${LINENO}": In function "${FUNCNAME}: "We dont have duplicated columns for this table: "${table}
		fi
	elif [[ ${#x_cols[@]} -gt 0 ]]; then
		echo "Line "${LINENO}": In function "${FUNCNAME}: "Seems that the _y suffix cols has been removed, need to remove _x suffix for those labels"
		awk -F '\t' '{ if(NR==1) {for (i=1;i<NF;i++) {gsub(/_x$/, "", $i); printf $i"\t";} printf $NF"\n";} else print; }' < ${table} > ${table}.tmp && mv ${table}.tmp ${table}
	else
		echo "Line "${LINENO}": In function "${FUNCNAME}: "We seems to detect duplicated rows but the number between _x suffix and _y suffix cols is not consistent. Check the table."
		display_table ${table}
	fi

    if [[ ${#extra_cols[@]} -gt 0 ]]; then
        local tobe_dropped_cols=$(echo ${extra_cols[@]} | awk '{for(i=1;i<NF;i++) {printf $i",";} printf $NF}')
        cut --complement -f ${tobe_dropped_cols} < ${table} > ${table}.tmp && mv ${table}.tmp ${table}
    fi

    if [[ ${#tmp_cols[@]} -gt 0 ]]; then
        local tobe_dropped_cols=$(echo ${tmp_cols[@]} | awk '{for(i=1;i<NF;i++) {printf $i",";} printf $NF}')
        cut --complement -f ${tobe_dropped_cols} < ${table} > ${table}.tmp && mv ${table}.tmp ${table}
    fi
}

function awk_merge_two_tables(){
	# We assume that two tables are with headers.
	# Better merge two tables on a common sharing column with the same header label.
	local left_table=$1
	local right_table=$2
	local left_on_column=$3
	local right_on_column=$4
	local merged_table=$5


	if [[ ${left_on_column} =~ ^[0-9]+$ ]]
	then
		local left_col_ind=${left_on_column}
		local right_col_ind=${right_on_column}
	elif [[ ${left_on_column} == ${right_on_column} ]]
	then
		local left_col_ind=$(head -1 ${left_table} | awk -F '\t' '{for(i=1;i<=NF;i++) {if ($i == "'${left_on_column}'") printf i;}}')
		local right_col_ind=$(head -1 ${right_table} | awk -F '\t' '{for(i=1;i<=NF;i++) {if ($i == "'${right_on_column}'") printf i;}}')
	elif [[ ${left_on_column} != ${right_on_column} ]]
	then
		awk -F '\t' 'BEGIN{OFS=FS="\t";} {if (NR == 1) {gsub("'${right_on_column}'", "'${left_on_column}'"); print;} else print;}' < ${right_table} > ${right_table}.tmp && \
		mv ${right_table}.tmp ${right_table}
		check_return_code
		check_result_size ${right_table}
		local left_col_ind=$(head -1 ${left_table} | awk -F '\t' '{for(i=1;i<=NF;i++) {if ($i == "'${left_on_column}'") printf i;}}')
		local right_col_ind=$(head -1 ${right_table} | awk -F '\t' '{for(i=1;i<=NF;i++) {if ($i == "'${left_on_column}'") printf i;}}')
	fi

	echo "Line "${LINENO}": In function "${FUNCNAME}: Left Table is $(ls -lh ${left_table})
	echo "Line "${LINENO}": In function "${FUNCNAME}: Right Table is $(ls -lh ${right_table})
	echo "Line "${LINENO}": In function "${FUNCNAME}: Left table target column is ${left_on_column} while the column index is ${left_col_ind}
	echo "Line "${LINENO}": In function "${FUNCNAME}: Right table target column is ${right_on_column} while the column index is ${right_col_ind}

	awk 'BEGIN {FS=OFS="\t";} \
	NR == FNR {aa[$'${right_col_ind}'] = $0; next;} \
	{print $0, aa[$'${left_col_ind}'];} \
	' ${right_table} ${left_table} > ${merged_table}

	check_return_code
	check_result_size ${merged_table}
}


function awk_filter_out_columns(){
	#$1 should be target column labels joined by ","
	table=$1
	columns=$2
	local -a col_arr=$(echo ${columns} | awk 'BEGIN{RS=",";} {printf $1" "}')
	echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "The to be removed columns are "${col_arr[@]}

	for col in ${col_arr[@]}
	do
		if [[ ${col} =~ ^[0-9]+$ ]]
		then
			local if_index=index
		else
			if_index="not indices" && break
		fi
	done
	check_return_code
	echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): "Are we inputting a list of column labels or a list of column indices? "${if_index}

	if [[ ${if_index} == "not indices" ]]
	then
		local -a col_arr=$(head -1 ${table} | awk -v cols="${col_arr[*]}" 'BEGIN{split(cols,col_ar," ");RS = "\t";} {for (n in col_ar) {if($1 == col_ar[n]) printf NR" ";}}')
		check_return_code
		echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): We have converted the column labels to column indices so we can use cut command to cut the column off:$'\n' ${col_arr[@]}
	else
		echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp): The to be filtered out column indices are:$'\n' ${col_arr[@]}.
	fi

	cut --complement -f $(join_by , "${col_arr[@]}") < ${table} > ${table}.tmp && \
	mv ${table}.tmp ${table}
	check_result_size ${table}
}


function return_array_intersection(){
    # Both arrays should not contain item values with space
    local -a array1=($(echo ${1}))
    local -a array2=($(echo ${2}))

    for item in ${array1[@]}; do
        if [[ ${2} =~ ${item} ]]; then
            result+=(${item})
        fi
    done
    echo ${result[@]}
}


function check_vcf_bam_index(){
    if [[ ${1} =~ \.bam$ ]]; then
        local bam=${1}
        if [[ -f ${bam}.bai || -f ${bam::-1}i ]]; then
            echo "Line "${LINENO}": In function "${FUNCNAME}: "${bam} file does have an index file"
        else
            module load samtools
            echo "Line "${LINENO}": In function "${FUNCNAME}: ${bam} file does not have a corresponding index file. Generating one with samtools
            samtools index ${bam}
            module unload samtools
        fi
    elif [[ ${1} =~ \.vcf$ ]]; then 
        local vcf=${1}
        if [[ -f ${vcf}.idx ]]; then
            echo "Line "${LINENO}": In function "${FUNCNAME}: ${vcf} file does have an index file generated by GATK IndexFeatureFile
        else
            echo "Line "${LINENO}": In function "${FUNCNAME}: ${vcf} file does not have a corresponding index file generated by GATK IndexFeatureFile, generating one now.
            local gatk=/home/yangyxt/software/gatk-4.1.8.1/gatk
            ${gatk} IndexFeatureFile -I ${vcf}
        fi
    elif [[ ${1} =~ \.vcf.gz$ ]]; then
        local gzvcf=${1}
        if [[ -f ${gzvcf}.tbi ]]; then
            echo "Line "${LINENO}": In function "${FUNCNAME}: ${gzvcf} is gzipped and does have an index file generated by tabix.
        else
            echo "Line "${LINENO}": In function "${FUNCNAME}: ${gzvcf} is gzipped and does not have a corresponding index file generated by tabix. Generating one now.
            tabix -f -p vcf ${gzvcf}
        fi
    fi 
}

function reset_contig_order_in_vcfhead(){
    local vcf=$1

    awk -F '\t' '{if ($1 ~ /^##contig=/) next; else print}' < ${vcf} > ${vcf}.without_contig
    first_line=$(awk -F '\t' '{if ($1 ~ /^##contig=/) printf NR"\n"; else next;}' < ${vcf} | sort -n - | awk 'NR==1{printf $1;}')
    awk -F '\t' '{if (NR < '${first_line}') print; else exit 0;}' < ${vcf} > ${vcf}.before_contig
    awk -F '\t' '{if (NR >= '${first_line}') print;}' < ${vcf}.without_contig > ${vcf}.after_contig
    cat ${vcf}.before_contig /paedyl01/disk1/yangyxt/indexed_genome/backup_sorted_vcf_contig_header.txt ${vcf}.after_contig > ${vcf}.new && \
    rm ${vcf}.before_contig && rm ${vcf}.after_contig && rm ${vcf}.without_contig && \
    mv ${vcf}.new ${vcf}
}

function return_array_substraction(){
    # Only return the items which array1 has and array2 does not have
    # Cannot return the items which array2 has and array1 does not have
    local -a array1=($(echo ${1}))
    local -a array2=($(echo ${2}))

    for item in ${array1[@]}; do
        if [[ ${2} =~ ${item} ]]; then
            result+=(${item})
        else
            subs+=(${item})
        fi
    done
    echo ${subs[@]}
}

function merge_columns(){
    local delimiter=":"
    local header="no"

    local OPTIND d l i o c
    while getopts d::l::i:o::c: args
    do 
        case ${args} in
            d) local delimiter=$OPTARG ;;
            i) local input=$OPTARG ;;
            o) local output=$OPTARG ;;
            l) local merged_label=$OPTARG ;;
            c) local to_be_merged_columns=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    if [[ ${#output} -le 0 ]]; then local output=${input}; fi

    local -a col_arr=($(echo ${to_be_merged_columns} | awk -F ',' '{for(i=1;i<=NF;i++) printf $i" ";}'))
    if [[ ${col_arr} =~ [0-9]+ ]]; then
        :
    else
        for col in ${col_arr}; do
            local col_inds=${col_inds}$(head -1 ${input} | awk -F '\t' '{for(i=1;i<=NF;i++) {if($i == "'${col}'") printf i" ";}}')
        done
        local -a col_arr=(${col_inds})
    fi
    
    # Note here that the col_arr we have is 1-indexed, not 0-indexed.
    tmp_bed_right=$(gnws ${input}).$(randomID).bed
    tmp_bed_left=$(gnws ${input}).$(randomID).bed
    
    cut -f $(join_by , "${col_arr[@]}") ${input} | awk -F '\t' '{for(i=1;i<NF;i++) {printf $i"'${delimiter}'";} printf $NF"\n";}' > ${tmp_bed_right}
    cut --complement -f $(join_by , "${col_arr[@]}") ${input} > ${tmp_bed_left}
    
    awk 'BEGIN {FS=OFS="\t";} \
	NR == FNR {aa[NR] = $0; next;} \
	NR > FNR {for(i=1;i<NF;i++) {if (i != '${col_arr}') printf $i"\t"; else printf aa[FNR]"\t"$i"\t";} printf $NF"\n";}' ${tmp_bed_right} ${tmp_bed_left} > ${output} && \
    rm ${tmp_bed_left} && rm ${tmp_bed_right}
    
    echo ${output}
}

function insert_column(){
    local input=${1}
    local col_ind=${2}
    local insert_value=${3}

    if [[ -f ${insert_value} ]]; then
        awk 'BEGIN {FS=OFS="\t";} \
        NR == FNR {aa = $0; next;} \
        {for(i=1;i<NF;i++) {if (i != '${col_ind}') printf $i"\t"; else printf aa"\t"$i"\t";} printf $NF"\n";}' ${insert_value} ${input} > ${input}.tmp && \
        mv ${input}.tmp ${input}
    else
        awk -F '\t' '{for(i=1;i<NF;i++) {if(i != '${col_ind}') printf $i"\t"; else printf "'${insert_value}'\t"$i"\t";} printf $NF"\n";}' < ${input} > ${input}.tmp && \
        mv ${input}.tmp ${input}
    fi   
}

function merge_interval_bed(){
    local bed=${1}
    # The input bed file should follow the standard format.
    local tmp_bed=$(gnws ${bed}).$(randomID).bed

    sort -k1,1 -k2,2n ${bed} | bedtools merge -s -c 4 -o distinct -i - > ${tmp_bed} && \
    mv ${tmp_bed} ${bed}
}

function fetch_interval_bed_using_gene_symbols(){
    local genelist=${1}
    local __output_bed=${2}
    local tmp_bed=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).bed

    bash ${central_scripts}/generate_bed_by_gene_names.sh -g ${genelist} -o ${tmp_bed}
    display_table ${tmp_bed}

    # Do this step is only to consider the fact that Gene symbol and exon are respectively stored in column4 and column 5
    merge_columns -i ${tmp_bed} -c 4,5

    # Do this step is for inserting integer value in column5
    insert_column ${tmp_bed} 5 0

    eval $__output_bed="'${tmp_bed}'"
}

function remove_sec_align_rec(){
    module load samtools

    local input=${1}
    local tmp_header=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).sam.head.txt

    samtools view -H ${input} > ${tmp_header}
    samtools view ${input} | awk -F '\t' '{if ($2 < 256) print;}' > ${input}.tmp.sam

    if [[ ${input} =~ \.sam$ ]]; then
        cat ${tmp_header} ${input}.tmp.sam > ${input} && rm ${tmp_header} && rm ${input}.tmp.sam
    elif [[ ${input} =~ \.bam$ ]]; then
        cat ${tmp_header} ${input}.tmp.sam > ${input::-4}.sam && rm ${tmp_header} && rm ${input}.tmp.sam && \
        samtools view -hSb ${input::-4}.sam > ${input} && rm ${input::-4}.sam
    fi
}

function remove_error_line_bed(){
    local bed=$1

    awk -F '\t' '{if ($2 >= $3) next; else print;}' < ${bed} > ${bed}.tmp && \
    mv ${bed}.tmp ${bed}
}

function convert_bam_table(){
    local input=${1}
    local input_dir=$(get_current_dir ${input})

    if [[ -z ${2} ]]; then
        local output_table=${input_dir}/$(randomID).txt
    else
        local output_table=${2}
    fi

    samtools view ${input} > ${output_table}
    echo ${output_table}
}

function condense_on_col(){
    local input=${1}
    local col_name=${2}

    python3 /paedyl01/disk1/yangyxt/ngs_scripts/condense_on_certain_col.py ${input} ${col_name} ${3}
    drop_duplicated_rows ${input}
}

function get_RG_SM(){
    local alignment=${1}
    samtools view -H ${alignment} | grep @RG | awk -F '\t' '{for (i=1;i<=NF;i++) if($i ~ /^SM\:/) printf $i;}' | awk -F ':' '{printf $NF;}'
}

function get_RG_ID(){
    local alignment=${1}
    samtools view -H ${alignment} | grep @RG | awk -F '\t' '{for (i=1;i<=NF;i++) if($i ~ /^ID\:/) printf $i;}' | awk -F ':' '{printf $NF;}'
}

function insert_align_records_to_another(){
    local OPTIND i t 
    while getopts i:t: args
    do 
        case ${args} in
            i) local insert_part_alignment=$OPTARG ;;
            t) local tobe_inserted_alignment=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    local -a insert_queries=($(samtools view ${insert_part_alignment} | awk -F '\t' '{print $1}' | sort -k 1 | uniq - | awk '{printf $1" ";}')) && \
    samtools view -H ${tobe_inserted_alignment} > ${tobe_inserted_alignment}.tmp && \
    printf "@PG\tID:insert_align_records_to_another\tCL:records inserted from %s\n" "${insert_part_alignment}" >> ${tobe_inserted_alignment}.tmp && \
    samtools view ${tobe_inserted_alignment} | awk -F '\t' -v qs="${insert_queries[*]}" 'BEGIN{split(qs,qarr," "); for(i in qarr) q_arr[qarr[i]]=qarr[i];} {if($1 in q_arr) next; else print;}' >> ${tobe_inserted_alignment}.tmp && \
    samtools view ${insert_part_alignment} >> ${tobe_inserted_alignment}.tmp
    local insert_ID=$(get_RG_ID ${insert_part_alignment})
    local insert_SM=$(get_RG_SM ${insert_part_alignment})
    local tobe_insert_ID=$(get_RG_ID ${tobe_inserted_alignment})
    local tobe_insert_SM=$(get_RG_SM ${tobe_inserted_alignment})

    if [[ ${insert_ID} == ${tobe_insert_ID} ]] && [[ ${insert_SM} == ${tobe_insert_SM} ]]; then
        if [[ ${tobe_inserted_alignment} =~ \.bam$ ]]; then samtools sort -O bam -o ${tobe_inserted_alignment} < ${tobe_inserted_alignment}.tmp && samtools index ${tobe_inserted_alignment}; fi
        if [[ ${tobe_inserted_alignment} =~ \.sam$ ]]; then samtools sort -O sam -o ${tobe_inserted_alignment} < ${tobe_inserted_alignment}.tmp; fi
    fi
}

function identify_max_min_bash_array(){
    # When calling this function, input array must be enclosed in double quotes
    local OPTIND m n a
    while getopts m::n::a: args
    do 
        case ${args} in
            a) local array=$OPTARG ;; 
            m) local __max=$OPTARG ;;
            n) local __min=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    local -a array=($(echo ${array}))

    local max=${array[0]}
    local min=${array[0]}

    for i in "${array[@]}"; do
        (( i > max )) && local max=$i
        (( i < min )) && local min=$i
    done

    if [[ ${#__max} -gt 0 ]]; then eval $__max="'${max}'"; fi
    if [[ ${#__min} -gt 0 ]]; then eval $__min="'${min}'"; fi
}

function add_PG_to_bam(){
    local OPTIND a p c 
    while getopts a:p:c:: args
    do 
        case ${args} in
            a) local alignment_file=$OPTARG ;;
            p) local PG_ID=$OPTARG ;;
            c) local PG_CL=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    samtools view -H ${alignment_file} > ${alignment_file}.head
    if [[ ${#PG_CL} -gt 0 ]]; then
        printf "@PG\tID:%s\tCL:%s\n" "${PG_ID}" "${PG_CL}" >> ${alignment_file}.head
    else
        printf "@PG\tID:%s\n" "${PG_ID}" >> ${alignment_file}.head
    fi

    samtools view ${alignment_file} >> ${alignment_file}.head && \
    if [[ ${alignment_file} =~ \.bam$ ]]; then
        samtools view ${alignment_file}.head -b -o ${alignment_file} && \
        samtools index ${alignment_file}
    elif [[ ${alignment_file} =~ \.sam$ ]]; then
        samtools view ${alignment_file}.head -o ${alignment_file}
    fi
}


function merge_logs(){
    local -a logs=($(echo $1 | awk -F ',' '{for(i=1;i<=NF;i++) printf $i" ";}'))
    local merged_log=$2

    for log in ${logs[@]}; do
		awk '{print;} END{printf "\n\n";}' < ${log} > ${log}.txt && \
		mv ${log}.txt ${log}
	done

    cat ${logs[@]} > ${merged_log} && \
    rm ${logs[@]}
}


function gathervcfs(){
    # This function is allowing to gather vcf with different header info.
    # So the first thing to do is to merge the header info.
    local -a vcfs=($(echo ${1} | awk 'BEGIN{FS=",";} {for(i=1;i<=NF;i++) printf $i" ";}'))
    echo "Line "${LINENO}": In function "${FUNCNAME}: the vcfs are ${vcfs[@]}
    local merged_vcf=${2}
    local output_shell_script=${3}
    local gatk=/home/yangyxt/software/gatk-4.1.8.1/gatk

    local common_header=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).vcf.header
    echo "Line "${LINENO}": In function "${FUNCNAME}: The temporary temp header file is ${common_header}

    if [[ $(echo ${1} | awk -F ',' '{printf $1}') =~ \.vcf\.gz$ ]]; then
		zcat $(echo ${1} | awk -F ',' '{printf $1}') | head -1 - > ${common_header}.firstline
	else
		head -1 $(echo ${1} | awk -F ',' '{printf $1}') > ${common_header}.firstline
	fi

    for vcf in ${vcfs[@]}; do
        if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
            zcat ${vcf} | awk 'NR > 1{if($1 ~ /^##/) print; else next;}' >> ${common_header}
        else
            awk 'NR > 1{if($1 ~ /^##/) print; else next;}' < ${vcf} >> ${common_header}
        fi
    done

    cat ${common_header} | sort - | uniq - > ${common_header}.tmp
    cat ${common_header}.firstline ${common_header}.tmp > ${common_header} && rm ${common_header}.tmp && rm ${common_header}.firstline
    # Note here now the contig order is not consistent with the one in ucsc.hg19.fasta. We need to resort the order
    reset_contig_order_in_vcfhead ${common_header}

    # Substitude the common header to each vcf in the to_be merged list. Save their original head as a back up.
    for vcf in ${vcfs[@]}; do
        local tmp_vcf_body=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).vcf

        if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
            zcat ${vcf} | awk '{if($1 !~ /^##/) print;}' > ${tmp_vcf_body}
            display_table ${tmp_vcf_body}
            zcat ${vcf} | awk '{if($1 ~ /^##/) print;}' > ${vcf::-3}.backup.head
            cat ${common_header} ${tmp_vcf_body} > ${vcf::-3} && \
            bgzip -f ${vcf::-3} && tabix -f -p vcf ${vcf} && rm ${tmp_vcf_body}
        else
            cat ${vcf} | awk '{if($1 !~ /^##/) print;}' > ${tmp_vcf_body}
            display_table ${tmp_vcf_body}
            cat ${vcf} | awk '{if($1 ~ /^##/) print;}' > ${vcf}.backup.head
            cat ${common_header} ${tmp_vcf_body} > ${vcf} && \
            check_vcf_bam_index ${vcf} && \
            rm ${tmp_vcf_body}
        fi
    done
    
    # Run GATK's GatherVCF tool.
    echo ${vcfs[@]} | awk 'BEGIN{printf "\
    time '${gatk}' GatherVcfs \\\n";} {for (i=1;i<=NF;i++) {printf "\
    -I "$i" \\\n";} printf "\
    -O '${merged_vcf}'";}' | awk '{gsub(/^ {4}/, ""); print;}' > ${output_shell_script} && \
    bash ${output_shell_script} && \
    rm ${vcfs[@]}
}

function mergevcfs(){
    # This function is intended to merge variants info from multiple samples(each vcf contains a separate sample ID)
    # This function is allowing to gather vcf with different header info.
    # So the first thing to do is to merge the header info.
    local gatk=/home/yangyxt/software/gatk-4.1.8.1/gatk
    local -a vcfs=($(echo ${1} | awk 'BEGIN{FS=",";} {for(i=1;i<=NF;i++) printf $i" ";}'))
    echo "Line "${LINENO}": In function "${FUNCNAME}: the vcfs are ${vcfs[@]}
    local merged_vcf=${2}
    local output_shell_script=${3}

    local common_header=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).vcf.header
    echo "Line "${LINENO}": In function "${FUNCNAME}: The temporary temp header file is ${common_header}

    if [[ $(echo ${1} | awk -F ',' '{printf $1}') =~ \.vcf\.gz$ ]]; then
		zcat $(echo ${1} | awk -F ',' '{printf $1}') | head -1 > ${common_header}.firstline
	else
		head -1 $(echo ${1} | awk -F ',' '{printf $1}') > ${common_header}.firstline
	fi

    for vcf in ${vcfs[@]}; do
        if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
            zcat ${vcf} | awk 'NR > 1{if($1 ~ /^##/) print; else exit 0;}' >> ${common_header}
        else
            awk 'NR > 1{if($1 ~ /^##/) print; else exit 0;}' < ${vcf} >> ${common_header}
        fi
    done

    cat ${common_header} | sort - | uniq - > ${common_header}.tmp
    cat ${common_header}.firstline ${common_header}.tmp > ${common_header} && rm ${common_header}.tmp && rm ${common_header}.firstline
    # Note here now the contig order is not consistent with the one in ucsc.hg19.fasta. We need to resort the order
    reset_contig_order_in_vcfhead ${common_header}

    # Substitude the common header to each vcf in the to_be merged list. Save their original head as a back up.
    for vcf in ${vcfs[@]}; do
        local tmp_vcf_body=/paedyl01/disk1/yangyxt/test_tmp/$(randomID).vcf

        if [[ ${vcf} =~ \.vcf\.gz$ ]]; then
            zcat ${vcf} | awk '{if($1 !~ /^##/) print;}' > ${tmp_vcf_body}
            display_table ${tmp_vcf_body}
            zcat ${vcf} | awk '{if($1 ~ /^##/) print;}' > ${vcf::-3}.backup.head
            cat ${common_header} ${tmp_vcf_body} > ${vcf::-3} && \
            bgzip -f ${vcf::-3} && tabix -f -p vcf ${vcf} && rm ${tmp_vcf_body}
        else
            cat ${vcf} | awk '{if($1 !~ /^##/) print;}' > ${tmp_vcf_body}
            display_table ${tmp_vcf_body}
            cat ${vcf} | awk '{if($1 ~ /^##/) print;}' > ${vcf}.backup.head
            cat ${common_header} ${tmp_vcf_body} > ${vcf} && \
            check_vcf_bam_index ${vcf} && \
            rm ${tmp_vcf_body}
        fi
    done

    # check tabix index
    for vcf in ${vcfs[@]}; do
        if [[ ${vcf} =~ \.vcf$ ]]; then
            bgzip -c ${vcf} > ${vcf}.gz && tabix -f -p vcf ${vcf}.gz
        elif [[ ${vcf} =~ \.vcf\.gz$ ]]; then 
            tabix -f -p vcf ${vcf}
        fi
    done

    local -a gz_vcfs
    for vcf in ${vcfs[@]}; do
        if [[ ${vcf} =~ \.vcf$ ]]; then
            gz_vcfs+=("${vcf}.gz")
        elif [[ ${vcf} =~ \.vcf\.gz$ ]]; then 
            gz_vcfs+=("${vcf}")
        fi
    done

    # Run vcf-merge from  vcftools.
    if [[ ${merged_vcf} =~ \.gz$ ]]; then
        echo ${gz_vcfs[@]} | awk 'BEGIN{printf "\
        module load vcftools \n\
        time vcf-merge ";} {for (i=1;i<=NF;i++) {printf \
        $i" ";} printf "\
        | bgzip -c > '${merged_vcf}'\n";}' | awk '{gsub(/^ {4}/, ""); print;} END{printf "\nmodule unload vcftools"}' > ${output_shell_script} && \
        bash ${output_shell_script}
    elif [[ ${merged_vcf} =~ \.vcf$ ]]; then
        echo ${gz_vcfs[@]} | awk 'BEGIN{printf "\
        module load vcftools \n\
        time vcf-merge ";} {for (i=1;i<=NF;i++) {printf \
        $i" ";} printf "\
        > '${merged_vcf}'\n";}' | awk '{gsub(/^ {4}/, ""); print;} END{printf "\nmodule unload vcftools"}' > ${output_shell_script} && \
        bash ${output_shell_script}
    fi
}

function check_chaos_sorting(){
	local gene_col=$(awk -F '\t' '{if(NR==1) {for(i=1;i<=NF;i++) if($i == "Gene.refGene") printf i}}' < $1)
	local aachange_col=$(awk -F '\t' '{if(NR==1) {for(i=1;i<=NF;i++) if($i == "AAChange.refGene") printf i}}' < $1)
	echo "Line "${LINENO}": In function "${FUNCNAME}: The column index of Gene.refGene is ${gene_col}
	echo "Line "${LINENO}": In function "${FUNCNAME}: The column index of AAChange.refGene is ${aachange_col}
	local cant_match_rows=$(awk -F '\t' '{if (($'${aachange_col}' == ".")||($'${aachange_col}' == "UNKNOWN")) next; else printf $'${gene_col}'"\t"$'${aachange_col}'"\n";}' < $1 | awk -F '\t' '{regex = "^"$1":"; if (($2 !~ regex)&&($1 !~ /;/)&&($2 !~ /^[0-9]+-[0-9]+:/)) print NR}' | wc -l - | awk '{print $1}')
	local total_rows=$(wc -l $1 | awk '{print $1}')
	echo "Line "${LINENO}": In function "${FUNCNAME}: The number of rows that do not have a correct correspondence between AAChange.refGene and Gene.refGene is ${cant_match_rows}
	echo "Line "${LINENO}": In function "${FUNCNAME}: The number of total rows in the inspecting table is ${total_rows}
	if [[ ${cant_match_rows} -gt 1 ]]; then awk -F '\t' '{if (($'${aachange_col}' == ".")||($'${aachange_col}' == "UNKNOWN")) next; else printf $'${gene_col}'"\t"$'${aachange_col}'"\n";}' < $1 | awk -F '\t' '{regex = "^"$1":"; if (($2 !~ regex)&&($1 !~ /;/)&&($2 !~ /^[0-9]+-[0-9]+:/)) printf NR"\t"$1"\t"$2"\n"}' > ${1::-4}.chaos_order_rows.txt; fi
	local ratio=$(echo "${cant_match_rows}/${total_rows}" | bc -l)
	echo "Line "${LINENO}": In function "${FUNCNAME}: The ratio is ${ratio}
	local decision=$(awk 'BEGIN{ ratio='${cant_match_rows}'/'${total_rows}'; if ( ratio > 0.1 ) {printf 1;} else printf 0}')
	if [[ ${decision} -eq 1 ]]; then echo "The table $1 order is in a mess. Exiting now" && exit 1; else echo "The table $1 row order is normal. Continue"; fi
}


function filter_on_GT_from_sample(){
    local OPTIND v g s
    while getopts v:g:s: args
    do 
        case ${args} in
            v) local vcf=$OPTARG ;;
            g) local genotype=$OPTARG ;;
            s) local sample=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    local genotype=$(echo ${genotype} | awk '{gsub("/", "\\/"); gsub(/\|/, "\\|"); print;}')
    echo The genotype to be filtered out is ${genotype}
    if [[ ${vcf} =~ \.vcf$ ]]; then
        local sample_ind=$(awk -F '\t' '$1 == "#CHROM" {for(i=1;i<=NF;i++) {if($i == "'${sample}'") printf i;}}' < ${vcf}) && \
        awk -F '\t' '$1 ~ /^#/ {print;} $1 !~ /^#/ {if($'${sample_ind}' ~ /^'${genotype}':/) next; else print;}' < ${vcf} > ${vcf}.tmp && \
        mv ${vcf}.tmp ${vcf}
    elif [[ ${vcf} =~ \.vcf\.gz$ ]]; then 
        local sample_ind=$(zcat ${vcf} | head -1000 | awk -F '\t' '$1 == "#CHROM" {for(i=1;i<=NF;i++) {if($i == "'${sample}'") printf i;}}') && \
        zcat ${vcf} | awk -F '\t' '$1 ~ /^#/ {print;} $1 !~ /^#/ {if($'${sample_ind}' ~ /^'${genotype}':/) next; else print;}' > ${vcf}.tmp && \
        mv ${vcf}.tmp ${vcf::-3} && \
        bgzip -f ${vcf::-3} && tabix -f -p vcf ${vcf}
    fi
}


function extract_on_GT_from_sample(){
    local OPTIND v g s
    while getopts v:g:s: args
    do 
        case ${args} in
            v) local vcf=$OPTARG ;;
            g) local genotype=$OPTARG ;;
            s) local sample=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    local genotype=$(echo ${genotype} | awk '{gsub("/", "\\/"); gsub(/\|/, "\\|"); print;}')
    echo The genotype to be extracted is ${genotype}
    if [[ ${vcf} =~ \.vcf$ ]]; then
        local sample_ind=$(awk -F '\t' '$1 == "#CHROM" {for(i=1;i<=NF;i++) {if($i == "'${sample}'") printf i;}}' < ${vcf}) && \
        awk -F '\t' '$1 ~ /^#/ {print;} $1 !~ /^#/ {if($'${sample_ind}' !~ /^'${genotype}':/) next; else print;}' < ${vcf} > ${vcf}.tmp && \
        mv ${vcf}.tmp ${vcf}
    elif [[ ${vcf} =~ \.vcf\.gz$ ]]; then 
        local sample_ind=$(zcat ${vcf} | head -1000 | awk -F '\t' '$1 == "#CHROM" {for(i=1;i<=NF;i++) {if($i == "'${sample}'") printf i;}}') && \
        zcat ${vcf} | awk -F '\t' '$1 ~ /^#/ {print;} $1 !~ /^#/ {if($'${sample_ind}' !~ /^'${genotype}':/) next; else print;}' > ${vcf}.tmp && \
        mv ${vcf}.tmp ${vcf::-3} && \
        bgzip -f ${vcf::-3} && tabix -f -p vcf ${vcf}
    fi
}


function get_intersect_length(){
    local bed_a=${1}
    local bed_b=${2}

    # Noted here bed_b can be a string of multiple bed files path delimited by comma
    module load BEDTools/2.27.1 2> /dev/null
    bedtools intersect -a ${bed_a} -b ${bed_b} | awk -F '\t' '{len=$3-$2; sum += len;} END{printf sum;}'
}

function create_asso_array(){
    local OPTIND t k v a
    while getopts t:k:v:a: args
    do 
        case ${args} in
            t) local table_path=$OPTARG ;;
            k) local key_col_ind=$OPTARG ;;
            v) local val_col_ind=$OPTARG ;;
            a) local __output_array=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass bam sample path." ;;
        esac
    done

    # Make sure the key value is uniq.
    local uniq_len=$(awk -F '\t' '{printf $'${key_col_ind}'"\n";}' < ${table_path} | sort - | uniq - | wc -l | awk '{printf $1}')
    local table_len=$(wc -l ${table_path} | awk '{printf $1}')

    if [[ ${uniq_len} -lt ${table_len} ]]; then
        echo "Line "${LINENO}": In function "${FUNCNAME}: $(timestamp)": The key value column is not uniquified. Make sure the key values are unique before using this func."
        exit 1
    fi

    awk -F '\t' 'BEGIN{printf "( "} {printf "[\""$'${key_col_ind}'"\"]=\""$'${val_col_ind}'"\" ";} END{printf ")";}' < ${table_path}
}

function common_annovar(){
    local vcf=$1
    local output_dir=$2
    local andir=/paedyl01/disk1/yangyxt/annovar
    local annvar=$andir/annotate_variation.pl
    local tablan=$andir/table_annovar.pl

    module load Perl 

    time perl ${tablan} ${vcf} $andir/humandb/ \
    -buildver hg19 \
    -out ${output_dir} -remove \
    -protocol refGene,genomicSuperDups,gnomad_genome,1000g2015aug_all,exac03,esp6500siv2_all,cg69,dbscsnv11,dbnsfp33a,clinvar_20200316 \
    -operation g,r,f,f,f,f,f,f,f,f -otherinfo -nastring . -vcfinput
	
    check_return_code
    module unload Perl
}


if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]
then
    while getopts f::a:: args
    do 
        case ${args} in
            f) func_name=$OPTARG ;;
            a) arguments=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    echo Executing command: ${func_name} ${arguments}

    ${func_name} ${arguments}
fi

