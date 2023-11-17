#!/bin/bash

###Deal with all argument inputs
if [[ $# -lt 1 ]]; then
    echo "You must pass one of the following combinations of parameters (-h) or (-a, -o, -i, AND -n)" >&2
    exit 1
fi

help() {
   # Display Help
   echo "This script is designed to consolidate consensus SV calls across bedpe files"
   echo
   echo "Usage: optimizing_sv_consensus.sh [-a aliquot_id] [-o output_directory] [-i (bedpe_file_path_array)] [-n (institution_array)]"
   echo "options:"
   echo "a     Aliquot ID"
   echo "o     Output directory (should not already exist)"
   echo "i     Input a string of paths to bedpe files, with spaces between each path. Can accept any number of bedpe files"
   echo "n     Input an string of institutions that created corresponding bedpe files, with spaces between each institution. Length must be same as array in -i"
   echo "h     Print this help menu"
   echo
   exit 0;
}

usage() {
    echo "Usage: optimizing_sv_consensus.sh [-a aliquot_id] [-o output_directory] [-i 'bedpe_file_path_array'] [-n 'institution_array']" 1>&2;
    exit 1;
}

while getopts ":a:o:i:n:h" var; do  
    echo "$var"
    case "$var" in
        a) 
            aliquot_id=${OPTARG};; 
        o) 
            output_directory=${OPTARG};; 
        i) 
            input_string=${OPTARG} 
            IFS=' ' read -ra bedpe_files <<< "$input_string";;
        n) 
            caller_string=${OPTARG} 
            IFS=' ' read -ra caller_name_array <<< "$caller_string";;
        h)
            help;;
        *) 
            usage
            exit 1;;
    esac
done



# Create a new directory for the sample of interest and tmp subdirectory
mkdir -m 777 ${output_directory}
mkdir -m 777 ${output_directory}/tmp
tmp_directory="${output_directory}/tmp"
log_file=${output_directory}/${aliquot_id}.log


#Add more information to the log file
echo -e "\n## 1)\tSETUP FOR PROCESSING\n"
echo "SETUP FOR PROCESSING" >> $log_file
echo "Log file for running SV consensus on $aliquot_id" >> $log_file #Add the sample name to the log file
echo $(date '+%Y-%m-%d') >> $log_file #Add the date to the log file
echo "Files will be output to $output_directory" | tee -a $log_file #Add the output directory to the log file


# Print all of the input paths and add to the log file
# Error checking: remove any ^M instances in file to prevent issues downstream from Windows vs Linux use
counter1=0
for file in "${bedpe_files[@]}"; do
    center=${caller_name_array[counter1]} 
    echo "Bedpe file number $counter1 in the array is $file" | tee -a $log_file 
    sed -e "s/\r//g" $file > ${tmp_directory}/$aliquot_id.$center.m_removal
    sed '2,$s/chr//g' ${tmp_directory}/$aliquot_id.$center.m_removal > ${tmp_directory}/$aliquot_id.$center.chr_removal
    let counter1++
done


#Error checking: check the header in all of the bedpe files
bedpe_header_array=("chrom1" "start1" "end1" "chrom2" "start2" "end2" "sv_id" "tumreads" "strand1" "strand2" "svclass" "svmethod")
headercounter=0
for file in "${bedpe_files[@]}"; do
    IFS=$'\t' read -r -a bedpe_file_header <<< "$(head -n 1 "$file")" # Get the bedpe header and store as an array
    
    # Check if there are the expected number of columns
    expected_num_columns=${#bedpe_header_array[@]} # The expected number of columns in the bedpe header (length of bedpe_header_array)
    observed_num_columns=${#bedpe_file_header[@]} # The observed number of columns in the bedpe header
    if [ $expected_num_columns != $observed_num_columns ]; then
        echo "Number of columns in bedpe -- $observed_num_columns is not equal to the expected number of columns -- $expected_num_columns. Check column headers." >&2
        echo "The file that triggered this error is ${caller_name_array["$headercounter"]}" >&2
        exit 1
    fi

    # Check if each column header is as expected
    for ((col=0; col<$expected_num_columns; col++)); do
        expected_value=${bedpe_header_array[$col]}
        observed_value=${bedpe_file_header[$col]}

        if [ "$expected_value" != "$observed_value" ]; then
            echo "Warning: Expected header value of $expected_value but instead found $observed_value in ${caller_name_array["$headercounter"]}" >&2
        fi
    done

    let headercounter++
done



# Print all of the caller names and add to the log file
counter2=0
for caller in "${caller_name_array[@]}"; do
    echo "Caller number $counter2 is named $caller" | tee -a $log_file
    let counter2++
done

#Print an error if there are an unequal number of callers and caller names specified in the input arguments
if [ $counter1 -ne $counter2 ]; then
    echo "The number of components in -i and -n parameters must be equal" >&2
    exit 1
fi


#####################################################################
# Merge bedpe files into one large file (aka 'master file')
#####################################################################


#Print the status
echo -e "\n## 2)\tCOMBINE BEDPE SV CALLS\n"
echo "COMBINE BEDPE SV CALLS" >> $log_file


#Add a "center" column to the bedpe. This will be necessary for merge_reorder_pairs.py because this script counts the number of times each center is seen.
i=0
for file in "${bedpe_files[@]}"; do
    center=${caller_name_array[i]} 
    in_temp_file=${tmp_directory}/$aliquot_id.$center.chr_removal
    temp_file=${tmp_directory}/$aliquot_id.$center.tmp
    tmp_bedpe_files+=($temp_file)
    awk -v val="$center" 'BEGIN {OFS = "\t"} FNR==1{a="center"} FNR>1{a=val} {print $0"\t"a}' $in_temp_file | grep -v chrom1 > $temp_file
    let i++
done


#Variable definition
SV_ORIG=${output_directory}/origMasterSvList.bedpe # Master file path
SV_MASTER=${output_directory}/sortedMasterSvList.bedpe # Master file path
umask 002
touch ${SV_MASTER}


# Since we have a variable number of bedpe files, let's just append them all at once
i=1
for file in "${tmp_bedpe_files[@]}"; do
    echo $file

    if [ $i -ne 1 ]
    then
        tail -n +2 $file >> $SV_ORIG
        let i++
    else
        cat $file >> $SV_ORIG
        let i++
    fi
done


# Sort, save unique lines, remove header (which contains chrom1 string), and save updated file
sort -k1,1V $SV_ORIG | uniq | grep -v chrom1 >> $SV_MASTER


# Update log file
echo -e "\n$SV_MASTER" | tee -a $log_file


#####################################################################
# Pair SVs and find overlaps between bedpe files
# Dependency: bedtools pairToPair
#####################################################################


# Variable definition
PAIR2PAIR=${output_directory}/pair2pair_SV_merging.bed # Output file path
echo -n "" > $PAIR2PAIR # Create empty file


###Run pairToPair
# Flag information: 
## -slop = amount of extra space
## -rdn = require hits to have diff names
## -a and -b = input bedpe files; -a = master, -b = institution
SLOP=400 #todo, change so that we specify slop as an input argument
i=0
for file in "${tmp_bedpe_files[@]}"; do
    pairToPair -slop $SLOP -rdn -a ${SV_MASTER} -b <(awk -v aliquot_id=$aliquot_id '{print $0"\t"aliquot_id}' $file)  >> $PAIR2PAIR
    let i++
done



#Make SV overlap for each SV in pair2pair
#Prepare special data frame format to load into graph algorithm
inBEDPE=${output_directory}/${aliquot_id}_SV_overlap.txt
python3 scripts/pcawg_merge_reorder_pairs.py $PAIR2PAIR $aliquot_id > ${inBEDPE}



###Merge and get cliques of SVs
#Add information to the log file
echo -e "\n## 3)\tMERGE SVs\n"
echo "MERGE SVs" >> $log_file

#Define the outputs
outBEDPE=${output_directory}/${aliquot_id}.somatic.sv.bedpe
outSTAT=${outBEDPE/.bedpe/.stat}


#Merge the SVs using a similar approach to the one used in PCAWG
bedpe_python_input_array=$(IFS=" " ; echo "${tmp_bedpe_files[@]}") # Convert the array of bedpe files to a space-separated string to make it easier to use as a python argument
python3 scripts/pcawg6_sv_merge_graph.py -e ${inBEDPE} -z "$bedpe_python_input_array" -o ${outBEDPE} -s ${outSTAT} | tee -a $log


#Add information to the log file
echo "DONE" >> $log_file
exit 0