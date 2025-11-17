#!/bin/bash

#The purpose of this code is to compute the SV consensus using our SV consensus pipeline


#All samples in data freeze
included_samples="/data/consensus_samples.txt"
num_samples=$(cat $included_samples | wc -l) 
base_directory="/path/to/base/directory"


#Iterate through all of the files
for ((i=1; i<=num_samples; i++)); do

sample_name=$(cat $included_samples | tail -n +2 | head -n $i | tail -n 1 | awk '{print $1}')


###Build the arrays
file_path_array=() 
caller_name_array=() 

#Look for the files, if exists then store as a variable, if not exist then do not store as a variable
embl_files=($base_directory/embl/reformatted_${sample_name}*.bedpe)
if [ -e "${embl_files[0]}" ]; then
    embl_bedpe=$(ls "${embl_files[@]}" | head -n 1)
    file_path_array+=($embl_bedpe)
    caller_name_array+=("embl")
    unset embl_bedpe
else
    echo "No matching embl .bedpe files found for sample $sample_name"
fi

mskcc_files=($base_directory/mskcc/reformatted_${sample_name}*.bedpe)
if [ -e "${mskcc_files[0]}" ]; then
    mskcc_bedpe=$(ls "${mskcc_files[@]}" | head -n 1)
    file_path_array+=($mskcc_bedpe)
    caller_name_array+=("mskcc")
    unset mskcc_bedpe
else
    echo "No matching mskcc .bedpe files found for sample $sample_name"
fi

nygc_files=($base_directory/nygc/$sample_name*.bedpe)
if [ -e "${nygc_files[0]}" ]; then
    nygc_bedpe=$(ls "${nygc_files[@]}" | head -n 1)
    file_path_array+=($nygc_bedpe)
    caller_name_array+=("nygc")
    unset nygc_bedpe
else
    echo "No matching nygc .bedpe files found for sample $sample_name"
fi

washu_files=($base_directory/washu/reformatted_$sample_name*.bedpe)
if [ -e "${washu_files[0]}" ]; then
    washu_bedpe=$(ls "${washu_files[@]}" | head -n 1)
    file_path_array+=($washu_bedpe)
    caller_name_array+=("washu")
    unset washu_bedpe
else
    echo "No matching washu .bedpe files found for sample $sample_name"
fi

broad_files=($base_directory/broad/$sample_name*.bedpe)
if [ -e "${broad_files[0]}" ]; then
    broad_bedpe=$(ls "${broad_files[@]}" | head -n 1)
    file_path_array+=($broad_bedpe)
    caller_name_array+=("broad")
    unset broad_bedpe
else
    echo "No matching broad .bedpe files found for sample $sample_name"
fi


#Loop through all of the files and run the sv consensus pipeline :) 
/scripts/optimizing_sv_consensus.bash -a $sample_name -o "/path/to/consensus_calls"/$sample_name -i ${file_path_array[@]} -n ${caller_name_array[@]} 


#unset the caller name array
unset caller_name_array
unset file_path_array


done