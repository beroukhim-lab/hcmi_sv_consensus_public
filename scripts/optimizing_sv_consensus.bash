#!/bin/bash

# Original code: pcawg6_sv_merge_master.sh <run_id> <dranger_symlink> <snowman_symlink> <brass_symlink> <delly_symlink>

# Possible new code: optimizing_sv_consensus.sh -a <run_id> -o <output_directory> -i <[bedpe_file_path_array]>
# In our case, people are providing a varying number of bedpe files (for diff projects) and they could be created from various callers
# Let's change the code so that the bedpe files are passed as an array and the last argument so this is dynamic

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

echo -e "\n## 1)\tSETUP FOR PROCESSING\n" # Print the status to the terminal

# Create a new directory for the sample of interest
mkdir -m 777 ${output_directory}
log_file=${output_directory}/${aliquot_id}.log

#Make a log file
# touch $aliquot_id.log
echo "Log file for running SV consensus on $aliquot_id" >> $log_file #Add the sample name to the log file
echo $(date '+%Y-%m-%d') >> $log_file #Add the date to the log file

#Echo the information and add to the log file
echo "Files will be output to $output_directory" | tee -a $log_file #Add the output directory to the log file

#Print all of the input paths and add to the log file
counter1=1
for file in "${bedpe_files[@]}"; do
    echo "Bedpe file number $counter1 in the array is $file" | tee -a $log_file 
    let counter1++
done

#Print all of the caller names and add to the log file
counter2=1
for caller in "${caller_name_array[@]}"; do
    echo "Caller number $counter2 is named $caller" | tee -a $log_file
    let counter2++
done

if [ $counter1 -ne $counter2 ]; then
    echo "The number of components in -i and -n parameters must be equal" >&2
    exit 1
fi


# TODO: Figure out relevance of the below variables/files and edit to not hardcode later (we can create paths in github repo)
VARIANT_DIR="${output_directory}/sv_call_concordance" #This appears to be just an output directory that will be created later?
BLACKLIST_DIR="data/blacklist_files/"
CCDSGENE="${BLACKLIST_DIR}/ccdsGene.bed.gz" #A path to a blacklist file (to be used later in this script?)
PCAWG_QC="data/PCAWG-QC_Summary-of-Measures.tsv" #A data file (no clue what it is for)
pcawg_dataset="data/pcawg_release_mar2016.tsv" #A PCAWG dataset (should be in docker container? no idea what it is for)

# TODO: Are there blacklisted regions for our projects? We should figure out how these regions were determined/how to handle accordingly
# Emailed Seongmin and he suggested this may help remove false positives but A) groups may have internally done this B) not sure if this is commonly done anymore -- should we keep??
# back in black...
#PCAWG6_DATA_DIR="${ETC_PCAWG6_DIR}/data" #This is the base PCAWG data directory
BLACKLIST_BEDPE="${BLACKLIST_DIR}/pcawg6_blacklist.slop.bedpe.gz" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_BED="${BLACKLIST_DIR}/pcawg6_blacklist.slop.bed.gz" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_TE_BEDPE="${BLACKLIST_DIR}/pcawg6_blacklist_TE_pseudogene.bedpe.gz" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_TE_INS="${BLACKLIST_DIR}/pcawg6_blacklist_TE_pseudogene_insertion.txt.gz" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_ALIQUOT_ID="${BLACKLIST_DIR}/blacklist_aliquot_id.txt" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_FOLDBACK="${BLACKLIST_DIR}/pcawg6_blacklist_foldback_artefacts.slop.bedpe.gz" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_ALIQUOT_DIR="${output_directory}/blacklisted" #A path to a directory that will be created later?
#Ah, yes. Such premonition. If the blacklist directory doesn't exist, then make it.
if [[ ! -d ${BLACKLIST_ALIQUOT_DIR} ]]
then
mkdir $BLACKLIST_ALIQUOT_DIR
fi

# There used to be code here to convert the inputs from vcf -> bedpe. We already have everything in bedpe format, so it was deleted

#####################################################################
# Merge bedpe files into one large file (aka 'master file')
#####################################################################

echo -e "\n## 2)\tCOMBINE BEDPE SV CALLS\n" # Print the status to the terminal

# Variable definition
# It looks like the bedpe_header they were using versus ours will be different -_- (fingers crossed this works)
# bedpe_header="chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\tsvclass\tsvmethod"
# bedpe_header="chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tqual\tstrand1\tstrand2\tsv_class\tsv_type\tscore\tsupp_reads\tscna\tcenter\tread_id" # Header 
SV_ORIG=${output_directory}/origMasterSvList.bed # Master file path
SV_MASTER=${output_directory}/sortedMasterSvList.bed # Master file path
umask 002
touch ${SV_MASTER}

# Since we have a variable number of bedpe files, let's just append them all at once
i=1
for file in "${bedpe_files[@]}"; do
    echo $file
    # sed 's/chr//g' $file > $file
    if [ $i -ne 1 ]
    then
        tail -n +2 $file >> $SV_ORIG
        let i++
    else
        cat $file >> $SV_ORIG
        let i++
    fi
done

sed '2,$s/chr//g' $SV_ORIG > $SV_MASTER

# Sort, save unique lines, remove header (which contains chrom1 string), and save updated file
sort -k1,1V $SV_MASTER | uniq | grep -v chrom1 >> $SV_MASTER

# Update log file
echo -e "\n$SV_MASTER" | tee -a $log_file

# exit 0

#####################################################################
# Pair SVs and find overlaps between bedpe files
# Dependency: pybedtools pairToPair
#####################################################################

# Variable definition
SLOP = 400 # Number of base pairs of slop for combining SV breakpoint ends (e.g. callers may call the same sv breakpoint, but be a few bp off because of the algorithm)
PAIR2PAIR = ${output_directory}/pair2pair_SV_merging.bed # Output file path
echo -n "" > $PAIR2PAIR # Create empty file

# Since we have a variable number of bedpe files, let's just append them all at once
# Flag information: 
## -slop = amount of extra space
## -rdn = require hits to have diff names
## -a and -b = input bedpe files; -a = master, -b = institution
for file in "${bedpe_files[@]}"; do
    pairToPair -slop $SLOP -rdn -a <(cut -f -19 ${SV_MASTER}) -b <(cut -f -19 $file | awk -v aliquot_id=$aliquot_id)   '{print $0"\t"aliquot_id}' ) >> $PAIR2PAIR
done

#####################################################################
# Make SV overlap for each SV in pair2pair
# Prepare special data frame format to load into graph algorithm
#####################################################################

inBEDPE=$SAMPLE_DIR/${aliquot_id}_SV_overlap.txt
python ${CODE_DIR}/pcawg_merge_reorder_pairs.py $PAIR2PAIR $aliquot_id > ${inBEDPE}

# I'm finding that commenting the above python file is hard to do without understanding the input. 
# Let's run everything before this step to test and make sure logic is working but also to be able to comment on the script












# #####################################################################
# #############	  Merge and get cliques of SVs 	    #################
# #####################################################################
# echo -e "\n## 3)\tMERGE SVs\n"

# outVCF=$ANALYSIS_DIR/${aliquot_id}.${call_stamp}.${today}.somatic.sv.vcf

# outSTAT=${outVCF/.vcf/.stat}

# cd ${SAMPLE_DIR}

# python ${CODE_DIR}/pcawg6_sv_merge_graph.py -e ${inBEDPE} -o ${outVCF} -s ${outSTAT} -a $SANGER_ANNO_VCF -b $DELLY_ANNO_VCF -c $DRANGER_ANNO_VCF -d $SNOWMAN_ANNO_VCF | tee -a $log










# #####################################################################
# ##############    RNA/GERMLINE BLACKLIST  REMOVAL   #################
# #####################################################################
# echo -e "\n\tREMOVE BLACKLIST\n"

# BLACKLIST_SV_ID=${aliquot_id}.blacklist_svid.txt
# echo -en "" > ${BLACKLIST_SV_ID}
# ## remove SVs with both breaks overlapping blacklist regions - germline artefacts and TEs
# pairToPair -is -type both -a ${outVCF/.vcf/.bedpe} -b $BLACKLIST_BEDPE  > ${outVCF/.vcf/_blacklisted_tmp.bedpe}
# if [[ -s ${outVCF/.vcf/_blacklisted_tmp.bedpe}  ]];then
# cut -f 7,19 ${outVCF/.vcf/_blacklisted_tmp.bedpe} | sed 's/|.*//' | sort -u   >> ${BLACKLIST_SV_ID}
# fi

# ## remove SVs with both breaks overlapping blacklist regions and read orientation match - foldback inversions
# pairToPair -type both -a ${outVCF/.vcf/.bedpe} -b $BLACKLIST_FOLDBACK  > ${outVCF/.vcf/_blacklisted_tmp.bedpe}
# if [[ -s ${outVCF/.vcf/_blacklisted_tmp.bed}  ]];then
# cut -f 7,18 ${outVCF/.vcf/_blacklisted_tmp.bed} | sed 's/|.*//' | sort -u   >> ${BLACKLIST_SV_ID}
# fi

# ## remove SVs with one break overalpping blacklist bed file regions
# pairToBed -a ${outVCF/.vcf/.bedpe}  -b ${BLACKLIST_BED}  > ${outVCF/.vcf/_blacklisted_tmp.bed}
# if [[ -s ${outVCF/.vcf/_blacklisted_tmp.bed}  ]];then
# cut -f 7,18 ${outVCF/.vcf/_blacklisted_tmp.bed} | sed 's/|.*//' | sort -u   >> ${BLACKLIST_SV_ID}
# fi
# ## remove SVs with either break in pseudogene exon-exon region
# pairToPair -is -type either -a ${outVCF/.vcf/.bedpe}  -b ${BLACKLIST_TE_BEDPE} > ${outVCF/.vcf/_blacklisted_TE_pseudogene_tmp.bed}
# if [[ -s ${outVCF/.vcf/_blacklisted_TE_pseudogene_tmp.bed}  ]];then
# cut -f 7,19 ${outVCF/.vcf/_blacklisted_TE_pseudogene_tmp.bed} | sed 's/|.*//' | sort -u   >> ${BLACKLIST_SV_ID}
# fi
# echo -e "\n## 4)\tRNA CONTAMINATION FILTER\n"

# ## RNA
# bash ${CODE_DIR}/remove_splicing_type_svs_yl.sh ${outVCF/.vcf/.bedpe} ${CCDSGENE} ${BLACKLIST_TE_INS}  > splice_type_id_blacklist.txt
# if [[ -s splice_type_id_blacklist.txt  ]];then
# awk '{print $1"\tcdna_contamination"}' splice_type_id_blacklist.txt >> ${BLACKLIST_SV_ID}
# fi
# rm splice_type_id_blacklist.txt
# ## FILTER and ANNOTATE
# python ${CODE_DIR}/pcawg6_sv_merge_annotate.py ${outVCF} ${BLACKLIST_SV_ID}  ${aliquot_id} ${PCAWG_QC} | tee -a $log

# echo -e "\n\ngenerate BEDPE fil"
# python ${CODE_DIR}/pcawg6Vcf2Bedpe.py ${outVCF}
# python ${CODE_DIR}/pcawg6Vcf2Bedpe.py ${outVCF/.vcf/_full.vcf}

# bgzip -f $outVCF

# if [[ -f ${outVCF/.vcf/_tmp.vcf}  ]];then
# bgzip -f ${outVCF/.vcf/_tmp.vcf}
# fi

# if [[ -f ${outVCF/.vcf/_full.vcf}  ]];then
# bgzip -f ${outVCF/.vcf/_full.vcf}
# fi

# echo -e "\n\nnumber of merged SVs" | tee -a  $log
# zgrep -vP "^#"  $outVCF | awk '$3~/_1/' | wc -l | tee -a  $log
# tabix -f -p vcf ${outVCF}.gz
# tabix -f -p vcf ${outVCF/.vcf/_full.vcf}.gz
# cd -

# if [[ $(stat -c %s ${outVCF}.gz) -gt 100  ]]; then echo -e "\n\n>>\t" $aliquot_id vcf file successful"\n>>\t${outVCF}.gz\n"; else echo -e ">>\t"  $aliquot_id generation error;fi  | tee -a $log

# echo -e "output clique count: ${aliquot_id}_SV_overlap.clique.txt"
# echo -e "output clique stat: ${aliquot_id}.stats"

# mv -f ${ANALYSIS_DIR}/${aliquot_id}*tmp*  $SAMPLE_DIR
# mv -f ${ANALYSIS_DIR}/${aliquot_id}*txt  $SAMPLE_DIR
# mv -f ${ANALYSIS_DIR}/${aliquot_id}*no_exon_exon_artefacts* $SAMPLE_DIR

# if grep ${aliquot_id} ${BLACKLIST_ALIQUOT_ID} > /dev/null;then
# echo -e "\n>>> sample blacklisted. Will be excluded from the set <<<\n"
# mv ${ANALYSIS_DIR}/${aliquot_id}* ${BLACKLIST_ALIQUOT_DIR}
# fi

# echo -e "complete, log file\n\n$log"














# #####################################################################
# #####	Rearrange clique calls and make binary tree for plotting ####
# #####################################################################

# clique=$SAMPLE_DIR/${aliquot_id}_SV_overlap.clique.txt
# cliqueCenter=${clique/.txt/.center.txt}
# echo -ne "" > ${cliqueCenter};
# for i in $(cut -f1 ${clique});
# do

# if echo $i |grep -P '^BRASS'  > /dev/null;
# then center=brass;

# elif echo $i | grep -P '^dRANGER'  > /dev/null;
# then center=dranger;
# elif echo $i | grep -P '^SNOWMAN'  > /dev/null;
# then center=snowman;

# elif echo $i  | grep -P '^DELLY'  > /dev/null;
# then center=delly;
# fi; awk -v id=$i '$1==id' ${clique} | sed "s/$i/$center/"  >> ${cliqueCenter};
# done

# ### COUNT TOTAL SVs

# dranger_count=$(grep -v "chrom1" $DRANGER_BEDPE  | wc -l)
# delly_count=$(grep -v "chrom1" $DELLY_BEDPE |   wc -l)
# brass_count=$(grep -v "chrom1" $BRASS_BEDPE |  wc -l)
# snowman_count=$(grep -v "chrom1" $SNOWMAN_BEDPE |  wc -l)

# ## Plotting

# binaryOut=$SAMPLE_DIR/${aliquot_id}.binaryTree.txt
# echo $binaryOut
# cat <(echo -e "brass\ndelly\ndranger\nsnowman") <(sort -k2,2n ${cliqueCenter}) |\
# awk  'BEGIN{OFS="\t";print "clique\tbrass\tdelly\tdranger\tsnowman" ;sampleCount=0; clusterId=0; }{ if(NF==1) { sampleCount++; samples[$1]=sampleCount; sampleOccur[sampleCount]=0; } else { if ($2==clusterId) { sampleOccur[samples[$1]]=1; } else { printf clusterId; for (i=1; i<=sampleCount; i++) { printf "\t"sampleOccur[i]; sampleOccur[i]=0; } print ""; clusterId=$2; sampleOccur[samples[$1]]=1;}  }} END { printf clusterId; for (i=1; i<=sampleCount; i++) { printf "\t"sampleOccur[i] } print ""}' | awk '$1!=0' > $binaryOut

# echo -e "count\t$brass_count\t$delly_count\t$dranger_count\t$snowman_count" >> $binaryOut
# Rscript ${CODE_DIR}/pcawg_multiple_sv_merging_heatmap.R $binaryOut

# echo -e "merged SV file:\n${outVCF}.gz\n" | tee -a $log
