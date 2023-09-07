
#The first argument of this shell script is to input the name of the sample 
#aliquot_id=$1
aliquot_id="test"

#Define the sv bedpe files
broad_bedpe=

#Save today's date as a variable
today=$(date +"%y%m%d")

#The "etc" directory? No idea what this is for
#Potentially some code is stored here? So we can probably just delete this 
#ETC_PCAWG6_DIR=/etc/pcawg6_merge_sv
#source <(grep = ${ETC_PCAWG6_DIR}/merge.ini)

#Specifc the data directory. We can maybe just delete this part, depending on how we structure things?
DATA_DIR="/xchip/beroukhimlab/Sean/software/pcawg_sv_merge/data" #Hard coded, need to change in the future? Let's see how lazy I am.

#If the analysis directory exists, then store it as a variable.
#If the analysis directory does not exist, then make it
#We can also probably just delete this?
if [[ -z ${ANALYSIS_DIR} ]]
then
ANALYSIS_DIR="/xchip/beroukhimlab/Sean/hcmi/sv_calls/consensus_optimization/output" #Also a temp directory while I optimize things. Will fix once everything is working.
fi
if [[ ! -d ${ANALYSIS_DIR}  ]]
then
mkdir -p $ANALYSIS_DIR
fi


#CODE_DIR="/opt/scripts" #The directory where all of the code is stored
VARIANT_DIR="${ANALYSIS_DIR}/sv_call_concordance" #This appears to be just an output directory that will be created later?
#DRANGER_DIR=${DIR}/${DRANGER} #Where all of the input SV VCFs are stored, can probably just delete?
#SNOWMAN_DIR=${DIR}/${SNOWMAN} #Where all of the input SV VCFs are stored, can probably just delete?
#BRASS_DIR=${DIR}/${BRASS} #Where all of the input SV VCFs are stored, can probably just delete?
#DELLY_DIR=${DIR}/${DELLY} #Where all of the input SV VCFs are stored, can probably just delete?
# back in black...
#PCAWG6_DATA_DIR="${ETC_PCAWG6_DIR}/data" #This is the base PCAWG data directory
BLACKLIST_DIR="/xchip/beroukhimlab/Sean/hcmi/sv_calls/consensus_optimization/blacklist_files" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_BEDPE="${BLACKLIST_DIR}/pcawg6_blacklist.slop.bedpe.gz" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_BED="${BLACKLIST_DIR}/pcawg6_blacklist.slop.bed.gz" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_TE_BEDPE="${BLACKLIST_DIR}/pcawg6_blacklist_TE_pseudogene.bedpe.gz" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_TE_INS="${BLACKLIST_DIR}/pcawg6_blacklist_TE_pseudogene_insertion.txt.gz" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_ALIQUOT_ID="${BLACKLIST_DIR}/blacklist_aliquot_id.txt" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_FOLDBACK="${BLACKLIST_DIR}/pcawg6_blacklist_foldback_artefacts.slop.bedpe.gz" #A path to a blacklist file (to be used later in this script?)
BLACKLIST_ALIQUOT_DIR=${ANALYSIS_DIR}/blacklisted #A path to a directory that will be created later?

#Ah, yes. Such premonition. If the blacklist directory doesn't exist, then make it.
if [[ ! -d ${BLACKLIST_ALIQUOT_DIR} ]]
then
mkdir $BLACKLIST_ALIQUOT_DIR
fi


CCDSGENE="${BLACKLIST_DIR}/ccdsGene.bed.gz" #A path to a blacklist file (to be used later in this script?)
PCAWG_QC="/xchip/beroukhimlab/Sean/software/pcawg_sv_merge/data/PCAWG-QC_Summary-of-Measures.tsv" #A data file (no clue what it is for)
pcawg_dataset="/xchip/beroukhimlab/Sean/software/pcawg_sv_merge/data/pcawg_release_mar2016.tsv" #A PCAWG dataset (should be in docker container? no idea what it is for)
log=${ANALYSIS_DIR}/${aliquot_id}.${call_stamp}.${today}.log #Create the log file and store it as a variable

#Add project id information to the log file
#Maybe I will get around to fixing this? Not really necessary for the code to run.
echo -en "" > ${log}
PROJECT=$(grep $aliquot_id $pcawg_dataset | cut -f 2)
if [[  -z $PROJECT  ]];then
PROJECT="UNKOWN"
fi

#Create a new directory for the sample of interest
SAMPLE_DIR=${VARIANT_DIR}/${PROJECT}/${aliquot_id}
mkdir -p $SAMPLE_DIR
cd $SAMPLE_DIR

#Print things and add it to the log file
echo -e "########################################\n\n"$PROJECT"\naliquot id=" $aliquot_id "\nsv merge set" ${call_stamp} "\nmerged on" $(date)  | tee  $log

###There used to be code here to convert the inputs from vcf -> bedpe. We already have everything in bedpe format, so I deleted it.








#####################################################################
#######################      make the master file   #################
#####################################################################


###Make the 'master file'. What is a 'master file'?
#I think the master file is just all of the merged bedpe calls
echo -e "\n## 2)\tCOMBINE BEDPE SV CALLS\n" #Print the status to the terminal
bedpe_header="chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tqual\tstrand1\tstrand2\tsv_class\tsv_type\tscore\tsupp_reads\tscna\tcenter\tread_id" #Define the bedpe header 

SV_MASTER=${SAMPLE_DIR}/masterSvList.bed #Define a file to store the sv master list 
echo -e $bedpe_header > $SV_MASTER #Add the header that we defined into the master file
cat ${BRASS_BEDPE} ${DRANGER_BEDPE} ${SNOWMAN_BEDPE} ${SMUFIN_BEDPE} ${DELLY_BEDPE} | sort -k1,1V | uniq | grep -v chrom1 >> $SV_MASTER #Combine all of the bedpe files, sort them, get unique svs, remove the header (which contains the chrom1 string) then add to the master file
echo -e "\n$SV_MASTER" | tee -a $log #Add the status to the output log









#####################################################################
#######################      PAIR the SVs   #########################
#####################################################################

###Next, write some code to pair the SVs.

SLOP=400 #Define a number of base pairs of slop for combining SV breakpoint ends (e.g. callers may call the same sv breakpoint, but be a few bp off because of the algorithm)
PAIR2PAIR=${SAMPLE_DIR}/pair2pair_SV_merging.bed #Create an output file for storing all of hte pairs
echo -n "" > $PAIR2PAIR #Add a new line character to the pair2pair file? Why?

#Loop through all of the bed files
#pairToPair is part of bedtools, so I will need to download that!!!!!!!!!!!!!!!!!!!!!!!!
#-a and -b arguemnts specify the inputs. In this case, a is the master file and b is the specific bed file of interest. 
#Save the output in the pair2pair file
for bed in ${BRASS_BEDPE} ${DRANGER_BEDPE} ${SNOWMAN_BEDPE} ${SMUFIN_BEDPE} ${DELLY_BEDPE}
do
pairToPair  -slop $SLOP -rdn -a <(cut -f -19 ${SV_MASTER}) -b <(cut -f -19  $bed | awk -v aliquot_id=$aliquot_id   '{print $0"\t"aliquot_id}' ) >> $PAIR2PAIR
done









#####################################################################
############ Make SV overlap for each SV in pair2pair ###############
#####################################################################

## prepare special data frame format to load into graph algorithm

inBEDPE=$SAMPLE_DIR/${aliquot_id}_SV_overlap.txt
python ${CODE_DIR}/pcawg_merge_reorder_pairs.py $PAIR2PAIR $aliquot_id > ${inBEDPE}












#####################################################################
#############	  Merge and get cliques of SVs 	    #################
#####################################################################
echo -e "\n## 3)\tMERGE SVs\n"

outVCF=$ANALYSIS_DIR/${aliquot_id}.${call_stamp}.${today}.somatic.sv.vcf

outSTAT=${outVCF/.vcf/.stat}

cd ${SAMPLE_DIR}

python ${CODE_DIR}/pcawg6_sv_merge_graph.py -e ${inBEDPE} -o ${outVCF} -s ${outSTAT} -a $SANGER_ANNO_VCF -b $DELLY_ANNO_VCF -c $DRANGER_ANNO_VCF -d $SNOWMAN_ANNO_VCF | tee -a $log










#####################################################################
##############    RNA/GERMLINE BLACKLIST  REMOVAL   #################
#####################################################################
echo -e "\n\tREMOVE BLACKLIST\n"

BLACKLIST_SV_ID=${aliquot_id}.blacklist_svid.txt
echo -en "" > ${BLACKLIST_SV_ID}
## remove SVs with both breaks overlapping blacklist regions - germline artefacts and TEs
pairToPair -is -type both -a ${outVCF/.vcf/.bedpe} -b $BLACKLIST_BEDPE  > ${outVCF/.vcf/_blacklisted_tmp.bedpe}
if [[ -s ${outVCF/.vcf/_blacklisted_tmp.bedpe}  ]];then
cut -f 7,19 ${outVCF/.vcf/_blacklisted_tmp.bedpe} | sed 's/|.*//' | sort -u   >> ${BLACKLIST_SV_ID}
fi

## remove SVs with both breaks overlapping blacklist regions and read orientation match - foldback inversions
pairToPair -type both -a ${outVCF/.vcf/.bedpe} -b $BLACKLIST_FOLDBACK  > ${outVCF/.vcf/_blacklisted_tmp.bedpe}
if [[ -s ${outVCF/.vcf/_blacklisted_tmp.bed}  ]];then
cut -f 7,18 ${outVCF/.vcf/_blacklisted_tmp.bed} | sed 's/|.*//' | sort -u   >> ${BLACKLIST_SV_ID}
fi

## remove SVs with one break overalpping blacklist bed file regions
pairToBed -a ${outVCF/.vcf/.bedpe}  -b ${BLACKLIST_BED}  > ${outVCF/.vcf/_blacklisted_tmp.bed}
if [[ -s ${outVCF/.vcf/_blacklisted_tmp.bed}  ]];then
cut -f 7,18 ${outVCF/.vcf/_blacklisted_tmp.bed} | sed 's/|.*//' | sort -u   >> ${BLACKLIST_SV_ID}
fi
## remove SVs with either break in pseudogene exon-exon region
pairToPair -is -type either -a ${outVCF/.vcf/.bedpe}  -b ${BLACKLIST_TE_BEDPE} > ${outVCF/.vcf/_blacklisted_TE_pseudogene_tmp.bed}
if [[ -s ${outVCF/.vcf/_blacklisted_TE_pseudogene_tmp.bed}  ]];then
cut -f 7,19 ${outVCF/.vcf/_blacklisted_TE_pseudogene_tmp.bed} | sed 's/|.*//' | sort -u   >> ${BLACKLIST_SV_ID}
fi
echo -e "\n## 4)\tRNA CONTAMINATION FILTER\n"

## RNA
bash ${CODE_DIR}/remove_splicing_type_svs_yl.sh ${outVCF/.vcf/.bedpe} ${CCDSGENE} ${BLACKLIST_TE_INS}  > splice_type_id_blacklist.txt
if [[ -s splice_type_id_blacklist.txt  ]];then
awk '{print $1"\tcdna_contamination"}' splice_type_id_blacklist.txt >> ${BLACKLIST_SV_ID}
fi
rm splice_type_id_blacklist.txt
## FILTER and ANNOTATE
python ${CODE_DIR}/pcawg6_sv_merge_annotate.py ${outVCF} ${BLACKLIST_SV_ID}  ${aliquot_id} ${PCAWG_QC} | tee -a $log

echo -e "\n\ngenerate BEDPE fil"
python ${CODE_DIR}/pcawg6Vcf2Bedpe.py ${outVCF}
python ${CODE_DIR}/pcawg6Vcf2Bedpe.py ${outVCF/.vcf/_full.vcf}

bgzip -f $outVCF

if [[ -f ${outVCF/.vcf/_tmp.vcf}  ]];then
bgzip -f ${outVCF/.vcf/_tmp.vcf}
fi

if [[ -f ${outVCF/.vcf/_full.vcf}  ]];then
bgzip -f ${outVCF/.vcf/_full.vcf}
fi

echo -e "\n\nnumber of merged SVs" | tee -a  $log
zgrep -vP "^#"  $outVCF | awk '$3~/_1/' | wc -l | tee -a  $log
tabix -f -p vcf ${outVCF}.gz
tabix -f -p vcf ${outVCF/.vcf/_full.vcf}.gz
cd -

if [[ $(stat -c %s ${outVCF}.gz) -gt 100  ]]; then echo -e "\n\n>>\t" $aliquot_id vcf file successful"\n>>\t${outVCF}.gz\n"; else echo -e ">>\t"  $aliquot_id generation error;fi  | tee -a $log

echo -e "output clique count: ${aliquot_id}_SV_overlap.clique.txt"
echo -e "output clique stat: ${aliquot_id}.stats"

mv -f ${ANALYSIS_DIR}/${aliquot_id}*tmp*  $SAMPLE_DIR
mv -f ${ANALYSIS_DIR}/${aliquot_id}*txt  $SAMPLE_DIR
mv -f ${ANALYSIS_DIR}/${aliquot_id}*no_exon_exon_artefacts* $SAMPLE_DIR

if grep ${aliquot_id} ${BLACKLIST_ALIQUOT_ID} > /dev/null;then
echo -e "\n>>> sample blacklisted. Will be excluded from the set <<<\n"
mv ${ANALYSIS_DIR}/${aliquot_id}* ${BLACKLIST_ALIQUOT_DIR}
fi

echo -e "complete, log file\n\n$log"














#####################################################################
#####	Rearrange clique calls and make binary tree for plotting ####
#####################################################################

clique=$SAMPLE_DIR/${aliquot_id}_SV_overlap.clique.txt
cliqueCenter=${clique/.txt/.center.txt}
echo -ne "" > ${cliqueCenter};
for i in $(cut -f1 ${clique});
do

if echo $i |grep -P '^BRASS'  > /dev/null;
then center=brass;

elif echo $i | grep -P '^dRANGER'  > /dev/null;
then center=dranger;
elif echo $i | grep -P '^SNOWMAN'  > /dev/null;
then center=snowman;

elif echo $i  | grep -P '^DELLY'  > /dev/null;
then center=delly;
fi; awk -v id=$i '$1==id' ${clique} | sed "s/$i/$center/"  >> ${cliqueCenter};
done

### COUNT TOTAL SVs

dranger_count=$(grep -v "chrom1" $DRANGER_BEDPE  | wc -l)
delly_count=$(grep -v "chrom1" $DELLY_BEDPE |   wc -l)
brass_count=$(grep -v "chrom1" $BRASS_BEDPE |  wc -l)
snowman_count=$(grep -v "chrom1" $SNOWMAN_BEDPE |  wc -l)

## Plotting

binaryOut=$SAMPLE_DIR/${aliquot_id}.binaryTree.txt
echo $binaryOut
cat <(echo -e "brass\ndelly\ndranger\nsnowman") <(sort -k2,2n ${cliqueCenter}) |\
awk  'BEGIN{OFS="\t";print "clique\tbrass\tdelly\tdranger\tsnowman" ;sampleCount=0; clusterId=0; }{ if(NF==1) { sampleCount++; samples[$1]=sampleCount; sampleOccur[sampleCount]=0; } else { if ($2==clusterId) { sampleOccur[samples[$1]]=1; } else { printf clusterId; for (i=1; i<=sampleCount; i++) { printf "\t"sampleOccur[i]; sampleOccur[i]=0; } print ""; clusterId=$2; sampleOccur[samples[$1]]=1;}  }} END { printf clusterId; for (i=1; i<=sampleCount; i++) { printf "\t"sampleOccur[i] } print ""}' | awk '$1!=0' > $binaryOut

echo -e "count\t$brass_count\t$delly_count\t$dranger_count\t$snowman_count" >> $binaryOut
Rscript ${CODE_DIR}/pcawg_multiple_sv_merging_heatmap.R $binaryOut

echo -e "merged SV file:\n${outVCF}.gz\n" | tee -a $log
