#!/bin/bash
pid=$1
echo $pid
rd_script=/home/weischej/Dropbox/job_scripts/pawg_scripts/PAWG/read_depth_ratio_sv_bedpe_plusCovSigSeg_ACEseqformat.R
acedir=/icgc/pcawg/analysis/final_variant_calls/projects/
#sv_merge_set="pcawg6_merge_1.00.160601"
#sv_merge_set="pcawg6_merge_1.2.160728"
sv_merge_set="pcawg_consensus_1.5.1.160912"
sv_merge_folder="160912.sv_merge_TE_free"
mainDir="/icgc/pcawg/analysis/train2full/March_2016_dataset/plots/${sv_merge_set}"
mkdir $mainDir

mkdir -p $mainDir/$pid
cd $mainDir/$pid

cov=$(find $acedir/*-*/$pid/ACEseq/plots -iname "*all_seg.gc_corrected.txt")
#sv=$(find /icgc/pcawg/analysis/train2full/March_2016_dataset/160601.sv_merge// -maxdepth 1 -iname "${pid}.${sv_merge_set}.somatic.sv.bedpe")
#sv=$(find /icgc/pcawg/analysis/train2full/March_2016_dataset/${sv_merge_folder}/ -maxdepth 1 -iname "${pid}.${sv_merge_set}.somatic.sv.bedpe")
sv="/icgc/pcawg/analysis/train2full/March_2016_dataset/${sv_merge_folder}/${pid}.${sv_merge_set}.somatic.sv.bedpe"
ls -lh $sv
ls -lh $cov
if [[ -f ${sv} && -f $cov ]];then
echo plotting
echo "dataFile=c('"$cov"'); svFile='"$sv"';sampleNames=c('"${pid}"','"normal"')" | cat - $rd_script | R --vanilla --slave
python /home/weischej/job_scripts/pawg_scripts/PAWG/dnacopy_segmentation.py `pwd`/${pid}_segmentation.txt
#echo "bedpe.file='${sv}'; header='${pid} Circos plot'; out.file='${pid}_circos.pdf'; " | cat - /home/weischej/job_scripts/plots/CircosPlot.R | R --vanilla --slave
fi

