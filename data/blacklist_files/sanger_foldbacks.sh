DIR=/icgc/pcawg/analysis/train2full/March_2016_dataset/
BLACKLIST_DIR=${DIR}/pcawg_sv_merge/blacklist_files
BLACKLIST_FOLDBACK="${BLACKLIST_DIR}/pcawg6_blacklist_foldback_artefacts.bedpe.gz"
BLACKLIST_FOLDBACK_SLOP="${BLACKLIST_DIR}/pcawg6_blacklist_foldback_artefacts.slop.bedpe.gz"
BLACKLIST_SLOP=200
cd $BLACKLIST_DIR
mkdir tmp
cp $BLACKLIST_DIR/sanger_pipeline_foldback_artefacts.20161113.tar.gz ./tmp 
cd tmp 
tar -xvzf sanger_pipeline_foldback_artefacts.20161113.tar.gz
cat *bedpe | bgzip -c > $BLACKLIST_DIR/${BLACKLIST_FOLDBACK}
cd $BLACKLIST_DIR
bedtools slop -i <(zcat ${BLACKLIST_FOLDBACK} | awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$1,$2,$3,$8,$9,$10}' ) -g ${DIR}/hg19.genome -b $BLACKLIST_SLOP | bedtools slop -i <(awk 'BEGIN{OFS="\t"}{print $4,$5,$6,$1,$2,$3, "foldback_artefact", $7,$8,$9}' - ) -g ${DIR}/hg19.genome -b $BLACKLIST_SLOP | sort -k1,1 -k2,2n |  bgzip -c > $BLACKLIST_FOLDBACK_SLOP
rm -r ./tmp
rm ${BLACKLIST_FOLDBACK}

