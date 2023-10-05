#!/usr/bin/env python
## created: 160727
## by: Joachim Weischenfeldt
## joachim.weischenfeldt@gmail.com

from __future__ import print_function
import pandas as pd
import glob, vcf
from pybedtools import BedTool
import glob, vcf, os, sys, shutil
vcfFile = sys.argv[1]
blacklist_svid = sys.argv[2]
aliquot_id = sys.argv[3]
qc = sys.argv[4]

bedpeFile = vcfFile.split('.vcf')[0] + ".bedpe"

vcfFileTmp = vcfFile.replace(".vcf", "_tmp.vcf")
bedpeFileTmp = bedpeFile.replace('.bedpe', '_tmp.bedpe')
shutil.move(vcfFile, vcfFileTmp)
shutil.move(bedpeFile, bedpeFileTmp)
vcfFile = vcfFile.split('.gz')[0]
vcfFile_full = vcfFile.replace('.vcf', '_full.vcf')
print (vcfFile)

print ("Annotating with QC star and TEs. Saving previous vcf file and bedpe file as", vcfFileTmp, "\nand\n", bedpeFileTmp)

df_qc = pd.read_csv(qc, sep="\t")

if not os.path.exists(bedpeFileTmp):
    print ("BEDPE file missing. Cannot annotate")
    bedpeFileTmp = ""
    #df_te = ""
#else:
    #teFile = bedpeFile.replace('.bedpe', '_TE_overlap_tmp.bedpe')
    #tmp = BedTool(bedpeFileTmp).pair_to_bed(BedTool(bed)).saveas(teFile)
    #if os.stat(teFile).st_size == 0:
    #    df_te = pd.DataFrame()
    #else:
    #    df_te = pd.read_csv(teFile, sep="\t", header=None)
    #    df_te = df_te.ix[(df_te[17] == ".") | (df_te[17] == aliquot_id)]
    #print (df_te.shape[0], " TE elements found\n")

#qc = "/icgc/pcawg/analysis/weischej/qc/PCAWG-QC_Summary-of-Measures.tsv"
qc = "/etc/pcawg6_merge_sv/data/PCAWG-QC_Summary-of-Measures.tsv"

df_qc = pd.read_csv(qc, sep="\t")
star_qc = False
if aliquot_id in df_qc.Tumour_WGS_aliquot_ID.tolist():
    star_qc = df_qc.ix[df_qc["Tumour_WGS_aliquot_ID"] == aliquot_id, 'Stars'].tolist()[0]
    print ("star rating =", star_qc)

blacklist_set = set()
blacklist_dict = dict()
with open(blacklist_svid) as readin:
    for f in readin:
        row = f.rsplit()
        if len(row) > 1:
            blacklist_dict[row[0]]=row[1]
            blacklist_set.add(row[1])


vcf_reader=vcf.Reader(open(vcfFileTmp), 'r', compressed=True) if vcfFileTmp.endswith('.gz') else vcf.Reader(open(vcfFileTmp), 'r', compressed=False)
if 'REJECT_CALL' not in vcf_reader.infos.keys():
    vcf_reader.infos['REJECT_CALL'] = vcf.parser._Info('REJECT_CALL', '.', 'String', 'SV call rejected based on QC, germline artefact or transposable element type', '', '')
if star_qc and 'QC_rating' not in vcf_reader.infos.keys():
    vcf_reader.infos['QC_rating'] = vcf.parser._Info('QC_rating', '1', 'Float', 'QC star-rating from PCAWG-QC group [1-5]', '', '')
vcf_writer = vcf.Writer(open(vcfFile, 'w'), vcf_reader, lineterminator='\n')
vcf_full_writer = vcf.Writer(open(vcfFile_full, 'w'), vcf_reader, lineterminator='\n')
for record in vcf_reader:
    te_absent = True
    sv_id_short = record.ID.split("_")[0]
    if star_qc:
        record.INFO['QC_rating'] = star_qc
    #if not  df_te.empty:
    #    if sv_id_short in df_te[6].tolist():
    #        te_absent = False
    #        record.INFO['REJECT_CALL'] = "TE_origin"
    if sv_id_short in blacklist_dict:
        print (sv_id_short, blacklist_dict[sv_id_short])
        record.INFO["REJECT_CALL"] = blacklist_dict[sv_id_short]
    vcf_full_writer.write_record(record)
    if sv_id_short not in blacklist_dict and te_absent:
        vcf_writer.write_record(record)


print (("## VCF file annotated with QC star rating and Transposable Elements.\nFull VCF file \n{0}\nVCF file with removed blacklist SVs \n{1}\nRemoved {2} SVs").format(vcfFile_full, vcfFile, len(list(blacklist_set))))

print (("annotate the vcf file {0} with:\n\tQC file: {1}\nStored un-annotated vcf file as {2}\n").format(vcfFile,  qc, vcfFileTmp))