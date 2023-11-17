from __future__ import division
import os,sys,re,vcf,gzip,csv
from collections import defaultdict

#Define arguments
merge_file = sys.argv[1] #output from pairToPair ($PAIR2PAIR)
pid = sys.argv[2] #$aliquot_id
slop=400.0 #to do: update so slop is passed to the python script as sys.argv[3]

# Define the headers for the output file
header = ["chrom1_SV1", "pos1_SV1", "chrom2_SV1", "pos2_SV1", "name_SV1", "score_SV1", "strand1_SV1", "strand2_SV1", "svtype_SV1", "center_SV1", "idSV1", "chrom1_SV2", "pos1_SV2", "chrom2_SV2", "pos2_SV2", "name_SV2", "score_SV2", "strand1_SV2", "strand2_SV2", "svtype_SV2", "center_SV2", "idSV2", "pid", "interSect", "BpOffset", "svOverlapType", "callerPair", "rn_share_count", "rn_share_fract"]
print('\t'.join(map(str, header)))

# Line size limit might be compromised.
old_limit = csv.field_size_limit()
csv.field_size_limit(13107200)

#Loop through all of the rows in the PAIR2PAIR output file, merge them, then output a bedpe file.
with open(merge_file, "r") as rin:
    reader = csv.reader(rin, delimiter="\t", quoting=csv.QUOTE_NONE)
    center_seen = set()
    rowadj = 0

    for row in reader:
        if row[-1] == pid: # Check if the value in the last column is equal to our specified sample id.
            center1 = row[12] #Specify the sequencing center
            center2 = row[25] #Specify the sequencing center
            centerset = '__'.join(map(str, sorted([row[6],row[19]])))

            if center1 != center2 and centerset not in center_seen: 
                center_seen.add(centerset) # Keep track of all of the unique pairs we have seen so far

                if row[0] == row[13]: # if chrom1-caller_1 == chrom1-caller_2 
                    distPos1 = abs(int(row[1]) - int(row[14])) # Then calculate the abs of the difference of start1-caller1 and start1-caller2 
                elif row[0] == row[16]: # if chrom1-caller_1 == chrom2-caller2 
                    distPos1 = abs(int(row[1]) - int(row[17])) # Then calculate the abs of the difference of start1-caller1 and start2-caller2 

                if row[3] == row[16]: # If chrom2-caller1 == chrom2-caller2 
                    distPos2 = abs(int(row[4]) - int(row[17])) # Then calculate the abs of the difference in start2-caller1 and start2-caller2 
                elif row[3] == row[13]: # If chrom2-caller1 == chrom1-caller2 
                    distPos2 = abs(int(row[4]) - int(row[14])) # Then calculate the abs of the difference in start2-caller1 and start1-caller2 

                interSect = (4*slop - distPos1 - distPos2)/(4*slop)
                BpOffset = distPos1 + distPos2
                svType1 = row[10] # The original bedpe has sv_type and sv_class. I think in our case we can just pretend that both are the same?
                svType2 = row[23]

                svOverlapType = svType1 + "_" + svType2
                callerPair = center1 + "_" +  svType1 + "|" + center2 + "_" + svType2

                idSV1 = center1 + "_" + row[6] + "_" + row[0] + ":" + row[1] + "-" + row[3] + ":" + row[4] #Construct a SV ID string
                idSV2 = center2 + "_" + row[19] + "_" + row[13] + ":" + row[14] + "-" + row[16] + ":" + row[17] #Construct a SV ID string

                rn_SV1 = set(i for i in row[6])
                rn_SV2 = set(i for i in row[19]) 
                rn_union = rn_SV1.union(rn_SV2) 
                rn_shared = rn_SV1.intersection(rn_SV2) 
                rn_share_count = len(rn_shared) 

                if not 'NA' in rn_SV1 and not 'NA' in rn_SV2 and not rn_shared:
                    continue # kick it out if each caller has readid but none are shared
                try:
                    rn_share_fract = round(rn_share_count/len(rn_union),2)
                except Exception as e:
                    rn_share_fract = 0

                #Format the output bedpe file
                out = row[0:2] + row[3:5] + row[6:11] + [center1, idSV1] + row[13:15] + row[16:18] + row[19:24] + [center2, idSV2, pid] + [interSect, BpOffset, svOverlapType, callerPair, rn_share_count, rn_share_fract] # index updated
                
                print ('\t'.join(map(str, out)))