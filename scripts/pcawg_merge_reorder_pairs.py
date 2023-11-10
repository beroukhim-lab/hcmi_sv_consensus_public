from __future__ import division
import os,sys,re,vcf,gzip,csv
from collections import defaultdict

merge_file = sys.argv[1]
pid = sys.argv[2]
slop=400.0

'''
# New pair2pair output headers should be:
0: chrom1 caller_1
1: start1 caller_1
2: end1 caller_1
3: chrom2 caller_1
4: start2 caller_1
5: end2 caller_1
6: sv_id (name) caller_1
7: pe_support caller_1
8: strand1 caller_1
9: strand2 caller_1
10: sv_class caller_1
11: sv_method caller_1
12: center caller_1
13: chrom1 caller_2
14: start1 caller_2
15: end1 caller_2
16: chrom2 caller_2
17: start2 caller_2
18: end2 caller_2
19: sv_id (name) caller_2
20: pe_support caller_2
21: strand1 caller_2
22: strand2 caller_2
23: sv_class caller_2
24: sv_method caller_2
25: center caller_2
26: aliquot_id
'''

'''
# Original pair2pair output headers should be:
0: chrom1 caller_1
1: start1 caller_1
2: end1 caller_1
3: chrom2 caller_1
4: start2 caller_1
5: end2 caller_1
6: name caller_1
7: qual caller_1
8: strand1 caller_1
9: strand2 caller_1
10: sv_class caller_1
11: sv_type caller_1
12: score caller_1
13: supp_reads caller_1
14: scna caller_1
15: center caller_1
16: read_id caller_1
17: chrom1 caller_2
18: start1 caller_2
19: end1 caller_2
20: chrom2 caller_2
21: start2 caller_2
22: end2 caller_2
23: name caller_2
24: qual caller_2
25: strand1 caller_2
26: strand2 caller_2
27: sv_class caller_2
28: sv_type caller_2
29: score caller_2
30: supp_reads caller_2
31: scna caller_2
32: center caller_2
33: read_id caller_2
34: sample_id
'''

# Define the headers for the output file
header = ["chrom1_SV1", "pos1_SV1", "chrom2_SV1", "pos2_SV1", "name_SV1", "score_SV1", "strand1_SV1", "strand2_SV1", "svtype_SV1", "center_SV1", "idSV1", "chrom1_SV2", "pos1_SV2", "chrom2_SV2", "pos2_SV2", "name_SV2", "score_SV2", "strand1_SV2", "strand2_SV2", "svtype_SV2", "center_SV2", "idSV2", "pid", "interSect", "BpOffset", "svOverlapType", "callerPair", "rn_share_count", "rn_share_fract"]
print('\t'.join(map(str, header)))

# Line size limit might be compromised. Increase
old_limit = csv.field_size_limit()
csv.field_size_limit(13107200)

with open(merge_file, "r") as rin:
    reader = csv.reader(rin, delimiter="\t", quoting=csv.QUOTE_NONE)
    center_seen = set()
    rowadj = 0

    for row in reader:
        if row[-1] == pid: # Check if the value in the last column is equal to our specified sample id. When would this ever be false?
            center1 = row[12]
            center2 = row[25]
            centerset = '__'.join(map(str, sorted([row[6],row[19]]))) # Index updated

            if center1 != center2 and centerset not in center_seen: 
                center_seen.add(centerset) # Keep track of all of the unique pairs we have seen so far

                if row[0] == row[13]: # if chrom1-caller_1 == chrom1-caller_2 -- index updated
                    distPos1 = abs(int(row[1]) - int(row[14])) # Then calculate the abs of the difference of start1-caller1 and start1-caller2 -- index updated
                elif row[0] == row[16]: # if chrom1-caller_1 == chrom2-caller2 -- index updated
                    distPos1 = abs(int(row[1]) - int(row[17])) # Then calculate the abs of the difference of start1-caller1 and start2-caller2 -- index updated

                if row[3] == row[16]: # If chrom2-caller1 == chrom2-caller2 -- index updated
                    distPos2 = abs(int(row[4]) - int(row[17])) # Then calculate the abs of the difference in start2-caller1 and start2-caller2 -- index updated
                elif row[3] == row[13]: # If chrom2-caller1 == chrom1-caller2 -- index updated
                    distPos2 = abs(int(row[4]) - int(row[14])) # Then calculate the abs of the difference in start2-caller1 and start1-caller2 -- index updated

                interSect = (4*slop - distPos1 - distPos2)/(4*slop)
                BpOffset = distPos1 + distPos2
                svType1 = row[10] # The original bedpe has sv_type and sv_class. I think in our case we can just pretend that both are the same? -- index updated
                svType2 = row[23]

                svOverlapType = svType1 + "_" + svType2
                callerPair = center1 + "_" +  svType1 + "|" + center2 + "_" + svType2

                idSV1 = center1 + "_" + row[6] + "_" + row[0] + ":" + row[1] + "-" + row[3] + ":" + row[4] # index updated
                idSV2 = center2 + "_" + row[19] + "_" + row[13] + ":" + row[14] + "-" + row[16] + ":" + row[17] # index updated

                # Not entirely sure we need to do the next 2 lines, may be specific for the pcawg bedpe read IDs. Will keep and go over SV_ID just in case
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

                out = row[0:2] + row[3:5] + row[6:11] + [center1, idSV1] + row[13:15] + row[16:18] + row[19:24] + [center2, idSV2, pid] + [interSect, BpOffset, svOverlapType, callerPair, rn_share_count, rn_share_fract] # index updated
                
                print ('\t'.join(map(str, out)))