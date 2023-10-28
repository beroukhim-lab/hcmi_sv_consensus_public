from __future__ import division
import os,sys,re,vcf,gzip,csv
from collections import defaultdict

merge_file = sys.argv[1]
pid = sys.argv[2]

# TODO: Why are different slop values being used (400 vs 500?
#Maybe let's just put this as an argument for the original bash script, then keep it consistent throughout? This also does not make sense to me.
slop=500.0 

'''
This might be unnecessary but writing out all the column indices so I don't have to keep scrolling below
0: chrom1_SV1
1: pos1_SV1
2: chrom2_SV1
3: pos2_SV1
4: name_SV1
5: score_SV1
6: strand1_SV1
7: strand2_SV1
8: svtype_SV1
9: center_SV1
10: idSV1
11: chrom1_SV2
12: pos1_SV2
13: chrom2_SV2
14: pos2_SV2
15: name_SV2
16: score_SV2
17: strand1_SV2
18: strand2_SV2
19: svtype_SV2
20: center_SV2
21: idSV2
22: pid
23: interSect
24: BpOffset
25: svOverlapType
26: callerPair
27: rn_share_count
28: rn_share_fract
'''

###The bedpe files from the original code should have the following headers:
#bedpe_header="chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tqual\tstrand1\tstrand2\tsv_class\tsv_type\tscore\tsupp_reads\tscna\tcenter\tread_id"
#I think that the original pair2pair output headers should be:
'''
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

#define the headers for the output file?
header = ["chrom1_SV1", "pos1_SV1", "chrom2_SV1", "pos2_SV1", "name_SV1", "score_SV1", "strand1_SV1", "strand2_SV1", "svtype_SV1", "center_SV1", "idSV1", "chrom1_SV2", "pos1_SV2", "chrom2_SV2", "pos2_SV2", "name_SV2", "score_SV2", "strand1_SV2", "strand2_SV2", "svtype_SV2", "center_SV2", "idSV2", "pid", "interSect", "BpOffset", "svOverlapType", "callerPair", "rn_share_count", "rn_share_fract"]


print('\t'.join(map(str, header)))

## line size limit might be compromised. Increase
old_limit = csv.field_size_limit()
csv.field_size_limit(13107200)

# Annotating below is hard without having a simple input file

with open(merge_file, "r") as rin:
    reader = csv.reader(rin, delimiter="\t", quoting=csv.QUOTE_NONE)
    center_seen = set()
    rowadj = 0

    for row in reader:

        #if len(row) == 33: #I don't know what this is supposed to be checking for, but I don't think we need it.
            #rowadj = -1

        if row[-1] == pid: #Check if the value in the last column is equal to our specified sample id. When would this ever be false?
            #center1 = row[15] #OG -- I don't think we have a value for this in our bedpe, so let's just set it to the caller name for now
            center1 = row[12] #Added this for troubleshooting. The original code is the line above.
            #center2 = row[32 + rowadj] #OG -- I don't think we have a value for this in our bedpe, so let's just set it to the caller name for now
            center2 = row[24] #Added this for troubleshooting. The original code is line above.
            #centerset = '__'.join(map(str, sorted([row[6],row[23 + rowadj]]))) #OG
            centerset = '__'.join(map(str, sorted([row[12],row[24]]))) #Updated

            ####Stopped editing here

            if center1 != center2 and centerset not in center_seen: # not same caller
                center_seen.add(centerset)

                if row[0] == row[17 + rowadj]:
                    distPos1 = abs(int(row[1]) - int(row[18 + rowadj]))
                elif row[0] == row[20 + rowadj]:
                    distPos1 = abs(int(row[1]) - int(row[21 + rowadj]))

                if row[3] == row[20 + rowadj]:
                    distPos2 = abs(int(row[4]) - int(row[21 + rowadj]))
                elif row[3] == row[17]:
                    distPos2 = abs(int(row[4]) - int(row[18 + rowadj]))

                interSect = (4*slop - distPos1 - distPos2)/(4*slop)
                BpOffset = distPos1 + distPos2
                svType1 = row[10]

                if svType1 == "NA":
                    svType1 = row[11]
                    svType2 = row[27]

                if svType2 == "NA":
                    svType2 = row[28 + rowadj]

                svOverlapType = svType1 + "_" + svType2
                callerPair = center1 + "_" +  svType1 + "|" + center2 + "_" + svType2
                idSV1 = center1 + "_" + row[6] + "_" + row[0] + ":" + row[1] + "-" + row[3] + ":" + row[4]
                idSV2 = center2 + "_" + row[23+ rowadj] + "_" + row[17] + ":" + row[18 + rowadj] + "-" + row[20+ rowadj] + ":" + row[21+ rowadj]
                rn_SV1 = set([i.strip('/[12]').replace('-', ':') for i in row[16].split(',')])
                rn_SV2 = set([i.strip('/[12]').replace('-', ':') for i in row[33 + rowadj].split(',')])
                rn_union = rn_SV1.union(rn_SV2)
                rn_shared = rn_SV1.intersection(rn_SV2)
                rn_share_count = len(rn_shared)

                if not 'NA' in rn_SV1 and not 'NA' in rn_SV2 and not rn_shared:
                    continue
                    # kick it out if each caller has readid but none are shared
                try:
                    rn_share_fract = round(rn_share_count/len(rn_union),2)
                except Exception, e:
                    rn_share_fract = 0

                if rowadj == -1:
                    out = row[0:2] + row[3:5] + row[6:11] + [center1, idSV1] + row[16:18] + row[19:21] + row[22:27] + [center2, idSV2, pid] + [interSect, BpOffset, svOverlapType, callerPair, rn_share_count, rn_share_fract]
                else:
                    out = row[0:2] + row[3:5] + row[6:11] + [center1, idSV1] + row[17:19] + row[20:22] + row[23:28] + [center2, idSV2, pid] + [interSect, BpOffset, svOverlapType, callerPair, rn_share_count, rn_share_fract]
                
                print ('\t'.join(map(str, out)))