#! /usr/bin/env python

from __future__ import print_function
from __future__ import division
import matplotlib
matplotlib.use('Agg')
import networkx
import networkx as nx
import matplotlib.pyplot
import argparse
import datetime
import sys
import collections
import numpy
import random
import uuid
import copy
import csv
import os
import vcf
import re
from collections import defaultdict, Counter
import gzip
import numpy as np
import pandas as pd

today = datetime.date.today().strftime("%Y%m%d")


##########################################################
# Part 1: Setup 
# Parse arguments 
# Setup graph structure and order clique calls
##########################################################

#This function will grab the read ID (everything before the undscore)
def _getShortID(recordID):
    shortID=recordID.split('_')[0]
    return shortID


# Parse command line
parser=argparse.ArgumentParser(description='Graph builder.')
parser.add_argument('-e', '--inBEDPE_overlap', metavar='overlap.bedpe', required=True, dest='inBEDPE_overlap', help='input _SV_overlap.txt table from pairToPair (required)')
parser.add_argument('-z', '--inBEDPE_array', metavar='bedpe_array_in.bedpe', required=True, dest='inBEDPE_array', help='input array of bedpe files as string') #Load in the array of bedpe files
parser.add_argument('-o', '--outBEDPE', metavar='out.bedpe', required=False, dest='outBEDPE', help='output bedpe file (optional)')
parser.add_argument('-s', '--stat', metavar='stat.bed', required=False, dest='outStat', help='output statistics (optional)')
args = parser.parse_args()

inBEDPE_overlap = args.inBEDPE_overlap
inBEDPE_array = args.inBEDPE_array.split(' ')
outBEDPE = args.outBEDPE
outSTAT = args.outStat

minCarrier = 2
copyConc = 0
priority=collections.defaultdict(float)


#Overlap graph
pid = os.path.basename(os.path.dirname(inBEDPE_overlap))
G=networkx.Graph()
#Parse bed file SV padding
f_reader=csv.DictReader(open(inBEDPE_overlap), delimiter="\t")
for row in f_reader:
    pid = row['pid']
    if not G.has_edge(row['name_SV1'], row['name_SV2']):
        G.add_edge(row['name_SV1'], row['name_SV2'], weight=float(row['rn_share_fract']))
        G.add_edge(row['name_SV1'], row['name_SV2'], weight=float(row['rn_share_fract']))
        svType1 = row['svtype_SV1']
        if row['svtype_SV1'] == "NA":
            svType1 = "BND"
        svType2 = row['svtype_SV2']
        if row['svtype_SV2'] == "NA":
            svType2 = "BND"
        G.nodes[row['name_SV1']]['SVtype']=(svType1, row['pid'], row['center_SV1'],row['idSV1']) 
        G.nodes[row['name_SV2']]['SVtype']=(svType2, row['pid'], row['center_SV2'],row['idSV2'])

#Parse graph structure
compAssign=dict()
compMembers=collections.defaultdict(list)
conCompCount=0
cliqueCount=0
for H in (G.subgraph(c) for c in nx.connected_components(G)):
    conCompCount+=1
    typeSet=set()
    sample=set()
    center=list()
    for n, d in H.nodes(data=True):
        typeSet.add(d['SVtype'][0])
        sample.add(d['SVtype'][1])
        center.append(d['SVtype'][2])
    baseDir = 'plots/' + pid + '/'
    try:
        os.makedirs(baseDir)
    except Exception as e:
        pass
    baseName=  '_'.join(list(sample) + sorted(typeSet)) + '.' + str(conCompCount) + '_' + ';'.join(map(str, list(set(center)))) + '.png' 
    for n in H.nodes(data=False):
         compAssign[n]=conCompCount
         compMembers[conCompCount].append(n)

    # Clique?
    if (float(H.number_of_edges())==float(H.number_of_nodes()*(H.number_of_nodes() - 1)) / 2.0):
        cliqueCount+=1
    else:
        pass

# Print summary
if outSTAT:
    bedOut = open(outSTAT, 'w')
    print('copyConc', 'minCarrier', 'number_of_nodes', 'connected_components_Count', 'clique_Count', sep="\t", file=bedOut)
    print(copyConc, minCarrier, G.number_of_nodes(), conCompCount, cliqueCount, sep="\t", file=bedOut)
    bedOut.close()

cliqueFile = inBEDPE_overlap.replace('.txt','.clique.txt')
with open(cliqueFile, 'w') as w:
    writer = csv.writer(w, delimiter="\t", lineterminator="\n")
    for k,v in compAssign.items(): #python 3 replaced dict.iteritems() with dict.items()
        out = [k,v]
        writer.writerow(out)


# Order the calls in a clique by assumed genotyping accuracy
for compID in compMembers.keys():
    outOrder=list()
    for callX in compMembers[compID]:
        xShortID=_getShortID(callX)
        insertPoint=0
        for callIndex, callY in enumerate(outOrder):
            yShortID=_getShortID(callY)
            if (xShortID==yShortID):
                if (priority[callX]>priority[callY]):
                    print ("\t\tSAMPLE")
                    break
            insertPoint+=1
        outOrder=outOrder[:insertPoint] + [callX] + outOrder[insertPoint:]
    compMembers[compID]=outOrder




##########################################################
# Part 2: Merge Calls
##########################################################
#Define a function that will extract the bedpe rows from a bedpe input and create the 'mastermerge' output file.
def extract_bedpe_rows(single_bedpe, mastermerge, header_concat, compAssign):
    with open(single_bedpe, 'rb') as rin:
        # I removed headers in my bedpe tmp files for pairToPair --> manually assign header
        header_concat.append("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tqual\tstrand1\tstrand2\tsv_class\tsv_type\tscore\tsupp_reads\tscna\tcenter\tread_id\tcenter")
        for f in rin:
            row = f.decode("utf-8").split() # Byte string formatting was screwing up recognition below -- edited to resolve
            # print(row)
            rinfos = dict()
            svid = row[6]
            # print(svid)
            if svid in compAssign:
                cliqid = compAssign.get(svid)
                mastermerge[cliqid].append(row)
    return mastermerge, header_concat

mastermerge = defaultdict(list)
header_concat = list()
for single_bedpe in inBEDPE_array:
    # print(single_bedpe)
    mastermerge, header_concat = extract_bedpe_rows(single_bedpe, mastermerge, header_concat, compAssign)


#This function will take the input (mastermerge) file, will collapse SVs into a single SV per clique, and will output that information in bedpe format.
#Input: mastermerge
#Output: collapsed bedpe with consensus SVs
def generate_bedpe_new(input_cliq_bedpe):
    
    output_bedpe=pd.DataFrame(columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'sv_id', 'pe_support', 'strand1', 'strand2', 'sv_class', 'svmethod', 'center']) 

    for key, sublist in input_cliq_bedpe.items():
            # Init some objects
            new_chrom1_set=set()
            new_chrom2_set=set()
            new_strand1_set=set()
            new_strand2_set=set()
            new_svtype_set=set()
            new_start1_list=list()
            new_start2_list=list()
            new_caller_set=set() # caller = the algorithms making the call. Note, in some cases one group uses multiple callers
            new_center_set=set() # center = institute (e.g broad, nygc, etc)
            new_qual_list=list() # Every group reports qual in a different way. Does anyone care about qual? Let's just assume that all groups pre-filter for their own qual threshold. I do not want to parse these strings.
            new_id_list=list()

            final_chrom1=None # This is the final chrom1 value for the output bedpe
            final_chrom2=None # This is the final chrom2 value for the output bedpe
            final_strand1=None # This is the final strand1 value for the output bedpe
            final_strand2=None # This is the final strand2 value for the output bedpe
            final_svtype=None # This is the final svtype for the output bedpe
            final_start1=None # This is the final start1 position for the output bedpe
            final_start2=None # This is the final start2 position for the output bedpe
            final_end1=None # This is the final end1 value for the output bedpe
            final_qual_string=None # This is the final qual string for the output bedpe, will be ";" sep and order will be same as final_center_string
            final_caller_string=None # This is the final caller string for the output bedpe, will be ";" sep.
            final_center_string=None # This is the final center string for the output bedpe, will be ";" sep and order will be the same as final_qual_string
            final_id_string=None # This is a list of all SV IDs for the output bedpe, will be ";" sep.
            final_bedpe_line=None # This is the final bedpe line (everything glued together) that will be output for each cliq

            for item in sublist: # For each sv call (one from each institution) for that cliq

                # Add all of the chrom1 and chrom2 values to a set.
                # We have downstream ambitions of checking to make sure for each cliq these values are always the same across all callers
                # Do something similar for strand
                # Also do something similar for sv type
                new_chrom1_set.add(item[0]) # 0 is the index of chrom1 in bedpe
                new_chrom2_set.add(item[3]) # 3 is the index of chrom2 in bedpe
                new_strand1_set.add(item[8]) # 8 is the index for strand1 in bedpe
                new_strand2_set.add(item[9]) # 9 is the index for strand2 in bedpe
                new_svtype_set.add(item[10]) # 10 is the index for sv type in bedpe

                # Add all of the start1 and start2 values to a list
                # We don't care about end1 and end2 because end1/2 == start1/2 +1, so we can just regenerate later.
                new_start1_list.append(item[1]) # 1 is the index for start1 in bedpe
                new_start2_list.append(item[4]) # 4 is the index for start2 in bedpe

                # Add all of the callers and centers to a set
                new_caller_set.update(item[11].split(',')) # Sometimes a single group uses multiple callers. Add each indivudal caller to the set
                new_center_set.add(item[12]) # This the index of the center (e.g. broad or nygc)

                # Add all of the qual scores to a list
                # Also add all of the sv IDs to a list
                new_qual_list.append(item[7]) # 7 is the inde of the qual score in bedpe
                new_id_list.append(item[6]) # 6 is the index of the SV ID in bedpe


        
            # Now that everything is assembled, calculate the values that will go into the actual bedpe
            # We can write this in a more beautiful way in the future, for now writing it like this may help with debugging.
            # to do: in the future just have it error and break if any of the elif conditions are met. Right now it is written this way for debugging purposes.
            # First, chrom1
            if len(new_chrom1_set) == 1:
                final_chrom1=new_chrom1_set.pop()
            elif len(new_chrom1_set) > 1:
                final_chrom1="error_multiple_chrom1"
            else:
                final_chrom1="error_no_chrom1"

            # Now chrom2
            if len(new_chrom2_set) == 1:
                final_chrom2=new_chrom2_set.pop()
            elif len(new_chrom2_set) > 1:
                final_chrom2="error_multiple_chrom2"
            else:
                final_chrom2="error_no_chrom2"

            # Now for strand1
            if len(new_strand1_set) == 1:
                final_strand1=new_strand1_set.pop()
            elif len(new_strand1_set) > 1:
                final_strand1="error_multiple_strand1"
            elif len(new_strand1_set) == 0:
                final_strand1="error_no_strand1"

            # Now for strand2
            if len(new_strand2_set) == 1:
                final_strand2=new_strand2_set.pop()
            elif len(new_strand2_set) > 1:
                final_strand2="error_multiple_strand2"
            elif len(new_strand2_set) == 0:
                final_strand2="error_no_strand2"

            # Now for sv type
            if len(new_svtype_set) == 1:
                final_svtype=new_svtype_set.pop()
            elif len(new_svtype_set) > 1:
                final_svtype="error_multiple_svtype"
            elif len(new_svtype_set) == 0:
                final_svtype="error_no_svtype"

            # Calculate the position values using the same logic as the original generate_vcf function
            if final_strand1 == '+':
                final_start1=max(new_start1_list)
                final_end1=int(final_start1)+1
            elif final_strand1 == '-':
                final_start1=min(new_start1_list)
                final_end1=int(final_start1)+1
        
            if final_strand2 == "+":
                final_start2=max(new_start2_list)
                final_end2=int(final_start2)+1
            elif final_strand2 == '-':
                final_start2=min(new_start2_list)
                final_end2=int(final_start2)+1

            # Construct a string for qual, caller, center, and id
            final_qual_string=';'.join(new_qual_list)
            final_caller_string=';'.join(new_caller_set)
            final_center_string=';'.join(new_center_set)
            final_id_string=';'.join(new_id_list)

            # Generate a final bedpe line with all of the information above
            final_bedpe_line=[final_chrom1, final_start1, final_end1, final_chrom2, final_start2, final_end2, final_id_string, final_qual_string, final_strand1, final_strand2, final_svtype, final_caller_string, final_center_string]

            # Add it to a pandas dataframe
            output_bedpe.loc[len(output_bedpe)] = final_bedpe_line

    #Now return the bedpe :)
    return output_bedpe


#Run the above function above to convert from mastermerge -> collapsed_bedpe
collapsed_bedpe = generate_bedpe_new(mastermerge)


# Search for high number of small SVs - likely artifacts
collapsed_bedpe['sv_class'] = collapsed_bedpe['sv_class'].str.replace('h2h', '').str.replace('t2t', '') #Trim the h2h and t2t substrings in the sv_class column
statCounterTmp = Counter(collapsed_bedpe['sv_class']) #Init counter
statCounter = dict() 
bins=[0, 5e3, 1e12] #Bins for counting numbers of small SVs
sumSVs= sum(statCounterTmp.values()) #Count the number of times each type of sv class was observed


for k,v in statCounterTmp.items():
    statCounter[k + '_count'] = v
    statCounter[k + '_count_fract'] = round(v/sumSVs,3)
    statCounter[k + '_size_hist'] = np.nan

    if k != 'TRA':
        # Filter the DataFrame
        filtered_df = collapsed_bedpe[(collapsed_bedpe['sv_class'].str.contains(k)) & (collapsed_bedpe['chrom1'] == collapsed_bedpe['chrom2'])]

        # Calculate the difference
        filtered_df = filtered_df.dropna(subset=['start1']).dropna(subset=['start2'])
        filtered_df['start1'] = filtered_df['start1'].astype(int)
        filtered_df['start2'] = filtered_df['start2'].astype(int)
        differences = filtered_df['start2'] - filtered_df['start1']

        # Create a histogram
        statCounter[k + '_size_hist'] = np.histogram(differences, bins=bins)


print ('\nTotal number of SVs:',sumSVs,'\nDistribution of SV types with bins', ' - '.join(map(str,bins)) )
for k,v in statCounter.items():
    print (k, v)


###I don't think this code block will work, because I don't think it will work for a pandas df
###But this code isn't producing an error because sumSVs < 200 for our example.
###We need to fix this in the future
mergeIdRemove = set()
if sumSVs > 200:

    for k in statCounterTmp.keys():
        size_small_sv = statCounter[k + '_count']

        if k != 'TRA':
            size_small_sv = statCounter[k + '_size_hist'][0][0]

        if (statCounter[k  + '_count_fract']) > 0.75 and size_small_sv/statCounter[k + '_count'] > 0.8:
            print ("Uneven high number of {0}. {1} out of {2} with {3} having a size smaller than 5e3".format(k, statCounterTmp[k], sumSVs, size_small_sv))
            
            for i in collapsed_bedpe:
                if k in i[10] and ((i[0] != i[3]) or (i[0] == i[3] and i[4]-i[1] <  statCounter[k + '_size_hist'][1][1])):
                    mergeIdRemove.add( i[6].replace('_1', '').replace('_2', ''))


d = dict()
df_filter = pd.DataFrame.from_dict(d, orient='index')

if statCounter: 

    for k,v in statCounter.items():

        if 'count' in k:
            d[k] = v

        if 'size_hist' in k:
            ks = k.split('_')[0] + "_count_lt_5kb"
            try:
                d[ks] = v[0][0]
            except Exception as e:
                d[ks] = np.nan    
    

collapsed_bedpe_filter = list()
if mergeIdRemove:
    for bed in collapsed_bedpe:
        if not re.sub(r'_[12]$', '',bed[6]) in mergeIdRemove:
            collapsed_bedpe_filter.append(bed)


bedpeheader = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'sv_id', 'pe_support', 'strand1', 'strand2', 'svclass', 'svmethod']


#Remove rows where there is an error
#NEED TO CHECK THE CAUSE OF THESE ERRORS IN THE FUTURE
#These errors come from string or chr mismatches in the original cliq assembly
collapsed_bedpe = collapsed_bedpe[~collapsed_bedpe.apply(lambda x: x.astype(str).str.contains('error').any(), axis=1)]

#Write the file
collapsed_bedpe.to_csv(outBEDPE, sep='\t', index=False)
print ("merged SVs into BEDPE format\n", collapsed_bedpe, "\n\n")