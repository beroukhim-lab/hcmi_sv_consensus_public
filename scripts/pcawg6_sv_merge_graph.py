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

# Original script call: python2 scripts/pcawg6_sv_merge_graph.py -e ${inBEDPE} -o ${outVCF} -s ${outSTAT} -a $SANGER_ANNO_VCF -b $DELLY_ANNO_VCF -c $DRANGER_ANNO_VCF -d $SNOWMAN_ANNO_VCF | tee -a $log
# New script call: python scripts/pcawg6_sv_merge_graph.py -e ${inBEDPE} -z "$bedpe_python_input_array" -o ${outVCF} -s ${outSTAT} | tee -a $log

##########################################################
# Part 1: Setup 
# Parse arguments 
# Setup graph structure and order clique calls
##########################################################

def _getShortID(recordID):
    shortID=recordID.split('_')[0]
    return shortID

# Plot graph
def plotGraph(g, fileName):
    #networkx.draw_spring(g)
    #matplotlib.rcParams.update({'font.size': 5})
    networkx.draw_circular(g)
    matplotlib.pyplot.savefig(fileName)
    matplotlib.pyplot.clf()
    matplotlib.pyplot.close('all')


# Use triangle coordinates
def inTriangle(x, y, cutX, cutY):
    v0 = numpy.array([1.0, cutY]) - 1
    v1 = numpy.array([cutX, 1.0]) - 1
    v2 = numpy.array([x, y]) - 1
    dot00 = numpy.dot(v0, v0)
    dot01 = numpy.dot(v0, v1)
    dot02 = numpy.dot(v0, v2)
    dot11 = numpy.dot(v1, v1)
    dot12 = numpy.dot(v1, v2)
    invDenom = 1.0 / float(dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom
    return (u >= 0) and (v >= 0) and (u + v < 1)

# Debugging: List all of the bedpe files in the array
def process_file_paths():
    for path in inBEDPE_array:
        print("File path: %s" % path)

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

# Overlap graph
process_file_paths()
pid = os.path.basename(os.path.dirname(inBEDPE_overlap))
G=networkx.Graph()
# Parse bed file SV padding
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
        G.nodes[row['name_SV1']]['SVtype']=(svType1, row['pid'], row['center_SV1'],row['idSV1']) # networkx 2.4 requires G.nodes[] instead of G.node[]
        G.nodes[row['name_SV2']]['SVtype']=(svType2, row['pid'], row['center_SV2'],row['idSV2'])

# Parse graph structure
compAssign=dict()
compMembers=collections.defaultdict(list)
conCompCount=0
cliqueCount=0
for H in (G.subgraph(c) for c in nx.connected_components(G)): # connected_components_subgraphs(G) deprecated --> use replacement
    conCompCount+=1
    typeSet=set()
    sample=set()
    center=list()
    for n, d in H.nodes(data=True): # nodes() instead of nodes.iter()
        typeSet.add(d['SVtype'][0])
        sample.add(d['SVtype'][1])
        center.append(d['SVtype'][2])
    baseDir = 'plots/' + pid + '/'
    try:
        os.makedirs(baseDir)
    except Exception as e:
        pass
    baseName=  '_'.join(list(sample) + sorted(typeSet)) + '.' + str(conCompCount) + '_' + ';'.join(map(str, list(set(center)))) + '.png' 
    
    
    '''
    # This section was used in original code if plotSubgraph was provided as an arg, but it was not in bash, so I removed the arg
    if plotSubgraph:
        elarge=[(u,v) for (u,v,d) in H.edges(data=True) if d['weight'] >0.3]
        esmall=[(u,v) for (u,v,d) in H.edges(data=True) if d['weight'] <=0.3]
        pos=nx.circular_layout(H) # positions for all nodes
        # nodes
        centercolor = [i.replace('sanger', '#ADD8E6').replace('embl', '#90EE90').replace('dRanger', '#F0E68C').replace('smufin', '#800080').replace('snowman', '#FFA500') for i in center]

        nx.draw_networkx_nodes(H,pos,node_size=700 , node_color=centercolor, font_family="arial")
        # edges
        weight = H.edges(data=True)[-1][2]['weight']
        nx.draw_networkx_edges(H,pos,edgelist=elarge,
                            width=2+4*weight)
        nx.draw_networkx_edges(H,pos,edgelist=esmall,
                            width=2+4*weight,alpha=0.5,edge_color='b',style='dashed')
        plotGraph(H, baseName)
    '''
    for n in H.nodes(data=False): ###I updated nodes_iter to nodes because the version of the python package I'm using changed nodes_iter -> nodes
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
    #print(compID, compMembers[compID])

# ^ Things seem to be working above this point (or at least they're running, no clue if the output is as intended)

##########################################################
# Part 2: Merge Calls
# Original code read in VCF files into a hash table; did not merge due to formatting differences
##########################################################

def read_orient(record):
    strands = ['NA', 'NA']
    if re.search(r'^[ACTGN]+\]',str(record)):
        strands = ['+', '+']
    if re.search(r'^[ACTGN]+\[',str(record)):
        strands = ['+', '-']
    if re.search(r'\][ACTGN]+$',str(record)):
        strands = ['-', '+']
    if re.search(r'\[[ACTGN]+$',str(record)):
        strands = ['-', '-']
    return strands

def mergeinfo(d1, d2):
    return reduce(lambda a,b: dict (a, **b), (d1,d2))

def sv_class(chrom1, strand1, chrom2, strand2):
    if chrom1 != chrom2:
        return 'TRA'
    elif strand1 == "+":
        if strand2 == "+":
            return 'h2hINV'
        elif strand2 == "-":
            return "DEL"
    elif strand1 == "-":
        if strand2 == "-":
            return 't2tINV'
        elif strand2 == "+":
            return "DUP"

# Create a new function that should do the same thing as extract_vcf_rows, except by using a bedpe input
def extract_bedpe_rows(single_bedpe, mastermerge, header_concat, compAssign):
    # I had to remove headers in the .tmp files for pairToPair, but we'll want to make use of them here:
    # print(compAssign)
    with open(single_bedpe, 'rb') as rin:
        # I removed headers in my bedpe tmp files for pairToPair --> manually assign header
        header_concat.append("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tqual\tstrand1\tstrand2\tsv_class\tsv_type\tscore\tsupp_reads\tscna\tcenter\tread_id\tcenter")
        for f in rin:
            row = f.decode("utf-8").split() # Byte string formatting was screwing up recognition below -- edited to resolve
            rinfos = dict()
            svid = row[6]
            if svid in compAssign:
                cliqid = compAssign.get(svid)
                mastermerge[cliqid].append(row)
    return mastermerge, header_concat

'''
mastermerge = defaultdict(list)
header_concat = list()
for inVCF in inVCF_DELLY, inVCF_BRASS, inVCF_DRANGER, inVCF_SNOWMAN:
    mastermerge, header_concat = extract_vcf_rows(inVCF, mastermerge, header_concat, compAssign)
'''

mastermerge = defaultdict(list)
header_concat = list()
for single_bedpe in inBEDPE_array:
    mastermerge, header_concat = extract_bedpe_rows(single_bedpe, mastermerge, header_concat, compAssign)

'''
###Debugging
#Print mastermerge to see what it looks like
print("This is the structure of mastermerge")
num_entries_to_print = 5  # Number of entries you want to print

# Iterating over the dictionary
for i, (cliqid, bedpe_rows) in enumerate(mastermerge.items()):
    print(f"Clique ID: {cliqid}")
    for row in bedpe_rows:
        print(row)
    print("\n")  # Adding a new line for readability

    # Break the loop after printing the desired number of entries
    if i + 1 == num_entries_to_print:
        break
'''

##########################################################
# Part 3:

##########################################################
'''

def generate_vcf(cliqset, n):
    
    # cliqset:    input list of clique SVs. 
    # n:          increasing number to to ID field
    # Merge read1 and read2 of a pair.
    # combine all read ids and all INFO fields into one aggregate
    # return two lists with read1 and read2
    
    chrom1set = set()
    chrom2set  = set()
    strand1set  = set()
    strand2set = set()
    ref1array = np.array([])
    ref2array = np.array([])
    pass1set = set()
    pass2set = set()
    qual1list = list()
    qual2list = list()
    pos1array  = np.array([])
    pos2array  = np.array([])
    alt1array = np.array([])
    alt2array = np.array([])
    index_alt = 0
    info1list = list()
    info2list = list()
    readnamelist = list()
    svcaller = set() #I'm a little confused because this gets filled with info from the ID column
    svid_1 = "SVMERGE" + str(n) + "_1"
    svid_2 = "SVMERGE" + str(n) + "_2"

    for row in cliqset: #Cliqset is a list of lists where each 'sublist' is a bedpe row. Each k is a cliq id. Cliqset typically has >1 rows, each corresponding to a different caller.
        svcaller.add ( row[2].split('_')[0]) #This line will pull from the ID column (example entry = TRA00043279_2) and will extract everything before the _.
        matechrom, matepos = get_mate(row[4]) #This line will take a complex string (example entry = [21:38498961[A) and will subset out the chr and pos and store as two values.
        infofull  = row[7] #This line will extract the data in the VCF INFO field
        chrom = row[0] #This line will extract the VCF chromosome field
        pos = int(row[1]) #This line will extract the VCF position and will convert to int

        # info = infofull.split(';')

        #This loop will look for the listed strings in the INFO column
        #If it finds the string, it will subset out that information from the INFO column
        #Example output = 'READ_ID=ABC'
        #I don't think we need this since we don't have an INFO column to begin with.
        for rn in ['READNAMES', 'READ_ID', 'TRDS', 'TSRDS']:
            if rn in row[7]:
                readnamelist += re.sub(r'.*' + rn + '=([^;]*).*', r'\1', row[7]).split(',') 

        #This code will remove all elements from 'rn' from the INFO string 
        #This will replace the original INFO string (infofull) in memory
        infofull = re.sub(r'(.*);' + rn + '=[^;]*;(.*)', r'\1;\2', infofull)

        #I think this code will take a string like this ("0/0:0.0,-12.3422,-246.0:123:PASS:430:0:0:41:0") and will:
        #First, split on ':', then figure out if the final chunk contains a ','.
        #If the above condition is true, it will add that value to readnamelist
        #I don't believe we care about this, because we don't have this information in our bedpe
        if len(row[-1].split(':')[-1].split(','))>1:
            readnamelist += row[-1].split(':')[-1].split(',')


        # Always take the smallest chrom or pos as the first pair
        ####note: It seems like all this is doing is adding a bunch of VCF information to variables
        ####note: There are two arrays/lists/etc for each type of information
        ####note: Information gets added to each array/list/etc depending on the chrom vs matechrom value  
        #This code block will first check for two conditions
        #If "mate chrom" is larger than "chrom" (for example, check if chr'12' is larger than chr'2')
        #If "mate chrome" is the same as "chrome", then do a similar comparison for the position
        if chrom < matechrom or (chrom == matechrom and pos < matepos):
            strand1, strand2 = read_orient(row[4]) #This function will take a string (e.g. [21:38498961[A) and depending on the orientation of the brackets it will determine if strand1 and strand2 are the - or + DNA strands
            strand1set.add(strand1) #Add strand1 to a set. Why? Just error checking?
            strand2set.add(strand2) #Add strand2 to a set. Why? Just error checking?
            chrom1set.add(row[0]) #Add the chromosome to a set. Why? Just error checking?
            ref1array = np.append(ref1array, row[3]) #Add the ref allele (e.g. A, T, C, or G) to a np array

            try:
                qual1list.append(str(row[5])) #Add the quality score to a list

            except Exception as e:
                pass

            alt1array = np.append(alt1array, row[4]) #Add the ALT allele string (e.g. [21:38498961[A) to a numpy array
            pos1array = np.append(pos1array , int(row[1])) #Add the VCF POS value to a numpy array
            info1list += infofull.split(';') #Split the INFO string (e.g. 'CIEND=-196,196;CIPOS=-196,196;PE=5') then add to a list

        else: #In cases where chrom > matechrome or chrom == matechrome and pos > matepos
        # elif row[2].endswith("2"):
            chrom2set.add(chrom) #Add the chromosome to a set
            ref2array = np.append(ref2array, row[3]) #Add the VCF REF value (e.g. A, T, C, or G) to an array

            try:
                qual2list.append(str(row[5])) #Add the VCF QUAL score to a list.

            except Exception as e:
                pass

            alt2array = np.append(alt2array, row[4]) #Add the alt allele string (e.g. '[21:38498961[A') to a numpy array 
            pos2array = np.append(pos2array , int(row[1])) #Add the VCF POS value to a numpy array
            info2list += infofull.split(';') #Split the INFO string (e.g. 'CIEND=-196,196;CIPOS=-196,196;PE=5') then add to a list 


    #This should always be true if things are working well. I don't know what type of situation would make this false.
    if len(strand1set) * len(strand2set) * len(chrom1set) * len(chrom2set) * len(pos1array)/len(pos2array) !=1:
        print ("error") #If it is false, then alert that there is an error
        
        #If it is an error, then add the troublemaker rows from the vcf file to a .txt file named unresolved.txt
        with open('unresolved.txt', 'a') as rin:
            for row in v:
                rin.writelines(','.join(map(str, v)) + "\n")
        return (False, False) #Return False, False as the ultimate output of this function (instead of a VCF row)
    

    chrom1 = ''.join(map(str, list(chrom1set))) #Convert chrom1set to a list, then convert each element to a string. concat all elements directly next to each other
    chrom2 = ''.join(map(str, list(chrom2set))) #Do the same thing for chr2
    strand1 = list(strand1set)[0] #Convert strand1set to a list, then take only the first element of the list
    strand2 = list(strand2set)[0] #Do the same thing for strand2set
    svclass = sv_class(chrom1, strand1, chrom2, strand2) #This function will look at characteristics of the SV and define what class it is (e.g. 'TRA'). We already have this info in our bedpe
    svmethod = '_'.join(map(str, list(svcaller))) #This will convert the stripped VCF ID column value (e.g. TRA00043279), convert to a list, convert to a string, then join them all together in a single string with a '_' sep
    svmethod = svmethod.replace("SANGER", "BRASS") #If it says SANGER, then replace it with "BRASS"
    svinfo = ['SVCLASS='+ svclass, 'SVMETHOD=' + svmethod] #Create two strings (the structure is something like 'SVCLASS=TRA' 'SVMETHOD=BRASS')
    
    # check if all agree. Then take the call with the longest nt stretch in ALT
    #This function will take a numpy array as input and will output a single index value
    #The function will loop through the input numpy array, will calculate the 'itemsize' and will store as a list.
    #argmax() will output the index of the largest alt
    #I'm really not sure why we do this
    def get_longest_alt(altarray):
        return np.array([i.itemsize for i in altarray]).argmax()

    #Basically, I think what the next four code blocks are doing are calculating the most extreme (max or min) position for each sv
    #Then this value is set as the single merge pos because we want to output only a single line for each sv
    #This should be straightforward to implement from the bedpe
    if ''.join(map(str, strand1set)) == '+': #If strand1 is on the plus strand, then
        pos1_index = pos1array.argmax() #Get the index of the largest value in pos1array
        mergepos1 = int(max(pos1array)) #Calculate the maximum value in pos1array, then convert to int

    elif ''.join(map(str, strand1set)) == '-': #If strand1 is on the minus strand, then
        pos1_index = pos1array.argmin() #Get the index of the smallest value in pos1array
        mergepos1 = int(min(pos1array)) #Get the minimum value in pos1array, then convert to int

    if ''.join(map(str, strand2set)) == '+': #If strand2 is on the plus strand, then
        pos2_index = pos2array.argmax() #Get the index of the largest value in pos2array
        mergepos2 = int(max(pos2array)) #Calculate the maximum value in pos2array, then convert to int

    elif ''.join(map(str, strand2set)) == '-': #If strand2 is on the minus strand, then
        pos2_index = pos2array.argmin() #Get the index of the smallest value in pos2array
        mergepos2 = int(min(pos2array)) #Calculate the minimum value in pos2array

    
    if len(list(set(pos1array))) == 1: #If there is only one unique element in pos1array, then
        pos1_index = get_longest_alt(alt1array) #Get the index of the largest value (in bytes) in alt1array


    if len(list(set(pos2array))) == 1: #If there is only one unique element in pos2array, then
        pos2_index = get_longest_alt(alt2array) #Get the index of the largest value (in bytes) in alt2array

    #I have no idea why we're doing this
    mergeref1 = ref1array[pos1_index] #ref1array should be an array of bases (e.g A, T, C, G). Get the value at the index of the largest "ALT"
    mergealt1 = alt1array[pos2_index].replace('chr', '') #subset alt1array to the value with the largest index in alt2array. Remove the 'chr' substring
    mergeref2 = ref2array[pos2_index] #subset ref2array to only include the index with the largest value in alt2array
    mergealt2 = alt2array[pos2_index].replace('chr', '') #Do something similar here

    #Loop through all quality scores. Set the merge quality score to the max quality score
    #If there are no int quality scores, then set the merge quality score to be equal to "."
    qual =  [int(i) for i in qual1list if i !='.']
    if qual:
        mergequal = max([int(i) for i in qual1list if i !='.'])
    else:
        mergequal = "."
    

    readname_uniq = sorted(list(set(re.sub(r'/[12]', r'',i) for i in readnamelist))) #Remove the trailing '1' and '2' from the strings, then return a sorted list of unique strings
    mergereadname = 'READ_ID=' + ','.join(map(str, readname_uniq)) #Convert to a string that looks something like this "READ_ID=0,1,2"

    #Take all of the stuff in info1 list, and remove all of the junk (defined in key_remove)
    #Save this as a new list named merge1info_raw
    key_remove= ['MATE','IMPRECISE','SOMATIC','SVMETHOD','SVCLASS','STRAND', 'PE', 'MATEID', 'MATEPOS', 'MATECHROM', 'MATESTRAND', 'MATEPOS', 'MATEMAPQ', 'SCTG', 'SPAN', 'TSPLIT','MAPQ', 'SVLEN', 'INSLEN', 'CONTROL'] #Define some terms (to remove?)
    merge1info_raw = list(set([i for i in info1list if not i.split('=')[0] in key_remove])) 

    #Defining a function inside of a function. Why?
    def merge_field(mergeinfo_raw, svinfo):
        homseqs = list()
        mergeinfo_tmp = list()
        for i in mergeinfo_raw:
            if i.startswith('HOMSEQ'):
                homseqs.append( i.split('=')[1])
            else:
                mergeinfo_tmp.append(i)
        mergeinfo_tmp.append('HOMSEQ=' + '|'.join(map(str, homseqs)))
        mergeinfo_tmp = mergeinfo_tmp + svinfo
        return mergeinfo_tmp

    #Using the function above:
    #1) Look for something called HOMSEQ in the INFO string, and if it exists, take the value associated with that entry and add to the 'homeseqs' list
    #2) If you can't find HOMSEQ, then append the element from mergeinfo_raw to mergeinfo_tmp
    #3) take mergeinfo_tmp (which shouldn't have a HOMSEQ value in it) and add the HOMSEQ value to the end.
    #4) Add the sv info to the list. svinfo has this structure: [SVCLASS=TRA, SVMETHOD=BRASS]
    merge1info_raw = merge_field(merge1info_raw, svinfo)

    #Create a list named mate1info that has more information about the chr, svid, strand, etc
    #Then create a list with only one element. I believ this list is just a megastring for the INFO column
    mate1info = ['MATECHROM='+ chrom2, 'MATEPOS=' + str(mergepos2), 'MATEID=' + svid_2, 'STRAND=' + strand1, 'MATESTRAND='+strand2, 'PE=' + str(len(readname_uniq))]
    merge1info = [';'.join(map(str, sorted(merge1info_raw + mate1info) + [mergereadname]))]

    #See above when merge_field() was invoked
    merge2info_raw = list(set([i for i in info2list if not i.split('=')[0] in key_remove]))
    merge2info_raw = merge_field(merge2info_raw, svinfo)

    #Do something similar for mate2
    mate2info = ['MATECHROM='+ chrom1, 'MATEPOS=' + str(mergepos1), 'MATEID=' + svid_1, 'STRAND=' + strand2, 'MATESTRAND=' + strand1, 'PE=' + str(len(readname_uniq))]
    merge2info = [';'.join(map(str, sorted(merge2info_raw + mate2info) + [mergereadname]))]
    
    
    ##  CHROM         POS      REF        ALT
    ## mergechrom   mergepos mergeref   mergealt
    ## todo: extract all genotype info from all callers and merge
    genotypecol = ['GT', '0/0', '0/1'] #Construct a list with genotype information
    ##

    #Construct some lines in the "VCF" format, then return them outside of the function.
    pos1info = [chrom1, mergepos1, svid_1, mergeref1, mergealt1, mergequal, "PASS"]
    pos2info = [chrom2, mergepos2, svid_2, mergeref2, mergealt2, mergequal, "PASS"]
    vcf1line = pos1info + merge1info + genotypecol
    vcf2line = pos2info + merge2info + genotypecol
    #print (svid_1, svid_2, svcaller)
    return (vcf1line, vcf2line)


#This function will take a string (example entry = [21:38498961[A) and will extract the chromosome and the position
def get_mate(record):
    matepos = int(re.sub(r'.*[^\:]:([0-9]*).*', r'\1' ,str(record)))
    matechrom = re.sub(r'[A-Z]*[\[\]]*([^:]*):.*', r'\1', record)
    # re.sub(r'.*([^\:]):([A-Z0-9]*).*', r'\1' ,str(record))
    return [matechrom.replace('chr', ''), matepos]


#The generate_vcf function returns two lines (lists) for each SV. Each of the two lines corresponds to each half of the breakpoint.
#This function takes a single list (a) and uses it to make a bedpe file

def generate_bedpe(a):
    chrom1, start1, end1 = a[0], a[1], int(a[1])+1 #Straightforward, just extract the chr and pos from the list, then add 1 to make a new end1 value
    #chrom2, start2 = re.sub(r'.*[\]\[][a-z]*([0-9]*):([0-9]*).*', r'\1,\2' ,a[4]).split(',') #This was already commented out
    chrom2, start2 = get_mate(a[4]) #This function (defined above) looks at the complex input string and extracts chr and pos for the mate
    end2 = int(start2)+1 #Add 1 to generate a new 'end2' value
    name = a[2].split("_")[0] #Take a string that looks like this "TRA00043279_2" and only keep everything before the _
    score = re.sub(r'.*;PE=([0-9]*).*', r'\1', a[7]) #Extract the indicated information from the INFO string
    strand1 = re.sub(r'.*;STRAND=([+-]).*', r'\1', a[7]) #Extract the indicated information from the INFO string
    strand2 = re.sub(r'.*;MATESTRAND=([+-]).*', r'\1', a[7]) #Extract the indicated information from the INFO string
    svclass = re.sub(r'.*;SVCLASS=([A-Za-z0-9]*).*', r'\1', a[7]) #Extract the indicated information from the INFO string
    svmethod = re.sub(r'.*;SVMETHOD=([A-Za-z0-9_]*);.*', r'\1', a[7]) #Extract the indicated information from the INFO string
    return [chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, svclass, svmethod] #Return a list with all of the info in the bedpe format
'''

#This function should replicate everything useful that the generate_vcf function does, except output bedpe format (like generate_bedpe does)
#
def generate_bedpe_new(input_cliq_bedpe): #Use mastermerge as input and output a collapsed bedpe
    
    output_bedpe=pd.DataFrame(columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'sv_id', 'pe_support', 'strand1', 'strand2', 'sv_class', 'svmethod', 'center' ]) #This is the bedpe that will be output from this function

    for key, sublist in input_cliq_bedpe.items():

            #init some objects
            #can rename later. for now I am naming them this way to make debugging less confusing for me.
            new_chrom1_set=set()
            new_chrom2_set=set()
            new_strand1_set=set()
            new_strand2_set=set()
            new_svtype_set=set()
            new_start1_list=list()
            new_start2_list=list()
            new_caller_set=set() #caller = the algorithms making the call. Note, in some cases one group uses multiple callers
            new_center_set=set() #center = institute (e.g broad, nygc, etc)
            new_qual_list=list() #Every group reports qual in a different way. Does anyone care about qual? Let's just assume that all groups pre-filter for their own qual threshold. I do not want to parse these strings.
            new_id_list=list()

            final_chrom1=None #This is the final chrom1 value for the output bedpe
            final_chrom2=None #This is the final chrom2 value for the output bedpe
            final_strand1=None #This is the final strand1 value for the output bedpe
            final_strand2=None #This is the final strand2 value for the output bedpe
            final_svtype=None #This is the final svtype for the output bedpe
            final_start1=None #This is the final start1 position for the output bedpe
            final_start2=None #This is the final start2 position for the output bedpe
            final_end1=None #This is the final end1 value for the output bedpe
            final_qual_string=None #This is the final qual string for the output bedpe, will be ";" sep and order will be same as final_center_string
            final_caller_string=None #This is the final caller string for the output bedpe, will be ";" sep.
            final_center_string=None #This is the final center string for the output bedpe, will be ";" sep and order will be the same as final_qual_string
            final_id_string=None #This is a list of all SV IDs for the output bedpe, will be ";" sep.
            final_bedpe_line=None #This is the final bedpe line (everything glued together) that will be output for each cliq

            for item in sublist: #And for each sv call (one from each institution) for that cliq

                #Add all of the chrom1 and chrom2 values to a set.
                #We have downstream ambitions of checking to make sure for each cliq these values are always the same across all callers
                #Do something similar for strand
                #Also do something similar for sv type
                new_chrom1_set.add(item[0]) #0 is the index of chrom1 in bedpe
                new_chrom2_set.add(item[3]) #3 is the index of chrom2 in bedpe
                new_strand1_set.add(item[8]) #8 is the index for strand1 in bedpe
                new_strand2_set.add(item[9]) #9 is the index for strand2 in bedpe
                new_svtype_set.add(item[10]) #10 is the index for sv type in bedpe

                #Add all of the start1 and start2 values to a list
                #We don't care about end1 and end2 because end1/2 == start1/2 +1, so we can just regenerate later.
                new_start1_list.append(item[1]) #1 is the index for start1 in bedpe
                new_start2_list.append(item[4]) #4 is the index for start2 in bedpe

                #Add all of the callers and centers to a set
                new_caller_set.update(item[11].split(',')) #Sometimes a single group uses multiple callers. Add each indivudal caller to the set
                new_center_set.add(item[12]) #This the index of the center (e.g. broad or nygc)

                #Add all of the qual scores to a list
                #Also add all of the sv IDs to a list
                new_qual_list.append(item[7]) #7 is the inde of the qual score in bedpe
                new_id_list.append(item[6]) #6 is the index of the SV ID in bedpe


        
            #Now that everything is assembled, calculate the values that will go into the actual bedpe
            #We can write this in a more beautiful way in the future, for now writing it like this may help with debugging.
            #to do: in the future just have it error and break if any of the elif conditions are met. Right now it is written this way for debugging purposes.
            #First, chrom1
            if len(new_chrom1_set) == 1:
                final_chrom1=new_chrom1_set.pop()
            elif len(new_chrom1_set) > 1:
                final_chrom1="error_multiple_chrom1"
            else:
                final_chrom1="error_no_chrom1"

            #Now chrom2
            if len(new_chrom2_set) == 1:
                final_chrom2=new_chrom2_set.pop()
            elif len(new_chrom2_set) > 1:
                final_chrom2="error_multiple_chrom2"
            else:
                final_chrom2="error_no_chrom2"

            #Now for strand1
            if len(new_strand1_set) == 1:
                final_strand1=new_strand1_set.pop()
            elif len(new_strand1_set) > 1:
                final_strand1="error_multiple_strand1"
            elif len(new_strand1_set) == 0:
                final_strand1="error_no_strand1"

            #Now for strand2
            if len(new_strand2_set) == 1:
                final_strand2=new_strand2_set.pop()
            elif len(new_strand2_set) > 1:
                final_strand2="error_multiple_strand2"
            elif len(new_strand2_set) == 0:
                final_strand2="error_no_strand2"

            #Now for sv type
            if len(new_svtype_set) == 1:
                final_svtype=new_svtype_set.pop()
            elif len(new_svtype_set) > 1:
                final_svtype="error_multiple_svtype"
            elif len(new_svtype_set) == 0:
                final_svtype="error_no_svtype"

            #Calculate the position values using the same logic as the original generate_vcf function
            if final_strand1 == '+':
                final_start1=max(new_start1_list)
                final_end1=int(final_start1)+1
            elif final_strand1 == '-':
                final_start1=min(new_start1_list)
                final_end1=int(final_start1)-1
        
            if final_strand2 == "+":
                final_start2=max(new_start2_list)
                final_end2=int(final_start2)+1
            elif final_strand2 == '-':
                final_start2=min(new_start2_list)
                final_end2=int(final_start2)-1

            #Construct a string for qual, caller, center, and id
            final_qual_string=';'.join(new_qual_list)
            final_caller_string=';'.join(new_caller_set)
            final_center_string=';'.join(new_center_set)
            final_id_string=';'.join(new_id_list)

            #Generate a final bedpe line with all of the information above
            final_bedpe_line=[final_chrom1, final_start1, final_end1, final_chrom2, final_start2, final_end2, final_id_string, final_qual_string, final_strand1, final_strand2, final_svtype, final_caller_string, final_center_string]

            #Add it to a pandas dataframe
            output_bedpe.loc[len(output_bedpe)] = final_bedpe_line

    #Now return the bedpe :)
    return output_bedpe


collapsed_bedpe=generate_bedpe_new(mastermerge)
            

'''
# fopen = open('unresolved.txt', 'w')
# fopen.close()
# vcflines = list()
bedpelines = list()
n = 0
for k, v in mastermerge.items():
    #print (k)
    if k:
        bedpelines = [*bedpelines, *v] # TO SELF: COME BACK TO THIS SECTION <- Not sure yet if this works (idea is to just concatenate bedpe lines in cliques)
        # n +=1
        # print(v)
        # print(n)
        # a,b = generate_vcf(v, n)
        # # print(a)
        # if a and b:
        #     vcflines.append(a)
        #     vcflines.append(b)
        #     bedpe = generate_bedpe(a)
        #     bedpelines.append(bedpe)

print(bedpelines)

'''


# vcflines.sort(key = lambda x: (x[0], int(x[1])))
#collapsed_bedpe.sort(key = lambda x: (x[0], int(x[1]))) #Comment this out. Sorting is for cowards.

print("This is the structure of the final collapsed bedpe")
pd.set_option('display.max_columns', None)
print(collapsed_bedpe.head())



## Search for high number of small SVs - likely artefacts
collapsed_bedpe['sv_class'] = collapsed_bedpe['sv_class'].str.replace('h2h', '').str.replace('t2t', '') #Rewrite for pandas df input
statCounterTmp = Counter(collapsed_bedpe['sv_class']) #Rewrite for pandas df input
#statCounterTmp = Counter([i[10].replace('h2h', '').replace('t2t', '') for i in collapsed_bedpe]) #In a string like 'h2hINV' remove the h2h so it is only INV
statCounter = dict() 
bins=[0, 5e3, 1e12] #Where did these numbers come from?
sumSVs= sum(statCounterTmp.values()) #Count the number of times each type of sv class was observed


for k,v in statCounterTmp.items():
    statCounter[k + '_count'] = v
    statCounter[k + '_count_fract'] = round(v/sumSVs,3)
    statCounter[k + '_size_hist'] = np.nan
    #if k != 'TRA':        
        #statCounter[k + '_size_hist'] = np.histogram([(i[4]-i[1])  for i in collapsed_bedpe if k in i[10] and i[0] == i[3]], bins=bins)

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
    #df_filter.columns = [pid]


#outFILTER = outVCF.replace('.vcf', '_sizeStat.txt') #Who cares?
#df_filter.to_csv(outFILTER, sep="\t", na_rep="NA") #Who cares?


#vcflines_filter = list()
collapsed_bedpe_filter = list()
if mergeIdRemove:
    #for vcf in vcflines:
        #if not re.sub(r'_[12]$', '',vcf[2]) in mergeIdRemove:
           #vcflines_filter.append(vcf)
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

''' #I modified this code, so I don't trust it. This was the original method for writing the bedpe
def write_bedpe(bedpefileOut, collapsed_bedpe):
    with open(bedpefileOut, 'w') as bedpeout:
        writer = csv.writer(bedpeout, delimiter="\t", lineterminator="\n")
        writer.writerow(bedpeheader)
        for bed in collapsed_bedpe:
            writer.writerow(bed)


if collapsed_bedpe_filter:
    bedpefileOut = outBEDPE.replace('.bedpe', '_raw.bedpe')
    write_bedpe(bedpefileOut, collapsed_bedpe)
    write_bedpe(outBEDPE, collapsed_bedpe_filter)
else:
    bedpefileOut = outBEDPE
    write_bedpe(bedpefileOut, collapsed_bedpe)

'''

print ("merged SVs into BEDPE format\n", collapsed_bedpe, "\n\n")



####################################################################################################
## HEADER
####################################################################################################
# fileformat

#####For now just comment this out. Do we even need to write a vcf output?
'''
merge_header = list()
fileformat = ["##fileformat=VCFv4.1"] #Edit this
filedate = ["##fileDate=" + today] #Edit this
codesource = ["##pcawg6_sv_merge=brass(Sanger),delly(EMBL),dranger(Broad),snowman(Broad)"] #Edit this  
genomeref = ["##reference=hs37d5,ftp://ftp.sanger.ac.uk/pub/project/PanCancer/"] #Edit this
merge_header += [fileformat]
merge_header += [filedate]
merge_header += [genomeref]
merge_header += [codesource]
header_clean = [list(j) for j in set(tuple(i) for i in header_concat) if not '##fileformat' in j[0] and not '##fileDate' in j[0] and not '##reference' in j[0] and not '##source' in j[0]]
header_clean.sort()
formats = list()
filters = list()
infos = list()
alts = list()
chroms = list()
headermisc = list()
for h in header_clean:
    hstring = ' '.join(map(str, h))
    if h[0].startswith("##FORMAT"):
        formats.append(h)
    elif h[0].startswith("##FILTER"):
        filters.append(h)
    elif h[0].startswith("##INFO"):
        infos.append(h)
    elif h[0].startswith("##ALT"):
        alts.append(h)
    elif h[0].startswith("#CHROM"):
        chroms.append(h)
    else:
        headermisc.append(h)
chromline = chroms[0][:9] + ['NORMAL', 'TUMOUR']

info_dict = defaultdict(list)
infos_pruned = list()
for i in infos:
    k = re.sub(r'##INFO=<ID=([^,]*)(.*)', r'\1', i[0]) 
    info_dict[k].append(i)
for k, v in info_dict.iteritems():
    # always take the longest string
    vprune = sorted(v, key=len)[-1]
    infos_pruned.append(vprune)

filter_dict = defaultdict(list)
filters_pruned = list()
for i in filters:
    k = re.sub(r'##FILTER=<ID=([^,]*)(.*)', r'\1', i[0]) 
    filter_dict[k].append(i)
for k, v in filter_dict.iteritems():
    # always take the longest string
    vprune = sorted(v, key=len)[-1]
    filters_pruned.append(vprune)

format_dict = defaultdict(list)
formats_pruned = list()
for i in formats:
    k = re.sub(r'##FORMAT=([^,]*)(.*)', r'\1', i[0]) 
    format_dict[k].append(i)
for k, v in format_dict.iteritems():
    # always take the longest string
    vprune = sorted(v, key=len)[-1]
    formats_pruned.append(vprune)


headermisc_dict = defaultdict(list)
headermisc_pruned = list()
for i in headermisc:
    k = re.sub(r'##contig=<ID=([^,]*),(.*)', r'\1', i[0])
    headermisc_dict[k].append(i)
for k, v in headermisc_dict.iteritems():
    # always take the longest string
    vprune = sorted(v, key=len)[-1]
    headermisc_pruned.append(vprune)

headermisc_pruned.sort()

merge_header += headermisc_pruned
merge_header += infos_pruned
merge_header += filters_pruned
merge_header += formats_pruned
merge_header += alts

# merge_header.append(chromline)

def write_vcf(outVCF, vcflines):
    with open(outVCF, 'w') as vcfout:
        for row in merge_header:
            row = ' '.join(map(str, row)) + "\n"
            vcfout.writelines(row )
        row = '\t'.join(map(str, chromline)) + "\n"
        vcfout.writelines(row)
        for vcfline in vcflines:
            row = '\t'.join(map(str, vcfline)) + "\n"
            vcfout.writelines(row)
    print ("merged SV VCF file\n", outVCF, "\n\n")


if vcflines_filter:
    write_vcf(outVCF.replace('.vcf', '_raw.vcf'), vcflines)
    write_vcf(outVCF, vcflines_filter)
else:
    write_vcf(outVCF, vcflines)
'''