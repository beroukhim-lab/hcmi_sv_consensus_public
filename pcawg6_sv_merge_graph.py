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



# Parse command line
parser=argparse.ArgumentParser(description='Graph builder.')
parser.add_argument('-a', '--inVCF_BRASS', metavar='brass_in.vcf', required=False, dest='inVCF_BRASS', help='input vcf file (optional)')
parser.add_argument('-b', '--inVCF_DELLY', metavar='delly_in.vcf', required=False, dest='inVCF_DELLY', help='input vcf file (optional)')
parser.add_argument('-c', '--inVCF_DRANGER', metavar='dranger_in.vcf', required=False, dest='inVCF_DRANGER', help='input vcf file (optional)')
parser.add_argument('-d', '--inVCF_SNOWMAN', metavar='snowman_in.vcf', required=False, dest='inVCF_SNOWMAN', help='input vcf file (optional)')

parser.add_argument('-e', '--inBEDPE', metavar='overlap.bedpe', required=True, dest='inBEDPE', help='input overlap table (required)')
# parser.add_argument('-u', '--unionDEL', metavar='uniondel.bed', required=False, dest='uDEL', help='union deletions merged calls (optional)')
parser.add_argument('-s', '--stat', metavar='stat.bed', required=False, dest='outStat', help='output statistics (optional)')
#parser.add_argument('-c', '--copynumberConcordance', metavar='0.5', required=True, dest='copyConc', help='required copynumber concordance (required)')
#parser.add_argument('-r', '--reciprocalOverlap', metavar='0.5', required=True, dest='recOver', help='required reciprocal overlap (required)')
#parser.add_argument('-m', '--minCarrier', metavar='5', required=True, dest='minCarrier', help='required minimum number of carrier samples (required)')
parser.add_argument('-t', '--triangle', dest='useTriangle', action='store_true', default=False, help='use upper triangle')
parser.add_argument('-p', '--plot', dest='plotSubgraph', action='store_true', default=False, help='plot subgraphs')
parser.add_argument('-o', '--outVCF', metavar='out.vcf', required=False, dest='outVCF', help='output vcf file (optional)')
parser.add_argument('-y', '--priority', metavar='priority.txt', required=False, dest='inPrior', help='priority scores (optional)')
args = parser.parse_args()


# Parse command-line parameters

#recOver=0.5
minCarrier=2
#if args.recOver:
#    recOver=float(args.recOver)
udel=collections.defaultdict(set)

priority=collections.defaultdict(float)
if args.inPrior:
    f_reader=csv.reader(open(args.inPrior), delimiter="\t")
    for fields in f_reader:
        if (fields[1]!="NA"):
            priority[fields[0]] = float(fields[1])


inBEDPE = args.inBEDPE
useTriangle = args.useTriangle
plotSubgraph = args.plotSubgraph
outStat = args.outStat
outVCF = args.outVCF
copyConc = 0

# Overlap graph

inVCF_BRASS = args.inVCF_BRASS
inVCF_DELLY = args.inVCF_DELLY
inVCF_DRANGER = args.inVCF_DRANGER
inVCF_SNOWMAN = args.inVCF_SNOWMAN

pid = os.path.basename(os.path.dirname(inBEDPE))
G=networkx.Graph()
# Parse bed fileSV padding
f_reader=csv.DictReader(open(inBEDPE), delimiter="\t")
for row in f_reader:
    pid = row['pid']
    #if (float(row['interSect'])>0.5):
    if not G.has_edge(row['name_SV1'], row['name_SV2']):
        G.add_edge(row['name_SV1'], row['name_SV2'], weight=float(row['rn_share_fract']))
        G.add_edge(row['name_SV1'], row['name_SV2'], weight=float(row['rn_share_fract']))
        svType1 = row['svtype_SV1']
        if row['svtype_SV1'] == "NA":
            svType1 = "BND"
        svType2 = row['svtype_SV2']
        if row['svtype_SV2'] == "NA":
            svType2 = "BND"
        G.node[row['name_SV1']]['SVtype']=(svType1, row['pid'], row['center_SV1'],row['idSV1'])
        G.node[row['name_SV2']]['SVtype']=(svType2, row['pid'], row['center_SV2'],row['idSV2'])




# Parse graph structure
compAssign=dict()
compMembers=collections.defaultdict(list)
conCompCount=0
cliqueCount=0
for H in networkx.connected_component_subgraphs(G):
    conCompCount+=1
    typeSet=set()
    sample=set()
    center=list()
    for n, d in H.nodes_iter(data=True):
        typeSet.add(d['SVtype'][0])
        sample.add(d['SVtype'][1])
        center.append(d['SVtype'][2])
    baseDir = 'plots/' + pid + '/'
    try:
        os.makedirs(baseDir)
    except Exception, e:
        pass
    baseName=  '_'.join(list(sample) + sorted(typeSet)) + '.' + str(conCompCount) + '_' + ';'.join(map(str, list(set(center)))) + '.png' 
    #print (baseName, n,d)
    if plotSubgraph:
        #labels=nx.draw_networkx_labels(H,pos=nx.spring_layout(H))
        #nx.draw_networkx_nodes(H,pos=nx.spring_layout(H),node_size=700)
        elarge=[(u,v) for (u,v,d) in H.edges(data=True) if d['weight'] >0.3]
        esmall=[(u,v) for (u,v,d) in H.edges(data=True) if d['weight'] <=0.3]
        pos=nx.circular_layout(H) # positions for all nodes
        # nodes
        ##centercolor = [i.replace('sanger', '#ADD8E6').replace('broad-dRanger', '#F0E68C').replace('embl-delly', '#90EE90') for i in center]
        centercolor = [i.replace('sanger', '#ADD8E6').replace('embl', '#90EE90').replace('dRanger', '#F0E68C').replace('smufin', '#800080').replace('snowman', '#FFA500') for i in center]

        nx.draw_networkx_nodes(H,pos,node_size=700 , node_color=centercolor, font_family="arial")
        # edges
        weight = H.edges(data=True)[-1][2]['weight']
        nx.draw_networkx_edges(H,pos,edgelist=elarge,
                            width=2+4*weight)
        nx.draw_networkx_edges(H,pos,edgelist=esmall,
                            width=2+4*weight,alpha=0.5,edge_color='b',style='dashed')
        plotGraph(H, baseName)
    for n in H.nodes_iter(data=False):
         compAssign[n]=conCompCount
         compMembers[conCompCount].append(n)

    # Clique?
    if (float(H.number_of_edges())==float(H.number_of_nodes()*(H.number_of_nodes() - 1)) / 2.0):
        cliqueCount+=1
    else:
        pass

# Print summary
if outStat:
    bedOut = open(outStat, 'w')
    print('copyConc', 'minCarrier', 'number_of_nodes', 'connected_components_Count', 'clique_Count', sep="\t", file=bedOut)
    print(copyConc, minCarrier, G.number_of_nodes(), conCompCount, cliqueCount, sep="\t", file=bedOut)
    bedOut.close()


cliqueFile = inBEDPE.replace('.txt','.clique.txt')
with open(cliqueFile, 'w') as w:
    writer = csv.writer(w, delimiter="\t", lineterminator="\n")
    for k,v in compAssign.iteritems():
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

# Merge calls
# read in all the vcf files into a hash table. Cannot merge due to too many different formats.

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

def get_mate(record):
    matepos = int(re.sub(r'.*[^\:]:([0-9]*).*', r'\1' ,str(record)))
    matechrom = re.sub(r'[A-Z]*[\[\]]*([^:]*):.*', r'\1', record)
    # re.sub(r'.*([^\:]):([A-Z0-9]*).*', r'\1' ,str(record))
    return [matechrom.replace('chr', ''), matepos]


def mergeinfo(d1, d2):
    return reduce(lambda a,b: dict (a, **b), (d1,d2))


def extract_vcf_rows(inVCF, mastermerge, header_concat, compAssign):
    svcaller = inVCF.split('/')[-2]
    with gzip.open(inVCF, 'rb') as rin:
        for f in rin:
            row = f.rsplit()
            rinfos = dict()
            if row[0].startswith("#"):
                if row[-1].endswith('>'):
                    row[-1] = re.sub(r'">$', ',tool=' + svcaller + '">', row[-1])
                header_concat.append(row)
            else:
                if re.search(r'_[12]$', row[2]):
                    svid = re.sub(r'_[12]$', '',row[2])
                elif re.search(r':[12]$', row[2]):
                    svid = re.sub(r':[12]$', '',row[2])
                else:
                    break
                if svid in compAssign:
                    cliqid = compAssign.get(svid)
                    mastermerge[cliqid].append(row)
    return mastermerge, header_concat




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


def generate_vcf(cliqset, n):
    '''
    cliqset:    input list of clique SVs
    n:          increasing number to to ID field
    Merge read1 and read2 of a pair.
    combine all read ids and all INFO fields into one aggregate
    return two lists with read1 and read2
    '''
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
    svcaller = set()
    svid_1 = "SVMERGE" + str(n) + "_1"
    svid_2 = "SVMERGE" + str(n) + "_2"
    for row in cliqset:
        svcaller.add ( row[2].split('_')[0])
        matechrom, matepos = get_mate(row[4])
        infofull  = row[7]
        chrom = row[0]
        pos = int(row[1])
        #print (chrom, pos, matechrom, matepos)
        # info = infofull.split(';')
        for rn in ['READNAMES', 'READ_ID', 'TRDS', 'TSRDS']:
            if rn in row[7]:
                readnamelist += re.sub(r'.*' + rn + '=([^;]*).*', r'\1', row[7]).split(',')
        infofull = re.sub(r'(.*);' + rn + '=[^;]*;(.*)', r'\1;\2', infofull)
        if len(row[-1].split(':')[-1].split(','))>1:
            readnamelist += row[-1].split(':')[-1].split(',')
        # Always take the smallest chrom or pos as the first pair
        if chrom < matechrom or (chrom == matechrom and pos < matepos):
            strand1, strand2 = read_orient(row[4])
            strand1set.add(strand1)
            strand2set.add(strand2)
            chrom1set.add(row[0])
            ref1array = np.append(ref1array, row[3])
            try:
                qual1list.append(str(row[5]))
            except Exception, e:
                pass
            alt1array = np.append(alt1array, row[4])
            pos1array = np.append(pos1array , int(row[1]))
            info1list += infofull.split(';')
        else:
        # elif row[2].endswith("2"):
            chrom2set.add(chrom)
            ref2array = np.append(ref2array, row[3])
            try:
                qual2list.append(str(row[5]))
            except Exception, e:
                pass
            alt2array = np.append(alt2array, row[4])
            pos2array = np.append(pos2array , int(row[1]))
            info2list += infofull.split(';')
    if len(strand1set) * len(strand2set) * len(chrom1set) * len(chrom2set) * len(pos1array)/len(pos2array) !=1:
        print ("error")
        with open('unresolved.txt', 'wa') as rin:
            for row in v:
                rin.writelines(','.join(map(str, v)) + "\n")
        return (False, False)
    chrom1 = ''.join(map(str, list(chrom1set)))
    chrom2 = ''.join(map(str, list(chrom2set)))
    strand1 = list(strand1set)[0]
    strand2 = list(strand2set)[0]
    svclass = sv_class(chrom1, strand1, chrom2, strand2)
    svmethod = '_'.join(map(str, list(svcaller)))
    svmethod = svmethod.replace("SANGER", "BRASS")
    svinfo = ['SVCLASS='+ svclass, 'SVMETHOD=' + svmethod]
    # check if all agree. Then take the call with the longest nt stretch in ALT
    def get_longest_alt(altarray):
        return np.array([i.itemsize for i in altarray]).argmax()

    if ''.join(map(str, strand1set)) == '+':
        pos1_index = pos1array.argmax()
        mergepos1 = int(max(pos1array))
    elif ''.join(map(str, strand1set)) == '-':
        pos1_index = pos1array.argmin()
        mergepos1 = int(min(pos1array))
    if ''.join(map(str, strand2set)) == '+':
        pos2_index = pos2array.argmax()
        mergepos2 = int(max(pos2array))
    elif ''.join(map(str, strand2set)) == '-':
        pos2_index = pos2array.argmin()
        mergepos2 = int(min(pos2array))
    if len(list(set(pos1array))) == 1:
        #print (pos1array)
        pos1_index = get_longest_alt(alt1array)
    if len(list(set(pos2array))) == 1:
        #print (pos2array)
        pos2_index = get_longest_alt(alt2array)
    mergeref1 = ref1array[pos1_index]
    mergealt1 = alt1array[pos2_index].replace('chr', '')
    mergeref2 = ref2array[pos2_index]
    mergealt2 = alt2array[pos2_index].replace('chr', '')
    qual =  [int(i) for i in qual1list if i !='.']
    if qual:
        mergequal = max([int(i) for i in qual1list if i !='.'])
    else:
        mergequal = "."
    readname_uniq = sorted(list(set(re.sub(r'/[12]', r'',i) for i in readnamelist)))
    mergereadname = 'READ_ID=' + ','.join(map(str, readname_uniq))
    key_remove= ['MATE','IMPRECISE','SOMATIC','SVMETHOD','SVCLASS','STRAND', 'PE', 'MATEID', 'MATEPOS', 'MATECHROM', 'MATESTRAND', 'MATEPOS', 'MATEMAPQ', 'SCTG', 'SPAN', 'TSPLIT','MAPQ', 'SVLEN', 'INSLEN', 'CONTROL']
    merge1info_raw = list(set([i for i in info1list if not i.split('=')[0] in key_remove])) 

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

    merge1info_raw = merge_field(merge1info_raw, svinfo)
    mate1info = ['MATECHROM='+ chrom2, 'MATEPOS=' + str(mergepos2), 'MATEID=' + svid_2, 'STRAND=' + strand1, 'MATESTRAND='+strand2, 'PE=' + str(len(readname_uniq))]
    merge1info = [';'.join(map(str, sorted(merge1info_raw + mate1info) + [mergereadname]))]
    merge2info_raw = list(set([i for i in info2list if not i.split('=')[0] in key_remove]))
    merge2info_raw = merge_field(merge2info_raw, svinfo)
    mate2info = ['MATECHROM='+ chrom1, 'MATEPOS=' + str(mergepos1), 'MATEID=' + svid_1, 'STRAND=' + strand2, 'MATESTRAND=' + strand1, 'PE=' + str(len(readname_uniq))]
    merge2info = [';'.join(map(str, sorted(merge2info_raw + mate2info) + [mergereadname]))]
    ##  CHROM         POS      REF        ALT
    ## mergechrom   mergepos mergeref   mergealt
    ## todo: extract all genotype info from all callers and merge
    genotypecol = ['GT', '0/0', '0/1']
    ##
    pos1info = [chrom1, mergepos1, svid_1, mergeref1, mergealt1, mergequal, "PASS"]
    pos2info = [chrom2, mergepos2, svid_2, mergeref2, mergealt2, mergequal, "PASS"]
    vcf1line = pos1info + merge1info + genotypecol
    vcf2line = pos2info + merge2info + genotypecol
    #print (svid_1, svid_2, svcaller)
    return (vcf1line, vcf2line)


mastermerge = defaultdict(list)
header_concat = list()
for inVCF in inVCF_DELLY, inVCF_BRASS, inVCF_DRANGER, inVCF_SNOWMAN:
    mastermerge, header_concat = extract_vcf_rows(inVCF, mastermerge, header_concat, compAssign)




def generate_bedpe(a):
    chrom1, start1, end1 = a[0], a[1], int(a[1])+1
    #chrom2, start2 = re.sub(r'.*[\]\[][a-z]*([0-9]*):([0-9]*).*', r'\1,\2' ,a[4]).split(',')
    chrom2, start2 = get_mate(a[4])
    end2 = int(start2)+1
    name = a[2].split("_")[0]
    score = re.sub(r'.*;PE=([0-9]*).*', r'\1', a[7])
    strand1 = re.sub(r'.*;STRAND=([+-]).*', r'\1', a[7])
    strand2 = re.sub(r'.*;MATESTRAND=([+-]).*', r'\1', a[7])
    svclass = re.sub(r'.*;SVCLASS=([A-Za-z0-9]*).*', r'\1', a[7])
    svmethod = re.sub(r'.*;SVMETHOD=([A-Za-z0-9_]*);.*', r'\1', a[7])
    return [chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, svclass, svmethod]



fopen = open('unresolved.txt', 'w')
fopen.close()
vcflines = list()
bedpelines = list()
n= 0
for k, v in mastermerge.iteritems():
    #print (k)
    if k:
        n +=1
        a,b = generate_vcf(v, n)
        if a and b:
            vcflines.append(a)
            vcflines.append(b)
            bedpe = generate_bedpe(a)
            bedpelines.append(bedpe)

vcflines.sort(key = lambda x: (x[0], int(x[1])))
bedpelines.sort(key = lambda x: (x[0], int(x[1])))


## Search for high number of small SVs - likely artefacts
statCounterTmp = Counter([i[10].replace('h2h', '').replace('t2t', '') for i in bedpelines])
statCounter = dict()
bins=[0, 5e3, 1e12]
sumSVs= sum(statCounterTmp.values())
for k,v in statCounterTmp.iteritems():
    statCounter[k + '_count'] = v
    statCounter[k + '_count_fract'] = round(v/sumSVs,3)
    statCounter[k + '_size_hist'] = np.nan
    if k != 'TRA':        
        statCounter[k + '_size_hist'] = np.histogram([(i[4]-i[1])  for i in bedpelines if k in i[10] and i[0] == i[3]], bins=bins)

print ('\nTotal number of SVs:',sumSVs,'\nDistribution of SV types with bins', ' - '.join(map(str,bins)) )
for k,v in statCounter.iteritems():
    print (k, v)

mergeIdRemove = set()
if sumSVs > 200:
    for k in statCounterTmp.keys():
        size_small_sv = statCounter[k + '_count']
        if k != 'TRA':
            size_small_sv = statCounter[k + '_size_hist'][0][0]
        if (statCounter[k  + '_count_fract']) > 0.75 and size_small_sv/statCounter[k + '_count'] > 0.8:
            print ("Uneven high number of {0}. {1} out of {2} with {3} having a size smaller than 5e3".format(k, statCounterTmp[k], sumSVs, size_small_sv))
            for i in bedpelines:
                if k in i[10] and ((i[0] != i[3]) or (i[0] == i[3] and i[4]-i[1] <  statCounter[k + '_size_hist'][1][1])):
                    mergeIdRemove.add( i[6].replace('_1', '').replace('_2', ''))


d = dict()
df_filter = pd.DataFrame.from_dict(d, orient='index')
if statCounter: 
    for k,v in statCounter.iteritems():
        if 'count' in k:
            d[k] = v
        if 'size_hist' in k:
            ks = k.split('_')[0] + "_count_lt_5kb"
            try:
                d[ks] = v[0][0]
            except Exception, e:
                d[ks] = np.nan    
    #df_filter.columns = [pid]
outFILTER = outVCF.replace('.vcf', '_sizeStat.txt')
df_filter.to_csv(outFILTER, sep="\t", na_rep="NA")


vcflines_filter = list()
bedpelines_filter = list()
if mergeIdRemove:
    for vcf in vcflines:
        if not re.sub(r'_[12]$', '',vcf[2]) in mergeIdRemove:
            vcflines_filter.append(vcf)
    for bed in bedpelines:
        if not re.sub(r'_[12]$', '',bed[6]) in mergeIdRemove:
            bedpelines_filter.append(bed)


bedpeheader = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'sv_id', 'pe_support', 'strand1', 'strand2', 'svclass', 'svmethod']

def write_bedpe(bedpefileOut, bedpelines):
    with open(bedpefileOut, 'w') as bedpeout:
        writer = csv.writer(bedpeout, delimiter="\t", lineterminator="\n")
        writer.writerow(bedpeheader)
        for bed in bedpelines:
            writer.writerow(bed)


if bedpelines_filter:
    bedpefileOut = outVCF.replace('.vcf', '_raw.bedpe')
    write_bedpe(bedpefileOut, bedpelines)
    bedpefileOut = outVCF.replace('.vcf', '.bedpe')
    write_bedpe(bedpefileOut, bedpelines_filter)
else:
    bedpefileOut = outVCF.replace('.vcf', '.bedpe')
    write_bedpe(bedpefileOut, bedpelines)


print ("merged SVs into BEDPE format\n", bedpefileOut, "\n\n")

####################################################################################################
## HEADER
####################################################################################################
# fileformat

merge_header = list()
fileformat = ["##fileformat=VCFv4.1"]
filedate = ["##fileDate=" + today]
codesource = ["##pcawg6_sv_merge=brass(Sanger),delly(EMBL),dranger(Broad),snowman(Broad)"]
genomeref = ["##reference=hs37d5,ftp://ftp.sanger.ac.uk/pub/project/PanCancer/"]
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