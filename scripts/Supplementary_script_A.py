#############################
#                           #
#    Ribosome_Assignment    #
#                           #
#############################

'''
Created on Mon Sep 24 2018
@author: Ulrike Friedrich 

-------------------------------------------------------------------------------
DESCRIPTION

This script performs the general analysis of selective 
riboSeq data. Deep sequencing data must be (quality) trimmed and aligned first.
The input file is a sam file given by Tophat2 performing alignment to the 
genome. Based on the reference files used for genome alignment the function 
chromConversion must be adjusted. 

The user can choose between the modes 'center' for center-weighting and '5-end' 
for 5'-end assignment of the ribosomal A-site. In case of '5-end' mode, a table
(Asite_offset.txt) must be provided defining the offset for each footprint length.
Additionally, the upper and lower border of accepted footprint lengths 
(including limits, default: 23-40) and the minimal number of reads per 
transcript to be included in further analyses (threshold) must be given. 


The pickle file 'Reads.pkl' is a combination of raw and normalized reads
with the following structure: 
key:      chrom --> position (e.g. '01\t0003467')
value1:   raw read + strand  [0]
value2:   '\t'               [1]
value3:   raw read - strand  [2]
value4:   '\t'               [3]
value5:   nucleotid          [4]
value6:   '\t'               [5]
value7:   norm read + strand [6]
value8:   '\t'               [7]
value9:   norm read - strand [8]

By setting the argument text to 'yes' (default: no), the content of this pickle
file is additionally written in a plain text file. 

The text file 'GeneExpression.txt' contains all transcripts, the sum of all 
reads per transcript [RPM], the transcript length, and the normalized sum of 
all reads [RPKM]. The tag 'included' or 'excluded' is added if the RPM-value 
is above or below the defined threshold, respectively. 



'''

import argparse
import pickle
import os
import sys
from time import gmtime, strftime

def chromConversion(chrom): 
    
    chromList = ['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII',
                 'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII',
                 'chrXIV', 'chrXV', 'chrXVI']
    try:
        chrom_new = str(chromList.index(chrom) + 1).zfill(2)
    except ValueError:
        if chrom == 'chrmt': 
            chrom_new = 'mi'
        else:
            print ('This chromosome ' + chrom + ' is not known.\n\
            Please re-check the reference list for chromosome conversion included in the python script.')
    return (chrom_new)


def getReferenceFiles(path_current):
    ''' This module uploads the required reference files. Default: Yeast. '''

    path_ref = os.path.join(path_current, 'references_yeast', '')
    try: 
        dictReads = pickle.load(open(path_ref + 'yeast_sequence.pkl', 'rb'))
        dictGenes = pickle.load(open(path_ref + 'yeast_transcripts.pkl', 'rb'))
        dicttRNA = pickle.load(open(path_ref + 'yeast_tRNA.pkl', 'rb'))
    except FileNotFoundError:
        print ('ERROR: Reference files not found. The reference files must be \
located in the folder "references_yeast" in the same directory as the scripts \
and the samples.')
        sys.exit()

    return (dictReads, dictGenes, dicttRNA)

def getOffset(path_current):
    ''' This module uploads the given offsets. '''

    path_ref = os.path.join(path_current, 'references_yeast', '')
    dictOffset = {}
    with open(path_ref + 'Asite_offset.txt', 'r') as f:
        for line in f: 
            if not line.startswith('length'):
                fields = line.split('\t')
                dictOffset[int(fields[0])] = int(fields[1].strip()) 
    return (dictOffset)



def RibosomeAssignment(input_path, input_sample, output_path, output_sample, mode, limits, threshold, text):

    dictFragments = {key: 0 for key in range(1,101)}
    path_current = os.path.dirname(os.path.realpath(__file__))
    dictReads, dictGenes, dicttRNA = getReferenceFiles(path_current)
    if mode == '5-end': 
        dictOffset = getOffset(path_current)

    # open input file
    inSam = open(input_path + input_sample, 'r')

    # read and process input data
    for line in inSam:
        if line.startswith('@') or line.startswith('CL:'):
            continue
        fields = line.split()
        strand = '+' if str(fields[1]) == '0' else '-'

        chrom = chromConversion(fields[2].strip())
        leftPos = int(fields[3])
        length = len(fields[9])
        dictFragments[length] += 1
        if limits[0] <= length <= limits[1]:
            if mode == 'center':
                start = leftPos + 11
                stop = leftPos + length - 11
            elif mode == '5-end':
                offset = dictOffset[length]
                start = leftPos + offset if strand == '+' else leftPos + length - offset
                stop = start + 1
            else: 
                print ('Please define analysis mode - it can be \'center\' or \'5-end\'.\nCheck for tipos.\n')
                return
            value = 1 / (stop - start)
            leftOrRight = 0 if strand == '+' else 2
            for pos in range(start, stop):
                key = chrom + '\t' + str(pos).zfill(7)
                dictReads[key][leftOrRight] += value

    # count total reads (per chromosome and all together )
    dictTotalReads = {}
    total  = 0
    for i in range(1, 17):
        key = str(i).zfill(2)
        dictTotalReads[key] = 0.0
    dictTotalReads['mi'] = 0.0
    for key_elem, value_elem in dictReads.items():
        dictTotalReads[key_elem.split('\t')[0]] += (float(value_elem[0]) + float(value_elem[2]))
        total += (float(value_elem[0]) + float(value_elem[2]))
        
    # count total reads (all together) -- only tRNA
    total_tRNA = 0
    for key in dicttRNA.keys():
        if key.startswith('t'):
            chrom = ''.join(list(dicttRNA[key][1])[:2])
            start = int(''.join(list(dicttRNA[key][0])[3:]))
            stop = int(''.join(list(dicttRNA[key][1])[3:])) + 1
            ind = 0 if (dicttRNA[key][2] == '+') else 2
            for pos in range(start, stop):
                total_tRNA += dictReads[chrom + '\t' + str(pos).zfill(7)][ind]    

    # normalize reads (--> RPM)
    for key_elem, value_elem in dictReads.items():
        value_elem.append('\t')
        value_elem.append((value_elem[0]/(total-total_tRNA))*1000000)
        value_elem.append('\t')
        value_elem.append((value_elem[2]/(total-total_tRNA))*1000000)
        
    # calculate gene expression
#    for value in dictGenes.values():
#        value.append(0.0)
#        value.append('--')
    threshold_norm = (threshold/(total-total_tRNA))*1000000
    ge_data = {}

    for gene in dictGenes.keys():
        raw_value = 0
        chrom = str(dictGenes[gene][1]).zfill(2)
        strand = 6 if dictGenes[gene][2] == '+' else 8
        local = dictGenes[gene][3].split(',')
        for exon in local:
            if exon == '':
                continue
            start = int(exon.split('-')[0]) if dictGenes[gene][2] == '+' else int(exon.split('-')[1])            # always the smaller value
            end = int(exon.split('-')[1]) if dictGenes[gene][2] == '+' else int(exon.split('-')[0])              # always the greater value
            length = (end - start) + 1
            for pos in range(length):
                key = chrom + '\t' + str(start+pos).zfill(7)
                value = dictReads[key][strand]
                raw_value += value

        gene_length = int(dictGenes[gene][7])*3
        norm_value = raw_value / gene_length * 1000
        ge_data[gene] = [str(raw_value), str(gene_length), str(norm_value)]
        if raw_value < threshold_norm: 
            ge_data[gene].append('excluded')
        else: 
            ge_data[gene].append('included')

    print ('Calculations done:')
    print (strftime('%a, %d %b %Y %H:%M:%S +0000', gmtime()))

    # write output file -- reads.pkl
    file = open(output_path + output_sample + '_Reads.pkl', 'wb')
    pickle.dump(dictReads, file)

    # write output file -- footprint lengths.txt
    listFragments = list(dictFragments.keys())
    listFragments.sort()
    with open(output_path + output_sample + '_FootprintLength.txt', 'w') as f:
        for length in listFragments:
            f.write(str(length) + '\t' + str(dictFragments[length]) + '\n')

    # write output file -- total reads.txt
    listTotal = list(dictTotalReads.keys())
    listTotal.sort()
    with open(output_path + output_sample + '_TotalReads.txt', 'w') as f:
        for i in listTotal: 
            f.write(str(i) + '\t' + str(dictTotalReads[i]) + '\n')
        f.write('total\t' + str(total) + '\n')
        f.write('total without tRNA\t' + str(total-total_tRNA))

    # write output file -- GeneExpression*.txt
    listGenes = list(ge_data.keys())
    listGenes.sort()
    with open(output_path + output_sample + '_GeneExpression.txt', 'w') as f:
        f.write('systematic transcript name\texpression [RPM]\ttranscript length[nt]\texpression [RPKM]\tquality tag\n')
        for gene in listGenes:
            f.write(gene)
            for item in ge_data[gene]:
                f.write('\t' + str(item))
            f.write('\n')            

    # write text files if required
    if text != 'no':
        with open(output_path + output_sample + '_Reads.txt', 'w') as f:
            listNorm = list(dictReads.keys())
            listNorm.sort()
            for pos in listNorm:
                f.write(pos + '\t')
                for item in dictReads[pos]:
                    f.write(str(item))
                f.write('\n')

    print ('Output files written:')
    print (strftime('%a, %d %b %Y %H:%M:%S +0000', gmtime()))

    return
 
 
if __name__ == '__main__':

    p = argparse.ArgumentParser(description='ribosomal A-site assignment of ribosome profiling data')

    # non-optional arguments
    p.add_argument('input_name', type=str, help = 'full name of input sam file, with file extension')
    p.add_argument('output_name', type=str, help = 'output sample name, NO file extension')
    p.add_argument('threshold', type=float, help = 'minimal number of footprints per transcript to be included in further analyses')
    # optional arguments
    p.add_argument('-i', '--input-path', dest = 'input_path', type=str, help = 'input path, default: cwd', default = os.path.join(os.getcwd(), ''))
    p.add_argument('-o', '--output-path', dest = 'output_path', type=str, help = 'output path, default: cwd', default = os.path.join(os.getcwd(), ''))
    p.add_argument('-m', '--mode', dest = 'mode', type=str, help = 'mode of ribosomal A-site assignment, default: 5-end', default = '5-end')
    p.add_argument('-l', '--limits', dest ='limits', type=tuple, help = 'min and max length of ribosome footprints, default: (23,40)' , default = (23,40))
    p.add_argument('-x', '--text', dest = 'text', type=str, help = 'required conversion of *_Reads.pkl to text file?, default: no', default = 'no')


    args = p.parse_args()

    input_path = args.input_path
    input_sample = args.input_name
    threshold = args.threshold
    output_path = args.output_path
    output_sample = args.output_name
    mode = args.mode
    limits = args.limits
    text = args.text

    RibosomeAssignment(input_path, input_sample, output_path, output_sample, mode, limits, threshold, text)

'''
>> to run this script via IDLE or another environment replace lines 257-280 by 
the section given below and manually type in the respective arguments: 

    input_path = '/path/to/data/'               # e.g. 'C:/SeRP/sequencing/data/'
    input_sample = 'sample file.sam'            # e.g. 'Sample1_accepted_hits.sam'
    output_path = '/path/to/output/'            # e.g. 'C:/SeRP/analysis/
    output_sample = 'sample name'               # e.g. 'S1' (no filename extension)
    threshold = xxx                             # minimal number of reads per gene to be "included" (int or float)

    mode = '5-end'                              # '5-end' or 'center'
    limits = (23, 40)                           # define the lower and upper border of allowed footprint lengths
    text = 'no'                                 # whether or not a text file version of *_Reads.pkl is required
'''

