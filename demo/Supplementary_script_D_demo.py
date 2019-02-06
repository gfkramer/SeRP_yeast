############################
#                          #
#     Total_Enrichment     #
#                          #
############################

'''
Created on Tue Sep 25 2018
@author: Ulrike Friedrich 

-------------------------------------------------------------------------------
DESCRIPTION

This script calculates the total enrichment (TE) value for each transcript 
comparing two biological replicates of total and selective samples. 

The output is a text file with the following columns: 
systematic transcript name 
(trivial) transcript name
for each of the four samples: 
    (normalized) gene expression value [RPKM]
    'included' or 'excluded'
ratio - replicate 1
ratio - replicate 2
average of both
log2 transform of average 
number of samples with 'included' 

Note: If the ratio, the average or the log2 transform cannot be calculated 
due to e.g. zero footprints for a transcript, the respective value is given as 
'n.a.'.


'''

import argparse
import math
import os
import pickle


def calcTE(input_path, file_total1, file_total2, file_selec1, file_selec2, output_name): 

    # read input data    
    data = {}
    spec = {}
    proc = {}
    
    for file in [file_total1, file_total2, file_selec1, file_selec2]: 
        with open(input_path + file + '_GeneExpression.txt', 'r') as f:
            for line in f:
                if line.startswith('systematic transcript name'):
                    continue
                fields = line.split('\t')
                gene = fields[0]
                ge = float(fields[3])
                sp = fields[4].strip()
                data.setdefault(gene, []).append(ge)
                spec.setdefault(gene, []).append(sp)

    # reference file for gene names
    path_current = os.path.dirname(os.path.realpath(__file__))
    path_ref = os.path.join(path_current, 'references_yeast_demo', '')
    dictGenes = pickle.load(open(path_ref + 'yeast_transcripts_demo.pkl', 'rb'))

    # process data / calculate ratio and TE
    for (gene, values) in data.items():
        try:
            ratio1 = values[2] / values[0]
        except ZeroDivisionError: 
            ratio1 = 'n.a.'
        try:
            ratio2 = values[3] / values[1]
        except ZeroDivisionError: 
            ratio2 = 'n.a.'
        try:
            TE = (ratio1 + ratio2) / 2
        except (TypeError, ValueError):
            TE = 'n.a.'
        try:
            log2 = math.log2(TE)
        except (TypeError, ValueError):
            log2 = 'n.a.'
        specNO = str(spec[gene].count('included')) + ' of 4'
        proc[gene] = [ratio1, ratio2, TE, log2, specNO]

    # write output file
    listGenes = list(data.keys())
    listGenes.sort()
    header = 'systematic transcript name\ttranscript name\ttotal 1\t\ttotal2\t\tselective 1\
    \t\tselective 2\t\tratio 1\tratio 2\taverage\taverage [log2]\tno. of included values\n'

    with open(input_path + output_name + '_TE.txt', 'w') as f: 
        f.write(header)
        for gene in listGenes:
            f.write(gene + '\t' + dictGenes[gene][0] + '\t')
            for n in range(4):
                f.write(str(data[gene][n]) + '\t' + spec[gene][n] + '\t')
            for value in proc[gene]:
                f.write(str(value) + '\t')
            f.write('\n')

    return


if __name__ == '__main__':

    p = argparse.ArgumentParser(description='calculation of total enrichment of SeRP data sets (selective vs. total translatomes)')

    # non-optional arguments
    p.add_argument('file_total_1', type=str, help = 'sample name total translatome 1, without file extension')
    p.add_argument('file_total_2', type=str, help = 'sample name total translatome 2, without file extension')
    p.add_argument('file_selec_1', type=str, help = 'sample name selective translatome 1, without file extension')
    p.add_argument('file_selec_2', type=str, help = 'sample name selective translatome 2, without file extension')
    p.add_argument('output_name', type=str, help = 'unique experiment name, without file extension')
    # optional arguments
    p.add_argument('-i', '--input-path', dest = 'input_path', type=str, help = 'input path, default: cwd', default = os.path.join(os.getcwd(), ''))

    args = p.parse_args()

    input_path = args.input_path
    file_total1 = args.file_total_1
    file_total2 = args.file_total_2
    file_selec1 = args.file_selec_1
    file_selec2 = args.file_selec_2
    output_name = args.output_name

    calcTE(input_path, file_total1, file_total2, file_selec1, file_selec2, output_name)

'''
>> to run this script via IDLE or another environment replace lines 107-125 by 
the section given below and manually type in the respective arguments: 

    input_path = '/path/to/data/'               # e.g. 'C:/SeRP/sequencing/data/'
    file_total1 = 'sample name'                 # sample name total translatome 1 (no file extension)
    file_total2 = 'sample name'                 # sample name total translatome 2 (no file extension)
    file_selec1 = 'sample name'                 # sample name selective translatome 1 (no file extension)
    file_selec2 = 'sample name'                 # sample name selective translatome 2 (no file extension)
    output_name = 'experiment name'             # experiment name (no file extension)

'''








