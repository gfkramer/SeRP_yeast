###########################
#                         #
#    Metagene_Profiles    #
#                         #
###########################

'''
Created on Mon Sep 24 2018
@author: Ulrike Friedrich 

-------------------------------------------------------------------------------
DESCRIPTION

This script performes the metagene analysis for a set of samples comprising 
2 biological replicates of selective and total translatomes, each. 
The metagene profiles are plotted separately for each data set or showing the 
enrichment between selective and total translatome, the latter one in log2 
scale. 
Furthermore, the earlier defined threshold of the minimal number of reads per 
transcripts (default: 64 reads) must be set again, if not 64.  

The graph is saved as png file and as pdf file. 


'''

import argparse
import os
import re
import pickle
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl


def removeZero(input_list): 
    ''' This module replaces all zeros in a given list by np.nan to allow the 
    meaningful plotting on a log-scale. np.nan values produce gaps in the
    plotted line while zeros would decrease the line to the lower plotting area.
    ''' 
    
    output_list = []
    for elem in input_list: 
        if elem == 0:
            output_list.append(np.nan)
        else: 
            output_list.append(elem)
    
    output_list = np.array(output_list)    
    return (output_list)


def metageneProfiles(input_path, file_total1, file_total2, file_selec1, file_selec2, output_name, threshold):

    # upload input data
    total_1 = pickle.load(open(input_path + file_total1 + '_Reads.pkl', 'rb'))
    total_2 = pickle.load(open(input_path + file_total2 + '_Reads.pkl', 'rb'))
    selec_1 = pickle.load(open(input_path + file_selec1 + '_Reads.pkl', 'rb'))
    selec_2 = pickle.load(open(input_path + file_selec2 + '_Reads.pkl', 'rb'))

    # reference files
    path_current = os.path.dirname(os.path.realpath(__file__))
    path_ref = os.path.join(path_current, 'references_yeast_demo', '')
    dictGenes = pickle.load(open(path_ref + 'yeast_transcripts_demo.pkl', 'rb'))
    dictIntrons = pickle.load(open(path_ref + 'yeast_introns_demo.pkl', 'rb'))

    # process data
    dictMeta_1 = dict([i, []] for i in range(-50, 1501))
    dictMeta_2 = dict([i, []] for i in range(-50, 1501))
    dictMeta_3 = dict([i, []] for i in range(-50, 1501))
    dictMeta_4 = dict([i, []] for i in range(-50, 1501))
    
    dictMeta_5 = dict([i, []] for i in range(-50, 1501))
    dictMeta_6 = dict([i, []] for i in range(-50, 1501))

    for gene in dictGenes.keys():

        chrom = str(dictGenes[gene][1]).zfill(2)
        strand = dictGenes[gene][2]
        pos_list = dictIntrons[gene].copy()
        maxi = int(max(re.split('\W+', dictGenes[gene][3].strip(','))))
        mini = int(min(re.split('\W+', dictGenes[gene][3].strip(','))))
        length = len(pos_list)
        ind_raw = 0 if strand == '+' else 2
        ind_norm = 6 if strand == '+' else 8

        sum_raw_1 = 0
        sum_raw_2 = 0
        sum_norm_1 = 0
        sum_norm_2 = 0
        sum_norm_3 = 0
        sum_norm_4 = 0
        sum_norm_5 = 0
        sum_norm_6 = 0
        
        for pos in pos_list:

            sum_raw_1 += total_1[pos][ind_raw]                                  # for threshold
            sum_raw_2 += total_2[pos][ind_raw]                                  # for threshold

            sum_norm_1 += total_1[pos][ind_norm]                                # for normalization
            sum_norm_2 += total_2[pos][ind_norm]                                # for normalization
            sum_norm_3 += selec_1[pos][ind_norm]                                # for normalization
            sum_norm_4 += selec_2[pos][ind_norm]                                # for normalization

            try:
                sum_norm_5 += selec_1[pos][ind_norm] / total_1[pos][ind_norm]
            except: 
                pass
            try:
                sum_norm_6 += selec_2[pos][ind_norm] / total_2[pos][ind_norm]
            except: 
                pass

        norm_val_1 = sum_norm_1 / length
        norm_val_2 = sum_norm_2 / length
        norm_val_3 = sum_norm_3 / length
        norm_val_4 = sum_norm_4 / length
        norm_val_5 = sum_norm_5 / length
        norm_val_6 = sum_norm_5 / length

        if sum_raw_1 < threshold or sum_raw_2 < threshold:                      # exclude genes with total1/2 below threshold
            continue 

        if strand == '+':                                                       # add 50nt before each ORF
            start_pos = mini
            for i in range(-1, -51, -1): 
                pos_list.insert(0, str(chrom) + '\t' + str(start_pos + i).zfill(7))
        elif strand == '-': 
            start_pos = maxi
            for i in range(1, 51, 1): 
                pos_list.insert(0, str(chrom) + '\t' + str(start_pos + i).zfill(7))

        pos_counter = -50
        for pos in pos_list:
            if pos_counter > 1500:
                break
            dictMeta_1[pos_counter].append(total_1[pos][ind_norm] / norm_val_1)
            dictMeta_2[pos_counter].append(total_2[pos][ind_norm] / norm_val_2)
            dictMeta_3[pos_counter].append(selec_1[pos][ind_norm] / norm_val_3)
            dictMeta_4[pos_counter].append(selec_2[pos][ind_norm] / norm_val_4)
            try: 
                ratio1 = selec_1[pos][ind_norm] / total_1[pos][ind_norm]
                dictMeta_5[pos_counter].append(ratio1 / norm_val_5)
            except: 
                pass
            try: 
                ratio2 = selec_2[pos][ind_norm] / total_2[pos][ind_norm]
                dictMeta_6[pos_counter].append(ratio2 / norm_val_6)
            except: 
                pass
            pos_counter += 1

    y_1 = []
    y_2 = []
    y_3 = []
    y_4 = []
    y_5 = []
    y_6 = []
    for n in range(-48,1498,3):
        y_1.append(np.mean([np.mean(dictMeta_1[n]), np.mean(dictMeta_1[n+1]), np.mean(dictMeta_1[n+2])]))
        y_2.append(np.mean([np.mean(dictMeta_2[n]), np.mean(dictMeta_2[n+1]), np.mean(dictMeta_2[n+2])]))
        y_3.append(np.mean([np.mean(dictMeta_3[n]), np.mean(dictMeta_3[n+1]), np.mean(dictMeta_3[n+2])]))
        y_4.append(np.mean([np.mean(dictMeta_4[n]), np.mean(dictMeta_4[n+1]), np.mean(dictMeta_4[n+2])]))
        y_5.append(np.mean([np.mean(dictMeta_5[n]), np.mean(dictMeta_5[n+1]), np.mean(dictMeta_5[n+2])]))
        y_6.append(np.mean([np.mean(dictMeta_6[n]), np.mean(dictMeta_6[n+1]), np.mean(dictMeta_6[n+2])]))

    y_1 = np.array(y_1)
    y_2 = np.array(y_2)
    y_3 = np.array(y_3)
    y_4 = np.array(y_4)
    y_5 = removeZero(y_5)
    y_6 = removeZero(y_6)

    x = np.arange(-16,500)

    # generate graph
    mpl.rcParams['axes.linewidth'] = 0.64
    mpl.rcParams['font.sans-serif'] = "Arial"
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6
    mpl.rcParams['legend.fontsize'] = 8

    fig = plt.figure(1, figsize = (3.2, 3.2))
    plt.subplots_adjust(0,0,1,1)
    ax = fig.add_subplot(211)
        
    plt.tick_params(axis='y', which='both',
                    left='on', right = 'off',
                    labelleft='on', labelright='off',
                    direction = 'out', width = 0.64, length = 2.0)
    plt.tick_params(axis='x', which='both',
                    bottom='on', top='off',
                    labelbottom='on', 
                    direction = 'out', width = 0.64, length = 2.0)
    x_label = ('position along average transcript [nt]')
    y_label = ('ribosome density [a.u.]')
    ax.set_ylabel(y_label, fontsize = 8)
    ax.set_xlim(-20, 500)   
    plt.plot(x, y_1, color = 'black', linewidth = 0.64)
    plt.plot(x, y_2, color = 'black', linewidth = 0.64)
    plt.plot(x, y_3, color = 'orange', linewidth = 0.64)
    plt.plot(x, y_4, color = 'orange', linewidth = 0.64)
    plt.legend(('Total 1', 'Total 2', 'Selective 1', 'Selective 2'), loc = 'upper right')


    ax = fig.add_subplot(212)

    plt.tick_params(axis='y', which='both',
                    left='on', right = 'off',
                    labelleft='on', labelright='off',
                    direction = 'out', width = 0.64, length = 2.0)
    plt.tick_params(axis='x', which='both',
                    bottom='on', top='off',
                    labelbottom='on', 
                    direction = 'out', width = 0.64, length = 2.0)
    y_label = ('enrichment [a.u.]')
    ax.set_xlabel(x_label, fontsize = 8)
    ax.set_ylabel(y_label, fontsize = 8)
    ax.set_xlim(-20, 500)
    plt.plot(x, y_5, color = 'brown', linewidth = 0.64)
    plt.plot(x, y_6, color = 'brown', linewidth = 0.64)
    plt.yscale('log', basey=2)
    plt.legend(('Enrichment 1', 'Enrichment 2'), loc = 'upper right')
    

    # save graph
    plt.savefig(input_path + output_name + '_metagene.png', bbox_inches='tight', dpi = 600)
    plt.savefig(input_path + output_name + '_metagene.pdf', bbox_inches='tight')
    plt.close()


    return
       


if __name__ == '__main__':

    p = argparse.ArgumentParser(description='metagene analysis of SeRP data sets')

    # non-optional arguments
    p.add_argument('file_total_1', type=str, help = 'sample name total translatome 1, without file extension')
    p.add_argument('file_total_2', type=str, help = 'sample name total translatome 2, without file extension')
    p.add_argument('file_selec_1', type=str, help = 'sample name selective translatome 1, without file extension')
    p.add_argument('file_selec_2', type=str, help = 'sample name selective translatome 2, without file extension')
    p.add_argument('output_name', type=str, help = 'unique experiment name, without file extension')
    # optional arguments
    p.add_argument('-i', '--input-path', dest = 'input_path', type=str, help = 'input path, default: cwd', default = os.path.join(os.getcwd(), ''))
    p.add_argument('-t', '--threshold', dest = 'threshold', type=float, help = 'minimal number of footprints per gene to be included in further analysis, default: 64.0', default = 64.0)

    args = p.parse_args()

    input_path = args.input_path
    file_total1 = args.file_total_1
    file_total2 = args.file_total_2
    file_selec1 = args.file_selec_1
    file_selec2 = args.file_selec_2
    output_name = args.output_name
    threshold = args.threshold

    metageneProfiles(input_path, file_total1, file_total2, file_selec1, file_selec2, output_name, threshold)

'''
>> to run this script via IDLE or another environment replace lines 237-257 by 
the section given below and manually type in the respective arguments: 

    input_path = '/path/to/data/'               # e.g. 'C:/SeRP/sequencing/data/'
    file_total1 = 'sample name'                 # sample name total translatome 1 (no file extension)
    file_total2 = 'sample name'                 # sample name total translatome 2 (no file extension)
    file_selec1 = 'sample name'                 # sample name selective translatome 1 (no file extension)
    file_selec2 = 'sample name'                 # sample name selective translatome 2 (no file extension)
    output_name = 'experiment name'             # experiment name (no file extension)
    threshold = xxx                             # minimal number of reads per gene to be "included" (int or float)

'''