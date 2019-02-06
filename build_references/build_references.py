#############################
#                           #
#      Reference_Files      #
#                           #
#############################

"""
Created on Wed Jan 16 2019
@author: Ulrike Friedrich


This script generates the required reference files for SeRP data analysis. 
Input files are: 

fasta file containing the sequence of each chromosome of the organism of interest
--> S. cerevisiae data is downloaded from SGD ('S288C_reference_sequence_R64-2-1_20150113.fsa')

text file containing all genes with additional information in tab separated format
--> S. cerevisiae data from SGD, present as excel file ('yeast_transcripts.xlsx')


"""



import pickle 



def genYeastSequence(input_path, output_path):
        
    fasta_file = 'S288C_reference_sequence_R64-2-1_20150113.fsa'
    seqref_file = 'yeast_sequence.pkl'
    dictSeq = {}
    
    with open(input_path + fasta_file, 'r') as f: 
        for line in f:
            if line.startswith('>'): 
                nc = int(line.split('\t')[0].split('|')[1][-2:]) - 32
                if nc >= 1: 
                    chrom = str(nc).zfill(2)
                else:
                    chrom = 'mi'
                pos_counter = 0
            else:
                line = line.strip()
                for nt in line:
                    pos_counter += 1
                    pos = str(pos_counter).zfill(7)
                    key = chrom + '\t' + pos
                    value = [0.0, '\t', 0.0, '\t', nt]
                    dictSeq[key] = value

    print (len(dictSeq.keys()))
    pickle.dump(dictSeq, open(output_path + seqref_file, 'wb'))

    return


def genYeastTranscripts(input_path, output_path): 
    
    genes_file = 'yeast_transcripts.txt'
    generef_file = 'yeast_transcripts.pkl'
    dictGenes = {}
    
    with open(input_path + genes_file, 'r') as f: 
        for line in f:
            if line.startswith('systematic'):
                continue
            else:
                fields = line.split('\t')
                gene = fields[0]
                value = fields[1:]
                dictGenes[gene] = value

    pickle.dump(dictGenes, open(output_path + generef_file, 'wb'))
        
    return


def genYeasttRNA(input_path, output_path):
    
    tRNA_file = 'rna_coding_R64-2-1_20150113.fasta'
    tRNA_out_file= 'yeast_tRNA.pkl'
    dicttRNA = {}

    dictChroms = {'I':'01', 'II':'02', 'III':'03', 'IV':'04', 'V':'05', 'VI':
        '06', 'VII':'07', 'VIII':'08', 'IX':'09', 'X':'10', 'XI':'11', 'XII':
            '12', 'XIII':'13', 'XIV':'14', 'XV':'15', 'XVI':'16', 'Mito':'mi'}
    
    with open(input_path + tRNA_file, 'r') as f: 
        for line in f: 
            if line.startswith('>t'):
                fields = line.split(' ')
                gene = fields[0].strip('>')
                chrom = dictChroms[fields[4]]

                pos = fields[6].strip(',')
                strand_spec = fields[10]
                if strand_spec == 'reverse':
                    strand = '-'
                    if ',' in pos:
                        pos_list = pos.strip(',').split('-')[1].split(',')
                        start = pos_list[0]
                        stop = pos_list[1]
                    else:
                        pos_list = pos.split('-')
                        start = pos_list[-1]
                        stop = pos_list[0]
                else:
                    strand = '+'
                    pos_list = pos.split('-')
                    start = pos_list[0]
                    stop = pos_list[-1].strip(',')

                dicttRNA[gene] = [chrom + '\t' + str(start).zfill(7), chrom + '\t' + str(stop).zfill(7), strand]
    
    pickle.dump(dicttRNA, open(output_path + tRNA_out_file, 'wb'))
    
    return


def genYeastIntrons(input_path, output_path):
    ''' This module generated a pickle file containing all chromosomal positions
    as a list for each transcript. Meaning that introns are excluded and the complete 
    transcript can be printen in graphs. 
    '''

    dictGenes = pickle.load(open(output_path + 'yeast_transcripts.pkl', 'rb'))
    dictIntrons = {}
    fileIntrons = 'yeast_introns.pkl'
    
    for gene in dictGenes:
        dictIntrons[gene] = []
        chroms = str(dictGenes[gene][1])
        strand = dictGenes[gene][2]
        splits = dictGenes[gene][3].split(',')
        splits = splits[:-1]

        if strand == '+': 
            for region in splits:
                region = region.split('-')
                start = int(region[0])
                stop = int(region[1])
                for pos in range(start, stop+1):
                    key = chroms.zfill(2) + '\t' + str(pos).zfill(7)
                    dictIntrons[gene].append(key)

        if strand == '-':
            for region in reversed(splits):
                region = region.split('-')
                start = int(region[0])
                stop = int(region[1])
                for pos in range(start, stop-1, -1):
                    key = chroms.zfill(2) + '\t' + str(pos).zfill(7)
                    dictIntrons[gene].append(key)

    pickle.dump(dictIntrons, open(output_path + fileIntrons, 'wb'))

    return
    


if __name__ == '__main__':
    
    input_path = '/path/to/folder/build_references/'
    output_path = '/path/to/folder/references_yeast/'

    genYeastSequence(input_path, output_path)
    genYeastTranscripts(input_path, output_path)
    genYeasttRNA(input_path, output_path)
    genYeastIntrons(input_path, output_path)
