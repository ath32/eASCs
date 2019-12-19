### IMPORTS ###

import os
import numpy as np
import csv
import pandas as pd

### CHOOSE SOURCE FILES - unhash genome of interest ###

# #A. Niger
# source_genes = '8_Conidia/CAGs_niger.txt'
# source_fasta = '8_Conidia/Filtered_FASTA/Aspergillus_niger.fa'

# A. Oryzae
source_genes = '8_Conidia/CAGs_oryzae.txt'
source_fasta = '8_Conidia/Filtered_FASTA/Aspergillus_oryzae.fa'

### MAIN ###

def main():

    csv_total = []

    #Get gene ids
    df = raw_data = open(source_genes).read()
    split = df.split('\n')

    CAGs = []

    for line in split[1:]:
        line_split = line.split('\t')
        id = line_split[0]
        fc = line_split[7]
        if fc != '#DIV/0!' and float(fc) > 4:
            CAGs.append(id)
        elif line_split[4] == '0':
            CAGs.append(id)

    #Get sequences for both sets
    raw_fasta = open(source_fasta).read()
    fasta_split = raw_fasta.split('>')

    CAG_seqs = []
    nonCAG_seqs = []

    for k in fasta_split:
        if k != '':
            chunks = k.split('\n')
            line = chunks[0]
            seq = chunks[1]
            line_split = line.split(';')
            fasta_id = line_split[0].replace('ID=gene:', '')
            if fasta_id in CAGs:
                CAG_seqs.append(seq)
            elif fasta_id not in CAGs:
                nonCAG_seqs.append(seq)

    #Calculate ASC frequencies for each position
    overall_CAG, positions_CAG, a, b = get_ASCs(CAG_seqs)
    overall_nonCAG, positions_nonCAG, a, b = get_ASCs(nonCAG_seqs)

    #Primary stop analysis
    print('MULTICELL')
    primary_stops(nonCAG_seqs)
    print('UNICELL')
    primary_stops(CAG_seqs)

### FUNCTIONS ###

def get_ASCs(sequences):

    overall_stops = 0
    stops = [0, 0, 0, 0, 0, 0]
    n = len(sequences)
    overall_codons = n * 6

    for sequence in sequences:
        codon_seq = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        if codon_seq[1] == 'TAA' or codon_seq[1] == 'TGA' or codon_seq[1] == 'TAG':
             overall_stops += 1
             stops[0] += 1
        if codon_seq[2] == 'TAA' or codon_seq[2] == 'TGA' or codon_seq[2] == 'TAG':
             overall_stops += 1
             stops[1] += 1
        if codon_seq[3] == 'TAA' or codon_seq[3] == 'TGA' or codon_seq[3] == 'TAG':
             overall_stops += 1
             stops[2] += 1
        if codon_seq[4] == 'TAA' or codon_seq[4] == 'TGA' or codon_seq[4] == 'TAG':
             overall_stops += 1
             stops[3] += 1
        if codon_seq[5] == 'TAA' or codon_seq[5] == 'TGA' or codon_seq[5] == 'TAG':
             overall_stops += 1
             stops[4] += 1
        if codon_seq[6] == 'TAA' or codon_seq[6] == 'TGA' or codon_seq[6] == 'TAG':
             overall_stops += 1
             stops[5] += 1

    overall_f = overall_stops / overall_codons
    f_array = np.array(stops) / n
    f_list = f_array.tolist()

    return overall_f, f_list, overall_stops, overall_codons


def primary_stops(genes):

    taa_seqs = []
    tga_seqs = []
    tag_seqs = []

    for sequence in genes:
        primary = sequence[:3]
        if primary == 'TAA':
            taa_seqs.append(sequence)
        if primary == 'TGA':
            tga_seqs.append(sequence)
        if primary == 'TAG':
            tag_seqs.append(sequence)

    x, y, taa_n, taa_codons = get_ASCs(taa_seqs)
    x, y, tga_n, tga_codons = get_ASCs(tga_seqs)
    x, y, tag_n, tag_codons = get_ASCs(tag_seqs)

    print ('TAA:', taa_n, taa_codons)
    print ('TGA:', tga_n, tga_codons)
    print ('TAG:', tag_n, tag_codons)


### RUN ###
if __name__ == '__main__':
    main()
