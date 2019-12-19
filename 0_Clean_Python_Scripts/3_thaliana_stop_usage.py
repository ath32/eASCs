### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FILES ###

source_pol = '4_Pollen/plntphys_pp.104.057935_57935Supplementary_Table_1.csv'
source_fasta = '4_Pollen/Arabidopsis_thaliana.fa'

### MAIN ###

def main():

    csv_total = []

    #Get pollen specific and depleted genes
    raw_pol = open(source_pol).read()
    pol_split = raw_pol.split('\n')

    pollen_genes = []
    nonpollen_genes = []

    for i in pol_split:
        items = i.split(',')
        pollen_selective = items[0]
        pollen_depleted = items[3]
        pollen_specific = items[4]
        id = items[6].upper()
        if pollen_selective == 'X':
            pollen_genes.append(id)
        if pollen_depleted == 'X':
            nonpollen_genes.append(id)

    #Get sequences for both sets
    raw_fasta = open(source_fasta).read()
    fasta_split = raw_fasta.split('>')

    pol_seqs = []
    nonpol_seqs = []

    for k in fasta_split:
        if k != '':
            chunks = k.split('\n')
            line = chunks[0]
            seq = chunks[1]
            line_split = line.split(';')
            fasta_id = line_split[0].replace('ID=gene:', '')
            if fasta_id in pollen_genes:
                pol_seqs.append(seq)
            elif fasta_id in nonpollen_genes:
                nonpol_seqs.append(seq)
            else:
                continue

    #Get stop usage for each set of genes
    print ('MULTI:')
    get_stops(nonpol_seqs)
    print ('SINGLE:')
    get_stops(pol_seqs)

### FUNCTIONS ###

def get_stops(genes):

    taa = 0
    tga = 0
    tag = 0

    for sequence in genes:
        primary = sequence[:3]
        if primary == 'TAA':
            taa += 1
        elif primary == 'TGA':
            tga += 1
        elif primary == 'TAG':
            tag += 1

    print ('TAA:', taa)
    print ('TGA:', tga)
    print ('TAG:', tag)

### RUN ###

if __name__ == '__main__':
    main()
