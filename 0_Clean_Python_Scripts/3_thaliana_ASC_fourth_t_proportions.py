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

    #Get position +1 +4T freq for both sets
    pollen_t_freq, pollen_t_hits, pollen_t_total = get_fourth_t(pol_seqs)
    nonpol_t_freq, nonpol_t_hits, nonpol_t_total= get_fourth_t(nonpol_seqs)

    #Get position +1 ASC freq for both sets
    pollen_asc_freq, pollen_asc_hits, pollen_asc_total = get_ASC_pos1(pol_seqs)
    nonpol_asc_freq, nonpol_asc_hits, nonpol_asc_total = get_ASC_pos1(nonpol_seqs)

    print (pollen_asc_hits, pollen_t_hits - pollen_asc_hits, pollen_t_hits)
    print (nonpol_asc_hits, nonpol_t_hits - nonpol_asc_hits, nonpol_t_hits)

### FUNCTIONS ###

def get_fourth_t(genes):

    n = len(genes)
    t = 0

    for gene in genes:
        fourth = gene[3]
        if fourth == 'T':
            t += 1

    overall_t = t / n

    return overall_t, t, n


def get_ASC_pos1(genes):

    n = len(genes)
    t = 0

    for gene in genes:
        first_position = gene[3:6]
        if first_position == 'TAA' or first_position == 'TGA' or first_position == 'TAG':
            t += 1

    overall_t = t / n

    return overall_t, t, n

### RUN ###

if __name__ == '__main__':
    main()
