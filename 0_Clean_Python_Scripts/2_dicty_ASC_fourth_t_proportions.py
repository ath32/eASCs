### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FILES ###

source_dicty = '2_Development/sociality_genes.csv'
source_fasta = '2_Development/Dictyostelium_ac.fa'

### MAIN ###

def main():

    csv_total = []

    #Get gene ids
    raw_data = open(source_dicty).read()
    dicty_split = raw_data.split('\n')

    social_genes = []
    single_genes = []

    for i in dicty_split:
        items = i.split(',')
        identifier = items[1]
        gene_id = items[0]

        if identifier == '1':
            social_genes.append(gene_id)
        if identifier == '0':
            single_genes.append(gene_id)

    #Get sequences for both sets
    raw_fasta = open(source_fasta).read()
    fasta_split = raw_fasta.split('>')

    social_seqs = []
    single_seqs = []

    for k in fasta_split:
        if k != '':
            chunks = k.split('\n')
            line = chunks[0]
            seq = chunks[1]
            line_split = line.split(';')
            fasta_id = line_split[0].replace('ID=gene:', '')
            if fasta_id in social_genes:
                social_seqs.append(seq)
            elif fasta_id in single_genes:
                single_seqs.append(seq)
            else:
                continue

    #Get position +1 +4T freq for both sets
    single_t_freq, single_t_hits, single_t_total = get_fourth_t(single_seqs)
    social_t_freq, social_t_hits, social_t_total = get_fourth_t(social_seqs)

    #Get position +1 ASC freq for both sets
    single_asc_freq, single_asc_hits, single_asc_total = get_ASC_pos1(single_seqs)
    social_asc_freq, social_asc_hits, social_asc_total = get_ASC_pos1(social_seqs)

    print (single_asc_hits, single_t_hits - single_asc_hits, single_t_hits)
    print (social_asc_hits, social_t_hits - social_asc_hits, social_t_hits)

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
