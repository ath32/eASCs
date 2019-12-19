### IMPORTS ###

import os
import numpy as np
import csv
import scipy.stats as st

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

    #Calculate ASC frequencies for each position
    overall_social, positions_social, overall_stops_m, overall_codons_m = get_ASCs(social_seqs)
    overall_single, positions_single, overall_stops_s, overall_codons_s = get_ASCs(single_seqs)

    #Primary stop analysis
    print('MULTICELL')
    primary_stops(social_seqs)
    print('UNICELL')
    primary_stops(single_seqs)

    social_line = ['social_genes', overall_social, positions_social[0], positions_social[1], positions_social[2], positions_social[3], positions_social[4], positions_social[5]]
    single_line = ['single_genes', overall_single, positions_single[0], positions_single[1], positions_single[2], positions_single[3], positions_single[4], positions_single[5]]

    csv_total.append(social_line)
    csv_total.append(single_line)
    headers = ['Set', 'Overall_f', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6']

    #Write output file
    filename = "social_vs_single.csv"
    subdir = "2_Development"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

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


def get_fourth_t(genes):

    n = len(genes)
    t = 0

    for gene in genes:
        fourth = gene[3]
        if fourth == 'T':
            t += 1

    t_percent = t / n * 100

    return t_percent


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
