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

    #Calculate ASC frequencies for each position
    overall_social, positions_social = get_t_codons(social_seqs)
    overall_single, positions_single = get_t_codons(single_seqs)

    social_line = ['social_genes', overall_social, positions_social[0], positions_social[1], positions_social[2], positions_social[3], positions_social[4], positions_social[5]]
    single_line = ['single_genes', overall_single, positions_single[0], positions_single[1], positions_single[2], positions_single[3], positions_single[4], positions_single[5]]

    csv_total.append(social_line)
    csv_total.append(single_line)
    headers = ['Set', 'Overall_f', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6']

    #Write output file
    filename = "social_vs_single_t.csv"
    subdir = "2_Development"
    filepath = os.path.join(subdir, filename)
    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### FUNCTIONS ###

def get_t_codons(sequences):

    overall_count = 0
    count = [0, 0, 0, 0, 0, 0]
    n = len(sequences)
    overall_codons = n * 6

    for sequence in sequences:
        codon_seq = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

        codon1 = codon_seq[1]
        codon2 = codon_seq[2]
        codon3 = codon_seq[3]
        codon4 = codon_seq[4]
        codon5 = codon_seq[5]
        codon6 = codon_seq[6]

        if codon1[0] == 'T' and codon1 != 'TAA' and codon1 != 'TGA' and codon1 != 'TAG':
             overall_count += 1
             count[0] += 1
        if codon2[0] == 'T' and codon2 != 'TAA' and codon2 != 'TGA' and codon2 != 'TAG':
             overall_count += 1
             count[1] += 1
        if codon3[0] == 'T' and codon3 != 'TAA' and codon3 != 'TGA' and codon3 != 'TAG':
             overall_count += 1
             count[2] += 1
        if codon4[0] == 'T' and codon4 != 'TAA' and codon4 != 'TGA' and codon4 != 'TAG':
             overall_count += 1
             count[3] += 1
        if codon5[0] == 'T' and codon5 != 'TAA' and codon5 != 'TGA' and codon5 != 'TAG':
             overall_count += 1
             count[4] += 1
        if codon6[0] == 'T' and codon6 != 'TAA' and codon6 != 'TGA' and codon6 != 'TAG':
             overall_count += 1
             count[5] += 1

    overall_f = overall_count / overall_codons
    f_array = np.array(count) / n
    f_list = f_array.tolist()

    return overall_f, f_list


def get_fourth_t(genes):

    n = len(genes)
    t = 0

    for gene in genes:
        fourth = gene[3]
        if fourth == 'T':
            t += 1

    t_percent = t / n * 100

    return t_percent

### RUN ###

if __name__ == '__main__':
    main()
