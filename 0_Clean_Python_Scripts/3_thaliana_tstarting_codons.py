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

    #Calculate T-starting codon frequencies for each position
    overall_pol, positions_pol = get_t_codons(pol_seqs)
    overall_nonpol, positions_nonpol = get_t_codons(nonpol_seqs)
    pol_line = ['pollen_genes', overall_pol, positions_pol[0], positions_pol[1], positions_pol[2], positions_pol[3], positions_pol[4], positions_pol[5]]
    nonpol_line = ['nonpollen_genes', overall_nonpol, positions_nonpol[0], positions_nonpol[1], positions_nonpol[2], positions_nonpol[3], positions_nonpol[4], positions_nonpol[5]]
    csv_total.append(pol_line)
    csv_total.append(nonpol_line)

    #Write output files
    headers = ['Set', 'Overall_f', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6']
    filename = "pollen_genes_tstarting_2005.csv"
    subdir = "4_Pollen"
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

    print (overall_count, overall_codons)

    return count, n


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
