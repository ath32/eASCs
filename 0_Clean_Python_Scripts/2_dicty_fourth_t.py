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
                if seq[3:6] != 'TAA' and seq[3:6] != 'TGA' and seq[3:6] != 'TAG':
                    social_seqs.append(seq)
            elif fasta_id in single_genes:
                if seq[3:6] != 'TAA' and seq[3:6] != 'TGA' and seq[3:6] != 'TAG':
                    single_seqs.append(seq)
            else:
                continue

    #Calculate +4T for each position
    T_social = get_fourth_t(social_seqs)
    T_single = get_fourth_t(single_seqs)

    #Get stop groups
    social_taa, social_tga, social_tag = get_stop_groups(social_seqs)
    single_taa, single_tga, single_tag = get_stop_groups(single_seqs)

    #Test each stop group
    soc_taa_a, soc_taa_t, soc_taa_g, soc_taa_c = get_fourth_freqs(social_taa)
    soc_tga_a, soc_tga_t, soc_tga_g, soc_tga_c = get_fourth_freqs(social_tga)
    soc_tag_a, soc_tag_t, soc_tag_g, soc_tag_c = get_fourth_freqs(social_tag)

    sin_taa_a, sin_taa_t, sin_taa_g, sin_taa_c = get_fourth_freqs(single_taa)
    sin_tga_a, sin_tga_t, sin_tga_g, sin_tga_c = get_fourth_freqs(single_tga)
    sin_tag_a, sin_tag_t, sin_tag_g, sin_tag_c = get_fourth_freqs(single_tag)

    social_line = ['social_genes', T_social, soc_taa_a, soc_taa_t, soc_taa_g, soc_taa_c, soc_tga_a, soc_tga_t, soc_tga_g, soc_tga_c, soc_tag_a, soc_tag_t, soc_tag_g, soc_tag_c]
    single_line = ['single_genes', T_single, sin_taa_a, sin_taa_t, sin_taa_g, sin_taa_c, sin_tga_a, sin_tga_t, sin_tga_g, sin_tga_c, sin_tag_a, sin_tag_t, sin_tag_g, sin_tag_c]
    csv_total.append(social_line)
    csv_total.append(single_line)

    #Write output file
    headers = ['Set', 'overall_t_freq', 'TAA_A', 'TAA_T', 'TAA_G', 'TAA_C', 'TGA_A', 'TGA_T', 'TGA_G', 'TGA_C', 'TAG_A', 'TAG_T', 'TAG_G', 'TAG_C']
    filename = "fourth.csv"
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

    return overall_f, f_list


def get_fourth_t(genes):

    taa = []
    tga = []
    tag = []

    n = len(genes)
    t = 0

    for gene in genes:
        fourth = gene[3]
        if fourth == 'T':
            t += 1
        if gene[0:3] == 'TAA':
            taa.append(gene)
        if gene[0:3] == 'TGA':
            tga.append(gene)
        if gene[0:3] == 'TAG':
            tag.append(gene)

    overall_t = t / n

    return overall_t


def get_stop_groups(genes):

    taa = []
    tga = []
    tag = []

    for gene in genes:
        if gene[0:3] == 'TAA':
            taa.append(gene)
        if gene[0:3] == 'TGA':
            tga.append(gene)
        if gene[0:3] == 'TAG':
            tag.append(gene)

    return taa, tga, tag


def get_fourth_freqs(genes):

    a = 0
    t = 0
    g = 0
    c = 0

    for gene in genes:
        fourth = gene[3]
        if fourth == 'A':
            a += 1
        if fourth == 'T':
            t += 1
        if fourth == 'G':
            g += 1
        if fourth == 'C':
            c += 1

    if len(genes) > 0:
        a_freq = a / len(genes)
        t_freq = t / len(genes)
        g_freq = g / len(genes)
        c_freq = c / len(genes)
    else:
        a_freq = 0
        t_freq = 0
        g_freq = 0
        c_freq = 0

    return a_freq, t_freq, g_freq, c_freq


### RUN ###

if __name__ == '__main__':
    main()
