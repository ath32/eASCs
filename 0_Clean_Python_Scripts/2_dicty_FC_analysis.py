### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FILES ###

source_fc = '2_Development/foldchanges_reverse.csv'
source_fasta = '2_Development/Dictyostelium_ac.fa'

### MAIN ###

def main():

    csv_total = []

    #Get raw fold change data
    raw_fc = open(source_fc).read()
    fc_split = raw_fc.split('\n')

    #Iterate through foldchange thresholds
    count = 0
    frequencies = []
    while count < 1000:
        social_genes = []
        for i in fc_split[1:]:
            if i != '':
                items = i.split(',')
                gene_id = items[0]
                fc = items[1]
                if fc != '':
                    if float(fc) > count or fc == 'inf':
                        social_genes.append(gene_id)
        raw_fasta = open(source_fasta).read()
        fasta_split = raw_fasta.split('>')

        #Get genes with fold-change above threshold
        social_seqs = []
        seqs = []
        for k in fasta_split:
            if k != '':
                chunks = k.split('\n')
                line = chunks[0]
                seq = chunks[1]
                line_split = line.split(';')
                fasta_id = line_split[0].replace('ID=gene:', '')
                if fasta_id in social_genes:
                    social_seqs.append(seq)
                    seqs.append(seq)
                else:
                    seqs.append(seq)
        if len(social_seqs) > 100:
            #Calculate ASC frequencies and TAA usage for the group of genes at this threshold
            overall_social, positions_social, overall_stops_social, overall_codons_social = get_ASCs(social_seqs)
            overall, positions, overall_stops, overall_codons = get_ASCs(seqs)
            taa_social = get_TAA(social_seqs)
            taa_overall = get_TAA(seqs)
            frequencies.append([count, overall_social, overall, taa_social, taa_overall])
            print (count, len(social_seqs))
            #Add to counter to move to next foldchange threshold
            count += 5

    #Set headers and write output files
    headers = ['FC_threshold', 'social_freq', 'overall', 'taa_social', 'taa_overall']
    filename = "dicty_thresh_samplesize_reverse_taa.csv"
    subdir = "2_Development"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in frequencies:
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


def get_TAA(sequences):

    taa = 0

    for sequence in sequences:
        if sequence[0:3] == 'TAA':
            taa += 1

    freq = taa / len(sequences)

    return freq


### RUN ###
if __name__ == '__main__':
    main()
