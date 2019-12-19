### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FILES ###

source_pol = '4_Pollen/plntphys_pp.104.057935_57935Supplementary_Table_1.txt'
source_fasta = '4_Pollen/Arabidopsis_thaliana.fa'

### MAIN ###

def main():

    csv_total = []

    #Get raw file
    raw_pol = open(source_pol).read()
    pol_split = raw_pol.split('\n')

    frequencies = []
    count = 0
    while count < 200:
        #Get pollen gene ids over F/C threshold
        pollen_genes = []
        for i in pol_split[1:]:
            items = i.split('\t')
            pollen_selective = items[0]
            fc = items[2]
            pollen_depleted = items[3]
            pollen_specific = items[4]
            id = items[6].upper()
            if fc != '':
                if float(fc) > count:
                    pollen_genes.append(id)

        #Get sequences
        raw_fasta = open(source_fasta).read()
        fasta_split = raw_fasta.split('>')
        pol_seqs = []
        seqs = []

        for k in fasta_split:
            if k != '':
                chunks = k.split('\n')
                line = chunks[0]
                seq = chunks[1]
                line_split = line.split(';')
                fasta_id = line_split[0].replace('ID=gene:', '')
                if fasta_id in pollen_genes:
                    pol_seqs.append(seq)
                    seqs.append(seq)
                else:
                    seqs.append(seq)

        #Ensure sample size is large enough
        if len(pol_seqs) > 100:
            #Calculate ASC frequencies for each position
            overall_pollen, positions_pollen, overall_stops_pollen, overall_codons_pollenG = get_ASCs(pol_seqs)
            overall, positions, overall_stops, overall_codons = get_ASCs(seqs)
            taa_pollen = get_TAA(pol_seqs)
            taa_overall = get_TAA(seqs)
            frequencies.append([count, overall_pollen, overall, taa_pollen, taa_overall])
            print (count, len(pol_seqs))
        count += 5

    #Write output file
    headers = ['FC_threshold', 'pollen_freq', 'overall', 'taa_social', 'taa_overall']
    filename = "pollen_thresh_samplesize_taa.csv"
    subdir = "4_Pollen"
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
