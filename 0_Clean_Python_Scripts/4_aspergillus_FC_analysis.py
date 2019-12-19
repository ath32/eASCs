### IMPORTS ###

import os
import numpy as np
import csv
import pandas as pd
import scipy.stats as st

### CHOOSE SOURCE FILES - unhash genome of interest ###

#A. Niger
source_genes = '8_Conidia/CAGs_niger.txt'
source_fasta = '8_Conidia/Filtered_FASTA/Aspergillus_niger.fa'

#A. Oryzae
# source_genes = '8_Conidia/CAGs_oryzae.txt'
# source_fasta = '8_Conidia/Filtered_FASTA/Aspergillus_oryzae.fa'

### MAIN ###

def main():

    csv_total = []

    #Get gene ids
    df = raw_data = open(source_genes).read()
    split = df.split('\n')

    #Get sequences for both sets
    raw_fasta = open(source_fasta).read()
    fasta_split = raw_fasta.split('>')

    frequencies = []
    count = 0
    while count < 1000:
        CAGs = []
        for line in split[1:]:
            line_split = line.split('\t')
            id = line_split[0]
            fc = line_split[6]
            if fc == '#DIV/0!' or float(fc) > count:
                CAGs.append(id)

        CAG_seqs = []
        nonCAG_seqs = []
        seqs = []
        for k in fasta_split:
            if k != '':
                chunks = k.split('\n')
                line = chunks[0]
                seq = chunks[1]
                line_split = line.split(';')
                fasta_id = line_split[0].replace('ID=gene:', '')
                if fasta_id in CAGs:
                    CAG_seqs.append(seq)
                    seqs.append(seq)
                elif fasta_id not in CAGs:
                    nonCAG_seqs.append(seq)
                    seqs.append(seq)

        if len(CAG_seqs) > 100:
            #Calculate ASC frequencies for each position
            overall_CAG, positions_CAG, overall_stops_CAG, overall_codons_CAG = get_ASCs(CAG_seqs)
            overall_nonCAG, positions_nonCAG, overall_stops_nonCAG, overall_codons_nonCAG = get_ASCs(nonCAG_seqs)
            overall, positions, overall_stops, overall_codons = get_ASCs(seqs)
            taa_CAG = get_TAA(CAG_seqs)
            taa_overall = get_TAA(seqs)
            frequencies.append([count, overall_CAG, overall, taa_CAG, taa_overall])

            #Do chisqs - two ways...
            observed = overall_stops_CAG
            expected = overall_codons_CAG * overall_nonCAG
            chisq = ((observed - expected) ** 2) / expected
            p = st.distributions.chi2.sf (chisq, 1)
            if p < 0.05:
                print (count, p)

            #Contingency table
            table = np.array([[overall_stops_CAG, overall_codons_CAG - overall_stops_CAG], [overall_stops_nonCAG, overall_codons_nonCAG - overall_stops_nonCAG]])
            p = st.chi2_contingency(table)[1]
            if p < 0.05:
                print (count, p)
        count += 5

    #Write output file
    headers = ['FC_threshold', 'CAG_freq', 'overall', 'taa_CAG', 'taa_overall']
    filename = "CAGs_thresh_Niger_taa.csv"
    subdir = "8_Conidia"
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
