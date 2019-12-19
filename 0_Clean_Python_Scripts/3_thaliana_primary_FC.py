### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FILES ###

source_pol = '4_Pollen/plntphys_pp.104.057935_57935Supplementary_Table_1.txt'
source_fasta = '4_Pollen/Arabidopsis_thaliana.fa'

### MAIN ###

def main():

    #Get pollen specific and depleted genes
    raw_pol = open(source_pol).read()
    pol_split = raw_pol.split('\n')

    #Get genes over and under 100 F/C expression
    under_100 = []
    over_100 = []

    for i in pol_split[1:]:
        items = i.split('\t')
        pollen_selective = items[0]
        fc = items[2]
        pollen_depleted = items[3]
        pollen_specific = items[4]
        id = items[6].upper()
        if fc != '':
            if float(fc) > 100:
                over_100.append(id)
            if 0 < float(fc) < 100:
                under_100.append(id)

    raw_fasta = open(source_fasta).read()
    fasta_split = raw_fasta.split('>')
    over_seqs = []
    under_seqs = []

    for k in fasta_split:
        if k != '':
            chunks = k.split('\n')
            line = chunks[0]
            seq = chunks[1]
            line_split = line.split(';')
            fasta_id = line_split[0].replace('ID=gene:', '')
            if fasta_id in over_100:
                over_seqs.append(seq)
            elif fasta_id in under_100:
                under_seqs.append(seq)

    #Calculate usage of TAA, TGA and TAG as stop codons in each set of genes
    oTAA, oTGA, oTAG = get_stop_usage(over_seqs)
    uTAA, uTGA, uTAG = get_stop_usage(under_seqs)

    #Write output file
    csv_total = [['over', oTAA, oTGA, oTAG], ['under', uTAA, uTGA, uTAG]]
    headers = ['Group', 'TAA', 'TGA', 'TAG']
    filename = "pollen_primarygroups.csv"
    subdir = "4_Pollen"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### FUNCTIONS ###

def get_stop_usage(sequences):


    taa = 0
    tga = 0
    tag = 0

    for sequence in sequences:
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

    return taa/len(sequences), tga/len(sequences), tag/len(sequences)

### RUN ###
if __name__ == '__main__':
    main()
