### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '7_Expression/LEGs_FASTA'

### MAIN ###

def main(source):

    csv_total = []
    headers = ['Genome', 'TAA', 'TGA', 'TAG', 'Total']

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            name = f.replace('.csv', '')
            path = os.path.join(source, f)
            raw = open(path).read()
            split = raw.split('>')

            #Count usage of each stop codon
            taa = 0
            tga = 0
            tag = 0
            total = 0
            sequences = []
            for section in split:
                if section != '':
                    line_split = section.split('\n')
                    sequence = line_split[1]
                    sequences.append(sequence)
            for seq in sequences:
                if seq[:3] == 'TAA':
                    taa += 1
                if seq[:3] == 'TGA':
                    tga += 1
                if seq[:3] == 'TAG':
                    tag += 1
            total = len(sequences)
            csv_total.append([name, taa, tga, tag, total])

    #Write output file
    filename = "pooled_leg_usage.csv"
    subdir = "7_Expression"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###
if __name__ == '__main__':
  main(source)
