### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FILES ###
source = '2_Development/Dev_RNAseq_Rosengarten_plus_veg_DESeqnorm_counts.txt'

### MAIN ###

def main():

    csv_total = []

    #Get raw file
    raw = open(source).read()
    split = raw.split('\n')

    #Get expression data of each gene during single and social growth stage
    for i in split[1:]:

        singles = []
        socials = []

        items = i.split('\t')
        singles.append(float(items[1]))
        singles.append(float(items[20]))
        for item in items[39:52]:
            singles.append(float(item))
        for social in items[2:20]:
            socials.append(float(social))
        for social2 in items[21:39]:
            socials.append(float(social2))

        #Get mean expression during both stages, calculate fold-change
        mean_single = np.mean(singles)
        mean_social = np.mean(socials)
        fc = mean_social / mean_single
        csv_total.append([items[0], fc])

    #Set headers and output file
    headers = ['Gene_ID', 'FC']
    filename = "foldchanges_reverse.csv"
    subdir = "2_Development"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main()
