### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FOLDER ###

source = '7_Expression/Chi1_LEGs'

### MAIN ###

def main(source):

    counts = [0, 0, 0, 0, 0]
    total_count = 0
    total_no = 0

    csv_total = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            #Get raw information
            path = os.path.join(source, f)
            raw = open(path).read()
            rows = raw.split('\n')

            #Count number of positions containing ASC enrichment
            local_count = 0

            for i,n in enumerate(rows[3:len(rows)-1]):
                split = n.split(',')
                if float(split[2]) > 0 and float(split[4]) < (0.05/5):
                    counts[i] += 1
                    local_count += 1

            #If enriched, add to enriched genomes
            if local_count > 0:
                total_count += 1
                print (f, 'ok')
            #If not enriched, add to not-enriched genomes
            else:
                total_no += 1

    #Write output files
    csv_total = [counts[0], counts[1], counts[2], counts[3], counts[4], total_count, total_no]
    headers = ['P2', 'P3', 'P4', 'P5', 'P6', 'Total_enrich', 'Total_no_enrich']
    filename = "LEGs_Ch1.2.csv"
    subdir = "7_Expression"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        writer.writerow(k for k in csv_total)

### RUN ###

if __name__ == '__main__':
    main(source)
