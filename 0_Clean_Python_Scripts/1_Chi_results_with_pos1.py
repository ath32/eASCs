### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FOLDER ###

source = '7_Expression/Chi1_LEGs'

### MAIN ###

def main(source):

    counts = [0, 0, 0, 0, 0, 0]
    total_count = 0
    total_no = 0

    csv_total = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            #Get raw file
            path = os.path.join(source, f)
            raw = open(path).read()
            #Split into rows to get each position
            rows = raw.split('\n')

            #Count how many positions have significant enrichment
            local_count = 0

            for i,n in enumerate(rows[2:len(rows)-1]):
                split = n.split(',')
                if float(split[2]) > 0 and float(split[4]) < (0.05/6):
                    counts[i] += 1
                    local_count += 1

            #If more than one, count the genome as enriched
            if local_count > 0:
                total_count += 1
                print (f, 'ok')
            else:
            #Else count the genome as not enriched
                total_no += 1

    #Write output file
    csv_total = [counts[0], counts[1], counts[2], counts[3], counts[4], counts[5], total_count, total_no]
    headers = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'Total_enrich', 'Total_no_enrich']
    filename = "LEGs_Ch1.1.csv"
    subdir = "7_Expression"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        writer.writerow(k for k in csv_total)

### RUN ###

if __name__ == '__main__':
    main(source)
