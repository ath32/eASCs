### IMPORTS ###

import os
import numpy as np
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '6_Ne/Extra_genomes_2/Chi2'

### MAIN ###

def main(source):

    csv_total = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            name = f.replace('.csv', '')
            path = os.path.join(source, f)
            raw = open(path).read()
            lines = raw.split('\n')

            scores = []
            for line in lines[2:]:
                if line != '':
                    chunks = line.split(',')
                    observed = chunks[0]
                    expected = chunks[1]
                    excess = float(chunks[2]) / 100
                    scores.append(excess)

            mean_score = np.mean(np.array(scores))
            csv_line = [name, mean_score]
            csv_total.append(csv_line)

    #Write output file
    headers = ['Species', 'Enrichment_score']
    filename = "scores_w6.csv"
    subdir = "6_Ne/Extra_genomes_2"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main(source)
