### IMPORTS ###

import os
import csv

### CHOOSE SOURCE FOLDER (EMBL) ###

source = '5_Protists/Protist_UTR'

### MAIN ###

def main(source):

    total_csv = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            #Get name
            path = os.path.join(source, f)
            raw = open(path).read()
            name = f.replace('.fa', '')

            #Get genes, split into codons
            split = raw.split('>')
            sequences = []
            for i in split:
                if i != '':
                    line_split = i.split('\n')
                    info = line_split[0]
                    utr = line_split[1]
                    codons = [utr[i:i+3] for i in range(0, len(utr), 3)]
                    sequences.append(codons)

            #Calculate frequency at position +1
            count1 = 0
            count_null = 0

            for sequence in sequences:
                position_1 = sequence[1]
                position_3 = sequence[3]
                position_4 = sequence[4]
                position_5 = sequence[5]
                position_6 = sequence[6]
                if position_1[1] == 'T':
                    count1 += 1
                if position_3[1] == 'T':
                    count_null += 1
                if position_4[1] == 'T':
                    count_null += 1
                if position_5[1] == 'T':
                    count_null += 1
                if position_6[1] == 'T':
                    count_null += 1

            #Convert to frequencies
            one_freq = count1 / len(sequences)
            null_freq = count_null / (len(sequences) * 4)
            #Append to output file
            total_csv.append([name, one_freq, null_freq])

    #Write output file
    headers = ["Accession", "one", "null"]
    get_CSV(total_csv, headers)

### FUNCTIONS ###

def get_CSV(total, headers):

    filename = "+4T_Protist.csv"
    subdir = ""
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        for j in total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    main(source)
