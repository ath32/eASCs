### IMPORTS ###

import os
import csv
import numpy as np

### CHOOSE SOURCE FOLDER ###

fasta_source = '9_Modeling/UTR_FASTA'
gff_source = '9_Modeling/GFF'

### MAIN ###

def main():

    csv_total = []

    for root, dirs, filenames in os.walk(fasta_source):
        for f in filenames:
            name = f.replace('.fa', '')
            fasta_path = os.path.join(fasta_source, f)
            raw_fasta = open(fasta_path).read()
            gc = get_gc(raw_fasta)
            csv_total.append([name, gc])

    for root, dirs, filenames in os.walk(gff_source):
        for f in filenames:
            split = f.split('_')
            genus = split[0]
            print (genus)
            gff_path = os.path.join(gff_source, f)
            raw_gff = open(gff_path).read()
            median_gene, median_igr = get_lengths(raw_gff, f)
            for entry in csv_total:
                if genus in entry[0]:
                    entry.append(median_gene)
                    entry.append(median_igr)

    headers = ['Name', 'GC', 'Median_gene_length', 'Median_3_IGR']

    #Write outputs
    filename = "new_correlates.csv"
    subdir = "9_Modeling"
    filepath = os.path.join(subdir, filename)
    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### FUNCTION ###

def get_gc(fasta):

    null = ''
    split = fasta.split('>')
    for chunk in split:
        if chunk != '':
            sections = chunk.split('\n')
            line = sections[0]
            sequence = sections[1]
            null += sequence
    gc = (null.count('G') + null.count('C')) / len(null)

    return gc


def get_lengths(gff, f):

    split = gff.split('\tgene\t')
    gene_sizes = []
    coords= []
    for section in split:
        if 'ID=gene:' in section or 'gene_id' in section:
            if 'supercontig' not in section and 'chromosome' not in section and 'scaffold' not in section.replace('scaffold_32', ''):
                smaller_split = section.split('\t')
                coord1 = smaller_split[0]
                coord2 = smaller_split[1]
                gene_size = float(coord2) - float(coord1)
                gene_sizes.append(gene_size)
                coords.append([float(coord1), float(coord2), smaller_split[3]])

    igrs = []
    for i,n in enumerate(coords):
        if coords[i][2] == '+':
            if i != len(coords)-1:
                IGR = coords[i+1][0] - coords[i][1]
            elif coords[i][2] == '-':
                if i != 0:
                    IGR = coords[i][0] - coords[i-1][1]
            igrs.append(IGR)
    median_gene = np.median(gene_sizes)
    median_igr = np.median(igrs)

    return median_gene, median_igr


### RUN ###

if __name__ == '__main__':
    main()
