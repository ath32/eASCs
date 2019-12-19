## IMPORTS ###

import os
import re
import numpy as np
import csv

### SOURCE ###
fasta = './7_Expression/Final_Unicells'
expression = './7_Expression/Eukaryotes'

### MAIN ###

def main():

    genomes = []
    csv_total = []
    headers = ['hegs_hits', 'hegs_no_hit', 'hegs_total_genes', 'legs_hits', 'legs_no_hit', 'legs_total_genes', 'asc_hegs_hits', 'asc_hegs_no_hit', 'asc_hegs_total', 'asc_legs_hits', 'asc_legs_no_hit', 'asc_legs_total']
    pooled_HEGs = ''
    pooled_LEGs = ''

    for root, dirs, filenames in os.walk(fasta):
        for f in filenames:
            genomes.append(f.replace('.fa', ''))

    for genome in genomes:
        exp_source = expression + '/' + genome

        for root, dirs, filenames in os.walk(expression):
            for f in filenames:
                if exp_source.replace('.', '') in root:
                    if 'WHOLE_ORGANISM' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '7165-GPM_201408' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '5061-GPM_201408' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '9615-GPM_2012_09_Canis_familiaris' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '3055-GPM_201408' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '44689-GPM_201408' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '9796-PA_201304' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '39947-GPM_2012_09_Oryza_sativa' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '9598-Chimp_iBAQ_Khan_2013' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '4081-Solanum_lycopersicum_SC_biomart_18197_E__1reps' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '9823-Pig_PeptideAtlas_2011-08_Ens62' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '5691-GPM_201408' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '8364-GPM_2012_09_Xenopus_tropicalis' in f:
                        file_path = os.path.join(exp_source, f)
                    elif '4577-LEAF-integrated' in f:
                        file_path = os.path.join(exp_source, f)

        for root, dirs, filenames in os.walk(fasta):
            for f in filenames:
                if genome in f:
                    genome_path = os.path.join(fasta, f)
                    raw_fasta = open(genome_path).read()

                    #Pool HEGs and LEGs together into one 'FASTA' file each
                    raw = open(file_path).read()
                    HEGs, LEGs = get_expression(raw)
                    HEGs_fasta, LEGs_fasta = get_new_txt(raw_fasta, HEGs, LEGs)
                    pooled_HEGs += HEGs_fasta
                    pooled_LEGs += LEGs_fasta

    #Calculate +4T and ASC frequency for each group of genes
    hegs_hits, hegs_no_hit, hegs_total_genes = get_t(pooled_HEGs)
    legs_hits, legs_no_hit, legs_total_genes = get_t(pooled_LEGs)
    asc_hegs_hits, asc_hegs_no_hit, asc_hegs_total = get_ASCs(pooled_HEGs)
    asc_legs_hits, asc_legs_no_hit, asc_legs_total = get_ASCs(pooled_LEGs)
    csv_total.append([hegs_hits, hegs_no_hit, hegs_total_genes,legs_hits, legs_no_hit, legs_total_genes, asc_hegs_hits, asc_hegs_no_hit, asc_hegs_total, asc_legs_hits, asc_legs_no_hit, asc_legs_total])

    #Write output file
    filename = "expression_pooled_uni_tag.csv"
    subdir = "7_Expression"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### FUNCTIONS ###

def get_expression(raw_paxdb):

    ''' Obtain the IDs of HEGs and LEGs '''

    #Create nested list for pax data, incl locus tag and protein expression
    raw = []

    #Create list, through which we can access the information for each gene
    if 'raw_spectral_count' in raw_paxdb:
        i = raw_paxdb.split('raw_spectral_count\n')
        i1 = i[1].split('\n')
        del i1[-1]

        for i in i1:
            chunked = i.split('\t')
            internal_id = chunked[0]
            external_id = chunked[1]
            abundance = chunked[2]
            raw_count = chunked[3]
            raw.append([internal_id, external_id, float(abundance)])

    else:
        i = raw_paxdb.split('abundance\n')
        i1 = i[1].split('\n')
        del i1[-1]

        for i in i1:
            chunked = i.split('\t')
            internal_id = chunked[0]
            external_id = chunked[1]
            abundance = chunked[2]
            raw.append([internal_id, external_id, float(abundance)])

    srted = sorted(raw, key = lambda x: float(x[2]))

    #Calculate top and botton quartiles
    abundances = []

    for a in srted:
        abundances.append(a[2])

    array = np.array(abundances)
    top = np.percentile(array, 75)
    bottom = np.percentile(array, 25)

    #Create two lists - HEGs (top 25%) and LEGs (bottom 25%)
    HEGs_values = []
    HEGs_IDs = []

    for b in srted:
        if b[2] >= top:
            HEGs_values.append(b[2])

    for c in srted:
        if c[2] in HEGs_values:
            HEGs_IDs.append(c[1])

    LEGs_values = []
    LEGs_IDs = []

    for d in srted:
        if d[2] <= bottom:
            LEGs_values.append(d[2])

    for e in srted:
        if e[2] in LEGs_values:
            LEGs_IDs.append(e[1])

    return HEGs_IDs, LEGs_IDs


def get_new_txt(raw_fasta, HEGs, LEGs):

    ''' Build new FASTA files - one for HEGs and one for LEGs for each genome '''

    HEGs_fasta = ''
    LEGs_fasta = ''

    refined_HEGs = []
    refined_LEGs = []

    for heg in HEGs:
        split = heg.split('.')
        id1 = ".".join(split[1:])
        refined_HEGs.append(id1)

    for leg in LEGs:
        split = leg.split('.')
        id2 = ".".join(split[1:])
        refined_LEGs.append(id2)

    #Get chunks to iterate through
    chunks = raw_fasta.split('>')
    chunks.pop(0)

    for i in chunks:
        a = i.split('\n')
        a1 = a[0].split(';')
        refined_id = "".join(a1[0].replace('ID=gene:', '')).strip().replace(':pep', '')
        sequence = a[1]
        primary_stop = sequence[0:3]
        if primary_stop == 'TAG':
            for j in refined_HEGs:
                if refined_id == j:
                    HEGs_fasta += ">" + i

            for k in refined_LEGs:
                if refined_id == k:
                    LEGs_fasta += ">" + i

    return HEGs_fasta, LEGs_fasta


def get_t(genes):

    t = 0
    total = 0
    sections = genes.split('>')
    for section in sections:
        if section != '':
            lines = section.split('\n')
            sequence = lines[1]
            if sequence[3] == 'T':
                t += 1
                total += 1
            else:
                total += 1
    freq = t / total
    return t, total-t, total


def get_ASCs(genes):
    asc = 0
    total = 0
    sections = genes.split('>')
    for section in sections:
        if section != '':
            lines = section.split('\n')
            sequence = lines[1]
            codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
            if codons[1] == 'TGA' or codons[1] == 'TAA' or codons[1] == 'TAG':
                asc += 1
                total += 1
            else:
                total += 1
    freq = asc / total
    return asc, total-asc, total


### RUN ###

if __name__ == '__main__':
    main()
