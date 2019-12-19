### IMPORTS ###
import os
import numpy as np
import csv
import time
import scipy.stats as st
import multiprocessing

### CHOOSE SOURCE FOLDER (EMBL) ###

fasta_source = os.path.abspath("../2_Project_FSM/2_Development/Dictyostelium_ac.fa")
csv_source = os.path.abspath("../2_Project_FSM/2_Development/sociality_genes.csv")

### MAIN ###

def main():

    csv_total = []

    #First, get list of sociality and vegetative genes
    social_genes, veg_genes = get_seqs(csv_source, fasta_source)

    #Get stop frequencies (O)
    sOF_hits, sOF_array, sOF_list = get_stop_freq(social_genes)
    vOF_hits, vOF_array, vOF_list = get_stop_freq(veg_genes)

    #Get whole genome stuff
    path = os.path.join(fasta_source)
    raw = open(path).read()
    genes = raw.strip().split('>')
    genes_c = list(filter(None, genes))
    glist = get_list(genes_c)

    #Obtain UTR list and generate total UTR string, which is useful for later functions
    utr_list, total_utr = get_utr_stuff(glist)

    #Get simulated stop frequencies (E)
    sEF_hits, sEF_array, sEF_list = get_null(utr_list, total_utr, len(social_genes))
    vEF_hits, vEF_array, vEF_list = get_null(utr_list, total_utr, len(veg_genes))

    #Set hits arrays
    sOF_a = np.array(sOF_hits)
    vOF_a = np.array(vOF_hits)
    sEF_a = np.array(sEF_hits)
    vEF_a = np.array(vEF_hits)

    #Calculate excess
    sexcess = ((sOF_a- sEF_a) / sEF_a) * 100
    sexcess_list = sexcess.tolist()
    vexcess = ((vOF_a- vEF_a) / vEF_a) * 100
    vexcess_list = vexcess.tolist()

    #Chisq
    schisq = ((sOF_a - sEF_a) ** 2) / sEF_a
    schisq_list = schisq.tolist()
    vchisq = ((vOF_a - vEF_a) ** 2) / vEF_a
    vchisq_list = vchisq.tolist()

    #p-vals
    sp_list = []
    for i in schisq_list:
        sp = st.distributions.chi2.sf (i, 1)
        sp_list.append(sp)

    vp_list = []
    for j in vchisq_list:
        vp = st.distributions.chi2.sf (j, 1)
        vp_list.append(vp)

    #Write outputs
    sline = [list(a) for a in zip(sOF_hits, sEF_hits, sexcess_list, schisq_list, sp_list)]
    vline = [list(b) for b in zip(vOF_hits, vEF_hits, vexcess_list, vchisq_list, vp_list)]
    headers = ["Observed", "Expected", "Excess%", "Chisq", "p-value"]
    create_csv(headers, sline, 'sociality')
    create_csv(headers, vline, 'vegetative')

### FUNCTIONS ###

def get_list(genes_c):

    ''' Split UTR sequence into codons '''

    utr_list = []

    for i in genes_c:
        b = i.split('\n')
        sequence = b[1].upper()

        if sequence[0:3] == 'TGA' or sequence[0:3] == 'TAA' or sequence[0:3] == 'TAG':
            utr_sequence = sequence.upper()
            utr_seq = [utr_sequence[i:i+3] for i in range(0, len(utr_sequence), 3)]
            utr_list.append(utr_seq)

    return utr_list


def get_seqs(source_dicty, source_fasta):

    #Get gene ids
    raw_data = open(source_dicty).read()
    dicty_split = raw_data.split('\n')

    social_genes = []
    single_genes = []

    for i in dicty_split:
        items = i.split(',')
        identifier = items[1]
        gene_id = items[0]

        if identifier == '1':
            social_genes.append(gene_id)
        if identifier == '0':
            single_genes.append(gene_id)

    #Get sequences for both sets
    raw_fasta = open(source_fasta).read()
    fasta_split = raw_fasta.split('>')

    social_seqs = []
    single_seqs = []

    for k in fasta_split:
        if k != '':
            chunks = k.split('\n')
            line = chunks[0]
            seq = chunks[1].upper()
            codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
            line_split = line.split(';')
            fasta_id = line_split[0].replace('ID=gene:', '')
            if fasta_id in social_genes:
                if codons[0] == 'TGA' or codons[0] == 'TAA' or codons[0] == 'TAG':
                    social_seqs.append(codons)
            elif fasta_id in single_genes:
                if codons[0] == 'TGA' or codons[0] == 'TAA' or codons[0] == 'TAG':
                    single_seqs.append(codons)
            else:
                continue

    return social_seqs, single_seqs


def get_stop_freq(utr_list):

    ''' Get stop frequencies at each position '''

    codons = [0, 0, 0, 0, 0, 0, 0]
    stops = ['TGA', 'TAA', 'TAG']

    for i in utr_list:

        if i[0] in stops:
            codons[0] += 1
        if i[1] in stops:
            codons[1] += 1
        if i[2] in stops and i[1] not in stops:
            codons[2] += 1
        if i[3] in stops and i[1] not in stops and i[2] not in stops:
            codons[3] += 1
        if i[4] in stops and i[1] not in stops and i[2] not in stops and i[3] not in stops:
            codons[4] += 1
        if i[5] in stops and i[1] not in stops and i[2] not in stops and i[3] not in stops and i[4] not in stops:
            codons[5] += 1
        if i[6] in stops and i[1] not in stops and i[2] not in stops and i[3] not in stops and i[4] not in stops and i[5] not in stops:
            codons[6] += 1

    frequencies_array = np.array(codons) / len(utr_list)
    frequencies_list = frequencies_array.tolist()

    return codons, frequencies_array, frequencies_list


def get_utr_stuff(list):

    ''' Obtain UTR sequences from my FASTA format '''

    utr_list = []
    total_utr = ''

    #For each gene in the FASTA file...
    for i in list:
        if i[0] == 'TGA' or i[0] == 'TAA' or i[0] == 'TAG':
            utr_list.append(i)
            seq = "".join(i)
            #Now create total UTR list which will be useful later
            total_utr += seq[1:]

    return utr_list, total_utr


def get_null(utr_list, total_utr, n_genes):

    #Split total_utr into codons, to assess ASCs
    codons = [total_utr[i:i+3] for i in range(0, len(total_utr), 3)]

    #Set counters, n codons and ASCs
    n_codons = 0
    ASCs = 0

    for codon in codons:
        if codon == 'TGA' or codon == 'TAA' or codon == 'TAG':
            ASCs += 1
            n_codons += 1
        else:
            n_codons += 1

    p = ASCs / n_codons
    EF_frequencies = []

    position = 0

    while position < 11:

        frequency = p * ((1-p) ** (position - 1))
        EF_frequencies.append(frequency)
        position += 1

    position_zero = 1
    EF_frequencies.insert(0,position_zero)
    EF_f_array = np.array(EF_frequencies[:7])
    EF_hits = EF_f_array * n_genes
    print (EF_hits, n_genes)
    EF_hits_list = EF_hits.tolist()

    return EF_hits_list, EF_f_array, EF_frequencies


def create_csv(headers, csv_total, title):

    ''' Output results to CSV file '''

    filename = title + ".csv"
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
