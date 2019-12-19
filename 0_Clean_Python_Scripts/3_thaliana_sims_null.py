### IMPORTS ###
import os
import numpy as np
import csv
import time
import scipy.stats as st
import multiprocessing

### CHOOSE SOURCE FOLDER ###

fasta_source = os.path.abspath("../2_Project_FSM/4_Pollen/Arabidopsis_thaliana.fa")
csv_source = os.path.abspath("../2_Project_FSM/4_Pollen/plntphys_pp.104.057935_57935Supplementary_Table_1.csv")

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
    EF_hits, EF_array, EF_list = get_simulations(utr_list, total_utr)

    prob = np.mean(EF_array)
    s_hits = [prob * len(social_genes), prob * len(social_genes), prob * len(social_genes), prob * len(social_genes), prob * len(social_genes), prob * len(social_genes), prob * len(social_genes), prob * len(social_genes), prob * len(social_genes), prob * len(social_genes), prob * len(social_genes)]
    v_hits = [prob * len(veg_genes), prob * len(veg_genes), prob * len(veg_genes), prob * len(veg_genes), prob * len(veg_genes), prob * len(veg_genes), prob * len(veg_genes), prob * len(veg_genes), prob * len(veg_genes), prob * len(veg_genes), prob * len(veg_genes)]

    #Set hits arrays
    sOF_a = np.array(sOF_hits)
    vOF_a = np.array(vOF_hits)

    #Calculate excess
    sexcess = ((sOF_a- s_hits) / s_hits) * 100
    sexcess_list = sexcess.tolist()
    vexcess = ((vOF_a- v_hits) / v_hits) * 100
    vexcess_list = vexcess.tolist()

    #Chisq
    schisq = ((sOF_a - s_hits) ** 2) / s_hits
    schisq_list = schisq.tolist()
    vchisq = ((vOF_a - v_hits) ** 2) / v_hits
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

    #Write output files
    sline = [list(a) for a in zip(sOF_hits, s_hits, sexcess_list, schisq_list, sp_list)]
    vline = [list(b) for b in zip(vOF_hits, v_hits, vexcess_list, vchisq_list, vp_list)]
    headers = ["Observed", "Expected", "Excess%", "Chisq", "p-value"]
    create_csv(headers, sline, 'specific')
    create_csv(headers, vline, 'absent')

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


def get_seqs(source_pol, source_fasta):

    #Get pollen specific genes
    raw_pol = open(source_pol).read()
    pol_split = raw_pol.split('\n')

    pollen_genes = []
    nonpollen_genes = []

    for i in pol_split:
        items = i.split(',')
        pollen_selective = items[0]
        pollen_depleted = items[3]
        pollen_specific = items[4]
        id = items[6].upper()
        if pollen_selective == 'X':
            pollen_genes.append(id)
        if pollen_depleted == 'X':
            nonpollen_genes.append(id)

    #Get sequences for both sets
    raw_fasta = open(source_fasta).read()
    fasta_split = raw_fasta.split('>')

    pol_seqs = []
    nonpol_seqs = []

    for k in fasta_split:
        if k != '':
            chunks = k.split('\n')
            line = chunks[0]
            seq = chunks[1]
            codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
            line_split = line.split(';')
            fasta_id = line_split[0].replace('ID=gene:', '')
            if fasta_id in pollen_genes:
                pol_seqs.append(codons)
            elif fasta_id in nonpollen_genes:
                nonpol_seqs.append(codons)
            else:
                continue

    return pol_seqs, nonpol_seqs


def get_stop_freq(utr_list):

    ''' Get stop frequencies at each position '''

    codons = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    for i in utr_list:

        if i[0] == 'TAG' or i[0] == 'TGA' or i[0] == 'TAA':
            codons[0] += 1
        if i[1] == 'TAG' or i[1] == 'TGA' or i[1] == 'TAA':
            codons[1] += 1
        if i[2] == 'TAG' or i[2] == 'TGA' or i[2] == 'TAA':
            codons[2] += 1
        if i[3] == 'TAG' or i[3] == 'TGA' or i[3] == 'TAA':
            codons[3] += 1
        if i[4] == 'TAG' or i[4] == 'TGA' or i[4] == 'TAA':
            codons[4] += 1
        if i[5] == 'TAG' or i[5] == 'TGA' or i[5] == 'TAA':
            codons[5] += 1
        if i[6] == 'TAG' or i[6] == 'TGA' or i[6] == 'TAA':
            codons[6] += 1
        if i[7] == 'TAG' or i[7] == 'TGA' or i[7] == 'TAA':
            codons[7] += 1
        if i[8] == 'TAG' or i[8] == 'TGA' or i[8] == 'TAA':
            codons[8] += 1
        if i[9] == 'TAG' or i[9] == 'TGA' or i[9] == 'TAA':
            codons[9] += 1
        if i[10] == 'TAG' or i[10] == 'TGA' or i[10] == 'TAA':
            codons[10] += 1

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


def get_simulations(utr_list, total_utr):

    ''' Simulate n UTR sequences for each genome based upon dinucleotide content of total UTR regions '''

    #Probability of first base - dictionary
    nt_dict = {}
    length = len(total_utr)

    for i in range(length):
        base = total_utr[i]
        if base not in nt_dict:
            nt_dict[base] = 1
        else:
            nt_dict[base] += 1

    nt_freq = {k: v / length for k, v in nt_dict.items()}

    #Second base / Next base
    trans = {}
    for i in range(len(total_utr)-1):

        dinuc = total_utr[i:i+2]
        first_dinuc = dinuc[0]
        sec_dinuc = dinuc[1]

        if first_dinuc not in trans:
            trans[first_dinuc] = [sec_dinuc]
        else:
            trans[first_dinuc] += sec_dinuc

    #Generate simulations
    total_list = []
    count = 0

    #This while loop determines how many simulations will be completed
    while count < 5000:

        sim_n = []
        codon_list = []

        sim = []

        np.random.seed()

        #Calculate next base
        first_base = np.random.choice(['A', 'C', 'G', 'T'], p=[nt_freq['A'], nt_freq['C'], nt_freq['G'], (1-nt_freq['A']-nt_freq['C']-nt_freq['G'])])
        sim.append(first_base)

        #For one gene simulation
        while (len(sim) <= 30):

            #Generate next base
            prev_base = sim[-1]
            next_base = np.random.choice(trans[prev_base], replace=True)
            sim.append(next_base)

        sim = "".join(sim)
        codon_seq = [sim[i:i+3] for i in range(0, len(sim), 3)]
        total_list.append(codon_seq)
        count += 1

    x, n_array, n_list = get_stop_freq(total_list)

    return x, n_array, n_list


def create_csv(headers, csv_total, title):

    ''' Output results to CSV file '''

    filename = title + "_sims.csv"
    subdir = "4_Pollen"
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)


### RUN ###
if __name__ == '__main__':
    main()
