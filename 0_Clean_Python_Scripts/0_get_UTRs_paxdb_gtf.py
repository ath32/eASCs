### IMPORTS ###

import os
import generic as gen
import re
import collections
import copy
import numpy as np
import random
import shutil
from pathlib import Path


### CHOOSE SOURCE FOLDER (EMBL) ###

gff_source = '7_Expression/Release-17/GTF'
wgs_source = '7_Expression/Release-17/WGS'
fasta_source = '7_Expression/Release-17/FASTA'

### MAIN ###

def main():

    #Get one genome per genus
    independent_genomes = one_per_genus(gff_source)

    #Get BED file and FASTA file for each genome
    for root, dirs, filenames in os.walk(gff_source):
        for f in filenames:

            #First get name
            name_split = f.split('.')
            name = name_split[0]

            if name in independent_genomes:
                bed_path = get_utr_bed(name, gff_source, f)
                get_utr_fasta(name, bed_path, wgs_source)

    #Filter resultant FASTA files to ensure (i) the genome is valid (>70% of UTRs start with a stop) (ii) all non-stop starting UTRs are removed
    print ('filtering fasta files...')
    for root, dirs, filenames in os.walk(fasta_source):
        for f in filenames:
            path = os.path.join(fasta_source, f)
            raw = open(path).read()

            valid_seq = 0
            total = 0
            new_fasta = ''

            file_split = raw.split('>')
            for i in file_split:
                if i != '':
                    if ';' in i:
                        total += 1
                        section_split = i.split('\n')
                        line = section_split[0]
                        sequence = section_split[1]
                        if sequence[0:3] == 'TGA' or sequence[0:3] == 'TAA' or sequence[0:3] == 'TAG':
                            if NucleotideCheck(sequence) == True:
                                valid_seq += 1
                                new_line = '>' + line + '\n' + sequence + '\n'
                                new_fasta += new_line

            percent_valid = valid_seq / total
            new_source = "7_Expression/Release-17/Filtered_FASTA/"
            new_fasta_path = os.path.join(new_source, f)

            if percent_valid > 0.70:
                f = open(new_fasta_path, 'a')
                f.write(new_fasta)
                f.close()

### FUNCTIONS ###

def NucleotideCheck(sequence):
    return all(base.upper() in ('A', 'C', 'T', 'G') for base in sequence)


def one_per_genus(gff_source):

    genus_list = []
    names_list = []

    for root, dirs, filenames in os.walk(gff_source):
        for f in filenames:
            print ('filtering for phylogenetic independence...')
            name_split = f.split('.')
            name = name_split[0]
            name_split2 = name.split('_')
            genus = name_split2[0]
            species = name_split2[1]
            if genus not in genus_list:
                names_list.append(name)
                genus_list.append(genus)

    return names_list


def get_utr_bed(name, source, file):

    print ('getting bed file...')

    path = os.path.join(source, file)
    raw = open(path).read()

    #Split into sections & IGR filter to remove overlapping genes
    split = raw.split('start_codon')
    filtered_ids = IGR_filter(split)

    if filtered_ids == '':
        print (name)

    #For each gene in the filtered list, get exon annotations
    bed = get_utr(split, filtered_ids)

    #Output bed
    bed_path = os.path.abspath("../2_Project_FSM/7_Expression/Release-17/BED/" + name + ".bed")

    f = open(bed_path, 'a')
    f.write(bed)
    f.close()

    return bed_path


def IGR_filter(list):

    print ('filtering genes for sufficient IGR...')

    ids = []
    temp_ids = []

    #For each set of coordinates...
    for i in list:
        starting_list = i.split('\t')
        start_start = starting_list[1]
        start_stop = starting_list[2]
        strand = starting_list[4]

        stopping_list = i.split('stop_codon')
        for section in stopping_list[1:]:
            tab_split = section.split('\t')
            stop_start = tab_split[1]
            stop_stop = tab_split[2]
            gene_id_list = tab_split[6].strip().split('"')
            gene_id = gene_id_list[1]
            if strand == '+':
                ids.append([gene_id, int(start_start)-1, int(stop_stop), strand])
            elif strand == '-':
                ids.append([gene_id, int(stop_stop)-1, int(start_start), strand])

    #Set empty list of filtered ids
    gene_samples = []

    #Calculate 3' UTRs and filter gene ids
    for i,n in enumerate(ids):
        if ids[i][3] == '+':
            if i != len(ids)-1:
                IGR = ids[i+1][1] - ids[i][2]
                if 100 < IGR:
                    gene_samples.append(ids[i][0])

        elif ids[i][3] == '-':
            if i != 0:
                IGR = ids[i][1] - ids[i-1][2]
                if 100 < IGR:
                    gene_samples.append(ids[i][0])

    return gene_samples


def get_utr(list, filtered_ids):

    total = ''
    useful_lines = []

    for i in list[1:]:
        starting_list = i.split('\t')
        start_start = starting_list[1]
        start_stop = starting_list[2]
        strand = starting_list[4]
        stopping_list = i.split('stop_codon')
        chromosome = i.split('\n')[1].split('\t')[0]
        for section in stopping_list[1:]:
            if 'protein_id' in section:
                tab_split = section.split('\t')
                stop_start = tab_split[1]
                stop_stop = tab_split[2]
                gene_id_list = tab_split[6].strip().split('"')
                gene_id = gene_id_list[1]
                protein_list = tab_split[-3].strip().split('"')
                protein_id = protein_list[-2]
                if gene_id in filtered_ids:
                    if strand == '+':
                        bed_start = str(int(stop_start.strip())-1)
                        bed_stop = str(int(bed_start)+100)
                    elif strand == '-':
                        bed_start = str(int(stop_start.strip())-98)
                        bed_stop = str(int(stop_start)+2)
                    line = "\t".join([chromosome, bed_start, bed_stop, protein_id + "; " + 'unknown_exon', '.', strand, "\n"])
                    total += line

    return total


def get_utr_fasta(name, bed, source):

    print ('getting fasta...')

    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            if name in f:
                utr_fasta = bed.replace('/BED/', '/FASTA/').replace('.bed', '.fa')
                genome_fasta = os.path.abspath("../2_Project_FSM/7_Expression/Release-17/WGS/" + f)
                fasta_from_intervals(bed, utr_fasta, genome_fasta, names = True)


def fasta_from_intervals(bed_file, fasta_file, genome_fasta, force_strand = True, names = False):
    '''
    Takes a bed file and creates a fasta file with the corresponding sequences.
    If names == False, the fasta record names will be generated from the sequence coordinates.
    If names == True, the fasta name will correspond to whatever is in the 'name' field of the bed file
    '''

    #if the index file exists, check whether the expected features are present
    genome_fasta_index = genome_fasta + '.fai'
    if(os.path.exists(genome_fasta_index)):
        bed_chrs = sorted(list(set([entry[0] for entry in gen.read_many_fields(bed_file, "\t")])))
        index_chrs = sorted(list(set([entry[0] for entry in gen.read_many_fields(genome_fasta_index, "\t")])))
        if(not set(bed_chrs).issubset(set(index_chrs))):
            gen.remove_file(genome_fasta_index)

    bedtools_args = ["bedtools", "getfasta", "-s", "-fi", genome_fasta, "-bed", bed_file, "-fo", fasta_file]
    if not force_strand:
        del bedtools_args[2]
    if names:
        bedtools_args.append("-name")
    gen.run_process(bedtools_args)
    names, seqs = gen.read_fasta(fasta_file)
    seqs = [i.upper() for i in seqs]
    gen.write_to_fasta(names, seqs, fasta_file)


### RUN ###

if __name__ == '__main__':
    main()
