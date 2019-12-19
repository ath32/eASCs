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
### DO A FIND AND REPLACE TO CHANGE FOLDER OF INTEREST ###

gff_source = 'x_Repeat/1_Vertebrates/GFF'
wgs_source = 'x_Repeat/1_Vertebrates/WGS'
fasta_source = 'x_Repeat/1_Vertebrates/FASTA'

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

            print (f)

            valid_seq = 0
            total = 0
            new_fasta = ''

            file_split = raw.split('>')
            for i in file_split:
                if i != '':
                    if 'ID=gene:' in i:
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
            print (percent_valid)
            new_source = "x_Repeat/1_Vertebrates/Filtered_FASTA/"
            new_fasta_path = os.path.join(new_source, f)
            print (new_fasta_path)

            if percent_valid > 0.75:
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
    split = raw.split('\n')
    filtered_ids = IGR_filter(split)

    if filtered_ids == '':
        print (name)

    #For each gene in the filtered list, get exon annotations
    special = raw.split('###')
    bed = get_utr(special, filtered_ids)

    #Output bed
    bed_path = os.path.abspath("../2_Project_FSM/x_Repeat/1_Vertebrates/BED/" + name + ".bed")

    f = open(bed_path, 'a')
    f.write(bed)
    f.close()

    return bed_path


def IGR_filter(list):

    print ('filtering genes for sufficient IGR...')

    ids = []

    #For each set of coordinates...
    for i in list:
        split_list = i.split(';')

        #Split to access gene info
        if 'gene\t' in split_list[0]:
            chunk1 = split_list[0].strip()
            chunk2 = chunk1.split('\t')

            #Gene info
            chromosome = chunk2[0]
            start = chunk2[3]
            stop = chunk2[4]
            strand = chunk2[6]
            gene_id = chunk2[8]

            #Split genes by strand - adjust coordinates to base 1 / reverse complement / to include start and stop codons
            if strand == '+':
                ids.append([gene_id, int(start)-1, int(stop)+3, strand])
            elif strand == '-':
                ids.append([gene_id, int(start)-4, int(stop), strand])

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
        if 'gene' in i and '\tCDS\t' in i:

            #First, determine strand
            split = i.split('\n')

            for bit in split:

                if '\tgene\t' in bit:
                    line = bit
                    line_split = line.split('\t')
                    strand = line_split[6]
                    chromosome = line_split[0]

                    #Next, get Gene ID
                    id_split = line_split[8].split(';')
                    gene_id = id_split[0]

                    #Next, get coordinates for + strand (last cds)
                    if strand == '+':
                        chunks = i.split('\tCDS\t')
                        last_CDS = "".join(chunks[len(chunks)-1:len(chunks)])
                        sections = last_CDS.split('\t')
                        start = str(int(sections[1].strip())-3)
                        stop = str(int(sections[1].strip())+97)

                    if strand == '-':
                        chunks = i.split('\tCDS\t')
                        first_CDS = "".join(chunks[1])
                        sections = first_CDS.split('\t')
                        start = str(int(sections[0].strip())-98)
                        stop = str(int(sections[0].strip())+2)

                    if gene_id in filtered_ids:
                        line = "\t".join([chromosome, start, stop, gene_id + "; " + 'unknown_exon', '.', strand, "\n"])
                        total += line

    return total


def get_utr_fasta(name, bed, source):

    print ('getting fasta...')

    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            if name in f:
                utr_fasta = bed.replace('/BED/', '/FASTA/').replace('.bed', '.fa')
                genome_fasta = os.path.abspath("../2_Project_FSM/x_Repeat/1_Vertebrates/WGS/" + f)
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
