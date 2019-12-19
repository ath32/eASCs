### IMPORTS ###

import os
import csv

### CHOOSE SOURCE FOLDER ###

# source = '1_Vertebrates/FASTA'
# source = '1_Metazoa/FASTA'
# source = '1_Fungi/FASTA'
source = '5_Protists/Protist_UTR'

### MAIN ###

def main(source):

    CSV_total = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:

            path = os.path.join(source, f)
            raw = open(path).read()
            name = f

            #Access sequences and split into primary stop groups
            a = raw.split('>')
            chunks = list(filter(None, a))
            taa_genes, tga_genes, tag_genes = get_lists(chunks)

            if len(taa_genes) > 0 and len(tga_genes) > 0 and len(tag_genes) > 0:

                #Get fourth site frequencies
                taa_a, taa_t, taa_g, taa_c = get_fourth_sites(taa_genes)
                tga_a, tga_t, tga_g, tga_c = get_fourth_sites(tga_genes)
                tag_a, tag_t, tag_g, tag_c = get_fourth_sites(tag_genes)

                #Create CSV line
                CSV_line = [name, taa_a, taa_t, taa_g, taa_c, tga_a, tga_t, tga_g, tga_c, tag_a, tag_t, tag_g, tag_c]
                CSV_total.append(CSV_line)

    #Write output files
    headers = ['Genome', 'TAA_FourthA', 'TAA_FourthT', 'TAA_FourthG', 'TAA_FourthC',
        'TGA_FourthA', 'TGA_FourthT', 'TGA_FourthG', 'TGA_FourthC',
        'TAG_FourthA', 'TAG_FourthT', 'TAG_FourthG', 'TAG_FourthC']

    get_CSV(CSV_total, headers)

### FUNCTIONS ###

def get_lists(chunks):

    taa = []
    tga = []
    tag = []

    for i in chunks:

        #Access utr
        b = i.split('\n')
        if len(b) > 1:
            utr = b[1]

            #Split into groups
            if utr[0:3] == 'TAA':
                taa.append(utr)
            elif utr[0:3] == 'TGA':
                tga.append(utr)
            elif utr[0:3] == 'TAG':
                tag.append(utr)

    return taa, tga, tag


def get_fourth_sites(utrs):

    total_genes = len(utrs)
    fourth_a = 0
    fourth_t = 0
    fourth_g = 0
    fourth_c = 0

    for i in utrs:

        #Get fourth site
        fourth = i[3]

        if fourth == 'A':
            fourth_a += 1
        if fourth == 'T':
            fourth_t += 1
        if fourth == 'G':
            fourth_g += 1
        if fourth == 'C':
            fourth_c += 1

    freq_a = fourth_a / total_genes
    freq_t = fourth_t / total_genes
    freq_g = fourth_g / total_genes
    freq_c = fourth_c / total_genes

    return freq_a, freq_t, freq_g, freq_c


def get_CSV(total, headers):

    filename = "fourth_site_uni.csv"
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
