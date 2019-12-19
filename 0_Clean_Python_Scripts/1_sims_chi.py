### IMPORTS ###
import os
import numpy as np
import csv
import time
import scipy.stats as st
import multiprocessing

### CHOOSE SOURCE FOLDER ###

source = '6_Ne/Extra_genomes_2/Filtered_FASTA'

### MAIN ###

def parallel():

    source = '6_Ne/Extra_genomes_2/Filtered_FASTA'
    filenames = get_files(source)
    workers = int(os.cpu_count()) - 1
    processes = run_in_parallel(filenames, ['foo', source], main, workers=workers)

    for process in processes:
        process.get()

def main(filenames, source):

    for i,f in enumerate(filenames):

        path = os.path.join(source, f)
        raw = open(path).read()

        print ('Doing genome {0}/{1}'.format(i+1, len(filenames)))

        #Split into chunks
        genes = raw.strip().split('>')
        genes_c = list(filter(None, genes))

        if len(genes_c) > 500:
            glist = get_list(genes_c)

            #Get stop frequencies (O)
            OF_hits, OF_array, OF_list = get_stop_freq(glist)

            #Obtain UTR list and generate total UTR string, which is useful for later functions
            utr_list, total_utr = get_utr_stuff(genes_c)

            #Get simulated stop frequencies (E)
            EF_hits, EF_array, EF_list = get_null(utr_list, total_utr)

            #Set hits arrays
            OF_a = np.array(OF_hits)
            EF_a = np.array(EF_hits)

            #Calculate excess
            excess = ((OF_a- EF_a) / EF_a) * 100
            excess_list = excess.tolist()

            #Chisq
            chisq = ((OF_a - EF_a) ** 2) / EF_a
            chisq_list = chisq.tolist()

            #p-vals
            p_list = []
            for i in chisq_list:
                p = st.distributions.chi2.sf (i, 1)
                p_list.append(p)

            #Output
            output = []
            output.extend(list(a) for a in zip(OF_hits, EF_hits, excess_list, chisq_list, p_list))
            headers = ["Observed", "Expected", "Excess%", "Chisq", "p-value"]
            create_csv(headers, output, f)

### FUNCTIONS ###

def get_files(source):

    ''' Obtain all file names from the target folder '''

    files = []

    for root, dirs, filenames in os.walk(source):
        for f in filenames:
            files.append(f)

    return files


def run_in_parallel(input_list, args, func, kwargs_dict = None, workers = None, onebyone = False):

    '''
    Take an input list, divide into chunks and then apply a function to each of the chunks in parallel.
    input_list: a list of the stuff you want to parallelize over (for example, a list of gene names)
    args: a list of arguments to the function. Put in "foo" in place of the argument you are parallelizing over.
    func: the function
    kwargs_dict: a dictionary of any keyword arguments the function might take
    workers: number of parallel processes to launch
    onebyone: if True, allocate one element from input_list to each process
    '''

    if not workers:
        #divide by two to get the number of physical cores
        #subtract one to leave one core free
        workers = int(os.cpu_count()/2 - 1)
    elif workers == "all":
        workers = os.cpu_count()
    #in the list of arguments, I put in "foo" for the argument that corresponds to whatever is in the input_list because I couldn't be bothered to do something less stupid
    arg_to_parallelize = args.index("foo")
    if not onebyone:
        #divide input_list into as many chunks as you're going to have processes
        chunk_list = [input_list[i::workers] for i in range(workers)]
    else:
        #each element in the input list will constitute a chunk of its own.
        chunk_list = input_list
    pool = multiprocessing.Pool(workers)
    results = []
    #go over the chunks you made and laucnh a process for each
    for elem in chunk_list:
        current_args = args.copy()
        current_args[arg_to_parallelize] = elem
        if kwargs_dict:
            process = pool.apply_async(func, tuple(current_args), kwargs_dict)
        else:
            process = pool.apply_async(func, tuple(current_args))
        results.append(process)
    pool.close()
    pool.join()
    return(results)


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

        #First generate UTR list, subnested for each codon
        a = i.split("\n")
        sequence = a[1].upper()
        utr_seq = sequence
        codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
        if codon_seq[0] == 'TGA' or codon_seq[0] == 'TAA' or codon_seq[0] == 'TAG':
            utr_list.append(codon_seq)
            #Now create total UTR list which will be useful later
            total_utr += utr_seq[1:]

    return utr_list, total_utr


def get_null(utr_list, total_utr):

    #No. of genes
    n_genes = len(utr_list)

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
    EF_hits_list = EF_hits.tolist()

    return EF_hits_list, EF_f_array, EF_frequencies


def create_csv(headers, csv_total, accession):

    ''' Output results to CSV file '''

    filename = accession.replace('.txt', '') + ".csv"
    subdir = '6_Ne/Extra_genomes_2/Chi2'
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)


### RUN ###

if __name__ == '__main__':
    parallel()
