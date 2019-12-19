### IMPORTS ###
import os
import numpy as np
import csv
import time
import scipy.stats as st
import multiprocessing

### MAIN FUNCTIONS ### - CHANGE SOURCE FOLDER, OUTPUT FOLDERS, OUTPUT FILENAMES FOR THE GROUP OF GENOMES OF INTEREST

def parallel():

    source = '1_Fungi/Filtered_FASTA'
    filenames = get_files(source)
    workers = int(os.cpu_count()) - 1
    processes = run_in_parallel(filenames, ['foo', source], main, workers=workers)

    csv_total = []

    for process in processes:
        output = process.get()
        csv_total.extend(output)

    write_outputs(csv_total)

def main(filenames, source):

    csv_total = []

    for i,f in enumerate(filenames):

        path = os.path.join(source, f)
        raw = open(path).read()

        print ('Doing genome {0}/{1}'.format(i+1, len(filenames)))

        #Split into chunks
        genes = raw.strip().split('>')
        genes_c = list(filter(None, genes))
        glist = get_list(genes_c)

        #Get +4T frequency
        real_f, real_hits = get_t(glist)

        #Get sims
        utr_list, total_utr = get_utr_stuff(genes_c)
        sim_f, sim_n = get_simulations(utr_list, total_utr)

        #Z-scores
        o = real_f * sim_n
        mean = sim_f * sim_n
        std = (sim_n * sim_f * (1-sim_f)) ** 0.5
        z = (o - mean) / std

        csv_line = [f, real_f, sim_f, z]
        csv_total.append(csv_line)

    return csv_total

### FUNCTIONS ###

def get_simulations(utr_list, total_utr):

    ''' Simulate n UTR sequences for each genome based upon dinucleotide content of total UTR regions '''

    #Set no of sims
    no_sims = 1000

    #Probability of stop codon
    stop_dict = {}
    length = len(utr_list)

    for i in range(length):
        codon = utr_list[i][0]
        if codon not in stop_dict:
            stop_dict[codon] = 1
        else:
            stop_dict[codon] += 1

    stop_freq = {k: v / length for k, v in stop_dict.items()}

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
    while count < no_sims:

        sim_n = []
        codon_list = []

        sim = []

        np.random.seed()
        taa = stop_freq['TAA']
        tag = stop_freq['TAG']
        tga = stop_freq['TGA']

        #Calculate next base
        stop = np.random.choice(['TAA', 'TGA', 'TAG'], p=[taa, tga, (1-taa-tga)])
        sim.append(stop)

        #For one gene simulation
        while (len(sim) < 2):

            #Generate next base
            prev_base = sim[0][2]
            next_base = np.random.choice(trans[prev_base], replace=True)
            sim.append(next_base)

        sim = "".join(sim)
        total_list.append(sim)
        count += 1

    t_total = 0
    for seq in total_list:
        if seq[3] == 'T':
            t_total += 1

    x = t_total / no_sims

    return x, no_sims


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

        if utr_seq[0:3] == 'TAA' or utr_seq[0:3] == 'TGA' or utr_seq[0:3] == 'TAG':
            if sequence[3:6] != 'TGA' or sequence[3:6] != 'TAA' or sequence[3:6] != 'TAG':
                codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
                utr_list.append(codon_seq)

                #Now create total UTR list (only up to 100 nts) which will be useful later
                total_utr += utr_seq[1:]

    return utr_list, total_utr


def get_t(genes):

    t_count = 0

    for gene in genes:
        sequence = "".join(gene)
        fourth = sequence[3]
        if fourth == 'T':
            t_count += 1

    freq = t_count / len(genes)

    return freq, t_count


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
            if sequence[3:6] != 'TGA' or sequence[3:6] != 'TAA' or sequence[3:6] != 'TAG':
                utr_sequence = sequence.upper()
                utr_seq = [utr_sequence[i:i+3] for i in range(0, len(utr_sequence), 3)]
                utr_list.append(utr_seq)

    return utr_list


def write_outputs(csv_total):

    headers = ["File", "Observed", "Expected", 'Z-score']
    create_csv(headers, csv_total)


def create_csv(headers, csv_total):

    filename = '+4T_sims_fungi.csv'
    subdir = ""
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(headers)
        for j in csv_total:
            writer.writerow(j)


### RUN ###

if __name__ == '__main__':
    parallel()
