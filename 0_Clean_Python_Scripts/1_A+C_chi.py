### IMPORTS ###
import os
import numpy as np
import csv
import time
import scipy.stats as st
import multiprocessing

### CHOOSE SOURCE FOLDER (EMBL)

source = '7_Expression/HEGs_FASTA'

### MAIN ###

def parallel():

    source = '7_Expression/HEGs_FASTA'
    filenames = get_files(source)
    workers = int(os.cpu_count()) - 3
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
            EF_hits, EF_array, EF_list = get_simulations(utr_list, total_utr)
            prob = np.mean(EF_array)
            prob_hits = [prob * len(utr_list), prob * len(utr_list), prob * len(utr_list), prob * len(utr_list), prob * len(utr_list), prob * len(utr_list), prob * len(utr_list)]

            #Set hits arrays
            OF_a = np.array(OF_hits)
            # EF_a = np.array(EF_hits)

            #Calculate excess
            excess = ((OF_a - prob_hits) / prob_hits) * 100
            excess_list = excess.tolist()

            #Chisq
            chisq = ((OF_a - prob_hits) ** 2) / prob_hits
            chisq_list = chisq.tolist()

            #p-vals
            p_list = []
            for i in chisq_list:
                p = st.distributions.chi2.sf (i, 1)
                p_list.append(p)

            #Output
            output = []
            output.extend(list(a) for a in zip(OF_hits, prob_hits, excess_list, chisq_list, p_list))
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

        if utr_seq[0:3] == 'TAA' or utr_seq[0:3] == 'TGA' or utr_seq[0:3] == 'TAG':
            codon_seq = [utr_seq[i:i+3] for i in range(0, len(utr_seq), 3)]
            utr_list.append(codon_seq)

            #Now create total UTR list (only up to 100 nts) which will be useful later
            total_utr += utr_seq[1:]

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
        a = nt_freq['A']
        c = nt_freq['C']
        g = nt_freq['G']
        t = nt_freq['T']

        #Calculate next base
        # first_base = np.random.choice(['A', 'C', 'G', 'T'], p=[nt_freq['A'], nt_freq['C'], nt_freq['G'], nt_freq['T']])
        first_base = np.random.choice(['A', 'C', 'G', 'T'], p=[a, c, g, (1-a-c-g)])
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


def create_csv(headers, csv_total, accession):

    ''' Output results to CSV file '''

    filename = accession.replace('.txt', '') + ".csv"
    subdir = '7_Expression/Chi1_HEGs'
    filepath = os.path.join(subdir, filename)

    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### RUN ###

if __name__ == '__main__':
    parallel()
