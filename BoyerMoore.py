from concurrent.futures import ProcessPoolExecutor

import threading

''' 
Returns a list containing a value for each position i in P and character x,
the position of the closest occurence of x in P to the left of i.

When a mismatch occurs at position i of P 
and the mismatched character in T is x, 
then shift P to the right 
so that the closest x to the left of position i in P 
is below the mismatched x in T.
'''
def bad_character_rule(kmer,seq_type):
    bcr_shifts = {}
    if seq_type.upper() == "DNA":
        alphabet = "AGTCN"
    elif seq_type.upper() == "RNA":
        alphabet = "AGUCN"

    for char in alphabet:
        bcr_shifts[char] = -1
    for i in range(len(kmer)):
        char = kmer[i]
        bcr_shifts[char] = i
    return bcr_shifts

'''
Example: in the case of "CAGAG" k-mer, wil return [5,5,5,2,5,1], because if
there's a mismatch at he second A at index 3, we can shift by two, so we
end up at the second instance of "AG".
'''
def good_character_rule(kmer):
    aux_table = [0] * (len(kmer)+1)
    gcr_shifts = [0] * (len(kmer)+1)
    i = len(kmer)
    j = len(kmer)+1
    aux_table[i] = j

    while i>0:
        while j <= len(kmer)  and kmer[i-1] != kmer[j-1]:
            if gcr_shifts[j] == 0:
                gcr_shifts[j] = j-i
            j = aux_table[j]
        i -= 1
        j -= 1
        aux_table[i] = j
    j = aux_table[0]
    for i in range(len(kmer)):
        if gcr_shifts[i] == 0:
            gcr_shifts[i] = j
        if i == j:
            j = aux_table[j]
    
    return gcr_shifts

'''
Returns a list containing the leftmost indexes of each occurence of the kmer 
in the sequence.
'''
def start_indexes(sequence, kmer, seq_type="DNA"):
    start_indexes = []
    bcr = bad_character_rule(kmer,seq_type)
    gcr = good_character_rule(kmer)

    shift = 0
    sequence_length = len(sequence)
    kmer_length = len(kmer)
    while shift <= sequence_length - kmer_length:

        i = kmer_length - 1
        while i >=0 and kmer[i] == sequence[i+shift]:
            i = i - 1

        if i == -1:
            start_indexes.append(shift)
            shift = shift + 1
        else:
            j = sequence[shift + i]
            gcr_shift = gcr[i+1]
            bcr_shift = i - bcr[j]
            if gcr_shift >= bcr_shift:
                shift = shift + gcr_shift
            else:
                shift = shift + bcr_shift

    return start_indexes

def edge_indexes(sequence, kmer, seq_type="DNA"):
    start_indexes = []
    end_indexes = []
    bcr = bad_character_rule(kmer,seq_type)
    gcr = good_character_rule(kmer)

    shift = 0
    sequence_length = len(sequence)
    kmer_length = len(kmer)
    while shift <= sequence_length - kmer_length:

        i = kmer_length - 1
        while i >=0 and kmer[i] == sequence[i+shift]:
            i = i - 1

        if i == -1:
            start_indexes.append(shift)
            end_indexes.append(shift+(kmer_length))
            shift = shift + 1
        else:
            j = sequence[shift + i]
            shift = shift + max(gcr[i+1], i - bcr[j])

    return start_indexes, end_indexes

def count(sequence, kmer, seq_type="DNA"):
    count = 0
    bcr = bad_character_rule(kmer,seq_type)
    gcr = good_character_rule(kmer)

    shift = 0
    sequence_length = len(sequence)
    kmer_length = len(kmer)
    while shift <= sequence_length - kmer_length:

        i = kmer_length - 1
        while i >=0 and kmer[i] == sequence[i+shift]:
            i = i - 1

        if i == -1:
            count = count + 1
            shift = shift + 1
        else:
            j = sequence[shift + i]
            shift = shift + max(gcr[i+1], i - bcr[j])
    return count

def start_indexes_parallel(seq,kmers,seq_type="DNA",workers=12):
    pool = ProcessPoolExecutor(max_workers=workers)
    results = pool.map(start_indexes,[seq]*len(kmers),kmers)
    return list(results)

def start_indexes_threading(seq,kmers,seq_type="DNA",threads=12):
    results = list()
    def thread_func(_seq,_kmer):
        results.append(start_indexes(_seq,_kmer))
    threads = list()
    for i, kmer in enumerate(kmers):
        arguments = (seq, kmer)
        x = threading.Thread(target=thread_func,args=arguments)
        x.start()
        threads.append(x)
    return list(results)
