from concurrent.futures import ProcessPoolExecutor
import pandas as pd

import BoyerMoore

kmers = ["AT", "UG", "AU", "CUU", "GA", "G", "GG", "TG", "C", "AA"]
sequences = ["GATATATGCATATACTT", "ATAGACCATATGTCAGTGACTGTGTAA", "CTAAGGGATTCCGGTAATTAGACAG"] 

def matrix_count_parallel(kmers, sequences, seq_type="DNA", max_workers=12):
    seq_queue = []
    kmer_queue = []
    for column, seq in enumerate(sequences):
        for row, kmer in enumerate(kmers):
            seq_queue.append(seq)
            kmer_queue.append(kmer)

    pool = ProcessPoolExecutor(max_workers)
    results = list(pool.map(BoyerMoore.count, seq_queue, kmer_queue)) 

    rows = len(kmers)
    columns = len(sequences)
    results = np.array(results).reshape((columns,rows))
    
    df = pd.DataFrame(results)
    return df

