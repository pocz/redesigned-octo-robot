from cProfile import Profile
from pstats import Stats

from concurrent.futures import ProcessPoolExecutor

import AhoCorasick
import BoyerMoore
import FastQ

def last_read(filepath):
    seq = ''
    for read in FastQ.Sample(filepath).reads:
        seq = FastQ.Sample(filepath).reads[read].seq
    return seq

def print_AC_profile(kmers,seq_type,seq):
    AC_instance = AhoCorasick.Automaton(kmers,seq_type)

    test = lambda: AC_instance.start_indexes(seq)
    profiler = Profile()
    profiler.runcall(test)

    stats = Stats(profiler)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    stats.print_stats()

def print_AC_profile_1kmer(kmer,seq_type,seq):
    kmers = []
    kmers.append(kmer)
    AC_instance = AhoCorasick.Automaton(kmers,seq_type)

    test = lambda: AC_instance.start_indexes(seq)
    profiler = Profile()
    profiler.runcall(test)

    stats = Stats(profiler)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    stats.print_stats()

def BM_on_all_kmers(kmers,seq_type,seq):
    for kmer in kmers:
        BoyerMoore.start_indexes(seq, kmer, seq_type)

def print_BM_profile(kmers,seq_type,seq):
    test = lambda: BM_on_all_kmers(kmers,"DNA",seq)
    profiler = Profile()
    profiler.runcall(test)

    stats = Stats(profiler)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    stats.print_stats()

def print_BM_profile_individual_kmer(kmer,seq_type,seq):
    test = lambda: BoyerMoore.start_indexes(seq,kmer,seq_type)
    profiler = Profile()
    profiler.runcall(test)

    stats = Stats(profiler)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    stats.print_stats()

def print_BM_profiles_1kmer(seqs,kmers):
    for seq in seqs:
        for kmer in kmers:
            print(f'Searching {len(kmer)} long kmer in {len(seq)} sequence...')
            print_BM_profile_individual_kmer(kmer,"DNA",seq) 

def parallel(seq, kmers):
    seqs = [seq] * len(kmers)
    pool = ProcessPoolExecutor(max_workers=12)
    results = list(pool.map(BoyerMoore.start_indexes, seqs, kmers))

def profile_processpool(seq,kmers):
    test = lambda: parallel(seq, kmers)
    profiler = Profile()
    profiler.runcall(test)

    stats = Stats(profiler)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    stats.print_stats()


