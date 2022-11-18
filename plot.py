import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import matplotlib.ticker as ticker
mplstyle.use('fast')

import BoyerMoore
import profiling

def plot_subsequences(kmers, sequence):
    for i, kmer in enumerate(kmers):
        start_indexes, end_indexes = BoyerMoore.edge_indexes(sequence, kmer)

        for l, y in enumerate(start_indexes):
            start_indexes[l] -= 0.5
        for l, y in enumerate(end_indexes):
            end_indexes[l] -= 0.5

        for j, start_index in enumerate(start_indexes):
            plt.hlines(i, start_index, end_indexes[j], lw=3)
            plt.vlines(start_index, i+0.05, i-0.05, lw=3)
            plt.vlines(end_indexes[j], i+0.05, i-0.05, lw=3)
    ax = plt.gca()
    plt.grid(color='black', linestyle='-', linewidth=0.5)
    plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(1.00))
    plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(1.00))

    current_values = plt.gca().get_yticks()
    print(kmers)
    plt.gca().set_yticklabels([f'{kmers[int(y)]}' for y in current_values])
    
    current_values = plt.gca().get_xticks()
    plt.gca().set_xticklabels([f'{sequence[int(x)]}\n{int(x)}' for x in current_values])

    #plt.gca().yaxis.set_major_formatter('')
    plt.show()

def plot_subsequences_barh(kmers, sequence):
    # Horizontal bar plot with gaps
    fig, ax = plt.subplots()
    for i, kmer in enumerate(kmers):
        #print(kmer)
        x_tuples = []
        start_indexes, end_indexes = BoyerMoore.edge_indexes(sequence, kmer)
        for j, start_index in enumerate(start_indexes):
            x_tuples.append( (start_index, end_indexes[j]-start_index) )
            #print(x_tuples)
        ax.broken_barh(x_tuples, (i-.25, .5))

    #plt.gca().set(ylim=(-0.5, len(kmers)-0.5))
    current_values = plt.gca().get_yticks()
    def y_tick_values(y):
        if y >= 0 and y < len(kmers):
            return kmers[int(y)] 
        else:
            return y
    plt.gca().set_yticklabels([f'{y_tick_values(y)}' for y in current_values])

    current_values = plt.gca().get_xticks()
    def x_tick_values(x):
        if x >= 0 and x < len(sequence):
            return sequence[int(x)] 
        else:
            return x
     
    plt.gca().set_xticklabels([f'{x_tick_values(x)}\n{int(x)}' for x in current_values])

    #plt.savefig('big_bars.png')
    plt.show()

def plot_subsequences_autosize(kmers, sequence):
    # Target resolution is 500
    divisor = len(sequence) / 500
    print(divisor)
    for i, kmer in enumerate(kmers):
        start_indexes, end_indexes = BoyerMoore.edge_indexes(sequence, kmer)
        for i, seq_pos in enumerate(start_indexes): 
            start_indexes[i] = int(seq_pos / divisor)
        for i, seq_pos in enumerate(end_indexes): 
            end_indexes[i] = int(seq_pos / divisor)

        for l, y in enumerate(start_indexes):
            start_indexes[l] -= 0.5 
        for l, y in enumerate(end_indexes):
            end_indexes[l] -= 0.5

        for j, start_index in enumerate(start_indexes):
            plt.hlines(i, start_index, end_indexes[j], lw=3)
            plt.vlines(start_index, i+0.25, i-.25, lw=3)
            plt.vlines(end_indexes[j], i+0.25, i-.25, lw=3)
    ax = plt.gca()

    '''
    plt.grid(color='black', linestyle='-', linewidth=0.5)
    plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(1.00))
    plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(1.00))
    '''

    current_values = plt.gca().get_yticks()
    plt.gca().set(ylim=(-0.5, len(kmers)-0.5))
    plt.gca().set_yticklabels([f'{kmers[int(y)%len(kmers)]}' for y in current_values])
    
    current_values = plt.gca().get_xticks()
    plt.gca().set_xticklabels([f'{int(x*divisor)}' for x in current_values])

    #plt.gca().yaxis.set_major_formatter('')
    plt.show()
    print(kmers)
    plt.savefig('50m_test.png')

def plot_subsequences_bargraph(kmers, sequence):
    bars = 50
    counts = []
    divisor = len(sequence) / bars
    print(divisor)
    for i, kmer in enumerate(kmers):
        start_indexes, end_indexes = BoyerMoore.edge_indexes(sequence, kmer)
        for l, y in enumerate(start_indexes):
            start_indexes[l] -= 0.5 
        for l, y in enumerate(end_indexes):
            end_indexes[l] -= 0.5
    

    ax = plt.gca()

    current_values = plt.gca().get_yticks()
    plt.gca().set(ylim=(-0.5, len(kmers)-0.5))
    plt.gca().set_yticklabels([f'{kmers[int(y)%len(kmers)]}' for y in current_values])
    
    current_values = plt.gca().get_xticks()
    plt.gca().set_xticklabels([f'{int(x*divisor)}' for x in current_values])

    #plt.gca().yaxis.set_major_formatter('')
    plt.show()
    print(kmers)
    #plt.savefig('50m_test.png')
'''
def density(kmer, sequence, resolution):
    start_indexes = BoyerMoore.start_indexes(sequence, kmer)
    sequence_length = len(sequence)
    step = sequence_length/resolution
    counts = []
'''

kmer_paths = [  "data/plot/a.fq",
                "data/plot/b.fq",
                "data/plot/c.fq",
                "data/plot/d.fq",
                "data/plot/e.fq",
                "data/plot/f.fq"]

seq_paths = ["data/timings/18000.py",
"data/timings/50m.fq"]

kmers = []
for kmer_path in kmer_paths:
    kmers.append(profiling.last_read(kmer_path))

seqs = []
for seq_path in seq_paths:
    seqs.append(profiling.last_read(seq_path))

#plot_subsequences_autosize(kmers, seqs[0])
plot_subsequences_barh(kmers, seqs[1])


