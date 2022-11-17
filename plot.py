import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import matplotlib.ticker as ticker
mplstyle.use('fast')

import BoyerMoore

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

plot_subsequences(kmers, sequences[0])
