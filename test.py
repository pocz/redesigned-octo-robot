import AhoCorasick
import BoyerMoore
import FastBoyerMoore

seq = []
kmer = []
seq_type = []
answer = []
function = []

seq.append("AGCNNNNNNNNNNNAGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGC") 
kmer.append("AGC")
seq_type.append("DNA")
answer.append([0,14,107])
function.append(BoyerMoore.start_indexes)
#function.append(FastBoyerMoore.start_indexes)


for j, y in enumerate(function):
    for i, x in enumerate(answer):
        if answer[i] == function[j](seq[i],kmer[i],seq_type[i]):
            print(f'{y} passed check #{i}. {x}')
        else:
            print(f'{y} failed check #{i}. {x}')


for i, x in enumerate(answer):
    y = "The Aho-Corasick algorithm"
    AC_instance = AhoCorasick.AhoCorasick([kmer[i]],seq_type[i])
    if x == AC_instance.start_indexes(seq[i]):
        print(f'{y} passed check #{i}. {x}')
    else:
        print(f'{y} failed check #{i}. {x}')


