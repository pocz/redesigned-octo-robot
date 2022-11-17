import AhoCorasick
import BoyerMoore

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

for j, y in enumerate(function):
    for i, x in enumerate(answer):
        if answer[i] == function[j](seq[i],kmer[i],seq_type[i]):
            print(f'{y} passed check #{i}. {x}')
        else:
            print(f'{y} failed check #{i}. {x}')

kmers = ['AGC', 'AGCC', 'GAAC', 'AAC']
seq_type = "DNA"
seq = "NNNNNAGGAGAAGCCAGAACGAGCNNN"

AC_instance = AhoCorasick.Automaton(kmers,seq_type)
y = 'The Aho-Corasick algorithm'
x = AC_instance.start_indexes(seq)
answer = {'AGC': [11, 21], 'AGCC': [11], 'GAAC': [16], 'AAC': [17]}
if answer == x:
    print(f'{y} passed the check.\n{x}')
else:
    print(f'{y} failed the check.\n {x}')

