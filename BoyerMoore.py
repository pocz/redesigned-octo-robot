''' 
Returns a list containing a value for each position i in P and character x,
the position of the closest occurence of x in P to the left of i.

When a mismatch occurs at position i of P 
and the mismatched character in T is x, 
then shift P to the right 
so that the closest x to the left of position i in P 
is below the mismatched x in T.
'''
def bad_character_rule(kmer):
    bad_characters = {}
    for char in "ACTGN":
        bad_characters[char] = -1
    for i in range(len(kmer)):
        c = kmer[i]
        bad_characters[c] = i
    return bad_characters

def good_character_rule(kmer):
    f = [0] * (len(kmer)+1)
    s = [0] * (len(kmer)+1)
    i = len(kmer)
    j = len(kmer)+1
    f[i] = j

    while i>0:
        while j <= len(kmer)  and kmer[i-1] != kmer[j-1]:
            if s[j] == 0:
                s[j] = j-i
            j = f[j]
        i -= 1
        j -= 1
        f[i] = j
    j = f[0]
    for i in range(len(kmer)):
        if s[i] == 0:
            s[i] = j
        if i == j:
            j = f[j]
    
    return s

'''
Returns a list containing the leftmost indexes of each occurence of the kmer 
in the sequence.
'''
def positions(sequence, kmer):
    positions = []
    bcr = bad_character_rule(kmer)
    gcr = good_character_rule(kmer)

    shift = 0
    while shift <= len(sequence)-len(kmer):

        i = len(kmer) - 1
        while i >=0 and kmer[i] == sequence[i+shift]:
            i = i - 1

        # if whole pattern matches
        if i == -1:
            positions.append(shift)
            shift = shift + 1
        else:
            j = sequence[shift + i]
            shift = shift + max(gcr[i+1], i - bcr[j])

    return positions

'''
from FastQ import Sample

print(bad_character_rule("CAGAG"))
print(good_character_rule("CAGAG"))
for label in sample.reads:
    print(positions(sample.reads[label].seq, "CAGAG"))

sample = fq_sample("./ExampleSample5.fq")
for label in sample.reads:
    print(positions(sample.reads[label].seq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"))
    break

sample = Sample("./ExampleSample2.fq")

for label in sample.reads:
    print(sample.reads[label].seq)
    print(naive_string_search(sample.reads[label].seq, "CAG"))

'''

