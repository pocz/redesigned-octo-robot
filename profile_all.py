
import profiling

'''
kmer_paths = ["data/timings/32.py",
"data/timings/281.py",
"data/timings/2000.py",
"data/timings/18000.py"]
'''
kmer_paths = ["data/timings/8.fq"] * 20

seq_paths = ["data/timings/18000.py",
"data/timings/18000.py"]

kmers = []
for kmer_path in kmer_paths:
    kmers.append(profiling.last_read(kmer_path))

seqs = []
for seq_path in seq_paths:
    seqs.append(profiling.last_read(seq_path))

for seq in seqs:
    profiling.print_BM_profile(kmers,"DNA",seq)
print(f'Boyer-Moore all kmers on 1 sequence')

for i, seq in enumerate(seqs):
    profiling.print_AC_profile(kmers,"DNA",seq)
    print(f'Searching on {len(seq)} long sequence')
print(f'Aho-Corasick all kmers on 1 sequence')

for seq in seqs:
    profiling.profile_processpool(seq,kmers)
print(f'processpool B-M')
