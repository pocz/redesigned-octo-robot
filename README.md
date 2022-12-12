Nukleotid szekvencia ábécéjét használó, különböző mintaillesztési algoritmusok alapvető implementációjait tartalmazó program. FastQ bemeneti fájlokat véve grafikusan megjeleníti a pontos egyezések szakaszait egy indexeket jelző számtengelyen.


A BoyerMoore  osztályban lett lehetővé téve az algoritmus több szálon való futása (valamint GIL megkerülésére - több processban). Így távolról becsülve összehasonlítható a prefix-fás és a párhuzamosításos megoldások sok rövid string (k-mer) illesztésére. 
<pre>
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
</pre>    
    
## Requirements
Uses the following pip packages:
<pre>
PyQt5 pyqtgraph
</pre>

## Running the visualizer (using test data files provided):
<pre>
./main.py --kmers test_data/6mer-* --seq test_data/seq-short.fq --mode ProcessPool
</pre>

## Parameters

```
--kmers     FASTQ file or files containing the k-mer to be searched;
--seq       FASTQ file containing the sequence to search in;
--mode      method used for the search: ProcessPool Uses the ProcessPoolExecutor from concurrent.futures,
            Threading=Uses multiple threads via the Threading package, 
            AhoCorasick=Uses an implementation of the Aho-Corasick prefix-tree algorithm.
``` 
