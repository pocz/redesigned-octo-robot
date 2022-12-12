Nukleotid szekvencia ábécéjét használó, különböző mintaillesztési algoritmusok alapvető implementációjait tartalmazó program. FastQ bemeneti fájlokat véve grafikusan megjeleníti a pontos egyezések szakaszait egy indexeket jelző számtengelyen.

## Dependencies
<pre>
</pre>

## Running the visualizer (using test data files provided):
<pre>
./main.py --kmers test_data/6mer-* --seq test_data/seq-short.fq --mode ProcessPool
</pre>

## Parameters

--kmers     FASTQ file or files containing the k-mer to be searched;
--seq       FASTQ file containing the sequence to search in;
--mode      function used for the search: ProcessPool=Uses the ProcessPoolExecutor from concurrent.futures,
            Threading=Uses multiple threads via the Threading package, AhoCorasick=Uses an implementation of the Aho-Corasick prefix-tree algorithm.
 
