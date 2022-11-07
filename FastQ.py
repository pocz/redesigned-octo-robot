class Read:
    def __init__(self, label = '', seq = '', qlt = ''):
        self.seq = seq
        self.label = label
        self.qlt = qlt
    
        
    """Quality score treshold given as an integer between 1 and 94"""
    def drop_calls(self, treshold):
        def drop_chars(string, indexes):
            s = bytearray(string.encode())
            for i in sorted(indexes, reverse=True):
                del s[i]
            return s.decode()

        to_drop = []
        for i, c in enumerate(self.qlt):
            if ord(c) < (treshold+32): 
                to_drop.append(i)
        self.seq = drop_chars(self.seq, to_drop)
        self.qlt = drop_chars(self.qlt, to_drop)
    
    def validate(self, seq_type):
        A = self.seq.count("A")  
        G = self.seq.count("G")
        C = self.seq.count("C")
        T = self.seq.count("T")
        U = self.seq.count("U")
        N = self.seq.count("N")

        dna_valid = A+G+C+T+N == len(self.seq) and seq_type=="dna" 
        rna_valid = A+G+C+U+N == len(self.seq) and seq_type=="rna"
        if dna_valid == False and rna_valid == False:
            return False

        for c in self.qlt:
            if ord(c) < 33 or ord(c) > 128:
                return False

        return True 
       
class Sample:
    def __init__(self, filepath):
        self.reads = {}
        with open(filepath, 'r') as f:
            label = f.readline().strip()
            seq = f.readline().strip()
            f.readline()
            qlt = f.readline().strip()
            while label != "":
                self.reads[label] = fq_read(label,seq,qlt)
                label = f.readline().strip()
                seq = f.readline().strip()
                f.readline()
                qlt = f.readline().strip()

    def drop_calls(self, treshold):
        for read in self.reads:
            self.reads[read].drop_calls(treshold)

    def validate(self, seq_type):
        for read in self.reads:
            if self.reads[read].validate(seq_type) == False:
                return False
        return True
