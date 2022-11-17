class Read:
    def __init__(self, label = '', seq = ''):
        self.seq = seq
        self.label = label
    
        
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
            while label != "":
                self.reads[label] = Read(label,seq,qlt)
                label = f.readline().strip()
                seq = f.readline().strip()

    def validate(self, seq_type):
        for read in self.reads:
            if self.reads[read].validate(seq_type) == False:
                return False
        return True
