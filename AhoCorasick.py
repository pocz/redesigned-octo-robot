class AhoCorasick:
    def __init__(self, kmers, seq_type):
        if seq_type == "DNA":
            self.alphabet = "AGTCN"
        elif seq_type == "RNA":
            self.alphabet = "AGUCN"

        self.max_states = sum([len(kmer) for kmer in kmers]) 
        self.max_characters = len(self.alphabet)
        self.kmers = kmers
        
        self.goto = [[-1]*self.max_characters for _ in range(self.max_states+1)]
        #self.goto = [[-1] * self.max_characters] * (self.max_characters+1)
        self.fail = [-1] * (self.max_states+1)
        self.out = [0] * (self.max_states+1)

        self.states_count = self.__build()

    def __build(self):
        states = 1
        
        for i, kmer in enumerate(self.kmers):
            curr_state = 0
            for nucleotide in kmer:
                index_in_alphabet = self.alphabet.rfind(nucleotide)
                if self.goto[curr_state][index_in_alphabet] == -1:
                    self.goto[curr_state][index_in_alphabet] = states
                    states = states + 1
                curr_state = self.goto[curr_state][index_in_alphabet]

            self.out[curr_state] |= (1<<i)

        for i in range(self.max_characters):
            if self.goto[0][i] == -1:
                self.goto[0][i] = 0

        # For self.fail
        # All nodes of depth 1 have failure value as 0.
        queue = []
        for i in range(self.max_characters):
            if self.goto[0][i] != 0:
                self.fail[self.goto[0][i]] = 0
                queue.append(self.goto[0][i])
        while queue:
            state = queue.pop(0)
            for i in range(self.max_characters):
                if self.goto[state][i] != -1:
                    failure = self.fail[state]
                    while self.goto[failure][i] == -1:
                        failure = self.fail[failure]
                    failure=self.goto[failure][i]
                    self.fail[self.goto[state][i]] = failure
                    self.out[self.goto[state][i]] |= self.out[failure]
                    queue.append(self.goto[state][i])
        return states

    def __next_state(self, curr_state, next_nucleotide):
        state = curr_state
        index_in_alphabet = self.alphabet.rfind(next_nucleotide)

        while self.goto[state][index_in_alphabet] == -1:
            state = self.fail[state]

        return self.goto[state][index_in_alphabet]

    def start_indexes(self, sequence):
        start_indexes = []
        curr_state = 0
        
        for i, nucleotide in enumerate(sequence):
            curr_state = self.__next_state(curr_state, nucleotide)

            if self.out[curr_state] == 0:
                continue
        
            for j in range(len(self.kmers)):
                if self.out[curr_state] > 0 and (1<<j) > 0:
                    kmer = self.kmers[j]

                    start_index = i - len(kmer) + 1
                    start_indexes.append(start_index)
        print(start_indexes)
        return start_indexes
