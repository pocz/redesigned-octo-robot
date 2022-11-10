'''
Does preprocessing on a list of k-mers. Afterwards sequences can be input into
the search function, which will look for occurences of the preprocessed k-mers.
''' 
class Automaton:
    def __init__(self, kmers, seq_type):
        if seq_type == "DNA":
            self.alphabet = "AGTCN"
        elif seq_type == "RNA":
            self.alphabet = "AGUCN"

        self.kmers = kmers
        self.max_characters = len(self.alphabet)
    
        self.max_states = 1
        for kmer in kmers:
            self.max_states = self.max_states + len(kmer)
        
        ''' goto stores the trie. It's a list of lists. Every state has a list
        containing directions for every character of the kmer. Eg.: when the
        alphabet is "AGTCN" then if the first member in the list is [1,5,0,0,0]
        it tells us to got to state one in case of A and to state 5 in case of
        G and stay on state 0 for other characters. It is initialized with -1
        values.
            fail stores values for each state. When there is no direction for
        the character (eg. go to is -1) the fail list tells us if there's a
        direction for the same character in the root. If there is we can skip
        to that state.
            output stores for each state,
        the indexes of all words ending at the state as bitmaps '''
        self.goto = [[-1]*self.max_characters for i in range(self.max_states+1)]
        self.fail = [-1] * (self.max_states+1)
        self.out = [0] * (self.max_states+1)

        self.states_count = self.__build()
        #print(f'GO TO list:\n{self.goto}\nFAIL list:\n{self.fail}\nOUTPUT list:\n{self.out}\n')

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

        # Queue for self.fail
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
        start_indexes = {}
        for kmer in self.kmers:
            start_indexes[kmer] = []   
        curr_state = 0
        for i, nucleotide in enumerate(sequence):
            curr_state = self.__next_state(curr_state, nucleotide)

            if self.out[curr_state] != 0:
                for j, kmer in enumerate(self.kmers):
                    if (self.out[curr_state] & (1<<j)) > 0:
                        kmer = self.kmers[j]

                        start_index = i - len(kmer) + 1
                        start_indexes[kmer].append(start_index)

        return start_indexes
