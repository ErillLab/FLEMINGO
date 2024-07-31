import random
from Markov_DNA import MCM

class Generator:

    KMER = 0
    MUT = 1
    MCM = 2
    GC = 3
    RANDOM = 4

    NUCLEOTIDES = "ACGT"

    def __init__(self, method):
        self.method = method

    def generate(self, N, n=None, seq=None, k=None, p=0.005):
        if self.method == self.KMER:
            return self.generate_kmer(seq, N, k)
        elif self.method == self.MUT:
            return self.generate_mut(seq, N, p)
        elif self.method == self.MCM:
            return self.generate_mcm(seq, N, k)
        elif self.method == self.GC:
            return self.generate_gc(seq, N)
        elif self.method == self.RANDOM:
            return self.generate_random(N, n)
        else:
            print("ERROR: INVALID SEQUENCE GENERATION METHOD")
            exit(-1)

    def generate_kmer(self, seq, N, k):
        fragments = []
        seqs = []
        for i in range(0, len(seq), k):
            if i+k < len(seq):
                fragments.append(seq[i:i+k])
            else:
                fragments.append(seq[i:])

        for i in range(N):
            new_seq = ""
            aux = fragments[:]
            for j in range(len(fragments)):
                f = random.randint(0, len(aux)-1)
                new_seq += aux[f]
            seqs.append(new_seq)

        return seqs

    def generate_mut(self, seq, N, p):
        fragments = []
        seqs = []

        for i in range(N):
            new_seq = ""
            aux = fragments[:]
            for j in range(len(seq)):
                if (random.random() <= p):
                    pos = random.randint(0, 3)
                    new_seq += self.NUCLEOTIDES[pos]
                else:
                    new_seq += seq[j]
            seqs.append(new_seq)

        return seqs

    def generate_mcm(self, seq, N, k):
        mcm = MCM(k, mode=MCM.CIRCULAR)
        mcm.train(seq)
        seqs = mcm.generate(len(seq), N)
        return seqs

    def generate_gc(self, seq, N):
        gc = 0
        for nuc in seq:
            if nuc == "C" or nuc == "c" or nuc == "G" or nuc == "g":
                gc += 1

        gc /= len(seq)

        seqs = []
        for n in range(N):
            new_seq = ""
            for i in range(len(seq)):
                r = random.random()
                if r < gc:
                    new_seq += self.NUCLEOTIDES[random.randint(0,1)+1]
                else:
                    new_seq += self.NUCLEOTIDES[random.randint(0,1)*-1]
            seqs.append(seq)
        return seqs

    def generate_random(self, N, n):
        seqs = []
        for _ in range(N):
            new_seq = ""
            for _ in range(n):
                r = random.randint(0, 3)
                new_seq += self.NUCLEOTIDES[r]
            seqs.append(new_seq)
        return seqs