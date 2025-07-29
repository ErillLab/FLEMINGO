"""
Sequence object stores a DNA sequence, it has two strands and a length.
It can return the reverse complement of the DNA sequence.
"""


from Bio.Seq import Seq


class SeqObject():
    
    def __init__(self, sequence: str):
        
        # "forward" (f) and "reverse" (r) strands
        self.f = sequence
        self.r = self.reverse_complement(sequence)
        # Length (in base-pairs)
        self.len = len(sequence)
        
    def reverse_complement(self, sequence):
        return str(Seq(sequence).reverse_complement())
    





