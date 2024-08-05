
from Bio import SeqIO
from .genbank_seq_object import GenBankSeq

class GenBank:

    def __init__(self, filename):
        self.filename = filename
        self.seqs = []
        self._read_genbank_file(filename)

    def _read_genbank_file(self, filename):
        genbank_sequences = SeqIO.parse(open(filename, "r"), "genbank")
        for genbank in genbank_sequences:
            seq = GenBankSeq(genbank)
            self.seqs.append(seq)

