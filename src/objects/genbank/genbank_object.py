
from Bio import SeqIO
from .genbank_seq_object import GenBankSeq

class GenBank:

    def __init__(self, filename):
        """ Initializes a GenBank object reading the GenBank
            file and initializing its subobjects (sequences, 
            genes and condign regions)
        """
        self.filename = filename
        self.seqs = []
        self._read_genbank_file(filename)

    def _read_genbank_file(self, filename):
        """ Reads the genbank file and creates an GenBankSeq 
            object for each sequence contained on the file
        """
        genbank_sequences = SeqIO.parse(open(filename, "r"), "genbank")
        for genbank in genbank_sequences:
            seq = GenBankSeq(genbank)
            self.seqs.append(seq)

