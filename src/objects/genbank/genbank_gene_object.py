class GenBankGene:

    def __init__(self, info):
        """ Initializes a gene using the information obtained 
            from GenBank file (info)
        """
        self.location = info.location
        self.strand = "+" if (self.location.strand==1) else "-"
        if "gene" in info.qualifiers.keys():
            self.name = info.qualifiers["gene"][0]
        else:
            self.name = None
        if "locus_tag" in info.qualifiers.keys():
            self.tag = info.qualifiers["locus_tag"][0]
        else:
            self.tag = None
        if "gene_synonym" in info.qualifiers.keys():
            self.synonym = info.qualifiers["gene_synonym"][0]
        else:
            self.synonym = None
        self.product = None
        self.protein_id = None
        self.score = None
        self.p_value = None

    def print(self):
        """ Prints gene information
        """
        print("Name:", self.name)
        print("TAG:", self.tag)
        print("Location: ", self.location.start, "-", self.location.end)
        print("Strand: ", self.strand)
