class GenBankCDS:

    def __init__(self, info):
        """ Initializes a cds object that represents a coding region 
            using the information obtained from GenBank file (info)
        """
        self.info = info
        self.location = info.location
        self.codon_start = info.qualifiers["codon_start"]
        self.strand = "+" if (self.location.strand==1) else "-"
        if "gene" in info.qualifiers.keys():
            self.gene = info.qualifiers["gene"]
        else:
            self.gene = None
        if "locus_tag" in info.qualifiers.keys():
            self.tag = info.qualifiers["locus_tag"][0]
        else:
            self.tag = None
        if "protein_id" in info.qualifiers.keys():
            self.protein_id = info.qualifiers["protein_id"][0]
        else:
            self.protein_id = None
        if "translation" in info.qualifiers.keys():
            self.trans = info.qualifiers["translation"]
        else:
            self.trans = None
        if "product" in info.qualifiers.keys():
            self.product = info.qualifiers["product"][0]
        else:
            self.product = None
