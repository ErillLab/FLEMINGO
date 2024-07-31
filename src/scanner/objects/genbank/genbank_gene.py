

class GenBankGene:

    def __init__(self, feature):
        self.location = feature.location
        self.strand = "+" if (self.location.strand==1) else "-"
        if "gene" in feature.qualifiers.keys():
            self.name = feature.qualifiers["gene"][0]
        else:
            self.name = None
        if "locus_tag" in feature.qualifiers.keys():
            self.tag = feature.qualifiers["locus_tag"][0]
        else:
            self.tag = None
        if "gene_synonym" in feature.qualifiers.keys():
            self.synonym = feature.qualifiers["gene_synonym"][0]
        else:
            self.synonym = None
        self.product = None
        self.protein_id = None

    def print(self):
        print("Name:", self.name)
        print("TAG:", self.tag)
        print("Location: ", self.location.start, "-", self.location.end)
        print("Strand: ", self.strand)
